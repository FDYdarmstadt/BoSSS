/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Solution.Timestepping;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Platform;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.XNSECommon {
    

    /// <summary>
    /// Projection of XDG Velocity to DG Space
    /// </summary>
    public class XVelocityProjection {

        /// <summary>
        /// Options for projecting an XDG Field in a cut cell to a DG Field
        /// </summary>
        public enum CutCellVelocityProjectiontype {
            /// <summary>
            /// L2-projection weighted by factors for each phase
            /// </summary>
            L2_weighted,

            /// <summary>
            /// L2-Projection for the whole cell
            /// </summary>
            L2_plain,

            /// <summary>
            /// Exact projection at the interface and the boundaries of the cell
            /// </summary>
            BoundaryAndInterface
        }

        /// <summary>
        /// Implement this class in a Control File to handle the Velocity-Projection
        /// </summary>
        public class Configuration {

            /// <summary>
            /// Projection-Type - Default: Exact at Boundaries and Interface
            /// </summary>
            public CutCellVelocityProjectiontype CutCellVelocityProjectiontype = CutCellVelocityProjectiontype.BoundaryAndInterface;

            /// <summary>
            /// filtering the filtered velocity once more with patch recovery.
            /// </summary>
            public bool UsePatchRecoveryFiltering = false;

            internal IDictionary<SpeciesId, double> ScalingFactors;
        }

        /// <summary>
        /// configuration for the velocity filtering, see <see cref="Configuration"/>
        /// </summary>
        public Configuration Config = new Configuration();
                

        /// <summary>
        /// Do Projection from XDG to DG
        /// </summary>
        /// <param name="Velocity">Input XDG Velocity</param>
        /// <param name="FilteredVelocity">Output DG Velocity</param>
        public XVelocityProjection(LevelSetTracker _LsTrk, DGField[] Velocity, VectorField<SinglePhaseField> FilteredVelocity) {
            this.m_FilteredVelocity = FilteredVelocity;
            this.LsTrk = _LsTrk;
            
            Config.ScalingFactors = new Dictionary<SpeciesId, double>();
            foreach (var kv in LsTrk.SpeciesIdS.Select(spc => new KeyValuePair<SpeciesId, double>(spc, 0.5)))
                Config.ScalingFactors.Add(kv);
            
            m_CC0 = LsTrk.Regions.GetNearFieldMask(0);
            if (this.Config.UsePatchRecoveryFiltering) {
                UpdatePatchRecoveryFilter();
            }
            this.FilterVelocity(Velocity, FilteredVelocity);
        }

        LevelSetTracker LsTrk;

        /// <summary>
        /// Update NearFieldMask
        /// </summary>
        public void Push() {
            this.m_CC0 = this.LsTrk.Regions.GetNearFieldMask(0);
        }



        /// <summary>
        /// Projects the 
        /// </summary>
        /// <param name="XVelocity"> Velocity in XDG Space</param>
        /// onto the 
        /// <param name="m_FilteredVelocity">Velocity in DG Space</param>
        public void Perform(DGField[] XVelocity, VectorField<SinglePhaseField> m_FilteredVelocity) {

            // update patch-recovery filter
            // ----------------------------
            UpdatePatchRecoveryFilter();

            // filter velocity
            // ---------------
            this.FilterVelocity(XVelocity, this.m_FilteredVelocity);

        }

        private void UpdatePatchRecoveryFilter() {
            Basis B = m_FilteredVelocity[0].Basis;
            var CC1 = LsTrk.Regions.GetNearFieldMask(0);

            if (this.Config.UsePatchRecoveryFiltering) {
                var filterDom = m_CC0.Union(CC1);
                if (this.Filter == null)
                    this.Filter = new L2PatchRecovery(B, B, filterDom, RestrictToCellMask: false);
                else
                    this.Filter.UpdateDomain(filterDom, RestrictToCellMask: false);

            } else {
                this.Filter = null;
            }
        }

        CellMask m_CC0;


        VectorField<SinglePhaseField> m_FilteredVelocity;

        L2PatchRecovery Filter;


        int FindQuadRuleOrder() {
            int a = this.m_FilteredVelocity[0].Basis.Degree * 2;

            // a hack, to provide a sufficient quad rule order
            a++;
            if (a % 2 != 0)
                a++;
            if (a % 4 != 0)
                a += 2;
            Debug.Assert(a % 4 == 0);

            //a = 8;
            //Console.WriteLine("running on order " + a);

            return a;
        }


        void SpecialProjection(LevelSetTracker LsTrk, int order, DGField[] Uin, VectorField<SinglePhaseField> Uout) {
            var CC = LsTrk.Regions.GetCutCellMask();
            var gDat = LsTrk.GridDat;

            // get quadrature rules
            // ====================

            //XQuadSchemeHelper H = new XQuadSchemeHelper(LsTrk, MomentFittingVariant);
            var H = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = H.GetLevelSetquadScheme(0, CC);
            ICompositeQuadRule<QuadRule> surfRule = cqs.Compile(gDat, order);
            ICompositeQuadRule<CellBoundaryQuadRule> bndyRule = (new CellBoundaryQuadratureScheme(true, CC)).Compile(gDat, order);



            // Compute Mass matrix and RHS for the 'strange' projection
            // ========================================================
            int L = CC.NoOfItemsLocally;
            var Q = new QuadratureKernels(Uin, L);

            Q.BlockCnt = 0;
            CellQuadrature.GetQuadrature(
                new int[] { Q.Nnx + Q.D, Q.Nnx },
                gDat,
                surfRule,
                Q.Evaluate, Q.SaveIntegrationResults_surf).Execute();
            Debug.Assert(Q.BlockCnt == L);

            Q.BlockCnt = 0;
            CellBoundaryQuadrature<CellBoundaryQuadRule>.GetQuadrature(
                 new int[] { Q.Nnx + Q.D, Q.Nnx },
                gDat,
                bndyRule,
                Q.Evaluate, Q.SaveIntegrationResults_bndy).Execute();
            Debug.Assert(Q.BlockCnt == L);

            // solve the non-diagonal mass matrix systems
            // ==========================================

            int BlkCnt = 0;
            foreach (int jCell in CC.ItemEnum) {
                var MassMatrix = Q.MassMatrix.ExtractSubArrayShallow(BlkCnt, -1, -1);
                var RHS = Q.RHS.ExtractSubArrayShallow(BlkCnt, -1, -1);

                // Die "Massenmatrix" muss nicht unbedingt invbar sein, daher: Least-Squares solve
                MassMatrix.LeastSquareSolve(RHS);

                for (int d = 0; d < Q.D; d++) {
                    for (int n = 0; n < Q.Nnx; n++) {
                        Uout[d].Coordinates[jCell, n] = RHS[n, d];
                    }
                }

                BlkCnt++;
            }
        }

        class QuadratureKernels {
            public QuadratureKernels(DGField[] _Uin, int L) {
                this.Uin = _Uin;
                if (_Uin[0] is XDGField) {
                    this.b = ((XDGField)_Uin[0]).Basis.NonX_Basis;
                } else {
                    this.b = _Uin[0].Basis;
                }
                Nnx = b.Length;
                gDat = (GridData)b.GridDat;
                D = gDat.SpatialDimension;
                MassMatrix = MultidimensionalArray.Create(L, Nnx, Nnx);
                RHS = MultidimensionalArray.Create(L, Nnx, D);
                UevalBuf = null;
            }

            Basis b;
            public int Nnx;
            GridData gDat;
            public int D;
            public MultidimensionalArray MassMatrix;
            public MultidimensionalArray RHS;
            public MultidimensionalArray UevalBuf;
            public int BlockCnt = 0;
            DGField[] Uin;

            public void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                // Del_Evaluate
                // ~~~~~~~~~~~~~
                NodeSet NS = QR.Nodes;
                int NoOfNodes = NS.NoOfNodes;
                var BasisVal = b.CellEval(NS, i0, Length);
                EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, Nnx - 1, Nnx - 1 })
                    .Multiply(1.0, BasisVal, BasisVal, 0.0, "ikmn", "ikm", "ikn");

                if (UevalBuf == null)
                    UevalBuf = MultidimensionalArray.Create(Length, NoOfNodes);
                if (UevalBuf.GetLength(0) < Length || UevalBuf.GetLength(1) != NoOfNodes) {
                    UevalBuf.Allocate(Length, NoOfNodes);
                }
                MultidimensionalArray _UevalBuf;
                if (UevalBuf.GetLength(0) > Length) {
                    _UevalBuf = UevalBuf.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Length - 1, NoOfNodes - 1});
                } else {
                    Debug.Assert(UevalBuf.GetLength(0) == Length);
                    _UevalBuf = UevalBuf;
                }

                for (int d = 0; d < D; d++) {
                    _UevalBuf.Clear();
                    Uin[d].Evaluate(i0, Length, NS, _UevalBuf);

                    EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, Nnx + d, 0 }, new int[] { Length - 1, NoOfNodes - 1, Nnx + d - 1, Nnx - 1 })
                        .Multiply(1.0, BasisVal, _UevalBuf, 0.0, "ikm", "ikm", "ik");
                }

            }
            public void SaveIntegrationResults_surf(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                // Del_SaveIntegrationResults
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                var TempBlock = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Length - 1, Nnx - 1, Nnx - 1 });
                var StorBlock = MassMatrix.ExtractSubArrayShallow(new int[] { BlockCnt, 0, 0 }, new int[] { BlockCnt + Length - 1, Nnx - 1, Nnx - 1 });
                StorBlock.Acc(1.0, TempBlock);

                var TempRhs = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { 0, Nnx, 0 }, new int[] { Length - 1, Nnx + D - 1, Nnx - 1 });
                var StorRhs = RHS.ExtractSubArrayShallow(new int[] { BlockCnt, 0, 0 }, new int[] { BlockCnt + Length - 1, Nnx - 1, D - 1, });
                //StorRhs.Acc(1.0, TempRhs);
                for (int i = 0; i < Length; i++) {
                    for (int d = 0; d < D; d++) {
                        for (int n = 0; n < Nnx; n++) {
                            StorRhs[i, n, d] += TempRhs[i, d, n];
                        }
                    }
                }

                BlockCnt += Length;
            }
            public void SaveIntegrationResults_bndy(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                // Del_SaveIntegrationResults
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~

                int NoFaces = ResultsOfIntegration.GetLength(1);
                for (int iFace = 0; iFace < NoFaces; iFace++) {
                    var TempBlock = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { 0, iFace, 0, 0 }, new int[] { Length - 1, iFace - 1, Nnx - 1, Nnx - 1 });
                    var StorBlock = MassMatrix.ExtractSubArrayShallow(new int[] { BlockCnt, 0, 0 }, new int[] { BlockCnt + Length - 1, Nnx - 1, Nnx - 1 });
                    StorBlock.Acc(1.0, TempBlock);

                    var TempRhs = ResultsOfIntegration.ExtractSubArrayShallow(new int[] { 0, iFace, Nnx, 0 }, new int[] { Length - 1, iFace - 1, Nnx + D - 1, Nnx - 1 });
                    var StorRhs = RHS.ExtractSubArrayShallow(new int[] { BlockCnt, 0, 0 }, new int[] { BlockCnt + Length - 1, Nnx - 1, D - 1, });
                    //StorRhs.Acc(1.0, TempRhs);
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            for (int n = 0; n < Nnx; n++) {
                                StorRhs[i, n, d] += TempRhs[i, d, n];
                            }
                        }
                    }
                }

                BlockCnt += Length;
            }
        }


        void FilterVelocity(DGField[] Velocity, VectorField<SinglePhaseField> FilteredVelocity) {
            var gDat = this.LsTrk.GridDat;
            int D = gDat.SpatialDimension;

            FilteredVelocity.Clear();

            // Projection in non-cut cells
            // ===========================
            foreach (var species in LsTrk.SpeciesIdS) {
                var mask = LsTrk.Regions.GetSpeciesMask(species).Except(LsTrk.Regions.GetCutCellMask());

                for (int d = 0; d < D; d++) {
                    ConventionalDGField Vel_d;
                    if (Velocity[d] is XDGField) {
                        Vel_d = ((XDGField)(Velocity[d])).GetSpeciesShadowField(species);
                    } else {
                        Vel_d = ((ConventionalDGField)(Velocity[d]));
                    }

                    FilteredVelocity[d].AccLaidBack(1.0, Vel_d, mask);
                }
            }

            // Projection in cut cells
            // =======================

            switch (this.Config.CutCellVelocityProjectiontype) {
                case CutCellVelocityProjectiontype.L2_weighted:
                {
                    foreach (var species in LsTrk.SpeciesIdS) {
                        var mask = LsTrk.Regions.GetCutCellMask();
                        for (int d = 0; d < D; d++) {
                            ConventionalDGField Vel_d;
                            if (Velocity[d] is XDGField) {
                                Vel_d = ((XDGField)(Velocity[d])).GetSpeciesShadowField(species);
                            } else {
                                Vel_d = ((ConventionalDGField)(Velocity[d]));
                            }

                            FilteredVelocity[d].Acc(Config.ScalingFactors[species], Vel_d, mask);
                        }
                    }
                    break;
                }

                case CutCellVelocityProjectiontype.L2_plain:
                {

                    var mask = LsTrk.Regions.GetCutCellMask();
                    for (int d = 0; d < D; d++) {
                        FilteredVelocity[d].ProjectField(1.0, Velocity[d].Evaluate, new CellQuadratureScheme(domain: mask));
                    }

                    break;
                }

                case CutCellVelocityProjectiontype.BoundaryAndInterface:
                {

                    this.SpecialProjection(LsTrk, this.FindQuadRuleOrder(), Velocity, FilteredVelocity);

                    break;
                }

                default:
                throw new NotSupportedException();

            }

            // optional Patch-Recovery filtering
            // =================================

            if (this.Config.UsePatchRecoveryFiltering) {
                for (int d = 0; d < D; d++) {
                    var U_b4Filter = FilteredVelocity[d].CloneAs();
                    var U = FilteredVelocity[d];
                    this.Filter.Perform(U, U_b4Filter);
                }
            }
        }

    }

    
}
