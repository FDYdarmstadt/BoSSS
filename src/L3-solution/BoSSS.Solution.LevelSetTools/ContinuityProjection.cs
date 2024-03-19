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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.ConstrainedDGprojection;
using ilPSP.Tracing;

namespace BoSSS.Solution.LevelSetTools {

    /// <summary>
    /// Options for enforcing the continuity of the level-set function
    /// </summary>
    public enum ContinuityProjectionOption {

        /// <summary>
        /// Do not perform continuity projection
        /// </summary>
        None,

        /// <summary>
        /// projection on a spectral finite element field
        /// </summary>
        SpecFEM,

        /// <summary>
        /// L2-projection with continuity constraints at inner cell boundaries
        /// </summary>
        ConstrainedDG
    }


    /// <summary>
    /// Projects a DG Field onto another DGField of higher order, without any discontinuities at the boundary.
    /// </summary>
    public class ContinuityProjection {
        
        /// <summary>
        /// Selects algorithm variant based on the <paramref name="Option"/>
        /// </summary>
        /// <param name="ContBasis">DG basis for the continuous output data</param>
        /// <param name="DGBasis">DG basis for the discontinuous input data</param>
        /// <param name="gridData"></param>
        /// <param name="Option">Choice of algorithm</param>
        public ContinuityProjection(Basis ContBasis, Basis DGBasis, IGridData gridData, ContinuityProjectionOption Option) {
            myOption = Option;
            switch (Option) {
                case ContinuityProjectionOption.SpecFEM: {
                        var ContinuousLevelSetBasis = new SpecFemBasis((GridData)gridData, DGBasis.Degree + 1);
                        MyProjection = new ContinuityProjectionSpecFem(ContinuousLevelSetBasis);
                        break;
                    }
                case ContinuityProjectionOption.ConstrainedDG: {
                        MyProjection = new ContinuityProjectionCDG(ContBasis);
                        break;
                    }
                case ContinuityProjectionOption.None: {
                        MyProjection = new NoProjection();
                        break;
                    }
                default:
                    throw new ArgumentException();
            }

        }

        /// <summary>
        /// Creates the LevelSet field for the **continuous representation**,
        /// for the Tracker based on the <paramref name="Option"/> selected
        /// </summary>
        /// <param name="DGLevelSet">the unfiltered Level-set (might be discontinuous)</param>
        /// <param name="gridData"></param>
        /// <param name="Option"></param>
        /// <returns>The Level-Set Field for storing the filtered, i.e. continuous representation.</returns>
        public static LevelSet CreateField(SinglePhaseField DGLevelSet, Foundation.Grid.Classic.GridData gridData, ContinuityProjectionOption Option) {
            int k = DGLevelSet.Basis.Degree;
            switch (Option) {
                case ContinuityProjectionOption.SpecFEM: {
                        var ContinuousLevelSetBasis = new SpecFemBasis(gridData, k + 1);
                        return new LevelSet(ContinuousLevelSetBasis.ContainingDGBasis, VariableNames.LevelSetCG);
                    }
                case ContinuityProjectionOption.ConstrainedDG: {
                        var ContinuousLevelSetDGBasis = new Basis(gridData, k + 1);
                        return new LevelSet(ContinuousLevelSetDGBasis, VariableNames.LevelSetCG);
                    }
                case ContinuityProjectionOption.None: {
                        Console.WriteLine("WARNING: No additional enforcement of the level-set continuity!");
                        LevelSet SmoothedLevelSet = new LevelSet(DGLevelSet.Basis, VariableNames.LevelSetCG);
                        return SmoothedLevelSet;
                    }
                default:
                    throw new ArgumentException();
            }
        }


        IContinuityProjection MyProjection;

        ContinuityProjectionOption myOption;

        /// <summary>
        /// Makes <paramref name="DGLevelSet"/> a continuous function <paramref name="LevelSet"/>,
        /// according to the Option from the initialization
        /// </summary>
        /// <param name="DGLevelSet">input; may be discontinuous on cell boundaries</param>
        /// <param name="LevelSet">output; should be continuous</param>
        /// <param name="Domain"></param>
        /// <param name="PosMask">
        /// If the FAR-field is set to a constant value (<paramref name="setFarFieldConstant"/> is true), this determines which cells belong to the positive part of the FAR-field;
        /// The negative part is determined by a complement-operation (<see cref="ExecutionMask.Complement{T}"/>)
        /// </param>
        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain, CellMask PosMask, bool setFarFieldConstant = true) {
            using(var ft = new FuncTrace()) {
                //Console.WriteLine("calling ContinuityProjection.MakeContinuous() ...");

                MyProjection.MakeContinuous(DGLevelSet, LevelSet, Domain);

                double Jnorm = JumpNorm(LevelSet, Domain);
                ft.Info($"jump norm after continuity projection = {Jnorm}");

                if(myOption != ContinuityProjectionOption.None && setFarFieldConstant) {
                    SetFarField(LevelSet, Domain, PosMask);
                }
            }
        }

        /// <summary>
        /// Sets the positive Far-field of the level-set to +1 and the negative side to -1 
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="Domain"></param>
        /// <param name="PosMask"></param>
        public void SetFarField(SinglePhaseField LevelSet, CellMask Domain, CellMask PosMask) {
            // set cells outside narrow band to +/- 1
            var Pos = PosMask.Except(Domain);
            var Neg = PosMask.Complement().Except(Domain);

            LevelSet.Clear(Pos);
            LevelSet.AccConstant(1, Pos);

            LevelSet.Clear(Neg);
            LevelSet.AccConstant(-1, Neg);
        }

        /// <summary>
        /// computes the jump norm of field <paramref name="f"/> at inner edges on <paramref name="mask"/>
        /// </summary>
        /// <param name="f"></param>
        /// <param name="mask"> if null full mask is chosen </param>
        /// <returns></returns>
        static double JumpNorm(DGField f, CellMask mask = null) {
            GridData grd = (GridData)f.GridDat;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            if (mask == null) {
                mask = CellMask.GetFullMask(grd);
            }
            SubGrid maskSG = new SubGrid(mask);
            EdgeMask innerEM = maskSG.InnerEdgesMask;

            return JumpNorm(f, innerEM);
        }


        /// <summary>
        /// computes the jump norm of field <paramref name="f"/> at inner edges on <paramref name="mask"/>
        /// </summary>
        /// <param name="f"></param>
        /// <param name="mask"> if null full mask is chosen </param>
        /// <returns></innerEM>
        static double JumpNorm(DGField f, EdgeMask innerEM = null) {
            GridData grd = (GridData)f.GridDat;
            int D = grd.SpatialDimension;
            var e2cTrafo = grd.Edges.Edge2CellTrafos;

            if (innerEM == null) {
                innerEM = EdgeMask.GetFullMask(grd, MaskType.Geometrical);
            }
            

            f.MPIExchange();

            double Unorm = 0;

            EdgeQuadrature.GetQuadrature(
                new int[] { D + 1 }, grd,
                (new EdgeQuadratureScheme(true, innerEM)).Compile(grd, f.Basis.Degree * 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) { // Evaluate
                    NodeSet NS = QR.Nodes;
                    EvalResult.Clear();
                    int NoOfNodes = NS.NoOfNodes;
                    for (int j = 0; j < Length; j++) {
                        int iEdge = j + i0;
                        int jCell_IN = grd.Edges.CellIndices[iEdge, 0];
                        int jCell_OT = grd.Edges.CellIndices[iEdge, 1];
                        var uDiff = EvalResult.ExtractSubArrayShallow(new int[] { j, 0, 0 }, new int[] { j, NoOfNodes - 1, -1 });

                        if (jCell_OT >= 0) {

                            int iTrafo_IN = grd.Edges.Edge2CellTrafoIndex[iEdge, 0];
                            int iTrafo_OT = grd.Edges.Edge2CellTrafoIndex[iEdge, 1];

                            MultidimensionalArray uIN = MultidimensionalArray.Create(1, NoOfNodes);
                            MultidimensionalArray uOT = MultidimensionalArray.Create(1, NoOfNodes);

                            NodeSet NS_IN = NS.GetVolumeNodeSet(grd, iTrafo_IN, false);
                            NodeSet NS_OT = NS.GetVolumeNodeSet(grd, iTrafo_OT, false);

                            f.Evaluate(jCell_IN, 1, NS_IN, uIN);
                            f.Evaluate(jCell_OT, 1, NS_OT, uOT);

                            uDiff.Acc(+1.0, uIN);
                            uDiff.Acc(-1.0, uOT);

                            //if (uDiff.L2Norm() > 1e-10) {
                            //    Console.WriteLine("uDiff at edge {0} between cell {1} and cell {2}: {3}", iEdge, jCell_IN, jCell_OT, uDiff.L2Norm());
                            //}
                        } else {
                            uDiff.Clear();
                        }
                    }

                    EvalResult.ApplyAll(x => x * x);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    Unorm += ResultsOfIntegration.Sum();
                }).Execute();

            Unorm = Unorm.MPISum();

            return Unorm.Sqrt();
        }


    }

    /// <summary>
    /// Driver interface
    /// </summary>
    interface IContinuityProjection {
        void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain);
    }

    // ===============================
    // The actual implementations
    // ===============================

    ///<summary>
    /// Smoothing based on SpecFem
    /// => Actually ContinuousFunctionSpace
    ///</summary>
    class ContinuityProjectionSpecFem : IContinuityProjection {

        public ContinuityProjectionSpecFem(SpecFemBasis myBasis) {
            FEMLevSet = new SpecFemField(myBasis);
        }
        SpecFemField FEMLevSet;

        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {
            FEMLevSet.ProjectDGField(1.0, DGLevelSet, Domain);
            LevelSet.Clear();
            FEMLevSet.AccToDGField(1.0, LevelSet, Domain);
            LevelSet.AccLaidBack(1.0, DGLevelSet, Domain.Complement());
        }
    }

    /// <summary>
    /// Smoothing based on <see cref="ConstrainedDGFieldMk3"/> 
    /// => Lagrange-Multiplier Approach
    /// </summary>
    /// <remarks>
    /// - Developed by Martin Smuda, described in his PhD thesis
    /// - the main advantage of the continuity projection using <see cref="ContinuityProjectionCDG"/>
    ///   is that it works with hanging nodes and in 3D.
    /// </remarks>
    public class ContinuityProjectionCDG : IContinuityProjection {

        public ContinuityProjectionCDG(Basis myBasis) {
            m_myBasis = myBasis;
        }
        Basis m_myBasis;

        ConstrainedDGFieldMk3 m_projector;

        ConstrainedDGFieldMk3 Update(CellMask domain) {
            
            bool DomainChanged = m_projector != null ? (!m_projector.domainLimit.Equals(domain)).MPIOr() : false;
            if (m_projector == null || DomainChanged) {
                if(m_projector != null)
                    m_projector.Dispose();
                m_projector = new ConstrainedDGField_Global(m_myBasis, domain);
                //m_projector = new ConstrainedDgField_Patchwise(m_myBasis, domain);

            }
            //m_projector = new ConstrainedDGField_Global(m_myBasis, domain);

            return m_projector;
        }
        
        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {
            if(Domain.NoOfItemsLocally.MPISum() > 0) {
                var p = Update(Domain);
                p.ProjectDGField(DGLevelSet);
                LevelSet.Clear();
                p.AccToDGField(1.0, LevelSet);
            } else {
                LevelSet.Clear();
                LevelSet.AccLaidBack(1.0, DGLevelSet);
            }
        }
    }

    /// <summary>
    /// Does nothing => no enforcement of continuity  
    /// </summary>
    class NoProjection : IContinuityProjection {

        public void MakeContinuous(SinglePhaseField DGLevelSet, SinglePhaseField LevelSet, CellMask Domain) {

            if (ReferenceEquals(DGLevelSet, LevelSet)) {
                // Do nothing
            } else {
                LevelSet.Clear();
                LevelSet.AccLaidBack(1.0, DGLevelSet);
            }

        }
    }

}
