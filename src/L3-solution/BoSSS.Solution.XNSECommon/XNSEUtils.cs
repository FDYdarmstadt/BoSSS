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
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using BoSSS.Platform;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.XNSECommon {

    public static class XNSEUtils {

        /// <summary>
        /// modifies a matrix <paramref name="Mtx"/> and a right-hand-side <paramref name="rhs"/>
        /// in order to fix the pressure at some reference point
        /// </summary>
        /// <param name="map">row mapping for <paramref name="Mtx"/> as well as <paramref name="rhs"/></param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="LsTrk"></param>
        /// <param name="Mtx"></param>
        /// <param name="rhs"></param>
        static public void SetPressureReferencePoint<T>(UnsetteledCoordinateMapping map, int iVar, LevelSetTracker LsTrk, IMutableMatrixEx Mtx, T rhs)
            where T : IList<double> //
        {
            using (new FuncTrace()) {
                var GridDat = map.GridDat;


                if (rhs.Count != map.LocalLength)
                    throw new ArgumentException();
                if (!Mtx.RowPartitioning.Equals(map) || !Mtx.ColPartition.Equals(map))
                    throw new ArgumentException();

                XDGBasis PressureBasis = (XDGBasis)map.BasisS[iVar];
                var grd = GridDat;
                int D = GridDat.SpatialDimension;

                long GlobalID, GlobalIndex;
                bool IsInside, onthisProc;
                double[] pt = (D == 2) ? new double[] { -5, -5 } : new double[] { -5.0, -5.0, -5.0 };
                grd.LocatePoint(pt, out GlobalID, out GlobalIndex, out IsInside, out onthisProc, LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement());


                int iRowGl = -111;

                if (onthisProc) {
                    int jCell = (int)GlobalIndex - GridDat.CellPartitioning.i0;
                    NodeSet CenterNode = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell), new double[D]);

                    MultidimensionalArray LevSetValues = LsTrk.DataHistories[0].Current.GetLevSetValues(CenterNode, jCell, 1); ;

                    MultidimensionalArray CenterNodeGlobal = MultidimensionalArray.Create(1, D);
                    GridDat.TransformLocal2Global(CenterNode, CenterNodeGlobal, jCell);
                    //Console.WriteLine("Pressure Ref Point @( {0:0.###E-00} | {1:0.###E-00} )", CenterNodeGlobal[0,0], CenterNodeGlobal[0,1]);

                    LevelSetSignCode scode = LevelSetSignCode.ComputeLevelSetBytecode(LevSetValues[0, 0]);
                    ReducedRegionCode rrc;
                    int No = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);
                    int iSpc = LsTrk.GetSpeciesIndex(rrc, scode);

                    iRowGl = (int)map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalIndex, PressureBasis.DOFperSpeciesPerCell * iSpc);
                }

                unsafe {
                    int SndBuf = iRowGl, RcvBuf = -1231;
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&SndBuf), (IntPtr)(&RcvBuf), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.INT, MPI.Wrappers.csMPI.Raw._OP.MAX, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    iRowGl = RcvBuf;
                }

                // clear row
                // ---------
                if (onthisProc) {
                    // ref. cell is on local MPI process
                    int jCell = (int)GlobalIndex - GridDat.CellPartitioning.i0;

                    ReducedRegionCode rrc;
                    int NoOfSpc = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);

                    // set matrix row to identity
                    int[] ColIdx = null;
                    int L = Mtx.GetOccupiedColumnIndices(iRowGl, ref ColIdx);
                    Mtx.SetValues(iRowGl, ColIdx, new double[L]);
                    Mtx.SetDiagonalElement(iRowGl, 1.0);


                    // clear RHS
                    int iRow = iRowGl - Mtx.RowPartitioning.i0;
                    rhs[iRow] = 0;
                }

                // clear column
                // ------------
                {
                    for (int i = Mtx.RowPartitioning.i0; i < Mtx.RowPartitioning.iE; i++) {
                        if (i != iRowGl) {
                            Mtx[i, iRowGl] = 0;
                        }
                    }
                }
            }
        }



        /// <summary>
        /// modifies a matrix <paramref name="Mtx"/> and a right-hand-side <paramref name="rhs"/>
        /// in order to fix the pressure at some reference point
        /// </summary>
        /// <param name="map">row mapping for <paramref name="Mtx"/> as well as <paramref name="rhs"/></param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="LsTrk"></param>
        /// <param name="Mtx"></param>
        /// <param name="rhs"></param>
        static public void SetPressureReferencePointResidual<T>(CoordinateVector currentState, int iVar, LevelSetTracker LsTrk, T Residual)
            where T : IList<double> //
        {
            using (new FuncTrace()) {
                var map = currentState.Mapping;
                var GridDat = map.GridDat;


                if (Residual.Count != map.LocalLength)
                    throw new ArgumentException();


                XDGBasis PressureBasis = (XDGBasis)map.BasisS[iVar];
                var grd = GridDat;
                int D = GridDat.SpatialDimension;

                long GlobalID, GlobalIndex;
                bool IsInside, onthisProc;
                grd.LocatePoint(new double[] { 5, 0 }, out GlobalID, out GlobalIndex, out IsInside, out onthisProc, LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement());


                int iRowGl = -111;

                if (onthisProc) {
                    int jCell = (int)GlobalIndex - GridDat.CellPartitioning.i0;
                    NodeSet CenterNode = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell), new double[D]);

                    MultidimensionalArray LevSetValues = LsTrk.DataHistories[0].Current.GetLevSetValues(CenterNode, jCell, 1); ;

                    MultidimensionalArray CenterNodeGlobal = MultidimensionalArray.Create(1, D);
                    GridDat.TransformLocal2Global(CenterNode, CenterNodeGlobal, jCell);
                    //Console.WriteLine("Pressure Ref Point @( {0:0.###E-00} | {1:0.###E-00} )", CenterNodeGlobal[0,0], CenterNodeGlobal[0,1]);

                    LevelSetSignCode scode = LevelSetSignCode.ComputeLevelSetBytecode(LevSetValues[0, 0]);
                    ReducedRegionCode rrc;
                    int No = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);
                    int iSpc = LsTrk.GetSpeciesIndex(rrc, scode);

                    iRowGl = (int)map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalIndex, PressureBasis.DOFperSpeciesPerCell * iSpc);
                }

                unsafe {
                    int SndBuf = iRowGl, RcvBuf = -1231;
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&SndBuf), (IntPtr)(&RcvBuf), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.INT, MPI.Wrappers.csMPI.Raw._OP.MAX, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    iRowGl = RcvBuf;
                }

                // clear row
                // ---------
                if (onthisProc) {
                    // set entry in residual vector equal to corresponding value in domain vector
                    // (as if the corresponding matrix would have a 1 in the diagonal element and 0 everywhere else)


                    int iRow = iRowGl - map.i0;
                    Residual[iRow] = currentState[iRow];
                }


            }
        }



        #region velocity jump (mass balance at interface)

        static ScalarFunctionEx GetVelocityJumpErrFunc(VectorField<XDGField> U, bool OnlyNormalComponent, bool squared) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;
            var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UB = U.Select(u => u.GetSpeciesShadowField("B")).ToArray();

            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster
                MultidimensionalArray uARes = MultidimensionalArray.Create(D, Len, K);
                MultidimensionalArray uBRes = MultidimensionalArray.Create(D, Len, K);

                for (int d = 0; d < D; d++) {
                    UA[d].Evaluate(j0, Len, NS, uARes.ExtractSubArrayShallow(d, -1, -1));
                    UB[d].Evaluate(j0, Len, NS, uBRes.ExtractSubArrayShallow(d, -1, -1));
                }

                if (OnlyNormalComponent) {
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {

                            double uJ = 0;
                            for (int d = 0; d < _D; d++) {
                                double nx = Normals[j, k, d];

                                uJ += (uBRes[d, j, k] - uARes[d, j, k]) * nx;
                            }
                            result[j, k] = squared ? uJ.Pow2() : uJ;
                        }
                    }
                } else {
                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {

                            double uJ = 0;
                            for (int d = 0; d < _D; d++) {
                                double q = (uBRes[d, j, k] - uARes[d, j, k]);
                                uJ += squared ? q.Pow2() : q;
                            }
                            result[j, k] = uJ;
                        }
                    }
                }
            };

        }

        /// <summary>
        /// Projects a surface property \f$ f\f$  (see below) of some vector field <paramref name="U"/>,
        /// onto a single-phase field.
        /// </summary>
        /// <param name="OnlyNormalComponent">
        /// if true, the norm of \f$ \llbracket \vec{u} \cdot \vec{n} \rrbracket =: f\f$  <br/>
        /// if false, the norm of \f$ \llbracket \vec{u} \rrbracket =: f\f$  <br/>
        /// </param>
        /// <param name="U">
        /// some XDG vector field
        /// </param>
        /// <param name="err">
        /// output accumulator; 
        /// </param>
        /// <param name="alpha">
        /// the usual scaling factor for the accumulation.
        /// </param>
        /// <param name="quadScheme">
        /// optional specification of the quadrature scheme.
        /// </param>
        public static void ProjectVelocityJumpNorm(this SinglePhaseField err, double alpha, VectorField<XDGField> U, bool OnlyNormalComponent, CellQuadratureScheme quadScheme = null) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetVelocityJumpErrFunc(U, OnlyNormalComponent, false);

            //double Jump_NORM = 0.0;
            int order = (U[0].Basis.Degree + err.Basis.Degree + 2);
            if (quadScheme == null)
                quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

            err.ProjectField(alpha, ErrFunc,
                quadScheme);
        }


        public static void ProjectMeanVelocityJumpNorm(this SinglePhaseField err, double alpha, VectorField<XDGField> U, bool OnlyNormalComponent, int momentFittingOrder) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetVelocityJumpErrFunc(U, OnlyNormalComponent, false);


            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            double mini = double.MaxValue, maxi = double.MinValue;

            CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        int jCell = i0 + i;
                        double val = err.GetMeanValue(jCell);

                        double JumpErr = ResultsOfIntegration[i, 0];
                        double ArcLen = ResultsOfIntegration[i, 1];

                        mini = Math.Min(JumpErr / ArcLen, mini);
                        maxi = Math.Max(JumpErr / ArcLen, maxi);

                        val += alpha * (JumpErr / ArcLen);
                        err.SetMeanValue(jCell, JumpErr / ArcLen);
                    }
                }
            ).Execute();

            // Console.WriteLine("jmp err mini:{0}, maxi {1}", mini, maxi);

        }


        /// <summary>
        /// The surface norm, i.e. \f$ \oint_{\mathfrak{I}} f^2 \ \mathrm{dS}\f$  for the jump 
        /// \f$ f\f$  (see below) of some vector field <paramref name="U"/>.
        /// </summary>
        /// <param name="OnlyNormalComponent">
        /// if true, the norm of \f$ \llbracket \vec{u} \cdot \vec{n} \rrbracket =: f\f$  <br/>
        /// if false, the norm of \f$ \llbracket \vec{u} \rrbracket =: f\f$  <br/>
        /// </param>
        /// <param name="U">
        /// some XDG vector field
        /// </param>
        /// <param name="momentFittingOrder">
        /// optional specification of quadrature order; if not positive, some value will be picked.
        /// </param>
        /// <returns></returns>
        public static double VelocityJumpNorm(VectorField<XDGField> U, bool OnlyNormalComponent, int momentFittingOrder) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetVelocityJumpErrFunc(U, OnlyNormalComponent, true);

            double Jump_NORM = 0.0;

            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        Jump_NORM += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return Math.Sqrt(Jump_NORM);
        }


        public static void ProjectMassBalanceNorm(this SinglePhaseField err, double alpha, VectorField<XDGField> U, int momentFittingOrder) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetVelocityJumpErrFunc(U, true, true);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        err.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0].Sqrt());
                    }
                }
            ).Execute();

        }

        #endregion


        #region momentum balance at interface

        static public double[] MomentumJumpNorm(VectorField<XDGField> U, XDGField P, double muA, double muB, int momentFittingOrder) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;
            var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UB = U.Select(u => u.GetSpeciesShadowField("B")).ToArray();

            ConventionalDGField pA = null, pB = null;
            if (P != null) {
                pA = P.GetSpeciesShadowField("A");
                pB = P.GetSpeciesShadowField("B");
            }

            //if (this.physParams.Sigma != 0.0)
            //    throw new NotImplementedException();
            //var FIx = this.SurfaceForce[0];
            //var FIy = this.SurfaceForce[1];

            //double muA = this.physParams.mu_A;
            //double muB = this.physParams.mu_B;

            double[] Jump_NORM = new double[D];
            for (int d = 0; d < D; d++) {
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D);
                    MultidimensionalArray Grad_UBRes = MultidimensionalArray.Create(Len, K, D, D);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray pBRes = MultidimensionalArray.Create(Len, K);


                    for (int i = 0; i < D; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1));
                        UB[i].EvaluateGradient(j0, Len, Ns, Grad_UBRes.ExtractSubArrayShallow(-1, -1, i, -1));
                    }
                    bool UsePressure = P != null;
                    if (UsePressure) {
                        pA.Evaluate(j0, Len, Ns, pARes);
                        pB.Evaluate(j0, Len, Ns, pBRes);
                    } else {
                        pARes.Clear();
                        pBRes.Clear();
                    }
                    //FIx.Evaluate(j0, Len, NodeSetIndex, FIRes.ExtractSubArrayShallow(-1, -1, 0));
                    //FIy.Evaluate(j0, Len, NodeSetIndex, FIRes.ExtractSubArrayShallow(-1, -1, 1));

                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);

                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {


                            double acc = 0.0;


                            // druck
                            if (UsePressure)
                                acc += (pBRes[j, k] - pARes[j, k]) * Normals[j, k, d];

                            // Nabla U
                            for (int dd = 0; dd < D; dd++) {
                                acc -= (muB * Grad_UBRes[j, k, d, dd] - muA * Grad_UARes[j, k, d, dd]) * Normals[j, k, dd];
                                acc -= (muB * Grad_UBRes[j, k, dd, d] - muA * Grad_UARes[j, k, dd, d]) * Normals[j, k, dd];
                            }
                            /*
                             */
                            //acc -= FIRes[j, k, d]*Normals[j, k, d];


                            result[j, k] = acc.Pow2();
                        }
                    }
                };

                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            Jump_NORM[d] += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            return Jump_NORM.Select(x => x.Sqrt()).ToArray();
        }


        static public double[] SurfaceTensionForceNorm(LevelSetTracker LsTrk, SinglePhaseField curv, double sigma, int momentFittingOrder) {

            int D = LsTrk.GridDat.SpatialDimension;

            double[] STF_Norm = new double[D];

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(QR.Nodes, i0, Length);

                    for (int d = 0; d < D; d++) {
                        curv.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));

                        double acc = 0.0;

                        for (int j = 0; j < Length; j++) {
                            int K = EvalResult.GetLength(1);
                            for (int k = 0; k < K; k++) {

                                acc = sigma * EvalResult[j, k, d] * Normals[j, k, d];

                                EvalResult[j, k, d] = acc.Pow2();
                            }
                        }

                    }

                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int d = 0; d < D; d++) {
                        for (int i = 0; i < Length; i++)
                            STF_Norm[d] += ResultsOfIntegration[i, d];
                    }
                }
            ).Execute();

            return STF_Norm.Select(x => x.Sqrt()).ToArray();

        }


        static ScalarFunctionEx GetMomentumBalanceFunc(XDGField P, VectorField<XDGField> U, SinglePhaseField C,
            PhysicalParameters physParam, SurfaceSressTensor sst, int d, bool squared) {

            int D = P.Basis.GridDat.SpatialDimension;

            ConventionalDGField pA = P.GetSpeciesShadowField("A");
            ConventionalDGField pB = P.GetSpeciesShadowField("B");

            var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UB = U.Select(u => u.GetSpeciesShadowField("B")).ToArray();

            double muA = physParam.mu_A;
            double muB = physParam.mu_B;
            double sigma = physParam.Sigma;

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pBRes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray Grad_UBRes = MultidimensionalArray.Create(Len, K, D, D);
                //MultidimensionalArray U_Res = MultidimensionalArray.Create(Len, K, D);
                //MultidimensionalArray GradU_Res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray curvRes = MultidimensionalArray.Create(Len, K);

                pA.Evaluate(i0, Len, nds, pARes);
                pB.Evaluate(i0, Len, nds, pBRes);

                for(int _d = 0; _d < D; _d++) {
                    UA[_d].EvaluateGradient(i0, Len, nds, Grad_UARes.ExtractSubArrayShallow(-1, -1, _d, -1));
                    UB[_d].EvaluateGradient(i0, Len, nds, Grad_UBRes.ExtractSubArrayShallow(-1, -1, _d, -1));
                }

                //for (int i = 0; i < D; i++) {
                //    Umean[i].Evaluate(i0, Len, nds, U_Res.ExtractSubArrayShallow(-1, -1, i));
                //    Umean[i].EvaluateGradient(i0, Len, nds, GradU_Res.ExtractSubArrayShallow(-1, -1, i, -1));
                //}

                C.Evaluate(i0, Len, nds, curvRes);

                var Normals = P.Basis.Tracker.DataHistories[0].Current.GetLevelSetNormals(nds, i0, Len);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double acc = 0.0;

                        // druck
                        acc -= (pBRes[j, k] - pARes[j, k]) * Normals[j, k, d];

                        // Nabla U
                        for (int dd = 0; dd < D; dd++) {
                            acc += (muB * Grad_UBRes[j, k, d, dd] - muA * Grad_UARes[j, k, d, dd]) * Normals[j, k, dd];
                            acc += (muB * Grad_UBRes[j, k, dd, d] - muA * Grad_UARes[j, k, dd, d]) * Normals[j, k, dd];
                        }

                        // isotropic surface tension force
                        acc -= sigma * curvRes[j, k] * Normals[j, k, d];

                        #region dynamic parts
                        // dynamic part of surface tension tensor
                        //MultidimensionalArray Nsurf = Normals.ExtractSubArrayShallow(j, k, -1);
                        //double[,] Psurf = new double[D, D];

                        //if (sst != SurfaceSressTensor.Isotropic) {
                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            if (d2 == d1)
                        //                Psurf[d1, d2] = (1 - Nsurf[d1] * Nsurf[d2]);
                        //            else
                        //                Psurf[d1, d2] = (0 - Nsurf[d1] * Nsurf[d2]);
                        //        }
                        //    }
                        //}

                        // surface rate of deformation
                        //if (sst == SurfaceSressTensor.SurfaceRateOfDeformation || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                        //    double muI = physParam.mu_I;

                        //    // GradU
                        //    double[,] GradUsurf = new double[D, D];
                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            for (int dd = 0; dd < D; dd++) {
                        //                GradUsurf[d1, d2] += Psurf[d1, dd] * GradU_Res[j, k, dd, d2];
                        //            }
                        //        }
                        //    }

                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            acc -= muI * Psurf[d, d2] * GradUsurf[d2, d1];
                        //        }
                        //    }

                        //    // GradU transpose
                        //    double[,] Psurf2 = new double[D, D];
                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            if (d2 == d1)
                        //                Psurf2[d1, d2] = (1 - 2 * Nsurf[d1] * Nsurf[d2]);
                        //            else
                        //                Psurf2[d1, d2] = (0 - 2 * Nsurf[d1] * Nsurf[d2]);
                        //        }
                        //    }

                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            acc -= muI * GradU_Res[j, k, d2, d] * Psurf2[d2, d1];
                        //        }
                        //    }

                        //}

                        // surface velocity divergence
                        //if (sst == SurfaceSressTensor.SurfaceVelocityDivergence || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                        //    double muI = physParam.mu_I;
                        //    double lamI = physParam.lambda_I;

                        //    double divUsurf = 0.0;
                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        for (int d2 = 0; d2 < D; d2++) {
                        //            divUsurf += Psurf[d1, d2] * GradU_Res[j, k, d2, d1];
                        //        }
                        //    }

                        //    for (int d1 = 0; d1 < D; d1++) {
                        //        acc -= (lamI - muI) * divUsurf * Psurf[d, d1];
                        //    }
                        //}
                        #endregion

                        if (squared) {
                            result[j, k] = acc.Pow2();
                        } else {
                            result[j, k] = acc;
                        }
                    }
                }

            };

        }

        public static double[] MomentumBalanceNormAtInterface(XDGField P, VectorField<XDGField> U, SinglePhaseField C,
            PhysicalParameters physParam, SurfaceSressTensor sst, int momentFittingOrder) {

            LevelSetTracker LsTrk = P.Basis.Tracker;

            int D = LsTrk.GridDat.SpatialDimension;

            double[] momBal_Norm = new double[D];

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    for (int d = 0; d < D; d++) {
                        ScalarFunctionEx momBalFunc = GetMomentumBalanceFunc(P, U, C, physParam, sst, d, true);
                        momBalFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, d));
                    }

                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int d = 0; d < D; d++) {
                        for (int i = 0; i < Length; i++)
                            momBal_Norm[d] += ResultsOfIntegration[i, d];
                    }
                }
            ).Execute();

            return momBal_Norm.Select(x => x.Sqrt()).ToArray();

        }

        public static void ProjectMomentumBalanceNorm(this SinglePhaseField err, double alpha, XDGField P, VectorField<XDGField> U, SinglePhaseField C,
            PhysicalParameters physParam, SurfaceSressTensor sst, int d, int momentFittingOrder) {

            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetMomentumBalanceFunc(P, U, C, physParam, sst, d, true);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        err.SetMeanValue(i0 + i, ResultsOfIntegration[i, 0].Sqrt());
                    }
                }
            ).Execute();

        }

        #endregion


        #region interface related properties (area, length, material points, velocities)

        public static double GetSpeciesArea(LevelSetTracker LsTrk, SpeciesId spcId) {

            double spcArea = 0.0;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            //XQuadSchemeHelper SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.SpeciesIdS.ToArray());
            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            CellQuadratureScheme vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                vqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        spcArea += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return spcArea;
        }


        public static double GetInterfaceLength(LevelSetTracker LsTrk) {

            double interLength = 0.0;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        interLength += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return interLength;
        }


        public static MultidimensionalArray GetInterfacePoints(LevelSetTracker LsTrk, LevelSet LevSet, SubGrid sgrd = null) {

            int D = LsTrk.GridDat.SpatialDimension;
            int p = LevSet.Basis.Degree;
            if (sgrd == null)
                sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);

            NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(p * 2).Nodes).ToArray();
            int Jsub = sgrd.LocalNoOfCells;
            int K = Nodes.Max(nds => nds.NoOfNodes);
            int numP = Jsub * K;

            var cp = new BoSSS.Solution.LevelSetTools.ClosestPointFinder(LsTrk, 0, sgrd, Nodes);

            MultidimensionalArray ClosestPoints = cp.X0_global_Resorted;

            MultidimensionalArray interfaceP = new MultidimensionalArray(2);
            interfaceP.Allocate(numP, D);

            for (int d = 0; d < D; d++) {
                MultidimensionalArray cp_d = ClosestPoints.ExtractSubArrayShallow(-1, -1, d).CloneAs().ResizeShallow(numP);
                interfaceP.ExtractSubArrayShallow(-1, d).Acc(1.0, cp_d);
            }

            // sort the interface points (non-equidistant)
            int[] permutation = new int[numP];
            for (int ind = 0; ind < numP; ind++)
                permutation[ind] = ind;

            Array.Sort(permutation, (int i, int j) => Math.Sign(interfaceP[i, 0] - interfaceP[j, 0]));

            MultidimensionalArray interfaceP_temp = interfaceP.CloneAs();
            for (int ip = 0; ip < numP; ip++) {
                for (int d = 0; d < D; d++) {
                    interfaceP[ip, d] = interfaceP_temp[permutation[ip], d];
                }
            }

            return interfaceP;

        }


        public static IDictionary<string, double[]> GetVelocityAtLevSet(LevelSetTracker LsTrk, VectorField<XDGField> velocity) {

            int D = LsTrk.GridDat.SpatialDimension;
            int p = velocity[0].Basis.Degree;
            SubGrid sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);
            NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(p * 2).Nodes).ToArray();
            int Jsub = sgrd.LocalNoOfCells;
            int K = Nodes.Max(nds => nds.NoOfNodes);

            var cp = new BoSSS.Solution.LevelSetTools.ClosestPointFinder(LsTrk, 0, sgrd, Nodes);

            ConventionalDGField[] VelocityA = velocity.Select(xVel => xVel.GetSpeciesShadowField("A")).ToArray();
            ConventionalDGField[] VelocityB = velocity.Select(xVel => xVel.GetSpeciesShadowField("B")).ToArray();

            MultidimensionalArray[] VelocityAEval = VelocityA.Select(sf => cp.EvaluateAtCp(sf)).ToArray();
            MultidimensionalArray[] VelocityBEval = VelocityB.Select(sf => cp.EvaluateAtCp(sf)).ToArray();

            MultidimensionalArray ClosestPoints = cp.X0_global_Resorted;

            IDictionary<string, double[]> velAtInterface = new Dictionary<string, double[]>();

            for (int d = 0; d < D; d++) {
                double[] x_d = ClosestPoints.ExtractSubArrayShallow(-1, -1, d).CloneAs().ResizeShallow(Jsub * K).To1DArray();
                velAtInterface.Add("x" + d, x_d);
            }

            for (int iSpc = 0; iSpc < 2; iSpc++) {
                MultidimensionalArray[] VelocityEval = iSpc == 0 ? VelocityAEval : VelocityBEval;
                string Species = iSpc == 0 ? "A" : "B";
                for (int d = 0; d < D; d++) {
                    double[] Vel_Spc_d = VelocityEval[d].ResizeShallow(Jsub * K).To1DArray();
                    velAtInterface.Add("u" + d + "_" + Species, Vel_Spc_d);
                }
            }

            return velAtInterface;

        }


        public static ConventionalDGField[] GetMeanVelocity(IEnumerable<DGField> Velocity, LevelSetTracker LsTrk, double rho_A = 0.0, double rho_B = 0.0) {

            ConventionalDGField[] meanVelocity;

            int D = LsTrk.GridDat.SpatialDimension;

            if (Velocity.ElementAt(0) is XDGField) {

                // XDG velocity -- take density-weighted mean value in cut-cells

                meanVelocity = new ConventionalDGField[D];

                CellMask CC = LsTrk.Regions.GetCutCellMask4LevSet(0);
                CellMask Neg = LsTrk.Regions.GetLevelSetWing(0, -1).VolumeMask;
                CellMask Pos = LsTrk.Regions.GetLevelSetWing(0, +1).VolumeMask;
                CellMask posNear = LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
                CellMask negNear = LsTrk.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

                for (int d = 0; d < D; d++) {
                    Basis b = ((XDGField)Velocity.ElementAt(d)).Basis.NonX_Basis;
                    meanVelocity[d] = new SinglePhaseField(b);


                    foreach (string spc in LsTrk.SpeciesNames) {
                        double rhoSpc;
                        switch (spc) {
                            case "A": rhoSpc = rho_A; break;
                            case "B": rhoSpc = rho_B; break;
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }

                        double scale = rhoSpc / (rho_A + rho_B);

                        meanVelocity[d].Acc(scale, ((XDGField)Velocity.ElementAt(d)).GetSpeciesShadowField(spc), CC);
                        switch (spc) {
                            //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                            case "A": meanVelocity[d].Acc(1.0, ((XDGField)Velocity.ElementAt(d)).GetSpeciesShadowField(spc), negNear); break;
                            case "B": meanVelocity[d].Acc(1.0, ((XDGField)Velocity.ElementAt(d)).GetSpeciesShadowField(spc), posNear); break;
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }
                    }

                }
            } else if (Velocity.ElementAt(0) is ConventionalDGField) {

                // plain DG velocity 

                meanVelocity = Velocity.Select(v => ((ConventionalDGField)v)).ToArray();
            } else {
                throw new ApplicationException("Should not happen.");
            }

            return meanVelocity;

        }

        #endregion



        //public static void ProjectEnergyJumpAtInterface(double alpha, SinglePhaseField acc, LevelSetTracker LsTrk, IEnumerable<XDGField> Velocity, XDGField Pressure, double muA, double muB,
        //    bool mean, CellQuadratureScheme quadScheme = null) {

        //    ScalarFunctionEx EnergyJumpFunc = GetEnergyJumpFunc(LsTrk, Velocity, Pressure, muA, muB);

        //    if (!mean) {

        //        int order = (acc.Basis.Degree + 2);
        //        if (quadScheme == null) {
        //            quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

        //            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
        //            //quadScheme = SchemeHelper.GetLevelSetquadScheme(0, LsTrk._Regions.GetCutCellMask());
        //        }

        //        acc.ProjectField(alpha, EnergyJumpFunc,
        //            quadScheme);

        //    } else {

        //        int order = 0;
        //        if (LsTrk.GetCachedOrders().Count > 0) {
        //            order = LsTrk.GetCachedOrders().Max();
        //        } else {
        //            order = 1;
        //        }

        //        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
        //        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

        //        double mini = double.MaxValue, maxi = double.MinValue;

        //        CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
        //            cqs.Compile(LsTrk.GridDat, order),
        //            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
        //                EnergyJumpFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
        //                EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
        //            },
        //            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
        //                for (int i = 0; i < Length; i++) {
        //                    int jCell = i0 + i;
        //                    double val = acc.GetMeanValue(jCell);

        //                    double JumpErr = ResultsOfIntegration[i, 0];
        //                    double ArcLen = ResultsOfIntegration[i, 1];

        //                    mini = Math.Min(JumpErr / ArcLen, mini);
        //                    maxi = Math.Max(JumpErr / ArcLen, maxi);

        //                    val += alpha * (JumpErr / ArcLen);
        //                    acc.SetMeanValue(jCell, JumpErr / ArcLen);
        //                }
        //            }
        //        ).Execute();

        //    }
        //}


        //public static void ProjectSurfaceEnergyChangerate(double alpha, SinglePhaseField acc, LevelSetTracker LsTrk, IEnumerable<DGField> Velocity, double sigma,
        //    bool mean, CellQuadratureScheme quadScheme = null) {

        //    ScalarFunctionEx SurfaceChangerate = GetSurfaceChangerateFunc(LsTrk, Velocity);

        //    if (!mean) {

        //        int order = (acc.Basis.Degree + 2);
        //        if (quadScheme == null) {
        //            quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

        //            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
        //            //quadScheme = SchemeHelper.GetLevelSetquadScheme(0, LsTrk._Regions.GetCutCellMask());
        //        }

        //        acc.ProjectField(alpha, SurfaceChangerate,
        //            quadScheme);

        //        acc.Scale(sigma);

        //    } else {

        //        int order = 0;
        //        if (LsTrk.GetCachedOrders().Count > 0) {
        //            order = LsTrk.GetCachedOrders().Max();
        //        } else {
        //            order = 1;
        //        }

        //        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
        //        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

        //        double mini = double.MaxValue, maxi = double.MinValue;

        //        CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
        //            cqs.Compile(LsTrk.GridDat, order),
        //            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
        //                SurfaceChangerate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
        //                EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
        //            },
        //            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
        //                for (int i = 0; i < Length; i++) {
        //                    int jCell = i0 + i;
        //                    double val = acc.GetMeanValue(jCell);

        //                    double JumpErr = ResultsOfIntegration[i, 0];
        //                    double ArcLen = ResultsOfIntegration[i, 1];

        //                    mini = Math.Min(JumpErr / ArcLen, mini);
        //                    maxi = Math.Max(JumpErr / ArcLen, maxi);

        //                    val += alpha * (JumpErr / ArcLen);
        //                    acc.SetMeanValue(jCell, JumpErr / ArcLen);
        //                }
        //            }
        //        ).Execute();

        //        acc.Scale(sigma);
        //    }
        //}


        //public static void ProjectCurvatureEnergy(double alpha, SinglePhaseField acc, LevelSetTracker LsTrk, IEnumerable<DGField> Velocity, double sigma,
        //    DGField Curvature, bool ExtVel, bool mean, CellQuadratureScheme quadScheme = null) {

        //    ScalarFunctionEx CurvEnergyFunc = GetCurvatureEnergyFunc(LsTrk, Velocity, sigma, Curvature, ExtVel);

        //    if (!mean) {

        //        int order = (acc.Basis.Degree + 2);
        //        if (quadScheme == null) {
        //            quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

        //            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
        //            //quadScheme = SchemeHelper.GetLevelSetquadScheme(0, LsTrk._Regions.GetCutCellMask());
        //        }

        //        acc.ProjectField(alpha, CurvEnergyFunc,
        //            quadScheme);

        //    } else {

        //        int order = 0;
        //        if (LsTrk.GetCachedOrders().Count > 0) {
        //            order = LsTrk.GetCachedOrders().Max();
        //        } else {
        //            order = 1;
        //        }

        //        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
        //        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

        //        double mini = double.MaxValue, maxi = double.MinValue;

        //        CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
        //            cqs.Compile(LsTrk.GridDat, order),
        //            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
        //                CurvEnergyFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
        //                EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
        //            },
        //            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
        //                for (int i = 0; i < Length; i++) {
        //                    int jCell = i0 + i;
        //                    double val = acc.GetMeanValue(jCell);

        //                    double JumpErr = ResultsOfIntegration[i, 0];
        //                    double ArcLen = ResultsOfIntegration[i, 1];

        //                    mini = Math.Min(JumpErr / ArcLen, mini);
        //                    maxi = Math.Max(JumpErr / ArcLen, maxi);

        //                    val += alpha * (JumpErr / ArcLen);
        //                    acc.SetMeanValue(jCell, JumpErr / ArcLen);
        //                }
        //            }
        //        ).Execute();

        //    }
        //}


        //public static void ProjectInterfaceDivergence(double alpha, DGField acc, LevelSetTracker LsTrk, IEnumerable<DGField> Velocity,
        //    bool OnlyNormalComp, bool mean, CellQuadratureScheme quadScheme = null) {

        //    ScalarFunctionEx InterfaceDivergence = GetInterfaceDivergenceFunc(LsTrk, Velocity, OnlyNormalComp);

        //    if (!mean) {

        //        int order = (acc.Basis.Degree + 2);
        //        if (quadScheme == null) {
        //            quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

        //            //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
        //            //quadScheme = SchemeHelper.GetLevelSetquadScheme(0, LsTrk._Regions.GetCutCellMask());
        //        }

        //        acc.ProjectField(alpha, InterfaceDivergence,
        //            quadScheme);

        //    } else {

        //        int order = 0;
        //        if (LsTrk.GetCachedOrders().Count > 0) {
        //            order = LsTrk.GetCachedOrders().Max();
        //        } else {
        //            order = 1;
        //        }

        //        var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
        //        CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

        //        double mini = double.MaxValue, maxi = double.MinValue;

        //        CellQuadrature.GetQuadrature(new int[] { 2 }, LsTrk.GridDat,
        //            cqs.Compile(LsTrk.GridDat, order),
        //            delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
        //                InterfaceDivergence(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
        //                EvalResult.ExtractSubArrayShallow(-1, -1, 1).SetAll(1.0);
        //            },
        //            delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
        //                for (int i = 0; i < Length; i++) {
        //                    int jCell = i0 + i;
        //                    double val = acc.GetMeanValue(jCell);

        //                    double JumpErr = ResultsOfIntegration[i, 0];
        //                    double ArcLen = ResultsOfIntegration[i, 1];

        //                    mini = Math.Min(JumpErr / ArcLen, mini);
        //                    maxi = Math.Max(JumpErr / ArcLen, maxi);

        //                    val += alpha * (JumpErr / ArcLen);
        //                    acc.SetMeanValue(jCell, JumpErr / ArcLen);
        //                }
        //            }
        //        ).Execute();

        //    }

        //}




        static public void ComputeGradientForParam(DGField f, VectorField<DGField> fGrad, LevelSetTracker LsTrk) {
            using(FuncTrace ft = new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                for(int d = 0; d < D; d++) {

                    foreach(var Spc in LsTrk.SpeciesIdS) { // loop over species...
                        // shadow fields
                        DGField f_Spc = ((f as XDGField).GetSpeciesShadowField(Spc));

                        (fGrad[d] as XDGField).GetSpeciesShadowField(Spc).Derivative(1.0, f_Spc, d);
                    }
                }

                fGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));

            }
        }


    }
}
