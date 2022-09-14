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
using BoSSS.Solution.LevelSetTools;
using System.Collections;

namespace BoSSS.Solution.XNSECommon {

    public static class XNSEUtils {

        /*
         * pressure ref pt is now controlled directly in the Multigrid Operator; FK, 11sep20

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
                double[] pt = (D == 2) ? new double[] { -100, -100 } : new double[] { -5.0, -5.0, -5.0 };
                grd.LocatePoint(pt, out GlobalID, out GlobalIndex, out IsInside, out onthisProc, LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement());


                int iRowGl = -111;

                if (onthisProc) {
                    int jCell = (int)GlobalIndex - GridDat.CellPartitioning.i0;
                    NodeSet CenterNode = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell), new double[D]);

                    MultidimensionalArray LevSetValues = LsTrk.DataHistories[0].Current.GetLevSetValues(CenterNode, jCell, 1);

                    MultidimensionalArray CenterNodeGlobal = MultidimensionalArray.Create(1, D);
                    GridDat.TransformLocal2Global(CenterNode, CenterNodeGlobal, jCell);
                    //Console.WriteLine("Pressure Ref Point @( {0:0.###E-00} | {1:0.###E-00} )", CenterNodeGlobal[0, 0], CenterNodeGlobal[0, 1]);

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
        /// modifies a residual (i.e. an operator evaluation)
        /// in order to fix the pressure at some reference point
        /// </summary>
        /// <param name="currentState">current state of velocity and pressure</param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="LsTrk"></param>
        /// <param name="Residual"></param>
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

                // somehow does not work for LDG...
                //grd.LocatePoint(new double[] { 5, 0 }, out GlobalID, out GlobalIndex, out IsInside, out onthisProc, LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement());

                //trying this:
                double[] pt = (D == 2) ? new double[] { -100, -100 } : new double[] { -5.0, -5.0, -5.0 };
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
                    // set entry in residual vector equal to corresponding value in domain vector
                    // (as if the corresponding matrix would have a 1 in the diagonal element and 0 everywhere else)


                    int iRow = iRowGl - map.i0;
                    Residual[iRow] = currentState[iRow];
                }


            }
        }

        */

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



        static ScalarFunctionEx GetFunc_MassBalanceAtInterface(VectorField<XDGField> u) {
            var LsTrk = u[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            var uA = u.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var uB = u.Select(v => v.GetSpeciesShadowField("B")).ToArray();

            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster

                MultidimensionalArray uARes = MultidimensionalArray.Create(D, Len, K);
                MultidimensionalArray uBRes = MultidimensionalArray.Create(D, Len, K);


                int JE = LsTrk.GridDat.Cells.Count;
                BitArray sbArray = new BitArray(JE);
                for (int j = j0; j < j0 + Len; j++) {
                    sbArray[j] = true;
                }
                CellMask sbmask = new CellMask(LsTrk.GridDat, sbArray);
                SubGrid sbgrd = new SubGrid(sbmask);

                ClosestPointFinder cp = new ClosestPointFinder(LsTrk, 0, sbgrd, NS.ToEnumerable());

                for (int d = 0; d < D; d++) {
                    uARes.ExtractSubArrayShallow(d, -1, -1).Set(cp.EvaluateAtCp(uA[d]));
                    uBRes.ExtractSubArrayShallow(d, -1, -1).Set(cp.EvaluateAtCp(uB[d]));
                }

                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double uJ = 0;
                        for (int d = 0; d < _D; d++) {
                            uJ += (uBRes[d, j, k] - uARes[d, j, k]) * Normals[j, k, d];
                        }
                        result[j, k] = uJ;
                    }
                }

            };

        }


        public static void ProjectMassBalanceAtInterface(this SinglePhaseField err, double alpha, VectorField<XDGField> u, CellQuadratureScheme quadScheme = null) {

            ScalarFunctionEx ErrFunc = GetFunc_MassBalanceAtInterface(u);

            var LsTrk = u[0].Basis.Tracker;
            int order = (u[0].Basis.Degree + err.Basis.Degree + 2);
            if (quadScheme == null)
                quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

            err.ProjectField(alpha, ErrFunc,
                quadScheme);

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
            PhysicalParameters physParam, bool Icoord, int dir, bool squared) {

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

                for (int _d = 0; _d < D; _d++) {
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

                        if (!Icoord) {

                            // druck
                            acc -= (pBRes[j, k] - pARes[j, k]) * Normals[j, k, dir];

                            // Nabla U
                            for (int dd = 0; dd < D; dd++) {
                                acc += (muB * Grad_UBRes[j, k, dir, dd] - muA * Grad_UARes[j, k, dir, dd]) * Normals[j, k, dd];
                                acc += (muB * Grad_UBRes[j, k, dd, dir] - muA * Grad_UARes[j, k, dd, dir]) * Normals[j, k, dd];
                            }

                            // isotropic surface tension force
                            acc -= sigma * curvRes[j, k] * Normals[j, k, dir];

                        } else {

                            double[] Vdir = new double[2];   // normal(0) or tangential(1)
                            if (dir == 0) {
                                Vdir = Normals.ExtractSubArrayShallow(j, k, -1).To1DArray();
                            } else if (dir == 1) {
                                Vdir[0] = -Normals[j, k, 1];
                                Vdir[1] = Normals[j, k, 0];
                            } else {
                                throw new NotImplementedException("ToDo 3D");
                            }

                            for (int d = 0; d < D; d++) {

                                // pressure
                                acc -= (pBRes[j, k] - pARes[j, k]) * Normals[j, k, d] * Vdir[d];

                                // velocity gradients
                                for (int dd = 0; dd < D; dd++) {
                                    acc += (muB * Grad_UBRes[j, k, d, dd] - muA * Grad_UARes[j, k, d, dd]) * Normals[j, k, dd] * Vdir[d];
                                    // transposed
                                    acc += (muB * Grad_UBRes[j, k, dd, d] - muA * Grad_UARes[j, k, dd, d]) * Normals[j, k, dd] * Vdir[d];
                                }

                                // surface tension
                                acc -= sigma * curvRes[j, k] * Normals[j, k, d] * Vdir[d];

                            }

                        }

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
            PhysicalParameters physParam, bool Icoord, int momentFittingOrder) {

            LevelSetTracker LsTrk = P.Basis.Tracker;

            int D = LsTrk.GridDat.SpatialDimension;

            double[] momBal_Norm = new double[D];

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    for (int d = 0; d < D; d++) {
                        ScalarFunctionEx momBalFunc = GetMomentumBalanceFunc(P, U, C, physParam, Icoord, d, true);
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

        public static void ProjectMomentumBalanceNorm(this VectorField<SinglePhaseField> err, XDGField P, VectorField<XDGField> U, SinglePhaseField C,
            PhysicalParameters physParam, bool Icoord) {

            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            for (int d = 0; d < D; d++) {

                ScalarFunctionEx ErrFunc = GetMomentumBalanceFunc(P, U, C, physParam, Icoord, d, true);

                int order = ((U[d].Basis.Degree - 1) + (P.Basis.Degree) + err[d].Basis.Degree + 2);
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++) {
                            err[d].SetMeanValue(i0 + i, ResultsOfIntegration[i, 0].Sqrt());
                        }
                    }
                ).Execute();

            }

        }



        static ScalarFunctionEx GetFunc_MomentumBalanceAtInterface(int dir, XDGField p, VectorField<XDGField> gradUx, VectorField<XDGField> gradUy, SinglePhaseField curv, 
            bool Icoord, PhysicalParameters physParam) {

            var LsTrk = gradUx[dir].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            var pA = p.GetSpeciesShadowField("A");
            var pB = p.GetSpeciesShadowField("B");

            var gradUxA = gradUx.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var gradUxB = gradUx.Select(v => v.GetSpeciesShadowField("B")).ToArray();

            var gradUyA = gradUy.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var gradUyB = gradUy.Select(v => v.GetSpeciesShadowField("B")).ToArray();

            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster

                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pBRes = MultidimensionalArray.Create(Len, K);

                MultidimensionalArray gradUxARes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUxBRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray gradUyARes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUyBRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray curvRes = MultidimensionalArray.Create(Len, K);

                int JE = LsTrk.GridDat.Cells.Count;
                BitArray sbArray = new BitArray(JE);
                for (int j = j0; j < j0 + Len; j++) {
                    sbArray[j] = true;
                }
                CellMask sbmask = new CellMask(LsTrk.GridDat, sbArray);
                SubGrid sbgrd = new SubGrid(sbmask);

                ClosestPointFinder cp = new ClosestPointFinder(LsTrk, 0, sbgrd, NS.ToEnumerable());

                for (int dd = 0; dd < D; dd++) {
                    pARes.Set(cp.EvaluateAtCp(pA));
                    pBRes.Set(cp.EvaluateAtCp(pB));
                    gradUxARes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUxA[dd]));
                    gradUxBRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUxB[dd]));
                    gradUyARes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUyA[dd]));
                    gradUyBRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUyB[dd]));
                    curvRes.Set(cp.EvaluateAtCp(curv));
                }

                MultidimensionalArray gradUARes = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray gradUBRes = MultidimensionalArray.Create(Len, K, D, D);

                gradUARes.ExtractSubArrayShallow(-1, -1, 0, 0).Set(gradUxARes.ExtractSubArrayShallow(-1, -1, 0));
                gradUARes.ExtractSubArrayShallow(-1, -1, 0, 1).Set(gradUxARes.ExtractSubArrayShallow(-1, -1, 1));
                gradUARes.ExtractSubArrayShallow(-1, -1, 1, 0).Set(gradUyARes.ExtractSubArrayShallow(-1, -1, 0));
                gradUARes.ExtractSubArrayShallow(-1, -1, 1, 1).Set(gradUyARes.ExtractSubArrayShallow(-1, -1, 1));

                gradUBRes.ExtractSubArrayShallow(-1, -1, 0, 0).Set(gradUxBRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUBRes.ExtractSubArrayShallow(-1, -1, 0, 1).Set(gradUxBRes.ExtractSubArrayShallow(-1, -1, 1));
                gradUBRes.ExtractSubArrayShallow(-1, -1, 1, 0).Set(gradUyBRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUBRes.ExtractSubArrayShallow(-1, -1, 1, 1).Set(gradUyBRes.ExtractSubArrayShallow(-1, -1, 1));


                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                double muA = physParam.mu_A;
                double muB = physParam.mu_B;
                double sigma = physParam.Sigma;

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {
                        double momB = 0.0;

                        if (!Icoord) {

                            // pressure
                            momB -= (pBRes[j, k] - pARes[j, k]) * Normals[j, k, dir];

                            // velocity gradients
                            for (int dd = 0; dd < _D; dd++) {
                                momB += (muB * gradUBRes[j, k, dir, dd] - muA * gradUARes[j, k, dir, dd]) * Normals[j, k, dd];
                                // transposed
                                momB += (muB * gradUBRes[j, k, dd, dir] - muA * gradUARes[j, k, dd, dir]) * Normals[j, k, dd];
                            }

                            // surface tension
                            momB -= sigma * curvRes[j, k] * Normals[j, k, dir];

                        } else {

                            double[] Vdir = new double[2];   // normal(0) or tangential(1)
                            if (dir == 0) {
                                Vdir = Normals.ExtractSubArrayShallow(j, k, -1).To1DArray();
                            } else if (dir == 1) {
                                Vdir[0] = -Normals[j, k, 1];
                                Vdir[1] = Normals[j, k, 0];
                            } else {
                                throw new NotImplementedException("ToDo 3D");
                            }


                            for (int d = 0; d < D; d++) {

                                // pressure
                                momB -= (pBRes[j, k] - pARes[j, k]) * Normals[j, k, d] * Vdir[d];

                                // velocity gradients
                                for (int dd = 0; dd < _D; dd++) {
                                    momB += (muB * gradUBRes[j, k, d, dd] - muA * gradUARes[j, k, d, dd]) * Normals[j, k, dd] * Vdir[d];
                                    // transposed
                                    momB += (muB * gradUBRes[j, k, dd, d] - muA * gradUARes[j, k, dd, d]) * Normals[j, k, dd] * Vdir[d];
                                }

                                // surface tension
                                momB -= sigma * curvRes[j, k] * Normals[j, k, d] * Vdir[d];

                            }

                        }

                        result[j, k] = momB;
                    }
                }

            };

        }


        public static void ProjectMomentumBalanceAtInterface(this VectorField<SinglePhaseField> err, double alpha, XDGField p, VectorField<XDGField> u, SinglePhaseField curv, 
            PhysicalParameters physParam, bool Icoord = true, CellQuadratureScheme quadScheme = null) {

            var LsTrk = u[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            VectorField<XDGField> GuX = new VectorField<XDGField>(D, u[0].Basis, XDGField.Factory);
            VectorField<XDGField> GuY = new VectorField<XDGField>(D, u[1].Basis, XDGField.Factory);

            foreach (SpeciesId spcId in LsTrk.SpeciesIdS) {
                CellMask sf = LsTrk.Regions.GetSpeciesMask(spcId);
                DGField f0_Spc = u[0].GetSpeciesShadowField(spcId);
                DGField f1_Spc = u[1].GetSpeciesShadowField(spcId);
                for (int dd = 0; dd < D; dd++) {
                    GuX[dd].GetSpeciesShadowField(spcId).Derivative(1.0, f0_Spc, dd, sf);
                    GuY[dd].GetSpeciesShadowField(spcId).Derivative(1.0, f1_Spc, dd, sf);
                }
            }
            GuX.ForEach(F => F.CheckForNanOrInf(true, true, true));
            GuY.ForEach(F => F.CheckForNanOrInf(true, true, true));


            for (int d = 0; d < D; d++) {

                ScalarFunctionEx ErrFunc = GetFunc_MomentumBalanceAtInterface(d, p, GuX, GuY, curv, Icoord, physParam);

                int order = ((u[d].Basis.Degree - 1) + (p.Basis.Degree) + err[d].Basis.Degree + 2);
                if (quadScheme == null)
                    quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

                err[d].ProjectField(alpha, ErrFunc,
                    quadScheme);

            }

        }



        static ScalarFunctionEx GetFunc_SurfaceTensionForce(int dir, LevelSetTracker LsTrk, VectorField<SinglePhaseField> gradUIx, VectorField<SinglePhaseField> gradUIy, SinglePhaseField curv,
            SurfaceSressTensor sst, bool Icoord, PhysicalParameters physParam) {

            int D = LsTrk.GridDat.SpatialDimension;

            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster

                MultidimensionalArray gradUIxRes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUIyRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray curvRes = MultidimensionalArray.Create(Len, K);

                int JE = LsTrk.GridDat.Cells.Count;
                BitArray sbArray = new BitArray(JE);
                for (int j = j0; j < j0 + Len; j++) {
                    sbArray[j] = true;
                }
                CellMask sbmask = new CellMask(LsTrk.GridDat, sbArray);
                SubGrid sbgrd = new SubGrid(sbmask);

                ClosestPointFinder cp = new ClosestPointFinder(LsTrk, 0, sbgrd, NS.ToEnumerable());

                for (int dd = 0; dd < D; dd++) {
                    gradUIxRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUIx[dd]));
                    gradUIyRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUIy[dd]));
                    curvRes.Set(cp.EvaluateAtCp(curv));
                }

                MultidimensionalArray gradUIRes = MultidimensionalArray.Create(Len, K, D, D);
                gradUIRes.ExtractSubArrayShallow(-1, -1, 0, 0).Set(gradUIxRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 0, 1).Set(gradUIxRes.ExtractSubArrayShallow(-1, -1, 1));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 1, 0).Set(gradUIyRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 1, 1).Set(gradUIyRes.ExtractSubArrayShallow(-1, -1, 1));


                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                double muA = physParam.mu_A;
                double muB = physParam.mu_B;
                double sigma = physParam.Sigma;

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {
                        double surfF = 0.0;

                        if (!Icoord) {

                            // isotropic part
                            surfF += sigma * curvRes[j, k] * Normals[j, k, dir];

                            if (sst == SurfaceSressTensor.SurfaceDivergence || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                                double divIuI = 0.0;
                                for (int d1 = 0; d1 < D; d1++) {
                                    for (int d2 = 0; d2 < D; d2++) {
                                        if (d2 == d1) {
                                            divIuI += (1.0 - Normals[j, k, d1] * Normals[j, k, d2]) * gradUIRes[j, k, d1, d2];
                                        } else {
                                            divIuI += (0.0 - Normals[j, k, d1] * Normals[j, k, d2]) * gradUIRes[j, k, d1, d2]; ;
                                        }
                                    }
                                }

                                surfF += sigma * physParam.lambda_I * divIuI * curvRes[j, k] * Normals[j, k, dir];
                            }

                            if (sst == SurfaceSressTensor.SurfaceRateOfDeformation || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                                throw new NotImplementedException("ToDo");
                            }

                        } else {

                            double[] Vdir = new double[2];   // normal(0) or tangential(1)
                            if (dir == 0) {
                                Vdir = Normals.ExtractSubArrayShallow(j, k, -1).To1DArray();
                            } else if (dir == 1) {
                                Vdir[0] = -Normals[j, k, 1];
                                Vdir[1] = Normals[j, k, 0];
                            } else {
                                throw new NotImplementedException("ToDo 3D");
                            }


                            for (int d = 0; d < D; d++) {

                                // isotropic part
                                surfF += sigma * curvRes[j, k] * Normals[j, k, d] * Vdir[d];

                                if (sst == SurfaceSressTensor.SurfaceDivergence || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                                    double divIuI = 0.0;
                                    for (int d1 = 0; d1 < D; d1++) {
                                        for (int d2 = 0; d2 < D; d2++) {
                                            if (d2 == d1) {
                                                divIuI += (1.0 - Normals[j, k, d1] * Normals[j, k, d2]) * gradUIRes[j, k, d1, d2];
                                            } else {
                                                divIuI += (0.0 - Normals[j, k, d1] * Normals[j, k, d2]) * gradUIRes[j, k, d1, d2]; ;
                                            }
                                        }
                                    }

                                    surfF += physParam.lambda_I * divIuI * curvRes[j, k] * Normals[j, k, d] * Vdir[d];
                                }

                                if (sst == SurfaceSressTensor.SurfaceRateOfDeformation || sst == SurfaceSressTensor.FullBoussinesqScriven) {
                                    //throw new NotImplementedException("ToDo");
                                    Console.WriteLine("WARNING: projection of surface tension force without surface rate of deformation terms");
                                }
                            }

                        }

                        result[j, k] = surfF;
                    }
                }

            };

        }


        public static void ProjectSurfaceTensionForce(this VectorField<SinglePhaseField> err, double alpha, ConventionalDGField[] uI, SinglePhaseField curv, 
            LevelSetTracker LsTrk, PhysicalParameters physParam, SurfaceSressTensor sst, bool Icoord = true, CellQuadratureScheme quadScheme = null) {

 
            int D = LsTrk.GridDat.SpatialDimension;

            VectorField<SinglePhaseField> GuIx = new VectorField<SinglePhaseField>(D, uI[0].Basis, SinglePhaseField.Factory);
            VectorField<SinglePhaseField> GuIy = new VectorField<SinglePhaseField>(D, uI[1].Basis, SinglePhaseField.Factory);
            for (int dd = 0; dd < D; dd++) {
                GuIx[dd].DerivativeByFlux(1.0, uI[0], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
                GuIy[dd].DerivativeByFlux(1.0, uI[1], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
            }
            GuIx.ForEach(F => F.CheckForNanOrInf(true, true, true));
            GuIy.ForEach(F => F.CheckForNanOrInf(true, true, true));


            for (int d = 0; d < D; d++) {

                ScalarFunctionEx ErrFunc = GetFunc_SurfaceTensionForce(d, LsTrk, GuIx, GuIy, curv, sst, Icoord, physParam);

                int order = ((uI[d].Basis.Degree - 1) + (uI[d].Basis.Degree - 2) + err[d].Basis.Degree + 2);
                if (quadScheme == null)
                    quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

                err[d].ProjectField(alpha, ErrFunc,
                    quadScheme);

            }

        }



        #endregion


        #region interface related properties (area, length, material points, velocities)

        public static double GetSpeciesArea(LevelSetTracker LsTrk, SpeciesId spcId, int quadRuleOrder = -1) {

            double spcArea = 0.0;

            int order = (quadRuleOrder < 0) ? 1 : quadRuleOrder;
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


        public static double GetInterfaceLength(LevelSetTracker LsTrk, int quadRuleOrder = -1) {

            double interLength = 0.0;

            int order  = (quadRuleOrder < 0) ? 1 : quadRuleOrder;
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


        public static double[] GetSurfaceTensionNetForce(LevelSetTracker LsTrk, double surfTcoeff) {

            int D = LsTrk.GridDat.SpatialDimension;
            double[] SurfTnetF = new double[D];

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            EdgeQuadratureScheme eqs = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(LsTrk.GetSpeciesId("A"), 0); // GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            EdgeQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                eqs.Compile(LsTrk.GridDat, order),
                delegate (int e0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    var EdgeNormals = LsTrk.GridDat.Edges.NormalsForAffine;

                    for (int j = 0; j < Length; j++) {

                        var NormalsIN = LsTrk.DataHistories[0].Current.GetLevelSetNormals(
                            QR.Nodes.GetVolumeNodeSet(LsTrk.GridDat, LsTrk.GridDat.Edges.Edge2CellTrafoIndex[e0 + j, 0], false), LsTrk.GridDat.Edges.CellIndices[e0 + j, 0], 1);
                        var NormalsOT = LsTrk.DataHistories[0].Current.GetLevelSetNormals(
                            QR.Nodes.GetVolumeNodeSet(LsTrk.GridDat, LsTrk.GridDat.Edges.Edge2CellTrafoIndex[e0 + j, 1], false), LsTrk.GridDat.Edges.CellIndices[e0 + j, 1], 1);

                        for (int k = 0; k < QR.NoOfNodes; k++) {

                            double[] tauIN = new double[D];
                            double[] tauOT = new double[D];
                            for (int d1 = 0; d1 < D; d1++) {
                                for (int d2 = 0; d2 < D; d2++) {
                                    double nn = NormalsIN[0, k, d1] * NormalsIN[0, k, d2];
                                    if (d1 == d2) {
                                        tauIN[d1] += (1 - nn) * EdgeNormals[e0 + j, d2];
                                    } else {
                                        tauIN[d1] += -nn * EdgeNormals[e0 + j, d2];
                                    }
                                    nn = NormalsOT[0, k, d1] * NormalsOT[0, k, d2];
                                    if (d1 == d2) {
                                        tauOT[d1] += (1 - nn) * EdgeNormals[e0 + j, d2];
                                    } else {
                                        tauOT[d1] += -nn * EdgeNormals[e0 + j, d2];
                                    }
                                }
                            }
                            tauIN.Normalize();
                            tauOT.Normalize();

                            for (int d = 0; d < D; d++) {
                                EvalResult[j, k, d] = tauOT[d] - tauIN[d];
                            }
            
                        }
                    }

                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        for (int d = 0; d < D; d++) {
                            SurfTnetF[d] += surfTcoeff * ResultsOfIntegration[i, d];
                        }
                    }
                }
            ).Execute();

            return SurfTnetF;

        }


        public static Tuple<List<int>, List<double[]>, List<double[]>> GetSurfaceTensionForceAtEdges(LevelSetTracker LsTrk) {

            int D = LsTrk.GridDat.SpatialDimension;

            int order = 0;
            if (LsTrk.GetCachedOrders().Count > 0) {
                order = LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            List<int> edgeID = new List<int>();
            List<double[]> surfTvectrosIN = new List<double[]>();
            List<double[]> surfTvectrosOT = new List<double[]>();

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;
            EdgeQuadratureScheme eqs = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(LsTrk.GetSpeciesId("A"), 0); // GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            EdgeQuadrature.GetQuadrature(new int[] { D }, LsTrk.GridDat,
                eqs.Compile(LsTrk.GridDat, order),
                delegate (int e0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                    var EdgeNormals = LsTrk.GridDat.Edges.NormalsForAffine;

                    for (int j = 0; j < Length; j++) {

                        if (LsTrk.GridDat.BoundaryEdges.Contains(e0 + j))
                            continue;

                        NodeSet vnds_loc = QR.Nodes.GetVolumeNodeSet(LsTrk.GridDat, LsTrk.GridDat.Edges.Edge2CellTrafoIndex[e0 + j, 0], false);
                        NodeSet vnds_glob = vnds_loc.CloneAs();
                        LsTrk.GridDat.TransformLocal2Global(vnds_loc, vnds_glob, LsTrk.GridDat.Edges.CellIndices[e0 + j, 0]);
                        Console.WriteLine("volume node ({0}, {1})", vnds_glob[0, 0], vnds_glob[0, 1]);

                        var NormalsIN = LsTrk.DataHistories[0].Current.GetLevelSetNormals(
                            QR.Nodes.GetVolumeNodeSet(LsTrk.GridDat, LsTrk.GridDat.Edges.Edge2CellTrafoIndex[e0 + j, 0], false), LsTrk.GridDat.Edges.CellIndices[e0 + j, 0], 1);
                        var NormalsOT = LsTrk.DataHistories[0].Current.GetLevelSetNormals(
                            QR.Nodes.GetVolumeNodeSet(LsTrk.GridDat, LsTrk.GridDat.Edges.Edge2CellTrafoIndex[e0 + j, 1], false), LsTrk.GridDat.Edges.CellIndices[e0 + j, 1], 1);

                        for (int k = 0; k < QR.NoOfNodes; k++) {

                            double[] tauIN = new double[D];
                            double[] tauOT = new double[D];
                            for (int d1 = 0; d1 < D; d1++) {
                                for (int d2 = 0; d2 < D; d2++) {
                                    double nn = NormalsIN[0, k, d1] * NormalsIN[0, k, d2];
                                    if (d1 == d2) {
                                        tauIN[d1] += (1 - nn) * EdgeNormals[e0 + j, d2];
                                    } else {
                                        tauIN[d1] += -nn * EdgeNormals[e0 + j, d2];
                                    }
                                    nn = NormalsOT[0, k, d1] * NormalsOT[0, k, d2];
                                    if (d1 == d2) {
                                        tauOT[d1] += (1 - nn) * EdgeNormals[e0 + j, d2];
                                    } else {
                                        tauOT[d1] += -nn * EdgeNormals[e0 + j, d2];
                                    }
                                }
                            }
                            tauIN.Normalize();
                            tauOT.Normalize();

                            //tauIN.ScaleV(surfTcoeff);
                            //tauOT.ScaleV(surfTcoeff);

                            edgeID.Add(e0+j);
                            surfTvectrosIN.Add(tauIN);
                            surfTvectrosOT.Add(tauOT);

                        }
                    }

                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    // do nothing
                }
            ).Execute();

            return new Tuple<List<int>, List<double[]>, List<double[]>>(edgeID, surfTvectrosIN, surfTvectrosOT);
        }


        public static MultidimensionalArray GetInterfacePoints(LevelSetTracker LsTrk, SinglePhaseField LevSet, SubGrid sgrd = null, int quadRuleOrderForNodeSet = -1) {

            int D = LsTrk.GridDat.SpatialDimension;
            int p = LevSet.Basis.Degree;
            if (sgrd == null)
                sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);

            int quadRule = (quadRuleOrderForNodeSet < 0) ? p * 2 : quadRuleOrderForNodeSet;
            NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(quadRule).Nodes).ToArray();
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




        static public void ComputeGradientForParam(DGField f, DGField[] fGrad, LevelSetTracker LsTrk, Dictionary<string, double> a = null, SubGrid optionalSubGrid = null) {
            using(FuncTrace ft = new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                for(int d = 0; d < D; d++) {

                    foreach(var Spc in LsTrk.SpeciesNames) { // loop over species...
                        // shadow fields
                        DGField f_Spc = ((f as XDGField).GetSpeciesShadowField(Spc));

                        double aSpc = (a != null && a.ContainsKey(Spc)) ? a[Spc] : 1.0;
                        SubGrid optSgrd;
                        if (optionalSubGrid == null) {
                            optSgrd = LsTrk.Regions.GetSpeciesSubGrid(Spc);
                        } else {
                            optSgrd = optionalSubGrid;
                        }
                        (fGrad[d] as XDGField).GetSpeciesShadowField(Spc).DerivativeByFlux(aSpc, f_Spc, d, optSgrd);
                    }
                }

                fGrad.ForEach(F => F.CheckForNanOrInf(true, true, true));
            }
        }


    }
}
