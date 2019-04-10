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
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Foundation.Grid;
using MPI;

namespace BoSSS.Solution.NSECommon {
    public static class IBMSolverUtils {


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
            where T : IList<double> {
            using (new FuncTrace()) {
                var GridDat = map.GridDat;
                
                if (rhs.Count != map.LocalLength)
                    throw new ArgumentException();
                if (!Mtx.RowPartitioning.EqualsPartition(map) || !Mtx.ColPartition.EqualsPartition(map))
                    throw new ArgumentException();

                Basis PressureBasis = (Basis)map.BasisS[iVar];
                int D = GridDat.SpatialDimension;

                long GlobalID, GlobalCellIndex;
                bool IsInside, onthisProc;
                GridDat.LocatePoint(new double[D], out GlobalID, out GlobalCellIndex, out IsInside, out onthisProc, 
                    LsTrk != null ? LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement() : null);
                
                int iRowGl = -111;
                if (onthisProc) {
                    //int jCell = (int)GlobalCellIndex - GridDat.CellPartitioning.i0;


                    //NodeSet CenterNode = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell), new double[D]);
                    //MultidimensionalArray LevSetValues = LsTrk.DataHistories[0].Current.GetLevSetValues(CenterNode, jCell, 1); ;


                    //MultidimensionalArray CenterNodeGlobal = MultidimensionalArray.Create(1, D);
                    //GridDat.TransformLocal2Global(CenterNode, CenterNodeGlobal, jCell);
                    //Console.WriteLine("Pressure Ref Point @( {0:0.###E-00} | {1:0.###E-00} )", CenterNodeGlobal[0,0], CenterNodeGlobal[0,1]);


                    //LevelSetSignCode scode = LevelSetSignCode.ComputeLevelSetBytecode(LevSetValues[0, 0]);
                    //ReducedRegionCode rrc;
                    //int No = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);
                    //int iSpc = LsTrk.GetSpeciesIndex(rrc, scode);

                    iRowGl = (int)map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalCellIndex, 0);

                }
                iRowGl = iRowGl.MPIMax();


                // clear row
                // ---------
                if (onthisProc) {
                    // ref. cell is on local MPI process
                    int jCell = (int)GlobalCellIndex - GridDat.CellPartitioning.i0;

                    //ReducedRegionCode rrc;
                    //int NoOfSpc = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);

                    // set matrix row to identity
                    Mtx.ClearRow(iRowGl);
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
        /// <param name="currentState">current state of velocity & pressure</param>
        /// <param name="iVar">the index of the pressure variable in the mapping <paramref name="map"/>.</param>
        /// <param name="LsTrk"></param>
        /// <param name="Residual"></param>
        static public void SetPressureReferencePointResidual<T>(CoordinateVector currentState, int iVar, LevelSetTracker LsTrk, T Residual)
            where T : IList<double> {
            using (new FuncTrace()) {
                var map = currentState.Mapping;
                var GridDat = map.GridDat;
                
                if (Residual.Count != map.LocalLength)
                    throw new ArgumentException();


                Basis PressureBasis = (Basis)map.BasisS[iVar];
                int D = GridDat.SpatialDimension;

                long GlobalID, GlobalCellIndex;
                bool IsInside, onthisProc;
                GridDat.LocatePoint(new double[D], out GlobalID, out GlobalCellIndex, out IsInside, out onthisProc, LsTrk.Regions.GetCutCellSubGrid().VolumeMask.Complement());
                
                int iRowGl = -111;
                if (onthisProc) {
                    int jCell = (int)GlobalCellIndex - GridDat.CellPartitioning.i0;


                    NodeSet CenterNode = new NodeSet(GridDat.iGeomCells.GetRefElement(jCell), new double[D]);
                    MultidimensionalArray LevSetValues = LsTrk.DataHistories[0].Current.GetLevSetValues(CenterNode, jCell, 1); ;


                    MultidimensionalArray CenterNodeGlobal = MultidimensionalArray.Create(1, D);
                    GridDat.TransformLocal2Global(CenterNode, CenterNodeGlobal, jCell);
                    //Console.WriteLine("Pressure Ref Point @( {0:0.###E-00} | {1:0.###E-00} )", CenterNodeGlobal[0,0], CenterNodeGlobal[0,1]);


                    LevelSetSignCode scode = LevelSetSignCode.ComputeLevelSetBytecode(LevSetValues[0, 0]);
                    ReducedRegionCode rrc;
                    int No = LsTrk.Regions.GetNoOfSpecies(jCell, out rrc);
                    int iSpc = LsTrk.GetSpeciesIndex(rrc, scode);

                    iRowGl = (int)map.GlobalUniqueCoordinateIndex_FromGlobal(iVar, GlobalCellIndex, 0);

                }
                iRowGl = iRowGl.MPIMax();


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




        /// <summary>
        /// Computes the energy stored in the fluid interface of a two-phase flow.
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <param name="sigma"></param>
        public static double GetSurfaceEnergy(LevelSetTracker LsTrk, double sigma) {
            using (new FuncTrace()) {
                if (LsTrk.LevelSets.Count != 1)
                    throw new NotImplementedException();

                double totSurface = 0;

                
                int order = 0;
                if (LsTrk.GetCachedOrders().Count > 0) {
                    order = LsTrk.GetCachedOrders().Max();
                } else {
                    order = 1;
                }

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, 1).XQuadSchemeHelper;
                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            totSurface += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

                return totSurface * sigma;
            }
        }

        /// <summary>
        /// Calculates the drag (x-component) and lift (y-component) forces acting on a IB contour
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="muA"></param>
        /// <returns></returns>
        static public double[] GetForces(VectorField<SinglePhaseField> U, SinglePhaseField P,
            LevelSetTracker LsTrk,
            double muA) {
            int D = LsTrk.GridDat.SpatialDimension;
           // var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UA = U.ToArray();

            int RequiredOrder = U[0].Basis.Degree * 3;
            //int RequiredOrder = LsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedSurfaceOrders(0).Max();
            //Console.WriteLine("Order reduction: {0} -> {1}", _RequiredOrder, RequiredOrder);

            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            ConventionalDGField pA = null;

            //pA = P.GetSpeciesShadowField("A");
            pA = P;

            double[] forces = new double[D];
            for (int d = 0; d < D; d++) {
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                    // Evaluate tangential velocity to level-set surface
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);


                    for (int i = 0; i < D; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }

                    pA.Evaluate(j0, Len, Ns, pARes);

                    if (LsTrk.GridDat.SpatialDimension == 2) {

                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;

                                // pressure
                                switch (d) {
                                    case 0:
                                    acc += pARes[j, k] * Normals[j, k, 0];
                                    acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                    break;
                                    case 1:
                                    acc += pARes[j, k] * Normals[j, k, 1];
                                    acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                    break;
                                    default:
                                    throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    } else {
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;

                                // pressure
                                switch (d) {
                                    case 0:
                                    acc += pARes[j, k] * Normals[j, k, 0];
                                    acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j,k,0,2] * Normals[j, k, 2];
                                    acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 2];
                                    break;
                                    case 1:
                                    acc += pARes[j, k] * Normals[j, k, 1];
                                    acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 2];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 2];
                                    break;
                                    case 2:
                                    acc += pARes[j, k] * Normals[j, k, 2];
                                    acc -= (2 * muA) * Grad_UARes[j, k, 2, 2] * Normals[j, k, 2];
                                    acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 1];
                                    acc -= (muA) * Grad_UARes[j, k, 0, 2] * Normals[j, k, 0];
                                    acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 1];
                                    break;
                                    default:
                                    throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    }
                
                };

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            forces[d] += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            for (int i = 0; i < D; i++)
                forces[i] = MPI.Wrappers.MPIExtensions.MPISum(forces[i]);

            return forces;
        }

        /// <summary>
        /// Calculates the Torque around the center of mass 
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="muA"></param>
        /// <param name="particleRadius"></param>
        /// <returns></returns>
        static public double GetTorque(VectorField<SinglePhaseField> U, SinglePhaseField P,
            LevelSetTracker LsTrk,
            double muA, double particleRadius) {
            var _LsTrk = LsTrk;
            int D = _LsTrk.GridDat.SpatialDimension;
            var UA = U.ToArray();

            //if (D > 2) throw new NotImplementedException("Currently only 2D cases supported");

            int RequiredOrder = U[0].Basis.Degree * 3;
            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine("Torque coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            ConventionalDGField pA = null;
            double force = new double();

            pA = P;

            ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                // Evaluate tangential velocity to level-set surface
                var Normals = _LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);

                for (int i = 0; i < D; i++) {
                    UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                }

                pA.Evaluate(j0, Len, Ns, pARes);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {


                        double acc = 0.0;
                        double acc2 = 0.0;


                        // Calculate the torque around a circular particle with a given radius (Paper Wan and Turek 2005)

                        acc += pARes[j, k] * Normals[j, k, 0];
                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                        acc *= -Normals[j, k, 1] * particleRadius;


                        acc2 += pARes[j, k] * Normals[j, k, 1];
                        acc2 -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                        acc2 -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                        acc2 -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                        acc2 *= Normals[j, k, 0] * particleRadius;

                        result[j, k] = acc + acc2;


                    }


                }

            };

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
            //var SchemeHelper = new XQuadSchemeHelper(_LsTrk, momentFittingVariant, _LsTrk.GetSpeciesId("A"));
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, _LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, _LsTrk.GridDat,
                cqs.Compile(_LsTrk.GridDat, RequiredOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        force += ResultsOfIntegration[i, 0];
                    }
                }

            ).Execute();


            return force;
        }

        /// <summary>
        /// Calculates the Torque around the center of mass 
        /// </summary>
        /// <param name="U"></param>
        /// <param name="P"></param>
        /// <param name="momentFittingVariant"></param>
        /// <param name="muA"></param>
        /// <param name="particleRadius"></param>
        /// <returns></returns>
        static public void GetCellValues(VectorField<XDGField> U, XDGField P,
            double muA, double particleRadius, SinglePhaseField P_atIB, SinglePhaseField gradU_atIB, SinglePhaseField gradUT_atIB) {
            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;
            var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();

            if (D > 2) throw new NotImplementedException("Currently only 2D cases supported");

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine("Cell values calculated by: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            ConventionalDGField pA = null;
            double circumference = new double();

            pA = P.GetSpeciesShadowField("A");

            for (int n = 0; n < 4; n++) {
                ScalarFunctionEx ErrFunc_CellVal = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                    // Evaluate tangential velocity to level-set surface
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);

                    for (int i = 0; i < D; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1));
                    }

                    pA.Evaluate(j0, Len, Ns, pARes);

                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {


                            double acc = 0.0;
                            double acc2 = 0.0;
                            switch (n) {

                                case 0: // Pressure part


                                acc += pARes[j, k] * Normals[j, k, 0];
                                acc *= -Normals[j, k, 1] * particleRadius;


                                acc2 += pARes[j, k] * Normals[j, k, 1];
                                acc2 *= Normals[j, k, 0] * particleRadius;

                                result[j, k] = acc + acc2;
                                break;
                                case 1: // GradU part

                                acc -= (1 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0]; // Attention was 2 times
                                acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                acc *= -Normals[j, k, 1] * particleRadius;

                                acc2 -= (1 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                acc2 -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                acc2 *= Normals[j, k, 0] * particleRadius;

                                result[j, k] = acc + acc2;
                                break;

                                case 2: // GradU_T part

                                acc -= (1 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0]; // Attention was 2 times
                                acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                acc *= -Normals[j, k, 1] * particleRadius;


                                acc2 -= (1 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1]; // Attention was 2 times
                                acc2 -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                acc2 *= Normals[j, k, 0] * particleRadius;

                                result[j, k] = acc + acc2;
                                break;

                                case 3: // Standardization with radians

                                result[j, k] = 1;
                                break;

                                default:
                                throw new NotImplementedException();

                            }

                        }


                    }

                };





                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper; //   new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc_CellVal(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++) {
                            switch (n) {
                                case 0:
                                P_atIB.SetMeanValue(i0, ResultsOfIntegration[i, 0]);
                                break;

                                case 1:
                                gradU_atIB.SetMeanValue(i0, ResultsOfIntegration[i, 0]);
                                break;

                                case 2:
                                gradUT_atIB.SetMeanValue(i0, ResultsOfIntegration[i, 0]);
                                break;

                                case 3:
                                circumference += ResultsOfIntegration[i, 0];
                                P_atIB.SetMeanValue(i0, P_atIB.GetMeanValue(i0) / ResultsOfIntegration[i, 0]);
                                gradU_atIB.SetMeanValue(i0, gradU_atIB.GetMeanValue(i0) / ResultsOfIntegration[i, 0]);
                                gradUT_atIB.SetMeanValue(i0, gradUT_atIB.GetMeanValue(i0) / ResultsOfIntegration[i, 0]);
                                break;

                                default:
                                throw new NotImplementedException();
                            }
                        }
                    }

                ).Execute();

            }

            Console.WriteLine("Circle circumference: " + circumference);

        }

        static public double[] GetParticleForces(VectorField<SinglePhaseField> U, SinglePhaseField P,
            LevelSetTracker LsTrk,
            double muA) {
            int D = LsTrk.GridDat.SpatialDimension;
            // var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UA = U.ToArray();

            int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            //int RequiredOrder = LsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedSurfaceOrders(0).Max();
            //Console.WriteLine("Order reduction: {0} -> {1}", _RequiredOrder, RequiredOrder);

            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);


            ConventionalDGField pA = null;

            //pA = P.GetSpeciesShadowField("A");
            pA = P;

            double[] forces = new double[D];
            for (int d = 0; d < D; d++) {
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_UARes = MultidimensionalArray.Create(Len, K, D, D); ;
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);

                    // Evaluate tangential velocity to level-set surface
                    var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);


                    for (int i = 0; i < D; i++) {
                        UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }

                    pA.Evaluate(j0, Len, Ns, pARes);

                    if (LsTrk.GridDat.SpatialDimension == 2) {

                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;

                                // pressure
                                switch (d) {
                                    case 0:
                                        acc += pARes[j, k] * Normals[j, k, 0];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        break;
                                    case 1:
                                        acc += pARes[j, k] * Normals[j, k, 1];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        break;
                                    default:
                                        throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    } else {
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {
                                double acc = 0.0;

                                // pressure
                                switch (d) {
                                    case 0:
                                        acc += pARes[j, k] * Normals[j, k, 0];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 0, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 2];
                                        break;
                                    case 1:
                                        acc += pARes[j, k] * Normals[j, k, 1];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 1, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 1] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 2];
                                        break;
                                    case 2:
                                        acc += pARes[j, k] * Normals[j, k, 2];
                                        acc -= (2 * muA) * Grad_UARes[j, k, 2, 2] * Normals[j, k, 2];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 0] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 2, 1] * Normals[j, k, 1];
                                        acc -= (muA) * Grad_UARes[j, k, 0, 2] * Normals[j, k, 0];
                                        acc -= (muA) * Grad_UARes[j, k, 1, 2] * Normals[j, k, 1];
                                        break;
                                    default:
                                        throw new NotImplementedException();
                                }

                                result[j, k] = acc;
                            }
                        }
                    }

                };

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;
                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, );
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            forces[d] += ResultsOfIntegration[i, 0];
                    }
                ).Execute();
            }

            for (int i = 0; i < D; i++)
                forces[i] = MPI.Wrappers.MPIExtensions.MPISum(forces[i]);

            return forces;
        }

        /// <summary>
        /// Calculates the drag (x-component) and lift (y-component) forces acting on a wall of a boundary fitted grid
        /// </summary>
        static public double[] GetForces_BoundaryFitted(VectorField<SinglePhaseField> GradU, VectorField<SinglePhaseField> GradV, SinglePhaseField StressXX, 
            SinglePhaseField StressXY, SinglePhaseField StressYY, SinglePhaseField P, LevelSetTracker LsTrk, double muA, double beta) {
            int D = LsTrk.GridDat.SpatialDimension;

            if (D > 2)
            {
                throw new ArgumentException("Method GetForces_BoundaryFitted only implemented for 2D (viscoelastic)!");
            }
            // var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            //var UA = U.ToArray();
            MultidimensionalArray Grad_U = new MultidimensionalArray(D);
            var _GradU = GradU.ToArray();
            var _GradV = GradV.ToArray();


            int RequiredOrder = _GradU[0].Basis.Degree * 3 + 2;
            //int RequiredOrder = U[0].Basis.Degree * 3 + 2;
            //int RequiredOrder = LsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedSurfaceOrders(0).Max();
            //Console.WriteLine("Order reduction: {0} -> {1}", _RequiredOrder, RequiredOrder);

            //if (RequiredOrder > agg.HMForder)
            //    throw new ArgumentException();

            Console.WriteLine();
            Console.WriteLine("Forces coeff: {0}, order = {1}", LsTrk.CutCellQuadratureType, RequiredOrder);

            SinglePhaseField _StressXX = StressXX;
            SinglePhaseField _StressXY = StressXY;
            SinglePhaseField _StressYY = StressYY;

            SinglePhaseField pA = null;

            //pA = P.GetSpeciesShadowField("A");
            pA = P;


            

            double[] forces = new double[D];
            for (int d = 0; d < D; d++)
            {
                ScalarFunctionEx ErrFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes
                    MultidimensionalArray Grad_URes = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray Grad_VRes = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray StressXXRes = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray StressXYRes = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray StressYYRes = MultidimensionalArray.Create(Len, K);

                    var Normals = LsTrk.GridDat.Edges.NormalsCache.GetNormals_Edge(Ns, j0, Len);
                    //var Normals = MultidimensionalArray.Create(1, Ns.Length, 1);
                    //var Normals = LsTrk.GridDat.Edges.NormalsForAffine;


                    for (int i = 0; i < D; i++)
                    {
                            _GradU[i].EvaluateEdge(j0, Len, Ns, Grad_URes.ExtractSubArrayShallow(-1, -1, i),
                                Grad_URes.ExtractSubArrayShallow(-1, -1, i), ResultIndexOffset: 0, ResultPreScale: 1);

                            _GradV[i].EvaluateEdge(j0, Len, Ns, Grad_VRes.ExtractSubArrayShallow(-1, -1, i),
                                Grad_VRes.ExtractSubArrayShallow(-1, -1, i), ResultIndexOffset: 0, ResultPreScale: 1);

                        //UA[i].EvaluateGradient(j0, Len, Ns, Grad_UARes.ExtractSubArrayShallow(-1, -1, i, -1), 0, 1);
                    }

                    //pA.Evaluate(j0, Len, Ns, pARes);
                    pA.EvaluateEdge(j0, Len, Ns, pARes, pARes, ResultIndexOffset: 0, ResultPreScale: 1);
                    _StressXX.EvaluateEdge(j0, Len, Ns, StressXXRes, StressXXRes, ResultIndexOffset: 0, ResultPreScale: 1);
                    _StressXY.EvaluateEdge(j0, Len, Ns, StressXYRes, StressXYRes, ResultIndexOffset: 0, ResultPreScale: 1);
                    _StressYY.EvaluateEdge(j0, Len, Ns, StressYYRes, StressYYRes, ResultIndexOffset: 0, ResultPreScale: 1);


                    //if (LsTrk.GridDat.SpatialDimension == 2)
                    //{


                    for (int j = 0; j < Len; j++) {
                        for (int k = 0; k < K; k++) {
                            double acc = 0.0;

                            // pressure
                            switch (d) {
                                case 0:
                                    acc += pARes[j, k] * Normals[j, k, 0];
                                    acc -= (2 * muA * beta) * Grad_URes[j, k, 0] * Normals[j, k, 0];
                                    acc -= (muA * beta) * Grad_URes[j, k, 1] * Normals[j, k, 1];
                                    acc -= (muA * beta) * Grad_VRes[j, k, 0] * Normals[j, k, 1];
                                    acc -= (muA) * StressXXRes[j, k] * Normals[j, k, 0];
                                    acc -= (muA) * StressXYRes[j, k] * Normals[j, k, 1];

                                    break;

                                case 1:
                                    acc += pARes[j, k] * Normals[j, k, 1];
                                    acc -= (2 * muA * beta) * Grad_VRes[j, k, 1] * Normals[j, k, 1];
                                    acc -= (muA * beta) * Grad_VRes[j, k, 0] * Normals[j, k, 0];
                                    acc -= (muA * beta) * Grad_URes[j, k, 1] * Normals[j, k, 0];
                                    acc -= (muA) * StressXYRes[j, k] * Normals[j, k, 0];
                                    acc -= (muA) * StressYYRes[j, k] * Normals[j, k, 1];
                                    break;
                                default:
                                    throw new NotImplementedException();
                            }

                            result[j, k] = acc;
                        }
                    }

                };


                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, RequiredOrder, 1).XQuadSchemeHelper;

                //EdgeMask Mask = new EdgeMask(LsTrk.GridDat, "Wall_bottom");
                EdgeMask Mask = new EdgeMask(LsTrk.GridDat, "Wall_cylinder");

                EdgeQuadratureScheme eqs = SchemeHelper.GetEdgeQuadScheme(LsTrk.GetSpeciesId("A"), IntegrationDomain: Mask);

                EdgeQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    eqs.Compile(LsTrk.GridDat, RequiredOrder), //  agg.HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        ErrFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            forces[d] += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

            }

            for (int i = 0; i < D; i++)
                forces[i] = MPI.Wrappers.MPIExtensions.MPISum(forces[i]);

            return forces;
        }
    }
}
