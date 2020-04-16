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
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools;
using System.Collections;
using BoSSS.Platform;

namespace BoSSS.Solution.XNSECommon {

    public static class EnergyUtils {

        #region global energies (kinetic energy in bulk and surface energy in interface)

        /// <summary>
        /// Computes the kinetic energy stored in the bulk of a multi-phase flow field.
        /// </summary>
        /// <param name="Velocity"></param>
        /// <param name="rho">Density of the fluids, ordering in the array correlates with species 
        /// ordering in the level-set tracker, see <see cref="LevelSetTracker.SpeciesIdS"/>.
        /// </param>
        /// <param name="lsTrk"></param>
        /// <param name="momentFittingOrder"></param>
        public static double GetKineticEnergy<T>(LevelSetTracker lsTrk, IEnumerable<T> Velocity, double[] rho, int momentFittingOrder, int HistInd = 1)
            where T : DGField //
        {
            using(new FuncTrace()) {
                int D = lsTrk.GridDat.SpatialDimension;
                if(Velocity.Count() != D) {
                    throw new ArgumentException();
                }
                if(lsTrk.SpeciesIdS.Count != rho.Length)
                    throw new ArgumentException();

                //int order = 0;
                //if (lsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedVolumeOrders(0).Length > 0) {
                //    order = lsTrk.GetXQuadFactoryHelper(momentFittingVariant).GetCachedVolumeOrders(0).Max();
                //}
                //order = Math.Max(order, Velocity.ElementAt(0).Basis.Degree * 2);

                var SchemeHelper = lsTrk.GetXDGSpaceMetrics(lsTrk.SpeciesIdS.ToArray(), momentFittingOrder, HistInd).XQuadSchemeHelper;
                //var SchemeHelper = new XQuadSchemeHelper(lsTrk, momentFittingVariant, lsTrk.SpeciesIdS.ToArray());
                double kinE = 0.0;
                for(int iSpc = 0; iSpc < lsTrk.SpeciesIdS.Count; iSpc++) {
                    double _rho = rho[iSpc];
                    SpeciesId spId = lsTrk.SpeciesIdS[iSpc];

                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                    for(int d = 0; d < D; d++) {
                        DGField U = Velocity.ElementAt(d);
                        if(U is XDGField) {
                            if(!object.ReferenceEquals((U as XDGField).Basis.Tracker, lsTrk))
                                throw new ArgumentException();
                            U = (U as XDGField).GetSpeciesShadowField(spId);
                        }
                        kinE += U.L2Error(null, momentFittingOrder, scheme).Pow2() * _rho * 0.5;
                    }
                }

                return kinE;
            }
        }

        /// <summary>
        /// Computes the energy stored in the fluid interface of a two-phase flow.
        /// </summary>
        /// <param name="LsTrk"></param>
        /// <param name="sigma">Surface tension \f$ \sigma \f$.</param>
        /// <param name="momentFittingOrder"></param>
        public static double GetSurfaceEnergy(LevelSetTracker LsTrk, double sigma, int momentFittingOrder, int HistInd = 1) {
            using(new FuncTrace()) {
                if(LsTrk.LevelSets.Count != 1)
                    throw new NotImplementedException();

                double totSurface = 0;

                int order = 0;
                if(LsTrk.GetCachedOrders().Count > 0) {
                    order = LsTrk.GetCachedOrders().Max();
                } else {
                    order = 1;
                }

                //var SchemeHelper = new XQuadSchemeHelper(LsTrk, momentFittingVariant, LsTrk.GetSpeciesId("A"));
                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, order, HistInd).XQuadSchemeHelper;
                CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

                CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                    cqs.Compile(LsTrk.GridDat, order),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for(int i = 0; i < Length; i++)
                            totSurface += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

                return totSurface * sigma;
            }
        }


        public static double GetKineticDissipation(LevelSetTracker LsTrk, DGField[] Velocity, double[] mu, int momentFittingOrder, DGField[] extF = null, int HistInd = 1) {
            using(new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                if(Velocity.Count() != D) {
                    throw new ArgumentException();
                }
                if(LsTrk.SpeciesIdS.Count != mu.Length)
                    throw new ArgumentException();

                double dE = 0.0;

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), momentFittingOrder, HistInd).XQuadSchemeHelper;

                for(int iSpc = 0; iSpc < LsTrk.SpeciesIdS.Count; iSpc++) {
                    SpeciesId spcId = LsTrk.SpeciesIdS[iSpc];
                    double _mu = mu[iSpc];

                    var Uspc = Velocity.Select(u => (u as XDGField).GetSpeciesShadowField(spcId)).ToArray();
                    ScalarFunctionEx changerate_dEspc = GetSpeciesKineticDissipationFunc(Uspc, _mu, extF);

                    CellQuadratureScheme vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
                    CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                        vqs.Compile(LsTrk.GridDat, momentFittingOrder),
                        delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                            changerate_dEspc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                        },
                        delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                            for(int i = 0; i < Length; i++)
                                dE += ResultsOfIntegration[i, 0];
                        }
                    ).Execute();

                }

                return dE;
            }
        }


        public static double GetInterfaceDilatationalViscosityEnergyCR(LevelSetTracker LsTrk, ConventionalDGField[] uI, double lamI, int momentFittingOrder) {

            double dilViscEnergy = 0.0;

            ScalarFunctionEx dilViscEnergyFunc = GetInterfaceDivergenceFunc(LsTrk, uI, true);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    dilViscEnergyFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        dilViscEnergy += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return lamI * dilViscEnergy;

        }

        static ScalarFunctionEx GetInterfaceShearViscosityEnergyCRFunc(LevelSetTracker LsTrk, ConventionalDGField[] uI, bool squared) {

            int D = LsTrk.GridDat.SpatialDimension;

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray GradU_Res = MultidimensionalArray.Create(Len, K, D, D);

                for(int i = 0; i < D; i++) {
                    uI.ElementAt(i).EvaluateGradient(i0, Len, nds, GradU_Res.ExtractSubArrayShallow(-1, -1, i, -1));
                }

                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(nds, i0, Len);

                for(int j = 0; j < Len; j++) {
                    for(int k = 0; k < K; k++) {

                        MultidimensionalArray Nsurf = Normals.ExtractSubArrayShallow(j, k, -1);
                        double[,] Psurf = new double[D, D];
                        for(int d1 = 0; d1 < D; d1++) {
                            for(int d2 = 0; d2 < D; d2++) {
                                if(d2 == d1)
                                    Psurf[d1, d2] = (1 - Nsurf[d1] * Nsurf[d2]);
                                else
                                    Psurf[d1, d2] = (0 - Nsurf[d1] * Nsurf[d2]);
                            }
                        }

                        double[,] GradUsurf = new double[D, D];
                        for(int d1 = 0; d1 < D; d1++) {
                            for(int d2 = 0; d2 < D; d2++) {
                                for(int dd = 0; dd < D; dd++) {
                                    GradUsurf[d1, d2] += Psurf[d1, dd] * GradU_Res[j, k, dd, d2];
                                }
                            }
                        }

                        double acc = 0.0;

                        // GradU
                        for(int d1 = 0; d1 < D; d1++) {
                            for(int d2 = 0; d2 < D; d2++) {
                                for(int dd = 0; dd < D; dd++) {
                                    acc += Psurf[d1, dd] * GradUsurf[dd, d2] * GradUsurf[d1, d2];
                                }
                            }
                        }

                        // GradU transpose
                        double[,] Psurf2 = new double[D, D];
                        for(int d1 = 0; d1 < D; d1++) {
                            for(int d2 = 0; d2 < D; d2++) {
                                if(d2 == d1)
                                    Psurf2[d1, d2] = (1 - 2 * Nsurf[d1] * Nsurf[d2]);
                                else
                                    Psurf2[d1, d2] = (0 - 2 * Nsurf[d1] * Nsurf[d2]);
                            }
                        }

                        for(int d1 = 0; d1 < D; d1++) {
                            for(int d2 = 0; d2 < D; d2++) {
                                for(int dd = 0; dd < D; dd++) {
                                    acc -= GradU_Res[j, k, d1, dd] * Psurf2[dd, d2] * GradUsurf[d1, d2];
                                }
                            }
                        }

                        if(squared) {
                            result[j, k] = acc.Pow2();
                        } else {
                            result[j, k] = acc;
                        }
                    }
                }

            };

        }

        public static double GetInterfaceShearViscosityEnergyCR(LevelSetTracker LsTrk, ConventionalDGField[] uI, double muI, int momentFittingOrder) {

            double shearViscEnergy = 0.0;

            ScalarFunctionEx shearViscEnergyFunc = GetInterfaceShearViscosityEnergyCRFunc(LsTrk, uI, false);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    shearViscEnergyFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for(int i = 0; i < Length; i++)
                        shearViscEnergy += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return muI * shearViscEnergy;

        }

        #endregion


        #region local energy fields


        static ScalarFunctionEx GetSpeciesKineticDissipationFunc(DGField[] U, double mu, DGField[] F = null) {

            int D = U[0].Basis.GridDat.SpatialDimension;

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes

                MultidimensionalArray GradU_res = MultidimensionalArray.Create(Len, K, D, D);
                for (int i = 0; i < D; i++) {
                    U[i].EvaluateGradient(i0, Len, nds, GradU_res.ExtractSubArrayShallow(-1, -1, i, -1));
                }

                MultidimensionalArray U_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray F_res = MultidimensionalArray.Create(Len, K, D);
                if (F != null) {
                    for (int i = 0; i < D; i++) {
                        U.ElementAt(i).Evaluate(i0, Len, nds, U_res.ExtractSubArrayShallow(-1, -1, i));
                        F.ElementAt(i).Evaluate(i0, Len, nds, F_res.ExtractSubArrayShallow(-1, -1, i));
                    }
                }

                double acc;
                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        acc = 0.0;

                        for (int d = 0; d < D; d++) {
                            // viscous term
                            for (int dd = 0; dd < D; dd++) {
                                acc -= mu * (GradU_res[j, k, d, dd] * GradU_res[j, k, dd, d] + GradU_res[j, k, dd, d] * GradU_res[j, k, dd, d]);
                            }
                            // external forces (gravity)
                            if (F != null) {
                                acc += F_res[j, k, d] * U_res[j, k, d];
                            }
                        }
                        result[j, k] = acc;
                    }
                }

            };

        }


        public static void ProjectKineticDissipation(this XDGField proj, LevelSetTracker LsTrk, DGField[] Velocity, double[] mu, int momentFittingOrder, DGField[] extF = null, int HistInd = 1) {
            using (new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                if (Velocity.Count() != D) {
                    throw new ArgumentException();
                }
                if (LsTrk.SpeciesIdS.Count != mu.Length)
                    throw new ArgumentException();

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), momentFittingOrder, HistInd).XQuadSchemeHelper;

                for (int iSpc = 0; iSpc < LsTrk.SpeciesIdS.Count; iSpc++) {
                    SpeciesId spcId = LsTrk.SpeciesIdS[iSpc];
                    double muSpc = mu[iSpc];

                    var Uspc = Velocity.Select(u => (u as XDGField).GetSpeciesShadowField(spcId)).ToArray();
                    ScalarFunctionEx spcKinDissip = GetSpeciesKineticDissipationFunc(Uspc, muSpc, extF);

                    proj.GetSpeciesShadowField(spcId).ProjectField(spcKinDissip);

                }
            }
        }


        public static void ProjectKineticEnergy(this XDGField proj, LevelSetTracker LsTrk, XDGField[] Velocity, double[] rho, int momentFittingOrder, int HistInd = 1) {
            using (new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                if (Velocity.Count() != D) {
                    throw new ArgumentException();
                }
                if (LsTrk.SpeciesIdS.Count != rho.Length)
                    throw new ArgumentException();

                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), momentFittingOrder, HistInd).XQuadSchemeHelper;

                for (int iSpc = 0; iSpc < LsTrk.SpeciesIdS.Count; iSpc++) {
                    SpeciesId spcId = LsTrk.SpeciesIdS[iSpc];
                    double _rho = rho[iSpc];

                    var Uspc = Velocity.Select(u => (u as XDGField).GetSpeciesShadowField(spcId)).ToArray();
                    //ScalarFunctionEx spcKinDissip = GetSpeciesKineticDissipationFunc(Uspc, _rho);

                    ScalarFunctionEx spcKinEnergy = delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                        int K = result.GetLength(1); // No nof Nodes

                        MultidimensionalArray U_res = MultidimensionalArray.Create(Len, K, D);
                        for (int i = 0; i < D; i++) {
                            Uspc[i].Evaluate(i0, Len, nds, U_res.ExtractSubArrayShallow(-1, -1, i));
                        }

                        double acc;
                        for (int j = 0; j < Len; j++) {
                            for (int k = 0; k < K; k++) {

                                acc = 0.0;

                                for (int d = 0; d < D; d++) {
                                    acc += U_res[j, k, d] * U_res[j, k, d];
                                }
                                result[j, k] = _rho * acc / 2.0;
                            }
                        }

                    };

                    proj.GetSpeciesShadowField(spcId).ProjectField(spcKinEnergy);

                }
            }
        }


        static ScalarFunctionEx GetSpeciesPowerOfStressesFunc(DGField P, DGField[] U, DGField[,] GradU, double mu) {

            int D = U[0].Basis.GridDat.SpatialDimension;

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes

                //MultidimensionalArray P_res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray GradP_res = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray U_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray GradU_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray GGradU_res = MultidimensionalArray.Create(Len, K, D, D, D);

                //P.Evaluate(i0, Len, nds, P_res.ExtractSubArrayShallow(-1, -1));
                P.EvaluateGradient(i0, Len, nds, GradP_res.ExtractSubArrayShallow(-1, -1, -1));
                for (int i = 0; i < D; i++) {
                    U.ElementAt(i).Evaluate(i0, Len, nds, U_res.ExtractSubArrayShallow(-1, -1, i));
                    //U[i].EvaluateGradient(i0, Len, nds, GradU_res.ExtractSubArrayShallow(-1, -1, i, -1));
                    for (int j = 0; j < D; j++) {
                        GradU[i, j].Evaluate(i0, Len, nds, GradU_res.ExtractSubArrayShallow(-1, -1, i, j));
                        GradU[i, j].EvaluateGradient(i0, Len, nds, GGradU_res.ExtractSubArrayShallow(-1, -1, i, j, -1));
                    }
                }

                double acc;
                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        if (D > 2)
                            throw new NotImplementedException("TODO");

                        double[] DivGradUT = new double[D];
                        double[] LapU = new double[D];
                        for (int d = 0; d < D; d++) {
                            DivGradUT[d] += GGradU_res[j, k, 0, d, 0] + GGradU_res[j, k, 1, d, 1];
                            LapU[d] += GGradU_res[j, k, d, 0, 0] + GGradU_res[j, k, d, 1, 1];
                        }

                        acc = 0.0;
                        for (int d = 0; d < D; d++) {
                            // gradP * U
                            acc -= GradP_res[j, k, d] * U_res[j, k, d];
                            // div(gradU^T)
                            acc += mu * DivGradUT[d] * U_res[j, k, d];
                            // lapU
                            acc += mu * LapU[d] * U_res[j, k, d];
                            for (int dd = 0; dd < D; dd++) {
                                // gradU^T:gradU
                                acc += mu * GradU_res[j, k, dd, d] * GradU_res[j, k, d, dd];
                                // gradU:gradU
                                acc += mu * GradU_res[j, k, d, dd] * GradU_res[j, k, d, dd];
                            }
                        }
                        result[j, k] = acc;
                    }
                }

            };

        }


        public static void ProjectPowerOfStresses(this XDGField proj, LevelSetTracker LsTrk, DGField Pressure, DGField[] Velocity, double[] mu, int momentFittingOrder, int HistInd = 1) {
            using (new FuncTrace()) {

                int D = LsTrk.GridDat.SpatialDimension;
                if (Velocity.Count() != D) {
                    throw new ArgumentException();
                }
                if (LsTrk.SpeciesIdS.Count != mu.Length)
                    throw new ArgumentException();


                var SchemeHelper = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), momentFittingOrder, HistInd).XQuadSchemeHelper;

                for (int iSpc = 0; iSpc < LsTrk.SpeciesIdS.Count; iSpc++) {
                    SpeciesId spcId = LsTrk.SpeciesIdS[iSpc];
                    double muSpc = mu[iSpc];

                    SubGrid sf = LsTrk.Regions.GetSpeciesSubGrid(spcId);

                    var Pspc = (Pressure as XDGField).GetSpeciesShadowField(spcId);
                    var Uspc = Velocity.Select(u => (u as XDGField).GetSpeciesShadowField(spcId)).ToArray();

                    VectorField<DGField> GradVelX = new VectorField<DGField>(D, Velocity[0].Basis, XDGField.Factory);
                    VectorField<DGField> GradVelY = new VectorField<DGField>(D, Velocity[1].Basis, XDGField.Factory);
                    for (int d = 0; d < D; d++) {
                        DGField f_Spc;
                        f_Spc = ((Velocity[0] as XDGField).GetSpeciesShadowField(spcId));
                        (GradVelX[d] as XDGField).GetSpeciesShadowField(spcId).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                        f_Spc = ((Velocity[1] as XDGField).GetSpeciesShadowField(spcId));
                        (GradVelY[d] as XDGField).GetSpeciesShadowField(spcId).DerivativeByFlux(1.0, f_Spc, d, optionalSubGrid: sf);
                    }
                    GradVelX.ForEach(F => F.CheckForNanOrInf(true, true, true));
                    GradVelY.ForEach(F => F.CheckForNanOrInf(true, true, true));

                    DGField[,] GradU = new DGField[D, D];
                    GradU[0, 0] = GradVelX[0];
                    GradU[0, 1] = GradVelX[1];
                    GradU[1, 0] = GradVelY[0];
                    GradU[1, 1] = GradVelY[1];

                    ScalarFunctionEx spcPowStress = GetSpeciesPowerOfStressesFunc(Pspc, Uspc, GradU, muSpc);

                    proj.GetSpeciesShadowField(spcId).ProjectField(spcPowStress);

                }
            }
        }


        #endregion




        #region energy balance

        static ScalarFunctionEx GetFunc_EnergyBalanceAtInterface(XDGField p, VectorField<XDGField> u, VectorField<XDGField> gradUx, VectorField<XDGField> gradUy,
             ConventionalDGField[] uI, VectorField<SinglePhaseField> gradUIx, VectorField<SinglePhaseField> gradUIy, SinglePhaseField curv, PhysicalParameters physParam) {

            var LsTrk = p.Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            var pA = p.GetSpeciesShadowField("A");
            var pB = p.GetSpeciesShadowField("B");

            var uA = u.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var uB = u.Select(v => v.GetSpeciesShadowField("B")).ToArray();

            var gradUxA = gradUx.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var gradUxB = gradUx.Select(v => v.GetSpeciesShadowField("B")).ToArray();

            var gradUyA = gradUy.Select(v => v.GetSpeciesShadowField("A")).ToArray();
            var gradUyB = gradUy.Select(v => v.GetSpeciesShadowField("B")).ToArray();


            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster

                MultidimensionalArray pARes = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pBRes = MultidimensionalArray.Create(Len, K);

                MultidimensionalArray uARes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray uBRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray gradUxARes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUxBRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray gradUyARes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUyBRes = MultidimensionalArray.Create(Len, K, D);

                MultidimensionalArray uIRes = MultidimensionalArray.Create(Len, K, D);
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
                    pARes.Set(cp.EvaluateAtCp(pA));
                    pBRes.Set(cp.EvaluateAtCp(pB));
                    uARes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(uA[dd]));
                    uBRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(uB[dd]));
                    gradUxARes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUxA[dd]));
                    gradUxBRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUxB[dd]));
                    gradUyARes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUyA[dd]));
                    gradUyBRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUyB[dd]));
                    uIRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(uI[dd]));
                    gradUIxRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUIx[dd]));
                    gradUIyRes.ExtractSubArrayShallow(-1, -1, dd).Set(cp.EvaluateAtCp(gradUIy[dd]));
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

                        double Ej = 0;
                        for (int d = 0; d < _D; d++) {
                            // pressure
                            Ej -= (pBRes[j, k] * uBRes[j, k, d] - pARes[j, k] * uARes[j, k, d]) * Normals[j, k, d];
                            for (int dd = 0; dd < _D; dd++) {
                                // velocity gradients
                                Ej += (muB * gradUBRes[j, k, d, dd] * uBRes[j, k, dd] - muA * gradUARes[j, k, d, dd] * uARes[j, k, dd]) * Normals[j, k, d];
                                Ej += (muB * gradUBRes[j, k, dd, d] * uBRes[j, k, dd] - muA * gradUARes[j, k, dd, d] * uARes[j, k, dd]) * Normals[j, k, d];
                                // surface divergence
                                if (dd == d) {
                                    Ej += sigma * (1.0 - Normals[j, k, d] * Normals[j, k, dd]) * gradUIRes[j, k, d, dd];
                                } else {
                                    Ej += sigma * (0.0 - Normals[j, k, d] * Normals[j, k, dd]) * gradUIRes[j, k, d, dd];
                                }
                            }
                            // surface tension
                            Ej -= sigma * curvRes[j, k] * uIRes[j, k, d] * Normals[j, k, d];
                        }
                        result[j, k] = Ej;
                    }
                }

            };

        }


        public static void ProjectEnergyBalanceAtInterface(this SinglePhaseField err, double alpha, XDGField p, VectorField<XDGField> u, ConventionalDGField[] uI, SinglePhaseField curv, PhysicalParameters physParam, CellQuadratureScheme quadScheme = null) {

            var LsTrk = u[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            VectorField<XDGField> GuX = new VectorField<XDGField>(D, u[0].Basis, XDGField.Factory);
            VectorField<XDGField> GuY = new VectorField<XDGField>(D, u[1].Basis, XDGField.Factory);

            foreach (SpeciesId spcId in LsTrk.SpeciesIdS) {
                SubGrid sf = LsTrk.Regions.GetSpeciesSubGrid(spcId);
                DGField f0_Spc = u[0].GetSpeciesShadowField(spcId);
                DGField f1_Spc = u[1].GetSpeciesShadowField(spcId);
                for (int dd = 0; dd < D; dd++) {
                    GuX[dd].GetSpeciesShadowField(spcId).DerivativeByFlux(1.0, f0_Spc, dd, optionalSubGrid: sf);
                    GuY[dd].GetSpeciesShadowField(spcId).DerivativeByFlux(1.0, f1_Spc, dd, optionalSubGrid: sf);
                }
            }
            GuX.ForEach(F => F.CheckForNanOrInf(true, true, true));
            GuY.ForEach(F => F.CheckForNanOrInf(true, true, true));


            VectorField<SinglePhaseField> GuIx = new VectorField<SinglePhaseField>(D, uI[0].Basis, SinglePhaseField.Factory);
            VectorField<SinglePhaseField> GuIy = new VectorField<SinglePhaseField>(D, uI[1].Basis, SinglePhaseField.Factory);
            for (int dd = 0; dd < D; dd++) {
                GuIx[dd].DerivativeByFlux(1.0, uI[0], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
                GuIy[dd].DerivativeByFlux(1.0, uI[1], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
            }
            GuIx.ForEach(F => F.CheckForNanOrInf(true, true, true));
            GuIy.ForEach(F => F.CheckForNanOrInf(true, true, true));



            ScalarFunctionEx ErrFunc = GetFunc_EnergyBalanceAtInterface(p, u, GuX, GuY, uI, GuIx, GuIy, curv, physParam);

            int order = ((u[0].Basis.Degree - 1) + (p.Basis.Degree) + u[0].Basis.Degree + err.Basis.Degree + 2);
            if (quadScheme == null)
                quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

            err.ProjectField(alpha, ErrFunc,
                quadScheme);



        }


        #endregion


        #region interface terms

        static ScalarFunctionEx GetFunc_InterfaceDivergence(LevelSetTracker LsTrk, VectorField<SinglePhaseField> gradUIx, VectorField<SinglePhaseField> gradUIy, PhysicalParameters physParam) {

            int D = LsTrk.GridDat.SpatialDimension;

            return delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                int _D = D; // local var may be a bit faster

                MultidimensionalArray gradUIxRes = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray gradUIyRes = MultidimensionalArray.Create(Len, K, D);

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
                }

                MultidimensionalArray gradUIRes = MultidimensionalArray.Create(Len, K, D, D);
                gradUIRes.ExtractSubArrayShallow(-1, -1, 0, 0).Set(gradUIxRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 0, 1).Set(gradUIxRes.ExtractSubArrayShallow(-1, -1, 1));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 1, 0).Set(gradUIyRes.ExtractSubArrayShallow(-1, -1, 0));
                gradUIRes.ExtractSubArrayShallow(-1, -1, 1, 1).Set(gradUIyRes.ExtractSubArrayShallow(-1, -1, 1));


                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                double sigma = physParam.Sigma;

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double divI = 0;
                        for (int d = 0; d < _D; d++) {
                            for (int dd = 0; dd < _D; dd++) {
                                // surface divergence
                                if (dd == d) {
                                    divI += sigma * (1.0 - Normals[j, k, d] * Normals[j, k, dd]) * gradUIRes[j, k, d, dd];
                                } else {
                                    divI += sigma * (0.0 - Normals[j, k, d] * Normals[j, k, dd]) * gradUIRes[j, k, d, dd];
                                }
                            }
                        }
                        result[j, k] = divI;
                    }
                }

            };

        }

        public static void ProjectInterfaceDivergence(this SinglePhaseField divI, double alpha, ConventionalDGField[] uI, LevelSetTracker LsTrk, PhysicalParameters physParam, CellQuadratureScheme quadScheme = null) {

            int D = LsTrk.GridDat.SpatialDimension;

            VectorField<SinglePhaseField> GuIx = new VectorField<SinglePhaseField>(D, uI[0].Basis, SinglePhaseField.Factory);
            VectorField<SinglePhaseField> GuIy = new VectorField<SinglePhaseField>(D, uI[1].Basis, SinglePhaseField.Factory);
            for (int dd = 0; dd < D; dd++) {
                GuIx[dd].DerivativeByFlux(1.0, uI[0], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
                GuIy[dd].DerivativeByFlux(1.0, uI[1], dd, optionalSubGrid: LsTrk.Regions.GetCutCellSubGrid());
            }
            GuIx.ForEach(F => F.CheckForNanOrInf(true, true, true));
            GuIy.ForEach(F => F.CheckForNanOrInf(true, true, true));



            ScalarFunctionEx ErrFunc = GetFunc_InterfaceDivergence(LsTrk, GuIx, GuIy, physParam);

            int order = ((uI[0].Basis.Degree - 1) + divI.Basis.Degree + 2);
            if (quadScheme == null)
                quadScheme = (new CellQuadratureScheme(false, LsTrk.Regions.GetCutCellMask())).AddFixedOrderRules(LsTrk.GridDat, order);

            divI.ProjectField(alpha, ErrFunc,
                quadScheme);

        }

        #endregion


        #region 

        static ScalarFunctionEx GetEnergyBalanceFunc(XDGField P, VectorField<XDGField> U, ConventionalDGField[] Umean, SinglePhaseField C, double muA, double muB, double sigma, bool squared) {

            int D = P.Basis.GridDat.SpatialDimension;

            ConventionalDGField pA = P.GetSpeciesShadowField("A");
            ConventionalDGField pB = P.GetSpeciesShadowField("B");

            var UA = U.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UB = U.Select(u => u.GetSpeciesShadowField("B")).ToArray();

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray pA_res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pB_res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray UA_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray UB_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray GradUA_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray GradUB_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray U_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray GradU_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray Curv_res = MultidimensionalArray.Create(Len, K);

                pA.Evaluate(i0, Len, nds, pA_res);
                pB.Evaluate(i0, Len, nds, pB_res);

                for (int i = 0; i < D; i++) {
                    UA[i].Evaluate(i0, Len, nds, UA_res.ExtractSubArrayShallow(-1, -1, i));
                    UB[i].Evaluate(i0, Len, nds, UB_res.ExtractSubArrayShallow(-1, -1, i));
                    Umean[i].Evaluate(i0, Len, nds, U_res.ExtractSubArrayShallow(-1, -1, i));

                    UA[i].EvaluateGradient(i0, Len, nds, GradUA_res.ExtractSubArrayShallow(-1, -1, i, -1));
                    UB[i].EvaluateGradient(i0, Len, nds, GradUB_res.ExtractSubArrayShallow(-1, -1, i, -1));
                    Umean[i].EvaluateGradient(i0, Len, nds, GradU_res.ExtractSubArrayShallow(-1, -1, i, -1));
                }

                C.Evaluate(i0, Len, nds, Curv_res);

                var Normals = P.Basis.Tracker.DataHistories[0].Current.GetLevelSetNormals(nds, i0, Len);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double acc = 0.0;

                        for (int d = 0; d < D; d++) {

                            // enrgy jump at interface
                            acc -= (pB_res[j, k] * UB_res[j, k, d] - pA_res[j, k] * UA_res[j, k, d]) * Normals[j, k, d];

                            for (int dd = 0; dd < D; dd++) {
                                acc += (muB * GradUB_res[j, k, d, dd] * UB_res[j, k, dd] - muA * GradUA_res[j, k, d, dd] * UA_res[j, k, dd]) * Normals[j, k, d];
                                acc += (muB * GradUB_res[j, k, dd, d] * UB_res[j, k, dd] - muA * GradUA_res[j, k, dd, d] * UA_res[j, k, dd]) * Normals[j, k, d];     // Transposed Term
                            }

                            // surface energy changerate
                            //for (int dd = 0; dd < D; dd++) {
                            //    if (dd == d) {
                            //        acc += sigma * (1 - Normals[j, k, d] * Normals[j, k, dd]) * GradU_res[j, k, dd, d];
                            //    } else {
                            //        acc += sigma * (-Normals[j, k, d] * Normals[j, k, dd]) * GradU_res[j, k, dd, d];
                            //    }
                            //}

                            // curvature energy
                            acc -= sigma * Curv_res[j, k] * U_res[j, k, d] * Normals[j, k, d];
                        }

                        if (squared) {
                            result[j, k] = acc.Pow2();
                        } else {
                            result[j, k] = acc;
                        }
                    }
                }

            };

        }

        public static double EnergyBalanceNormAtInterface(XDGField P, VectorField<XDGField> U, ConventionalDGField[] Umean, SinglePhaseField C, double muA, double muB, double sigma, int momentFittingOrder) {

            LevelSetTracker LsTrk = P.Basis.Tracker;

            double energyBal_Norm = 0.0;

            ScalarFunctionEx energyBalFunc = GetEnergyBalanceFunc(P, U, Umean, C, muA, muB, sigma, true);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    energyBalFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        energyBal_Norm += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return energyBal_Norm.Sqrt();

        }

        public static void ProjectEnergyBalanceNorm(this SinglePhaseField err, double alpha, XDGField P, VectorField<XDGField> U, ConventionalDGField[] Umean, SinglePhaseField C,
            double muA, double muB, double sigma, int momentFittingOrder) {

            var LsTrk = U[0].Basis.Tracker;
            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx ErrFunc = GetEnergyBalanceFunc(P, U, Umean, C, muA, muB, sigma, true);

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



        static ScalarFunctionEx GetInterfaceDivergenceFunc(LevelSetTracker LsTrk, ConventionalDGField[] uI, bool squared) {

            int D = LsTrk.GridDat.SpatialDimension;

            return delegate (int i0, int Len, NodeSet nds, MultidimensionalArray result) {

                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray GradU_Res = MultidimensionalArray.Create(Len, K, D, D);

                for (int i = 0; i < D; i++) {
                    uI.ElementAt(i).EvaluateGradient(i0, Len, nds, GradU_Res.ExtractSubArrayShallow(-1, -1, i, -1));
                }

                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(nds, i0, Len);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        MultidimensionalArray Nsurf = Normals.ExtractSubArrayShallow(j, k, -1);
                        double[,] Psurf = new double[D, D];
                        for (int d1 = 0; d1 < D; d1++) {
                            for (int d2 = 0; d2 < D; d2++) {
                                if (d2 == d1)
                                    Psurf[d1, d2] = (1 - Nsurf[d1] * Nsurf[d2]);
                                else
                                    Psurf[d1, d2] = (0 - Nsurf[d1] * Nsurf[d2]);
                            }
                        }

                        double acc = 0.0;

                        for (int d1 = 0; d1 < D; d1++) {
                            for (int dd = 0; dd < D; dd++) {
                                acc += Psurf[d1, dd] * GradU_Res[j, k, dd, d1];
                            }
                        }

                        if (squared) {
                            result[j, k] = acc.Pow2();
                        } else {
                            result[j, k] = acc;
                        }
                    }
                }

            };

        }

        public static double GetSurfaceChangerate(LevelSetTracker LsTrk, ConventionalDGField[] uI, int momentFittingOrder) {

            double SurfChangerate = 0.0;

            ScalarFunctionEx surfChangerateFunc = GetInterfaceDivergenceFunc(LsTrk, uI, false);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingOrder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingOrder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    surfChangerateFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        SurfChangerate += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            return SurfChangerate;

        }



        public static ScalarFunctionEx GetEnergyJumpFunc(LevelSetTracker LsTrk, VectorField<XDGField> Velocity, XDGField Pressure, double muA, double muB, bool squared) {

            var UA = Velocity.Select(u => u.GetSpeciesShadowField("A")).ToArray();
            var UB = Velocity.Select(u => u.GetSpeciesShadowField("B")).ToArray();

            ConventionalDGField pA = null, pB = null;
            bool UsePressure = Pressure != null;
            if (UsePressure) {
                pA = Pressure.GetSpeciesShadowField("A");
                pB = Pressure.GetSpeciesShadowField("B");
            }

            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx EnergyJumpFunc = delegate (int j0, int Len, NodeSet Ns, MultidimensionalArray result) {
                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray UA_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray UB_res = MultidimensionalArray.Create(Len, K, D);
                MultidimensionalArray GradUA_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray GradUB_res = MultidimensionalArray.Create(Len, K, D, D);
                MultidimensionalArray pA_res = MultidimensionalArray.Create(Len, K);
                MultidimensionalArray pB_res = MultidimensionalArray.Create(Len, K);

                for (int i = 0; i < D; i++) {
                    UA[i].Evaluate(j0, Len, Ns, UA_res.ExtractSubArrayShallow(-1, -1, i));
                    UB[i].Evaluate(j0, Len, Ns, UB_res.ExtractSubArrayShallow(-1, -1, i));

                    UA[i].EvaluateGradient(j0, Len, Ns, GradUA_res.ExtractSubArrayShallow(-1, -1, i, -1));
                    UB[i].EvaluateGradient(j0, Len, Ns, GradUB_res.ExtractSubArrayShallow(-1, -1, i, -1));
                }
                if (UsePressure) {
                    pA.Evaluate(j0, Len, Ns, pA_res);
                    pB.Evaluate(j0, Len, Ns, pB_res);
                } else {
                    pA_res.Clear();
                    pB_res.Clear();
                }

                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(Ns, j0, Len);

                for (int j = 0; j < Len; j++) {
                    for (int k = 0; k < K; k++) {

                        double acc = 0.0;

                        for (int d = 0; d < D; d++) {
                            // pressure
                            if (UsePressure) {
                                acc += (pB_res[j, k] * UB_res[j, k, d] - pA_res[j, k] * UA_res[j, k, d]) * Normals[j, k, d];
                            }

                            // Nabla U + (Nabla U) ^T
                            for (int dd = 0; dd < D; dd++) {
                                acc -= (muB * GradUB_res[j, k, d, dd] * UB_res[j, k, dd] - muA * GradUA_res[j, k, d, dd] * UA_res[j, k, dd]) * Normals[j, k, d];
                                acc -= (muB * GradUB_res[j, k, dd, d] * UB_res[j, k, dd] - muA * GradUA_res[j, k, dd, d] * UA_res[j, k, dd]) * Normals[j, k, d];     // Transposed Term
                            }

                        }
                        if (squared) {
                            result[j, k] = acc.Pow2();
                        } else {
                            result[j, k] = acc;
                        }

                    }
                }
            };

            return EnergyJumpFunc;

        }

        public static double EnergyJumpAtInterface(LevelSetTracker LsTrk, VectorField<XDGField> Velocity, XDGField Pressure, double muA, double muB, bool Norm, int momentFittingorder) {

            double EnergyJump = 0.0;

            ScalarFunctionEx EnergyJumpFunc = GetEnergyJumpFunc(LsTrk, Velocity, Pressure, muA, muB, Norm);

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingorder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingorder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EnergyJumpFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        EnergyJump += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            if (Norm) {
                EnergyJump.Sqrt();
            }

            return EnergyJump;

        }


        public static double SurfaceEnergyChangerate(LevelSetTracker LsTrk, ConventionalDGField[] uI, double sigma, bool Norm, int momentFittingorder) {

            double Changerate_Surface = 0.0;

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingorder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            ScalarFunctionEx SurfaceChangerate = GetInterfaceDivergenceFunc(LsTrk, uI, Norm);

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingorder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    SurfaceChangerate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        Changerate_Surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            double Changerate_Esurf;
            if (Norm) {
                Changerate_Esurf = sigma * Changerate_Surface.Sqrt();
            } else {
                Changerate_Esurf = sigma * Changerate_Surface;
            }

            return Changerate_Esurf;

        }


        public static ScalarFunctionEx GetCurvatureEnergyFunc(LevelSetTracker LsTrk, DGField Curvature, double sigma, ConventionalDGField[] uI, bool ExtVel, bool squared) {

            int D = LsTrk.GridDat.SpatialDimension;

            ScalarFunctionEx CurvatureEnergyFunc = delegate (int i0, int Length, NodeSet nds, MultidimensionalArray result) {

                Curvature.Evaluate(i0, Length, nds, result);

                int K = result.GetLength(1); // No nof Nodes
                MultidimensionalArray U_res = MultidimensionalArray.Create(Length, K, D);

                for (int i = 0; i < D; i++) {
                    uI.ElementAt(i).Evaluate(i0, Length, nds, U_res.ExtractSubArrayShallow(-1, -1, i));
                }

                var Normals = LsTrk.DataHistories[0].Current.GetLevelSetNormals(nds, i0, Length);

                for (int j = 0; j < Length; j++) {
                    for (int k = 0; k < K; k++) {

                        double acc = result[j, k];

                        for (int d = 0; d < D; d++) {
                            // U * N
                            if (!ExtVel) {
                                acc *= U_res[j, k, d] * Normals[j, k, d];
                            } else {
                                acc *= U_res[j, k, d];
                            }

                        }

                        if (squared) {
                            result[j, k] = (sigma * acc).Pow2();
                        } else {
                            result[j, k] = sigma * acc;
                        }
                    }
                }

            };

            return CurvatureEnergyFunc;

        }

        public static double CurvatureEnergy(LevelSetTracker LsTrk, SinglePhaseField Curvature, double sigma, ConventionalDGField[] uI, bool ExtVel, bool Norm, int momentFittingorder) {

            double EnergyCurv = 0.0;

            var SchemeHelper = LsTrk.GetXDGSpaceMetrics(new[] { LsTrk.GetSpeciesId("A") }, momentFittingorder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());

            ScalarFunctionEx CurvEnergyFunc = GetCurvatureEnergyFunc(LsTrk, Curvature, sigma, uI, ExtVel, Norm);

            CellQuadrature.GetQuadrature(new int[] { 1 }, LsTrk.GridDat,
                cqs.Compile(LsTrk.GridDat, momentFittingorder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    CurvEnergyFunc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));

                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        EnergyCurv += ResultsOfIntegration[i, 0];
                }
            ).Execute();

            if (Norm) {
                EnergyCurv.Sqrt();
            }

            return EnergyCurv;

        }


        #endregion

    }
}
