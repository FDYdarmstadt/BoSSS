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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools;
using ilPSP.LinSolvers.PARDISO;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.XNSECommon {


    /// <summary>
    /// A collection of algorithms to evaluate the curvature of a level-set.
    /// </summary>
    public static class CurvatureAlgorithms {

        

        public enum GradientOption {
            //external = 0,
            LevSet = 1,
            FiltLevSet = 2,
            ByFluxLevSet = 3,
            ByFluxFiltLevSet = 4
        }

        public enum HessianOption {
            LevSetH = 1,
            FiltLevSetH = 2,
            LevSetGrad = 3,
            FiltLevSetGrad = 4
        }

        public enum LevelSetSource {
            fromDG,
            fromC0,
            fromReinit
        }

        [DataContract]
        public class FilterConfiguration {

            /// <summary>
            /// Returns a configuration with all filters turned off.
            /// </summary>
            public static FilterConfiguration NoFilter {
                get {
                    FilterConfiguration r = new FilterConfiguration();
                    r.gradOpt = GradientOption.LevSet;
                    r.hessOpt = HessianOption.LevSetGrad;
                    r.useFiltLevSetGrad = false;
                    r.useFiltLevSetHess = false;
                    r.FilterCurvatureCycles = 0;
                    r.LevelSetSource = LevelSetSource.fromC0;
                    r.PatchRecoveryDomWidth = 0;
                    r.NoOfPatchRecoverySweeps = 0;
                    r.CurvatureLimiting = false;

                    return r;
                }
            }

            /// <summary>
            /// Returns a recommended filter configuration.
            /// </summary>
            public static FilterConfiguration Default {
                get {
                    FilterConfiguration r = new FilterConfiguration();
                    return r;
                }
            }


            [DataMember]
            public GradientOption gradOpt = GradientOption.FiltLevSet;

            [DataMember]
            public HessianOption hessOpt = HessianOption.FiltLevSetGrad;

            [DataMember]
            public bool useFiltLevSetGrad = true;

            [DataMember]
            public bool useFiltLevSetHess = true;

            [DataMember]
            public int FilterCurvatureCycles = 0;

            [DataMember]
            public LevelSetSource LevelSetSource = LevelSetSource.fromDG;

            [DataMember]
            public int PatchRecoveryDomWidth = 0;

            [DataMember]
            public int NoOfPatchRecoverySweeps = 3;

            [DataMember]
            public bool CurvatureLimiting = false;
        }



        /// <summary>
        /// The main driver routine, which provides unified access to all implemented variants.
        /// </summary>
        /// <param name="momentFittingVariant"></param>
        /// <param name="surfTenM">
        /// Configuration of surface tension implementation;
        /// </param>
        /// <param name="Curvature">
        /// Output: on exit, an (hopefully good) approximation to the curvature of the level set.
        /// Only used for curvature-based versions of <paramref name="surfTenM"/>.
        /// </param>
        /// <param name="LevelSetGradient">
        /// output: the filtered level-set gradient .
        /// Only used for Laplace-Beltrami-based versions of <paramref name="surfTenM"/>.
        /// </param>
        /// <param name="LsTrk"></param>
        /// <param name="DG_LevSet">
        /// The DG-representation of the Level-Set; this can be discontinuous, and curvature based on 
        /// the DG-Level set is usually more accurate.
        /// </param>
        /// <param name="config"></param>
        /// <param name="HMForder">
        /// If HMF is used anywhere, thats the order that will be used.
        /// </param>
        /// <returns></returns>
        static public SurfaceStressTensor_IsotropicMode CurvatureDriver(
            SurfaceStressTensor_IsotropicMode surfTenM, FilterConfiguration config,
            SinglePhaseField Curvature, out VectorField<SinglePhaseField> LevelSetGradient,
            LevelSetTracker LsTrk,
            int HMForder,
            SinglePhaseField DG_LevSet) //
       {

            
            CellMask CC = LsTrk.Regions.GetNearFieldMask(config.PatchRecoveryDomWidth);

            CurvatureBasedSurfaceTension.hmin = LsTrk.GridDat.Cells.h_minGlobal;

            Curvature.Clear();

            // select which level-set will be used for curvature computation
            // =============================================================

            var levSet = SourceLevelSet(config, LsTrk, DG_LevSet);

            // compute gradients, Hessian and apply Filters
            // =============================================
            SinglePhaseField[] G; // gradient
            SinglePhaseField[,] H; // Hessian
            L2PatchRecovery l2pr;
            Filter(levSet, Curvature.Basis, CC, surfTenM, config, out l2pr, out G, out H);

            // now, do all the rest...
            // =======================


            if(surfTenM == SurfaceStressTensor_IsotropicMode.Curvature_Projected) {


                int[] lim;
                ProjectTotalCurvature(G, H, Curvature, CC, out lim);

                if(!config.CurvatureLimiting)
                    lim = null;

                if(config.CurvatureLimiting)
                    PatchRecMultipassFilter(Curvature, config.FilterCurvatureCycles, l2pr, lim);

                Curvature.CheckForNanOrInf();


                //Basis BasisForFiltered = Curvature.Basis;
                //PatchRecFilteredCurvature(levSet,
                //    Curvature,
                //    LsTrk.GetNearFieldMask(config.PatchRecoveryDomWidth),
                //    config.NoOfPatchRecoverySweeps,
                //    config.gradOpt,
                //    config.hessOpt,
                //    config.useFiltLevSetGrad,
                //    config.useFiltLevSetHess,
                //    config.FilterCurvatureCycles,
                //    config.CurvatureLimiting);
                //Curvature.CheckForNanOrInf(true, true, true);

                LevelSetGradient = new VectorField<SinglePhaseField>(G); 

                return SurfaceStressTensor_IsotropicMode.Curvature_Projected;

            } else if(surfTenM == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint) {


                int[] lim;
                ClosestPointTotalCurvature(LsTrk, G, H, Curvature, out lim);

                if(!config.CurvatureLimiting)
                    lim = null;

                if(config.CurvatureLimiting)
                    PatchRecMultipassFilter(Curvature, config.FilterCurvatureCycles, l2pr, lim);

                Curvature.CheckForNanOrInf();

                
                LevelSetGradient = null;// new VectorField<SinglePhaseField>(G); 

                return SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint;

            } else if(surfTenM == SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean) {

                LaplaceBeltramiFiltered(
                    LsTrk,
                    Curvature,
                    config.FilterCurvatureCycles, 
                    G,
                    HMForder);
                Curvature.CheckForNanOrInf(true, true, true);

                LevelSetGradient = null;
                return SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean;

            } 

            throw new ArgumentOutOfRangeException("unknown AlgorithmClass");

        }

        private static SinglePhaseField SourceLevelSet(FilterConfiguration config, LevelSetTracker LsTrk,
            SinglePhaseField DG_LevSet)
        {
            SinglePhaseField levSet;
            SinglePhaseField C0_LevSet = (SinglePhaseField) LsTrk.LevelSets[0];
            switch (config.LevelSetSource)
            {
                case LevelSetSource.fromDG:
                    levSet = DG_LevSet;
                    break;
                case LevelSetSource.fromC0:
                    levSet = C0_LevSet;
                    break;
                case LevelSetSource.fromReinit:
                {
                    levSet = new SinglePhaseField(DG_LevSet.Basis);

                    QuadRule[] rulez = LsTrk.GridDat.Grid.RefElements.Select(
                        Kref => Kref.GetQuadratureRule(levSet.Basis.Degree*2)).ToArray();

                    var cpf = new ClosestPointFinder(LsTrk, 0, LsTrk.Regions.GetCutCellSubGrid(), rulez.Select(R => R.Nodes));

                    var jCell_2_jsub = cpf.sgrd.LocalCellIndex2SubgridIndex;
                    var Dist = cpf.Distance;

                    levSet.Clear(cpf.sgrd.VolumeMask);
                    levSet.ProjectField(1.0,
                        delegate(int j0, int Len, NodeSet NS, MultidimensionalArray result)
                        {
                            int K = result.GetLength(1);

                            for (int j = j0; j < (j0 + Len); j++)
                            {
                                int jSub = jCell_2_jsub[j];

                                for (int k = 0; k < K; k++)
                                {
                                    double r = Dist[jSub, k];
                                    //result[j - j0, k] = r; // signed distance
                                    result[j - j0, k] = 0.25*(r + 2)*(r + 2) - 1; // quadratic distance
                                }
                            }
                        },
                        new CellQuadratureScheme(false, cpf.sgrd.VolumeMask).AddFixedRuleS(rulez));

                    break;
                }
                default:
                    throw new ArgumentException();
            }
            return levSet;
        }


        static public SurfaceStressTensor_IsotropicMode LaplaceBeltramiDriver(

            SurfaceStressTensor_IsotropicMode surfTenM, FilterConfiguration config,
            out VectorField<SinglePhaseField> LevelSetGradient,
            LevelSetTracker LsTrk,
            SinglePhaseField DG_LevSet
            ){

            var levSet = SourceLevelSet(config, LsTrk, DG_LevSet);
            CellMask CC = LsTrk.Regions.GetNearFieldMask(config.PatchRecoveryDomWidth);

            // compute gradients, Hessian and apply Filters
            // =============================================

            SinglePhaseField[] G; // gradient
            SinglePhaseField[,] H; // Hessian
            L2PatchRecovery l2pr;
            Filter(levSet, levSet.Basis, CC, surfTenM, config, out l2pr, out G, out H);

            LevelSetGradient = new VectorField<SinglePhaseField>(G);
            return surfTenM;

        }




        private static void Filter(
            SinglePhaseField LevSet,
            Basis patchRecoveryBasis,
            CellMask CC,
            SurfaceStressTensor_IsotropicMode surfTenM, FilterConfiguration config,
            out L2PatchRecovery _l2pr,
            out SinglePhaseField[] G,
            out SinglePhaseField[,] H
            )
        {
            using (new FuncTrace()) {

                // init
                // ====

                var gdat = LevSet.GridDat;
                int D = gdat.SpatialDimension;
                Basis BasisForUnfiltered = LevSet.Basis;
                    

                // determine which vars we have to filter
                // ======================================

                bool HessReqAtAll = (surfTenM == SurfaceStressTensor_IsotropicMode.Curvature_Projected || surfTenM == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint);

                bool FiltLevSetHess_Req = config.useFiltLevSetHess && HessReqAtAll;
                bool FiltLevSetGrad_Req = config.useFiltLevSetGrad || ((config.hessOpt == HessianOption.FiltLevSetGrad) && HessReqAtAll);
                bool FiltLevSet_Req = (config.gradOpt == GradientOption.FiltLevSet)
                                        || ((config.hessOpt == HessianOption.FiltLevSetH) && HessReqAtAll);



                // Set Up Patch Recovery
                L2PatchRecovery l2pr = null;
                if(FiltLevSet_Req || FiltLevSetGrad_Req || FiltLevSetHess_Req || config.FilterCurvatureCycles > 0) {
                    l2pr = new L2PatchRecovery(patchRecoveryBasis, patchRecoveryBasis, CC, true);
                }
                _l2pr = l2pr;
                

                // Compute Filtered LevelSet
                SinglePhaseField FiltLevSet = FiltLevSet_Req ? new SinglePhaseField(patchRecoveryBasis) : null;
                
                if (FiltLevSet_Req) {
                    FiltLevSet.Clear();
                    FiltLevSet.AccLaidBack(1.0, LevSet);
                    PatchRecMultipassFilter(FiltLevSet, config.NoOfPatchRecoverySweeps, l2pr);
                }


                // Compute Gradient
                SinglePhaseField[] LevSetGrad;
                SinglePhaseField[] FiltLevSetGrad;
                ComputeGradient(LevSet, CC, config, D, out LevSetGrad, FiltLevSet, FiltLevSetGrad_Req, out FiltLevSetGrad, l2pr);



                // Compute Hessian
                SinglePhaseField[,] LevSetHess;
                SinglePhaseField[,] FiltLevSetHess;
                if (HessReqAtAll) {
                    ComputeHessian(LevSet, CC, config, out LevSetHess, out FiltLevSetHess, FiltLevSet, LevSetGrad, FiltLevSetGrad, l2pr);
                } else {
                    LevSetHess = null;
                    FiltLevSetHess = null;
                }

                // 
                G = config.useFiltLevSetGrad ? FiltLevSetGrad : LevSetGrad;
                H = config.useFiltLevSetHess ? FiltLevSetHess : LevSetHess;
            }
        }

        private static void ComputeGradient(SinglePhaseField LevSet, CellMask CC, FilterConfiguration config, int D,
            out SinglePhaseField[] LevSetGrad, SinglePhaseField FiltLevSet, bool FiltLevSetGrad_Req,
            out SinglePhaseField[] FiltLevSetGrad, L2PatchRecovery l2pr) {
            LevSetGrad = new SinglePhaseField[D];
            FiltLevSetGrad = new SinglePhaseField[D];

            switch (config.gradOpt) {
                case GradientOption.LevSet:
                for (int d = 0; d < D; d++) {
                    LevSetGrad[d] = new SinglePhaseField(LevSet.Basis, string.Format("G_{0}", d));
                    LevSetGrad[d].Derivative(1.0, LevSet, d, CC);
                }
                break;
                case GradientOption.ByFluxLevSet:
                for (int d = 0; d < D; d++) {
                    LevSetGrad[d] = new SinglePhaseField(LevSet.Basis, string.Format("G_{0}", d));
                    LevSetGrad[d].DerivativeByFlux(1.0, LevSet, d, new SubGrid(CC));
                }
                break;
                case GradientOption.FiltLevSet:
                for (int d = 0; d < D; d++) {
                    LevSetGrad[d] = new SinglePhaseField(FiltLevSet.Basis, string.Format("G_{0}", d));
                    LevSetGrad[d].Derivative(1.0, FiltLevSet, d, CC);
                }
                break;
                case GradientOption.ByFluxFiltLevSet:
                for (int d = 0; d < D; d++) {
                    LevSetGrad[d] = new SinglePhaseField(FiltLevSet.Basis, string.Format("G_{0}", d));
                    LevSetGrad[d].DerivativeByFlux(1.0, FiltLevSet, d, new SubGrid(CC));
                }
                break;
                default:
                throw new NotImplementedException();
            }


            for (int d = 0; d < D && FiltLevSetGrad_Req; d++) {
                FiltLevSetGrad[d] = new SinglePhaseField(FiltLevSet.Basis, string.Format("~G({0})", d));
                FiltLevSetGrad[d].AccLaidBack(1.0, LevSetGrad[d]);
                PatchRecMultipassFilter(FiltLevSetGrad[d], config.NoOfPatchRecoverySweeps, l2pr);
            }
        }

        private static void ComputeHessian(SinglePhaseField LevSet, CellMask CC, FilterConfiguration config,
            out SinglePhaseField[,] LevSetHess, out SinglePhaseField[,] FiltLevSetHess, 
            SinglePhaseField FiltLevSet, SinglePhaseField[] LevSetGrad,
            SinglePhaseField[] FiltLevSetGrad, L2PatchRecovery l2pr)
        {
            int D = LevSet.GridDat.SpatialDimension;
            LevSetHess = new SinglePhaseField[D, D];
            FiltLevSetHess = new SinglePhaseField[D, D];

            switch (config.hessOpt) {
                case HessianOption.LevSetH:
                for (int d1 = 0; d1 < D; d1++) {
                    for (int d2 = 0; d2 < D; d2++) {
                        //LevSetHess[d1, d2].Clear();
                        LevSetHess[d1, d2] = new SinglePhaseField(LevSet.Basis, string.Format("H({0},{1})", d1, d2));
                        LevSetHess[d1, d2].ProjectField(1.0,
                            delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                                var bffr = MultidimensionalArray.Create(Len, result.GetLength(1), D, D);
                                LevSet.EvaluateHessian(j0, Len, NS, bffr);
                                result.Set(bffr.ExtractSubArrayShallow(-1, -1, d1, d2));
                            }, new CellQuadratureScheme(UseDefaultFactories: true, domain: CC));
                    }
                }
                break;
                case HessianOption.FiltLevSetH: // 
                for (int d1 = 0; d1 < D; d1++) {
                    for (int d2 = 0; d2 < D; d2++) {
                        //LevSetHess[d1, d2].Clear();
                        LevSetHess[d1, d2] = new SinglePhaseField(FiltLevSet.Basis, string.Format("H({0},{1})", d1, d2));
                        LevSetHess[d1, d2].ProjectField(1.0,
                            delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                                var bffr = MultidimensionalArray.Create(Len, result.GetLength(1), D, D);
                                FiltLevSet.EvaluateHessian(j0, Len, NS, bffr);
                                result.Set(bffr.ExtractSubArrayShallow(-1, -1, d1, d2));
                            }, new CellQuadratureScheme(UseDefaultFactories: true, domain: CC));
                    }
                }
                break;
                case HessianOption.LevSetGrad:
                for (int d1 = 0; d1 < D; d1++) {
                    for (int d2 = 0; d2 < D; d2++) {
                        //LevSetHess[d1, d2].Clear();
                        LevSetHess[d1, d2] = new SinglePhaseField(LevSetGrad[d1].Basis, string.Format("H({0},{1})", d1, d2));
                        LevSetHess[d1, d2].Derivative(1.0, LevSetGrad[d1], d2, CC);
                    }
                }
                break;
                case HessianOption.FiltLevSetGrad:
                for (int d1 = 0; d1 < D; d1++) {
                    for (int d2 = 0; d2 < D; d2++) {
                        //LevSetHess[d1, d2].Clear();
                        LevSetHess[d1, d2] = new SinglePhaseField(FiltLevSetGrad[d1].Basis, string.Format("H({0},{1})", d1, d2));
                        LevSetHess[d1, d2].Derivative(1.0, FiltLevSetGrad[d1], d2, CC);
                    }
                }
                break;

                default:
                throw new NotImplementedException();
            }
            

            for (int d1 = 0; d1 < D && config.useFiltLevSetHess; d1++) {
                for (int d2 = 0; d2 < D; d2++) {
                    FiltLevSetHess[d1, d2] = new SinglePhaseField(FiltLevSetGrad[d1].Basis, string.Format("~H({0},{1})", d1, d2));
                    FiltLevSetHess[d1, d2].AccLaidBack(1.0, LevSetHess[d1, d2]);
                    PatchRecMultipassFilter(FiltLevSetHess[d1, d2], config.NoOfPatchRecoverySweeps, l2pr);
                }
            }
        }


        /// <summary>
        /// Computes curvature based on the Laplace-Beltrami Ansatz (some would also call this Stokes theorem).
        /// </summary>
        static void LaplaceBeltramiFiltered(LevelSetTracker LsTrk, 
            SinglePhaseField Curvature,
            int NoOfCurvaturePCsweeps, 
            SinglePhaseField[] LevSetGradient,
            int HMForder) {
            GridData GridDat = LsTrk.GridDat;
            int D = GridDat.SpatialDimension;
            Curvature.Clear();
            var CC = LsTrk.Regions.GetCutCellMask();

            var CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D));
            var Params = (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D);
            var DomName = VariableNames.VelocityVector(D);

            var Bs = new Basis(GridDat, 0);
            var Bss = new Basis[D];
            Bss.SetAll(Bs);
            var map = new UnsetteledCoordinateMapping(Bss);

                       

            // \lineint \vec{s} \dl
            // ============================
            var IntgrlOf_CurvTimesNormal = new VectorField<SinglePhaseField>(new SinglePhaseField(Bs), new SinglePhaseField(Bs));

            {
                XSpatialOperatorMk2 op = new XSpatialOperatorMk2(DomName, Params, CodName, (int[] A, int[] B, int[] C) => HMForder, null);
                for(int d = 0; d < D; d++) {
                    var H = new SurfaceTension_LaplaceBeltrami_BndLine(d, 1.0, true);
                    op.SurfaceElementOperator.EquationComponents[CodName[d]].Add(H);
                }
                op.Commit();
                                
                //op.ComputeMatrixEx(LsTrk,
                //    map, LevSetGradient, map,
                //    default(MsrMatrix), IntgrlOf_CurvTimesNormal.CoordinateVector, true, 0.0,
                //    false, LsTrk.GetSpeciesId("A"));
                XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = op.GetMatrixBuilder(LsTrk, map, LevSetGradient, map, LsTrk.GetSpeciesId("A"));
                mtxBuilder.time = 0.0;
                mtxBuilder.ComputeAffine(IntgrlOf_CurvTimesNormal.CoordinateVector);
            }

            // \oint \vec{n} \dS
            // ==================================

            var IntgrlOf_Normal = new VectorField<SinglePhaseField>(new SinglePhaseField(Bs), new SinglePhaseField(Bs));
            {
                XSpatialOperatorMk2 qr = new XSpatialOperatorMk2(DomName, new string[] { "Curvature" }, CodName, (int[] A, int[] B, int[] C) => HMForder, null);
                for(int d = 0; d < D; d++) {
                    qr.EquationComponents[CodName[d]].Add(new CurvatureBasedSurfaceTension(d, D, LsTrk, 1.0));
                }
                qr.Commit();
                
                SinglePhaseField fakeCurvature = new SinglePhaseField(Bs);
                fakeCurvature.AccConstant(1.0);
                
                //qr.ComputeMatrixEx(LsTrk,
                //    map, new DGField[] { fakeCurvature }, map,
                //    default(MsrMatrix), IntgrlOf_Normal.CoordinateVector, true, 0.0,
                //    false, LsTrk.GetSpeciesId("A"));
                XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = qr.GetMatrixBuilder(LsTrk, map, LevSetGradient, map, LsTrk.GetSpeciesId("A"));
                mtxBuilder.time = 0.0;
                mtxBuilder.ComputeAffine(IntgrlOf_CurvTimesNormal.CoordinateVector);
            }

            // compute mean curvature in each cell
            // ===================================

            double kappaMin = double.MaxValue, kappaMax = double.MinValue; // maximum and minimum curvature over all cells
            foreach(var jCell in CC.ItemEnum) {
                double Int_kappaN1 = IntgrlOf_CurvTimesNormal[0].GetMeanValue(jCell);
                double Int_kappaN2 = IntgrlOf_CurvTimesNormal[1].GetMeanValue(jCell);

                double Int_N1 = IntgrlOf_Normal[0].GetMeanValue(jCell);
                double Int_N2 = IntgrlOf_Normal[1].GetMeanValue(jCell);

                double kappa;
                if(Math.Abs(Int_N1) > 1000 * Math.Abs(Int_N2)) {
                    kappa = Int_kappaN1 / Int_N1;
                } else if(Math.Abs(Int_N2) > 1000 * Math.Abs(Int_N1)) {
                    kappa = Int_kappaN2 / Int_N2;
                } else {
                    kappa = 0.5 * (Int_kappaN1 / Int_N1 + Int_kappaN2 / Int_N2);
                }


                Curvature.SetMeanValue(jCell, kappa);

                kappaMin = Math.Min(kappaMin, kappa);
                kappaMax = Math.Max(kappaMax, kappa);
            }


            // finally, do patch-recovery filtering
            // ====================================

            for(int deg = 1; deg <= 3; deg++) {
                if(NoOfCurvaturePCsweeps > 0) {
                    SinglePhaseField CurvatureRed = new SinglePhaseField(new Basis(Curvature.GridDat, deg), "X");
                    CurvatureRed.AccLaidBack(1.0, Curvature);

                    var l2pr = new L2PatchRecovery(CurvatureRed.Basis, CurvatureRed.Basis, CC, true);
                    PatchRecMultipassFilter(CurvatureRed, NoOfCurvaturePCsweeps, l2pr);

                    Curvature.Clear();
                    Curvature.AccLaidBack(1.0, CurvatureRed);
                }
            }
        }

        
        private static void PatchRecMultipassFilter(SinglePhaseField F, int NoOfSweeps, L2PatchRecovery l2pr, int[] lim = null) {

            SinglePhaseField F_org = F.CloneAs();


            if (lim != null) {
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                // This is the CUrvature limiter  --  one of my GREATEST hacks
                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                if (NoOfSweeps == 0)
                    throw new NotSupportedException("curvature limiter cannot be used when curvature filtering is turned off");


                for (int pass = 0; pass < 4; pass++) {
                    //lim.SetAll(F_org.Basis.Polynomials.GetRow(0).Where(p => p.Degree <= deg).Count());

                    int[] thisLim = lim.CloneAs();

                    if (pass > 0) {
                        int J = F.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
                        for (int j = 0; j < J; j++) {
                            if (thisLim[j] >= 0)
                                thisLim[j] += pass;
                        }

                    }

                    l2pr.notchangeunlim = true;
                    //for(int i = 0; i <= pass && l2pr.notchangeunlim; i++)
                    //    LimMod(F.GridDat, thisLim, F.Basis.Degree);


                    for (int subpass = 0; subpass < 2; subpass++) {
                        F_org.Clear();
                        F_org.Acc(1.0, F);
                        F.Clear();
                        l2pr.SetLimiter(thisLim);
                        l2pr.Perform(F, F_org);
                    }


                }//*/

                l2pr.SetLimiter(null);
                l2pr.notchangeunlim = false;
            }


            {
                for (int pass = 0; pass < NoOfSweeps; pass++) {
                    F_org.Clear();
                    F_org.Acc(1.0, F);
                    F.Clear();
                    l2pr.Perform(F, F_org);
                }
            }
        }

        
        /// <summary>
        /// Project Curvature From Hessian and Gradient
        /// </summary>
        /// <remarks>
        /// see for clarification equations 3.6 (2D) and 4.2 (3D)
        /// Goldman, Ron. “Curvature Formulas for Implicit Curves and Surfaces.”
        /// Computer Aided Geometric Design, Geometric Modelling and Differential Geometry, 22, no. 7 (October 2005): 632–58.
        /// doi:10.1016/j.cagd.2005.06.005
        /// 
        /// The formulas only differ by a factor of 1/2.
        /// </remarks>
        /// <param name="_GradPhi">Gradient Vector of the Level-Set</param>
        /// <param name="_HessPhi">Hessian Matrix of the Level-Set</param>
        /// <param name="Curvature">The Resulting Curvature Field</param>
        /// <param name="CC">Cell Mask, on which the Curvature is calculated</param>
        /// <param name="lim">Cells in which the Curvature was Cut by the limiter</param>
        static void ProjectTotalCurvature(
            SinglePhaseField[] _GradPhi,
            SinglePhaseField[,] _HessPhi,
            SinglePhaseField Curvature,
            CellMask CC,
            out int[] lim)//
        {
            using (new FuncTrace()) {

                double scale = 1.0;
                GridData GridDat = (GridData)(_GradPhi[0].GridDat);
                int D = GridDat.SpatialDimension;
                int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                double[] Rmin = new double[J], Rmax = new double[J], Kmin = new double[J], Kmax = new double[J];


                Rmin.SetAll(double.MaxValue);
                Rmax.SetAll(double.MinValue);
                Kmin.SetAll(double.MaxValue);
                Kmax.SetAll(double.MinValue);


                // compute and project 
                // ===================

                // buffers:
                MultidimensionalArray Phi = new MultidimensionalArray(2);
                MultidimensionalArray GradPhi = new MultidimensionalArray(3);
                MultidimensionalArray HessPhi = new MultidimensionalArray(4);

                MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
                MultidimensionalArray Laplace = new MultidimensionalArray(2);
                MultidimensionalArray Q = new MultidimensionalArray(3);

                // evaluate/project:
                //double Erracc = 0;
                Curvature.ProjectField(scale,
                    (ScalarFunctionEx)delegate(int j0, int Len, NodeSet NS, MultidimensionalArray result) { // ScalarFunction2
                    Debug.Assert(result.Dimension == 2);
                    Debug.Assert(Len == result.GetLength(0));
                    int K = result.GetLength(1); // number of nodes


                    // alloc buffers
                    // -------------

                    if (Phi.GetLength(0) != Len || Phi.GetLength(1) != K) {
                        Phi.Allocate(Len, K);
                        GradPhi.Allocate(Len, K, D);
                        HessPhi.Allocate(Len, K, D, D);
                        ooNormGrad.Allocate(Len, K);
                        Laplace.Allocate(Len, K);
                        Q.Allocate(Len, K, D);
                    } else {
                        Phi.Clear();
                        GradPhi.Clear();
                        HessPhi.Clear();
                        ooNormGrad.Clear();
                        Laplace.Clear();
                        Q.Clear();
                    }

                    // evaluate Gradient and Hessian
                    // -----------------------------

                    {
                        for(int d = 0; d < D; d++)
                            _GradPhi[d].Evaluate(j0, Len, NS, GradPhi.ExtractSubArrayShallow(-1, -1, d));
                    }
                    {
                        for(int d1 = 0; d1 < D; d1++)
                            for(int d2 = 0; d2 < D; d2++)
                                _HessPhi[d1, d2].Evaluate(j0, Len, NS, HessPhi.ExtractSubArrayShallow(-1, -1, d1, d2));
                    }

                    // compute the monstrous formula
                    // -----------------------------

                    // norm of Gradient:
                    for (int d = 0; d < D; d++) {
                        var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                        ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

                    // laplacian of phi:
                    for (int d = 0; d < D; d++) {
                        var HessPhi_d_d = HessPhi.ExtractSubArrayShallow(-1, -1, d, d);
                        Laplace.Acc(1.0, HessPhi_d_d);
                    }

                    // result = Laplacian(phi)/|Grad phi|
                    result.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");


                    // result = Grad(1/|Grad(phi)|)
                    for (int d1 = 0; d1 < D; d1++) {
                        var Qd = Q.ExtractSubArrayShallow(-1, -1, d1);

                        for (int d2 = 0; d2 < D; d2++) {
                            var Grad_d2 = GradPhi.ExtractSubArrayShallow(-1, -1, d2);
                            var Hess_d2_d1 = HessPhi.ExtractSubArrayShallow(-1, -1, d2, d1);

                            Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "ik", "ik", "ik");
                        }
                    }

                    ooNormGrad.ApplyAll(x => x * x * x);

                    result.Multiply(1.0, GradPhi, Q, ooNormGrad, 1.0, "ik", "ikd", "ikd", "ik");

                    //In 3D formula contains Factor of 1/2
                    if (D == 3)
                        result.Scale(0.5);

                    for (int i = 0; i < Len; i++) {
                        int jCell = i + j0;
                        for (int k = 0; k < K; k++) {
                            double kappa = result[i, k];
                            if (double.IsNaN(kappa))
                                throw new ArithmeticException();

                            double Radius;
                            if (Math.Abs(kappa) >= 1.0e64 || double.IsInfinity(kappa)) {
                                Radius = 0;
                            } else {
                                Radius = 1.0 / kappa;
                            }

                            Rmin[jCell] = Math.Min(Radius, Rmin[jCell]);
                            Rmax[jCell] = Math.Max(Radius, Rmax[jCell]);

                            Kmin[jCell] = Math.Min(kappa, Kmin[jCell]);
                            Kmax[jCell] = Math.Max(kappa, Kmax[jCell]);
                        }
                    }



                },
                (new CellQuadratureScheme(domain: CC)).SaveCompile(GridDat, Curvature.Basis.Degree * 3)
                    );


                //var _Rmax = new double[CC.NoOfItemsLocally];
                //var _Rmin = new double[CC.NoOfItemsLocally];
                //Curvature.GetCellwiseExtremalValues(_Rmin, _Rmax, CC, delegate(double kappa) {
                //    double RadiusAbs;
                //    if(Math.Abs(kappa) >= 1.0e64 || double.IsInfinity(kappa)) {
                //        RadiusAbs = 0;
                //    } else {
                //        RadiusAbs = Math.Abs(1.0 / kappa);
                //    }
                //    return RadiusAbs;
                //});

                //int jSub = 0;
                //foreach(int jCell in CC.ItemEnum) {
                //    Rmax[jCell] = _Rmax[jSub];
                //    Rmin[jCell] = _Rmin[jSub];
                //    jSub++;
                //}


                //Basis b = new Basis(Curvature.GridDat, 0);
                //var RadiusMin = new SinglePhaseField(b, "Rmin");
                //var RadiusMax = new SinglePhaseField(b, "Rmax");
                //var CurvMin = new SinglePhaseField(b, "Kmin");
                //var CurvMax = new SinglePhaseField(b, "Kmax");
                //var RadiusRatio = new SinglePhaseField(b, "RminOverRmax");
                //var CurvRatio = new SinglePhaseField(b, "Kdiff");
                var h_min = GridDat.Cells.h_min;
                lim = new int[J];
                lim.SetAll(-1);
                //int N1 = Curvature.Basis.Polynomials.GetRow(0).Where(p => p.Degree <= 1).Count();
                //int N2 = Curvature.Basis.Polynomials.GetRow(0).Where(p => p.Degree <= 2).Count();

                bool limused = false;
                foreach (int j in CC.ItemEnum) {
                    //RadiusMin.SetMeanValue(j, Rmin[j]);
                    //RadiusMax.SetMeanValue(j, Rmax[j]);
                    //RadiusRatio.SetMeanValue(j, Math.Abs(Rmax[j] / Rmin[j]));
                    double kappa_crit = 1.0 / (h_min[j] * 0.5);
                    //CurvMin.SetMeanValue(j, Kmin[j]/kappa_crit);
                    //CurvMax.SetMeanValue(j, Math.Abs(Kmax[j]/kappa_crit));
                    double KdiffNorm = Math.Abs((Kmax[j] - Kmin[j])) / kappa_crit;
                    //CurvRatio.SetMeanValue(j, KdiffNorm);
                    if (KdiffNorm >= 1.0) {
                        lim[j] = 2;
                        limused = true;
                    }
                    if (KdiffNorm >= 1.8) {
                        lim[j] = 1;
                        limused = true;
                    }

                }

                //Tecplot.PlotFields(new DGField[] { RadiusMin, RadiusMax, RadiusRatio, CurvMin, CurvMax, CurvRatio, Curvature }, b.GridDat, "radius", "radius", 0, 3);
                if (!limused) {
                    lim = null;
                }

            }
        }


        /// <summary>
        /// Project Curvature From Hessian and Gradient
        /// </summary>
        /// <remarks>
        /// see for clarification equations 3.6 (2D) and 4.2 (3D)
        /// Goldman, Ron. “Curvature Formulas for Implicit Curves and Surfaces.”
        /// Computer Aided Geometric Design, Geometric Modelling and Differential Geometry, 22, no. 7 (October 2005): 632–58.
        /// doi:10.1016/j.cagd.2005.06.005
        /// 
        /// The formulas only differ by a factor of 1/2.
        /// </remarks>
        /// <param name="_GradPhi">Gradient Vector of the Level-Set</param>
        /// <param name="_HessPhi">Hessian Matrix of the Level-Set</param>
        /// <param name="Curvature">The Resulting Curvature Field</param>
        /// <param name="lim">Cells in which the Curvature was Cut by the limiter</param>
        static void ClosestPointTotalCurvature(
            LevelSetTracker lsTrk,
            SinglePhaseField[] _GradPhi,
            SinglePhaseField[,] _HessPhi,
            SinglePhaseField Curvature,
            out int[] lim)//
        {
            using(new FuncTrace()) {

                double scale = 1.0;
                GridData GridDat = (GridData)(_GradPhi[0].GridDat);
                int D = GridDat.SpatialDimension;
                int J = GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                double[] Rmin = new double[J], Rmax = new double[J], Kmin = new double[J], Kmax = new double[J];

                Rmin.SetAll(double.MaxValue);
                Rmax.SetAll(double.MinValue);
                Kmin.SetAll(double.MaxValue);
                Kmax.SetAll(double.MinValue);

                // find closest points 
                // ===================

                int deg = Curvature.Basis.Degree;
                QuadRule[] quadRulez = GridDat.iGeomCells.RefElements.Select(Kref => Kref.GetQuadratureRule(deg*2)).ToArray();
                NodeSet[] nodes = quadRulez.Select(qr => qr.Nodes).ToArray();
                var sgrd = lsTrk.Regions.GetCutCellSubgrid4LevSet(0);
                ClosestPointFinder cpf = new ClosestPointFinder(lsTrk, 0, sgrd, nodes);

                CellQuadratureScheme quadScheme = new CellQuadratureScheme(UseDefaultFactories:false, domain: sgrd.VolumeMask);
                quadScheme.AddFixedRuleS(quadRulez);

                var CompQuadRule = quadScheme.SaveCompile(GridDat, deg*2);
                
                // compute and project 
                // ===================

                MultidimensionalArray[] GradPhi = new MultidimensionalArray[D];
                MultidimensionalArray[,] HessPhi = new MultidimensionalArray[D,D];
                for(int d = 0; d < D; d++) {
                    GradPhi[d] = cpf.EvaluateAtCp(_GradPhi[d]);
                }
                for(int d = 0; d < D; d++) {
                    for(int e = 0; e < D; e++) {
                        HessPhi[d, e] = cpf.EvaluateAtCp(_HessPhi[d, e]);
                    }
                }

                MultidimensionalArray ooNormGrad = new MultidimensionalArray(1);
                MultidimensionalArray Laplace = new MultidimensionalArray(1);
                MultidimensionalArray Q = new MultidimensionalArray(2);

                int[] jCell2jSub = sgrd.LocalCellIndex2SubgridIndex;

                // evaluate/project:
                //double Erracc = 0;
                Curvature.ProjectField(scale,
                    (ScalarFunctionEx)delegate(int j0, int Len, NodeSet NS, MultidimensionalArray _result) { // ScalarFunction2
                    Debug.Assert(_result.Dimension == 2);
                    Debug.Assert(Len == _result.GetLength(0));
                    int K = _result.GetLength(1); // number of nodes

                    for(int jCell = j0; jCell < (j0 + Len); jCell++) {
                        int jSub = jCell2jSub[jCell];
                        MultidimensionalArray result = _result.ExtractSubArrayShallow(jCell - j0, -1);

                        // alloc buffers
                        // -------------

                        if(ooNormGrad.GetLength(0) != NS.NoOfNodes) {
                            ooNormGrad.Allocate(K);
                            Laplace.Allocate(K);
                            Q.Allocate(K, D);
                        } else {
                            ooNormGrad.Clear();
                            Laplace.Clear();
                            Q.Clear();
                        }

                        

                        // compute the monstrous formula
                        // -----------------------------

                        // norm of Gradient:
                        for(int d = 0; d < D; d++) {
                            var GradPhi_d = GradPhi[d].ExtractSubArrayShallow(jSub, -1);
                            ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "k", "k", "k");
                        }
                        ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

                        // laplacian of phi:
                        for(int d = 0; d < D; d++) {
                            var HessPhi_d_d = HessPhi[d,d].ExtractSubArrayShallow(jSub, -1);
                            Laplace.Acc(1.0, HessPhi_d_d);
                        }

                        // result = Laplacian(phi)/|Grad phi|
                        result.Multiply(1.0, Laplace, ooNormGrad, 0.0, "k", "k", "k");


                        // result = Grad(1/|Grad(phi)|)
                        for(int d1 = 0; d1 < D; d1++) {
                            var Qd = Q.ExtractSubArrayShallow(-1, d1);

                            for(int d2 = 0; d2 < D; d2++) {
                                var Grad_d2 = GradPhi[d2].ExtractSubArrayShallow(jSub, -1);
                                var Hess_d2_d1 = HessPhi[d2,d1].ExtractSubArrayShallow(jSub, -1);

                                Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "k", "k", "k");
                            }
                        }

                        ooNormGrad.ApplyAll(x => x * x * x);
                        for(int d = 0; d < D; d++) {
                            result.Multiply(1.0, GradPhi[d].ExtractSubArrayShallow(jSub, -1), Q.ExtractSubArrayShallow(-1, d), ooNormGrad, 1.0, "k", "k", "k", "k");
                        }

                        //In 3D formula contains Factor of 1/2
                        if(D == 3)
                            result.Scale(0.5);

                        {
                            for(int k = 0; k < K; k++) {
                                double kappa = result[k];
                                if(double.IsNaN(kappa))
                                    throw new ArithmeticException();

                                double Radius;
                                if(Math.Abs(kappa) >= 1.0e64 || double.IsInfinity(kappa)) {
                                    Radius = 0;
                                } else {
                                    Radius = 1.0 / kappa;
                                }

                                Rmin[jCell] = Math.Min(Radius, Rmin[jCell]);
                                Rmax[jCell] = Math.Max(Radius, Rmax[jCell]);

                                Kmin[jCell] = Math.Min(kappa, Kmin[jCell]);
                                Kmax[jCell] = Math.Max(kappa, Kmax[jCell]);
                            }
                        }


                    }
                }, quadScheme);


                
                var h_min = GridDat.Cells.h_min;
                lim = new int[J];
                lim.SetAll(-1);
                //int N1 = Curvature.Basis.Polynomials.GetRow(0).Where(p => p.Degree <= 1).Count();
                //int N2 = Curvature.Basis.Polynomials.GetRow(0).Where(p => p.Degree <= 2).Count();

                bool limused = false;
                foreach(int j in sgrd.VolumeMask.ItemEnum) {
                    //RadiusMin.SetMeanValue(j, Rmin[j]);
                    //RadiusMax.SetMeanValue(j, Rmax[j]);
                    //RadiusRatio.SetMeanValue(j, Math.Abs(Rmax[j] / Rmin[j]));
                    double kappa_crit = 1.0 / (h_min[j] * 0.5);
                    //CurvMin.SetMeanValue(j, Kmin[j]/kappa_crit);
                    //CurvMax.SetMeanValue(j, Math.Abs(Kmax[j]/kappa_crit));
                    double KdiffNorm = Math.Abs((Kmax[j] - Kmin[j])) / kappa_crit;
                    //CurvRatio.SetMeanValue(j, KdiffNorm);
                    if(KdiffNorm >= 1.0) {
                        lim[j] = 2;
                        limused = true;
                    }
                    if(KdiffNorm >= 1.8) {
                        lim[j] = 1;
                        limused = true;
                    }

                }

                //Tecplot.PlotFields(new DGField[] { RadiusMin, RadiusMax, RadiusRatio, CurvMin, CurvMax, CurvRatio, Curvature }, b.GridDat, "radius", "radius", 0, 3);
                if(!limused) {
                    lim = null;
                }

            }
        }



        public static void ProjectSurfaceForce(SurfaceStressTensor_IsotropicMode surfaceTensionMode, double sigma, double cellAggTrsh,
            LevelSetTracker LsTrk, MultiphaseCellAgglomerator Agglom, MassMatrixFactory massMtxFac,
            VectorField<SinglePhaseField> levSetGradient, SinglePhaseField Curvature,
            VectorField<SinglePhaseField> surfaceForce ) //
        {
            
            //var surfaceTensionMode = config.dntParams.surfTensionMode;

            GridData GridDat = LsTrk.GridDat;
            int D = GridDat.SpatialDimension;

            var codName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D));
            var parameters = (new string[] {"Curvature", "NX", "NY", "NZ" }).GetSubVector(0, D+1);
            var domName = VariableNames.VelocityVector(D);

            surfaceForce.Clear();


            var ParamsList = new List<DGField>();
            var normals = new SinglePhaseField[D];
            


            //dummy Operator Matrix
            var OpMatrix = new MsrMatrix(surfaceForce.Mapping, surfaceForce.Mapping);

            int HMForder = Agglom.CutCellQuadratureOrder;
            XSpatialOperatorMk2 xOp = new XSpatialOperatorMk2(domName, parameters, codName, (A,B,C) => HMForder, null);

            if(surfaceTensionMode == SurfaceStressTensor_IsotropicMode.Curvature_Projected 
                || surfaceTensionMode == SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean 
                || surfaceTensionMode == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint) {

                
                for (int d = 0; d < D; d++) {
                    xOp.EquationComponents[codName[d]].Add(new CurvatureBasedSurfaceTension(d, D, LsTrk, sigma));
                }
                
            } else if(surfaceTensionMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local 
                || surfaceTensionMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux) {

                for (int d = 0; d < D; d++) {
                        normals[d] = levSetGradient[d];
                }
                
                for (int d = 0; d < D; d++) {
                    var line = new SurfaceTension_LaplaceBeltrami_Surface(d, sigma * 0.5);
                    var surface = new SurfaceTension_LaplaceBeltrami_BndLine(d, sigma * 0.5, surfaceTensionMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);

                    xOp.SurfaceElementOperator.EquationComponents[codName[d]].Add(line);
                    xOp.SurfaceElementOperator.EquationComponents[codName[d]].Add(surface);
                }
                
            } else {
                throw new NotImplementedException(surfaceTensionMode.ToString());
            }

            ParamsList.Add(Curvature);
            ParamsList.AddRange(normals);

            xOp.Commit();


            //xOp.ComputeMatrixEx(LsTrk,
            //    surfaceForce.Mapping, ParamsList, surfaceForce.Mapping, 
            //    OpMatrix, surfaceForce.CoordinateVector, true, 0.0, true,
            //    LsTrk.SpeciesIdS.ToArray());
            XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = xOp.GetMatrixBuilder(LsTrk, surfaceForce.Mapping, ParamsList, surfaceForce.Mapping, LsTrk.SpeciesIdS.ToArray());
            mtxBuilder.time = 0.0;
            mtxBuilder.MPITtransceive = true;
            mtxBuilder.ComputeMatrix(OpMatrix, surfaceForce.CoordinateVector);


            //    SpeciesDictionary, momentFittingVariant, cellAggTrsh, out Agglom);
            //xOp.Evaluate(1.0,1.0,surfaceForce.Mapping, normals ,surfaceForce.Mapping);

            double[] TotalForce = new double[D];
            int J = GridDat.Cells.NoOfLocalUpdatedCells;
            for(int j = 0; j < J; j++) {
                for(int d = 0; d < D; d++) {
                    TotalForce[d] += surfaceForce[d].GetMeanValue(j);
                }
            }

            Console.WriteLine("Total force: ({0}|{1}).", TotalForce[0], TotalForce[1]);

            Agglom.Extrapolate(surfaceForce.Mapping);

            
            // Project Surface Force from Volume Force to Surface Force -> Mutliply with inverse Surface Mass-Matrix

            var CC = LsTrk.Regions.GetCutCellMask();
            var gDat = LsTrk.GridDat;


            // get quadrature rules
            // ====================
            XQuadSchemeHelper H = LsTrk.GetXDGSpaceMetrics(LsTrk.SpeciesIdS.ToArray(), HMForder, 1).XQuadSchemeHelper;
            CellQuadratureScheme cqs = H.GetLevelSetquadScheme(0, CC);
            ICompositeQuadRule<QuadRule> surfRule = cqs.Compile(gDat, HMForder);
            
            // Compute Mass matrix and RHS for the 'strange' projection
            // ========================================================
            int L = CC.NoOfItemsLocally;
            VectorField<SinglePhaseField> buffer = new VectorField<SinglePhaseField>(D, surfaceForce[0].Basis, SinglePhaseField.Factory);
            buffer.Acc(1.0, surfaceForce);
            var Q = new QuadratureKernels(buffer, L);
            
            Q.BlockCnt = 0;
            CellQuadrature.GetQuadrature(
                new int[] { Q.Nnx + D, Q.Nnx },
                gDat,
                surfRule,
                Q.Evaluate, Q.SaveIntegrationResults_surf).Execute();
            Debug.Assert(Q.BlockCnt == L);
            

            // solve the non-diagonal mass matrix systems
            // ==========================================

            int BlkCnt = 0;
            foreach (int jCell in CC.ItemEnum) {
                var LocalMassMatrix = Q.MassMatrix.ExtractSubArrayShallow(BlkCnt, -1, -1);
                var LocalRHS = Q.RHS.ExtractSubArrayShallow(BlkCnt, -1, -1);

                // Die "Massenmatrix" muss nicht unbedingt invbar sein, daher: Least-Squares solve
                LocalMassMatrix.LeastSquareSolve(LocalRHS);

                for (int d = 0; d < Q.D; d++) {
                    for (int n = 0; n < Q.Nnx; n++) {
                        surfaceForce[d].Coordinates[jCell, n] = LocalRHS[n, d];
                    }
                }

                BlkCnt++;
            }
#if Debug
            buffer.Acc(-1.0, surfaceForce);
            double L2Change = buffer.L2Norm();
            Console.WriteLine("Multiplication with InverseMassMatrix changed SurfaceForce by {0}", L2Change);
#endif
            // */

        }


        class QuadratureKernels {
            public QuadratureKernels(VectorField<SinglePhaseField> surfaceVectorField, int L) {
                b = surfaceVectorField[0].Basis;
                Nnx = b.Length;
                gDat = (GridData)b.GridDat;
                D = gDat.SpatialDimension;
                MassMatrix = MultidimensionalArray.Create(L, Nnx, Nnx);
                RHS = MultidimensionalArray.Create(L, Nnx, D);
                UevalBuf = null;
                Uin = surfaceVectorField;
            }

            Basis b;
            public int Nnx;
            GridData gDat;
            public int D;
            public MultidimensionalArray MassMatrix;
            public MultidimensionalArray RHS;
            public MultidimensionalArray UevalBuf;
            public int BlockCnt = 0;
            VectorField<SinglePhaseField> Uin;

            public void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                // Del_Evaluate
                // ~~~~~~~~~~~~~
                NodeSet NS = QR.Nodes;
                int NoOfNodes = NS.NoOfNodes;
                var BasisVal = b.CellEval(NS, i0, Length);
                EvalResult.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, Nnx - 1, Nnx - 1 })
                    .Multiply(1.0, BasisVal, BasisVal, 0.0, "ikmn", "ikm", "ikn");

                if (UevalBuf == null || UevalBuf.GetLength(1) != NoOfNodes)
                    UevalBuf = MultidimensionalArray.Create(Length, NoOfNodes);
                if (UevalBuf.GetLength(0) < Length) {
                    UevalBuf.Allocate(Length, NoOfNodes);
                }
                MultidimensionalArray _UevalBuf;
                if (UevalBuf.GetLength(0) < Length)
                    _UevalBuf = UevalBuf.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { Length - 1, NoOfNodes - 1, Nnx - 1, Nnx - 1 });
                else
                    _UevalBuf = UevalBuf;

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
        }




        public static void MakeItConservative(LevelSetTracker lsTrk,
            SinglePhaseField curvature,
            double sigma,
            VectorField<SinglePhaseField> surfaceForce,
            VectorField<SinglePhaseField> levSetGradient,
            XQuadFactoryHelper.MomentFittingVariants momentFittingVariant, int HMForder) //
        {
            //XQuadSchemeHelper qsh = new XQuadSchemeHelper(lsTrk, 0.0, null, momentFittingVariant, HMForder);

            //var EdgeQuadScheme = qsh.Get_SurfaceElement_EdgeQuadScheme(lsTrk.GetSpeciesId("A"));

            var gDat = lsTrk.GridDat;
            int D = gDat.SpatialDimension;
            SubGrid CutCellsGrid = lsTrk.Regions.GetCutCellSubgrid4LevSet(0);
            int JSub = CutCellsGrid.LocalNoOfCells;
            int[,] E2C = gDat.Edges.CellIndices;

            var Bs = new Basis(gDat, 0);
            var Bss = new Basis[D];
            Bss.SetAll(Bs);
            var map = new UnsetteledCoordinateMapping(Bss);

            var CodName = ((new string[] { "momX", "momY", "momZ" }).GetSubVector(0, D));
            var Params = (new string[] { "NX", "NY", "NZ" }).GetSubVector(0, D);
            var DomName = VariableNames.VelocityVector(D);

            int[] jSub2jCell = CutCellsGrid.SubgridIndex2LocalCellIndex;
            int[] jCell2jSub = CutCellsGrid.LocalCellIndex2SubgridIndex;

            //var SchemeHelper = new XQuadSchemeHelper(lsTrk, momentFittingVariant, lsTrk.GetSpeciesId("A"));
            var SchemeHelper = lsTrk.GetXDGSpaceMetrics(new[] { lsTrk.GetSpeciesId("A") }, HMForder).XQuadSchemeHelper; // new XQuadSchemeHelper(lsTrk, momentFittingVariant, );
            

            /*
            MultidimensionalArray ArcLen = MultidimensionalArray.Create(JSub);
            {
                var SchemeHelper = new XQuadSchemeHelper(lsTrk, 0.0, null, momentFittingVariant, HMForder);
                CellQuadratureScheme cqs;
                SchemeHelper.GetLevelSetquadScheme(out cqs, 0, CutCellsGrid.VolumeMask);

                CellQuadrature.GetQuadrature(new int[] { 1 }, lsTrk.GridDat,
                    cqs.Compile(lsTrk.GridDat, HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        EvalResult.SetAll(1.0);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            ArcLen[jCell2jSub[i0 + i]] = ResultsOfIntegration[i,0];
                    }
                ).Execute();
            }
            */

            CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, CutCellsGrid.VolumeMask);
            EdgeQuadratureScheme eqs = SchemeHelper.Get_SurfaceElement_EdgeQuadScheme(lsTrk.GetSpeciesId("A"));
            
            // =============================
            // compute force vector in cell.
            // =============================
            MultidimensionalArray Force = MultidimensionalArray.Create(JSub, D);
            {
                // compute, in each cell 
                // \oint \vec{n} \kappa \dS
                                

                CellQuadrature.GetQuadrature(new int[] { D }, lsTrk.GridDat,
                    cqs.Compile(lsTrk.GridDat, HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        MultidimensionalArray Normals = lsTrk.DataHistories[0].Current.GetLevelSetNormals(QR.Nodes, i0, Length);
                        MultidimensionalArray Curv = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                        curvature.Evaluate(i0, Length, QR.Nodes, Curv);
                        EvalResult.Multiply(1.0, Normals, Curv, 0.0, "jkd", "jkd", "jk");
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            for (int d = 0; d < D; d++) {
                                Force[jCell2jSub[i0 + i], d] = ResultsOfIntegration[i, d];
                            }
                    },
                    cs: CoordinateSystem.Physical
                ).Execute();


            }
                        
            // =========
            // Tangenten
            // =========
            MultidimensionalArray Tangents = MultidimensionalArray.Create(gDat.Edges.Count, D);
            {
                // \lineint \vec{s} \dl
                

                EdgeQuadrature.GetQuadrature(new int[] { D }, lsTrk.GridDat,
                    eqs.Compile(lsTrk.GridDat, HMForder),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        var EdgeNormals = gDat.Edges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, Length);

                        MultidimensionalArray levSetNormalsIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);
                        MultidimensionalArray levSetNormalsOT = MultidimensionalArray.Create(Length, QR.NoOfNodes, D);

                        for (int d = 0; d < D; d++) {
                            levSetGradient[d].EvaluateEdge(i0, Length, QR.Nodes, levSetNormalsIN.ExtractSubArrayShallow(-1,-1, d), levSetNormalsOT.ExtractSubArrayShallow(-1, -1, d));
                        }

                        double[] TangIn = new double[D];
                        double[] TangOt = new double[D];
                        double[] Tang = new double[D];

                        int K = QR.NoOfNodes;
                        for (int i = 0; i < Length; i++) {
                            for (int k = 0; k < K; k++) {
                                int jCellOt = E2C[i + i0, 1];

                                tangente(levSetNormalsIN.ExtractSubArrayShallow(i, k, -1).To1DArray(),
                                    EdgeNormals.ExtractSubArrayShallow(i, k, -1).To1DArray(), TangIn);

                                if (jCellOt >= 0) {
                                    tangente(levSetNormalsOT.ExtractSubArrayShallow(i, k, -1).To1DArray(),
                                        EdgeNormals.ExtractSubArrayShallow(i, k, -1).To1DArray(), TangIn);

                                    Tang.SetV(TangIn);
                                    Tang.AccV(1.0, TangOt);
                                    Tang.Normalize();
                                } else {
                                    Tang.SetV(TangIn);
                                }

                                EvalResult.ExtractSubArrayShallow(i, k, -1).SetVector(Tang);
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int _i = 0; _i < Length; _i++) {
                            for (int d = 0; d < D; d++) {
                                Tangents[i0 + _i, d] = ResultsOfIntegration[_i, d];
                            }
                        }
                    }
                ).Execute();
            }
#if DEBUG
            foreach (int iEdge in CutCellsGrid.InnerEdgesMask.Complement().ItemEnum) {
                Debug.Assert(Tangents.ExtractSubArrayShallow(iEdge, -1).L2Norm() == 0);
            }
#endif

            MultidimensionalArray ForceDef = Force.CloneAs();
            {
                for (int jSub = 0; jSub < JSub; jSub++) {
                    int jCell = jSub2jCell[jSub];

                    foreach (int i in gDat.Cells.Cells2Edges[jCell]) {
                        double sign = Math.Sign(i);
                        int iEdge = Math.Abs(i) - 1;

                        Debug.Assert(((sign > 0) && (gDat.Edges.CellIndices[iEdge, 0] == jCell)) 
                                  || ((sign < 0) && (gDat.Edges.CellIndices[iEdge, 1] == jCell)));

                        for (int d = 0; d < D; d++) {
                            ForceDef[jSub, d] += sign * Tangents[iEdge, d];
                        }
                    }

                }
            }
            Console.WriteLine("Force defect L2: " + ForceDef.L2Norm());


            int[] iEdge2i = new int[gDat.Edges.Count];
            {
                iEdge2i.SetAll(int.MinValue);
                int _i = 0;
                foreach (int iEdge in CutCellsGrid.InnerEdgesMask.ItemEnum) {
                    iEdge2i[iEdge] = _i;
                    _i++;
                }
            }

            MultidimensionalArray n0 = MultidimensionalArray.Create(gDat.Edges.Count, D);
            foreach (int iEdge in CutCellsGrid.InnerEdgesMask.ItemEnum) {
                n0[iEdge, 0] = +Tangents[iEdge, 1];
                n0[iEdge, 1] = -Tangents[iEdge, 0];
            }


            // if(i = iTest) then
            //     # quasi also: loop über die beiden Zellen an der Kante 'iTest'
            //     # iTest -> (j,r)
            
            //    RHW := RHW - 2*vz[j,r]*n0[i,d]*(-F[j,d] + add(s0[jr2i[j,r1],d]*vz[j,r1] ,r1=1..2));

            //    for r1 from 1 to 2 by 1 do
            //        i1 := jr2i[j,r1];
            //        E[i1] := E[i1] + 2*vz[j,r]*n0[i,d]*n0[i1,d]*vz[j,r1];
            //    end do;
            // end if;



            int I = CutCellsGrid.InnerEdgesMask.NoOfItemsLocally;
            MsrMatrix QQ = new MsrMatrix(I, I, 1, 1);
            double[] RR = new double[I];
            foreach (int iEdge in CutCellsGrid.InnerEdgesMask.ItemEnum) { // loop over all edges with a tangent vector assigned to
                int iSub = iEdge2i[iEdge];

                for (int iInOt = 0; iInOt < 2; iInOt++) { // loop over cells bound to that edge
                    int jCell = E2C[iEdge,iInOt];
                    int jSub = jCell2jSub[jCell];
                    double vz_jr = iInOt == 0 ? 1 : -1;

                    for (int d = 0; d < D; d++) {
                        RR[iSub] = RR[iSub] - 2 * vz_jr * n0[iEdge, d] * (Force[jSub, d]); // + add(s0[jr2i[j, r1], d] * vz[j, r1], r1 = 1..2))
                    }
                    
                    foreach (int r1 in gDat.Cells.Cells2Edges[jCell]) {

                        double vz_jr1 = Math.Sign(r1);
                        int i1Edge = Math.Abs(r1) - 1;
                        int i1Sub = iEdge2i[i1Edge];
                        if (i1Sub < 0)
                            continue;
                        
                        for (int d = 0; d < D; d++) {
                            RR[iSub] = RR[iSub] - 2 * vz_jr * n0[iEdge, d] * (Tangents[i1Edge, d] * vz_jr1);

                            QQ[iSub, i1Sub] = QQ[iSub, i1Sub] + 2 * vz_jr * n0[iEdge, d] * n0[i1Edge, d] * vz_jr1;
                        }
                    }
                }


            }
            
            //double check = RR.L2Norm();
            //Console.WriteLine("check norm: " + check);


            //SimpleInterfaceSolvers
            double[] delta = new double[RR.Length];
            using (var solver = new PARDISOSolver()) {
                solver.DefineMatrix(QQ);
                var SolRes = solver.Solve(delta, RR.ToArray());
            }

            Console.WriteLine(" |delta|_2 = " + delta.L2Norm());



            // project curvature-induced force
            // ===============================

            for(int d = 0; d < D; d++) {
                surfaceForce[d].ProjectField(sigma,
                    (ScalarFunctionEx)delegate(int j0, int Len, NodeSet nodes, MultidimensionalArray result) //
                    {
                        MultidimensionalArray normals_d;
                        
                        normals_d = lsTrk.DataHistories[0].Current.GetLevelSetNormals(nodes, j0, Len).ExtractSubArrayShallow(-1, -1, d); // normals from continuous level-Set
                        {
                            MultidimensionalArray normals = MultidimensionalArray.Create(Len, nodes.NoOfNodes, D);

                            for(int d1 = 0; d1 < D; d1++)
                                levSetGradient[d1].Evaluate(j0, Len, nodes, normals.ExtractSubArrayShallow(-1, -1, d1));
                            for(int j = 0; j < Len; j++) {
                                for(int k = 0; k < nodes.NoOfNodes; k++) {
                                    double len = 0;
                                    
                                    for(int d1 = 0; d1 < D; d1++)
                                        len += normals[j, k, d1].Pow2();
                                    
                                    double ooNorm = 1.0 / Math.Sqrt(len);
                                    
                                    for(int d1 = 0; d1 < D; d1++)
                                        normals[j, k, d1] *= ooNorm;
                                }
                            }

                            normals_d = normals.ExtractSubArrayShallow(-1, -1, d);
                        }
                        MultidimensionalArray Curv = MultidimensionalArray.Create(Len, nodes.NoOfNodes);
                        curvature.Evaluate(j0, Len, nodes, Curv);
                        result.Multiply(1.0, normals_d, Curv, 0.0, "jk", "jk", "jk");
                    }, cqs);
            }

            

            // Force computed via correted tangents
            // =====================================

            //delta.Clear();

            MultidimensionalArray ForceByTangents = MultidimensionalArray.Create(Force.Lengths);
            {
                for (int jSub = 0; jSub < JSub; jSub++) {
                    int jCell = jSub2jCell[jSub];

                    foreach (int i in gDat.Cells.Cells2Edges[jCell]) { // loop over cut edges
                        double sign = Math.Sign(i);
                        int iEdge = Math.Abs(i) - 1;
                        int iSub = iEdge2i[iEdge];
                        if(iSub < 0)
                            continue;

                        Debug.Assert(((sign > 0) && (gDat.Edges.CellIndices[iEdge, 0] == jCell)) 
                                  || ((sign < 0) && (gDat.Edges.CellIndices[iEdge, 1] == jCell)));

                        for (int d = 0; d < D; d++) {
                            ForceByTangents[jSub, d] -= sign * (Tangents[iEdge, d] * Math.Cos(delta[iSub]) + n0[iEdge, d] * Math.Sin(delta[iSub]));
                        }
                    }

                }
            }

           
            


            // replace 0-th mode with force from tangent-computation
            // =====================================================

            foreach(int jCell in CutCellsGrid.VolumeMask.ItemEnum) {
                
                
                for(int d = 0; d < D; d++) {
                    MultidimensionalArray BasisValues = surfaceForce[d].Basis.CellEval(gDat.Cells.GetRefElement(jCell).Center, jCell, 1);

                    double a = sigma * BasisValues[0, 0, 0];

                    double Fd1 = surfaceForce[d].Coordinates[jCell, 0];
                    double Fd2 = ForceByTangents[jCell2jSub[jCell], d]*a;

                    surfaceForce[d].Coordinates[jCell, 0] = ForceByTangents[jCell2jSub[jCell], d] * a;

                    //Console.WriteLine("F[{3},{4}] {0}  {1}, {2}", Fd1 , Fd2, Fd1/Fd2, jCell, d);
                }
            }


        }

        static void tangente(double[] SurfN, double[] EdgeN, double[] tan) {
            Debug.Assert(SurfN.Length == EdgeN.Length);
            if (SurfN.Length != 2)
                throw new NotSupportedException();

            tan[0] = -SurfN[1];
            tan[1] = SurfN[0];

            if (GenericBlas.InnerProd(tan, EdgeN) < 0.0)
                tan.ScaleV(-1.0);
        }
    }
}

