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
using System.Threading.Tasks;

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;

using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;
using BoSSS.Solution.RheologyCommon;

namespace BoSSS.Solution.XNSECommon {

    /// <summary>
    /// 
    /// </summary>
    public static partial class XOperatorComponentsFactory {

        //=======================
        // Navier Stokes equation
        //=======================
        #region nse

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        /// <param name="U0meanrequired"></param>
        public static void AddSpeciesNSE_component(XSpatialOperatorMk2 XOp, INSE_Configuration config, int d, int D, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            // check input
            if(XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already committed. Adding of new components is not allowed");

            string CodName = EquationNames.MomentumEquationComponent(d);
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoSpc, LFFSpc, muSpc, betaS;
            switch(spcName) {
                case "A": { rhoSpc = physParams.rho_A; LFFSpc = dntParams.LFFA; muSpc = physParams.mu_A; betaS = physParams.betaS_A; break; }
                case "B": { rhoSpc = physParams.rho_B; LFFSpc = dntParams.LFFB; muSpc = physParams.mu_B; betaS = physParams.betaS_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective operator
            // ===================
            U0meanrequired = false;
            if(physParams.IncludeConvection && config.isTransport) {
                var conv = new Operator.Convection.ConvectionInSpeciesBulk_LLF(D, BcMap, spcName, d, rhoSpc, LFFSpc);
                comps.Add(conv);
                U0meanrequired = true;
            }

            // pressure gradient
            // =================
            if(config.isPressureGradient) {
                var pres = new Operator.Pressure.PressureInSpeciesBulk(d, BcMap, spcName);
                comps.Add(pres);

                //problably necessary for LDG Simulation. Only one species parameter reynoldsA!!!!!!
                //var presStab = new PressureStabilizationInBulk(dntParams.PresPenalty2, physParams.reynolds_A, spcName, spcId);
                //comps.Add(presStab);
            }

            // viscous operator
            // ================
            if(config.isViscous && !(muSpc == 0.0)) {

                double penalty = dntParams.PenaltySafety;
                switch(dntParams.ViscosityMode) {
                    case ViscosityMode.Standard:
                    case ViscosityMode.TransposeTermMissing: {
                        // Bulk operator:
                        var Visc1 = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                            penalty, //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                            1.0,
                            BcMap, spcName, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);
                        comps.Add(Visc1);

                        //if(dntParams.UseGhostPenalties) {
                        //    var Visc1Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                        //        penalty, 0.0,
                        //        BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);
                        //    XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc1Penalty);
                        //}

                        break;
                    }
                    case ViscosityMode.FullySymmetric: {
                        // Bulk operator
                        var Visc1 = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                            penalty, //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                            1.0,
                            BcMap, spcName, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);
                        comps.Add(Visc1);

                        var Visc2 = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUtranspTerm(
                            penalty,  //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                            1.0,
                            BcMap, spcName, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);
                        comps.Add(Visc2);


                        //if(dntParams.UseGhostPenalties) {
                        //    var Visc1Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                        //        penalty, 0.0,
                        //        BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);
                        //    var Visc2Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUtranspTerm(
                        //        penalty, 0.0,
                        //        BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B, _betaS: betaS);

                        //    XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc1Penalty);
                        //    XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc2Penalty);

                        //}

                        break;
                    }
                    case ViscosityMode.Viscoelastic: {

                        //set species arguments
                        double ReSpc;
                        //double betaSpc;
                        switch(spcName) {
                            case "A": {
                                ReSpc = physParams.reynolds_A;
                                //betaSpc = ((PhysicalParametersRheology)physParams).beta_a;
                                break;
                            }
                            case "B": {
                                ReSpc = physParams.reynolds_B;
                                //betaSpc = ((PhysicalParametersRheology)physParams).beta_b;
                                break;
                            }
                            default: throw new ArgumentException("Unknown species.");
                        }

                        // Bulk operator:
                        var Visc1 = new Operator.Viscosity.DimensionlessViscosityInSpeciesBulk_GradUTerm(
                            penalty, // dntParams.UseGhostPenalties ? 0.0 : penalty, 
                            1.0,
                            BcMap, spcName, d, D, physParams.reynolds_A / ((PhysicalParametersRheology)physParams).beta_a, physParams.reynolds_B / ((PhysicalParametersRheology)physParams).beta_b);
                        comps.Add(Visc1);

                        var Visc2 = new Operator.Viscosity.DimensionlessViscosityInSpeciesBulk_GradUtranspTerm(
                            penalty, //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                            1.0,
                            BcMap, spcName, d, D, physParams.reynolds_A / ((PhysicalParametersRheology)physParams).beta_a, physParams.reynolds_B / ((PhysicalParametersRheology)physParams).beta_b);
                        comps.Add(Visc2);

                        var div = new StressDivergenceInBulk(d, BcMap, ReSpc, dntParams.Penalty1, dntParams.Penalty2, spcName);
                        comps.Add(div);

                        break;
                    }


                    default:
                    throw new NotImplementedException();
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        /// <param name="U0meanrequired"></param>
        public static void AddSpeciesNSE(XSpatialOperatorMk2 XOp, INSE_Configuration config, int D, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            U0meanrequired = false;

            for (int d = 0; d < D; d++) {
                AddSpeciesNSE_component(XOp, config, d, D, spcName, spcId, BcMap, LsTrk, out U0meanrequired);
            }

            XOp.FreeMeanValue[VariableNames.Pressure] = !BcMap.DirichletPressureBoundary;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_component(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int d, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.MomentumEquationComponent(d);
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double LFFA = dntParams.LFFA;
            double LFFB = dntParams.LFFB;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;

            //viscoelastic



            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective operator
            // ===================
            if (physParams.IncludeConvection && config.isTransport) {
                var conv = new Operator.Convection.ConvectionAtLevelSet_LLF(d, D, rhoA, rhoB, LFFA, LFFB, physParams.Material, BcMap, config.isMovingMesh);
                comps.Add(conv);
            }

            // pressure gradient
            // =================
            if (config.isPressureGradient) {
                var presLs = new Operator.Pressure.PressureFormAtLevelSet(d, D, _freeSurface: dntParams.freeSurfaceFlow, _pFree: physParams.pFree);
                comps.Add(presLs);
            }

            // viscous operator
            // ================
            if (config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {

                double penalty = dntParams.PenaltySafety;
                switch(dntParams.ViscosityMode) {
                    case ViscosityMode.Standard:
                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(muA, muB, penalty * 1.0, D, d, true));
                    break;

                    case ViscosityMode.TransposeTermMissing:
                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(muA, muB, penalty * 1.0, D, d, false));
                    break;

                    case ViscosityMode.FullySymmetric:
                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(D, muA, muB, penalty, d,
                        _freeSurface: dntParams.freeSurfaceFlow));
                    break;

                    case ViscosityMode.Viscoelastic:
                    //comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, 1 / reynoldsA, 1 / reynoldsB, penalty * 1.0, d, false));
                    double reynoldsA = physParams.reynolds_A;
                    double reynoldsB = physParams.reynolds_B;
                    double betaA = ((PhysicalParametersRheology)physParams).beta_a;
                    double betaB = ((PhysicalParametersRheology)physParams).beta_b;
                    double[] penalty1 = dntParams.Penalty1;
                    double penalty2 = dntParams.Penalty2;
                    comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(D, betaA / reynoldsA, betaB / reynoldsB, penalty, d));
                    comps.Add(new Operator.Viscosity.StressDivergenceAtLevelSet(reynoldsA, reynoldsB, penalty1, penalty2, D, d));
                    break;

                    default:
                    throw new NotImplementedException();
                }
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

            for (int d = 0; d < D; d++) {
                AddInterfaceNSE_component(XOp, config, d, D, BcMap, LsTrk);
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        /// <param name="degU"></param>
        /// <param name="NormalsRequired"></param>
        /// <param name="CurvatureRequired"></param>
        public static void AddSurfaceTensionForce_component(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int d, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, int degU,
            out bool NormalsRequired, out bool CurvatureRequired) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.MomentumEquationComponent(d);
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set arguments
            double sigma = physParams.Sigma;

            // set components
            //var comps = XOp.SurfaceElementOperator.EquationComponents[CodName];

            // surface stress tensor
            // =====================
            NormalsRequired = false;
            CurvatureRequired = false;
            if (config.isPressureGradient && physParams.Sigma != 0.0) {

                // isotropic part of the surface stress tensor
                // ===========================================

                if (dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux
                 || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local
                 || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {

                    if (dntParams.SST_isotropicMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                        IEquationComponent G = new SurfaceTension_LaplaceBeltrami_Surface(d, sigma * 0.5);
                        XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(G);
                        IEquationComponent H = new SurfaceTension_LaplaceBeltrami_BndLine(d, sigma * 0.5, dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
                        XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(H);
                    } else {
                        IEquationComponent isoSurfT = new IsotropicSurfaceTension_LaplaceBeltrami(d, D, sigma * 0.5, BcMap.EdgeTag2Type, BcMap, physParams.theta_e, physParams.betaL);
                        XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(isoSurfT);
                    }
                    NormalsRequired = true;

                } else if (dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Projected
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier) {

                    XOp.EquationComponents[CodName].Add(new CurvatureBasedSurfaceTension(d, D, sigma));

                    CurvatureRequired = true;

                } else {
                    throw new NotImplementedException("Not implemented.");
                }


                // dynamic part of the surface stress tensor
                // =========================================

                if (dntParams.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                    double muI = physParams.mu_I;
                    double lamI = physParams.lambda_I;
                    double lamI_t = (physParams.lambdaI_tilde < 0) ? (lamI - muI) : physParams.lambdaI_tilde;

                    double penalty_base = (degU + 1) * (degU + (D-1)) / (D-1);
                    double penalty = penalty_base * dntParams.PenaltySafety;

                    // surface dilatational viscosity
                    if (dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceDivergence ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        var surfDiv = new BoussinesqScriven_SurfaceVelocityDivergence(d, D, lamI_t * 0.5, penalty, BcMap.EdgeTag2Type, false);
                        XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(surfDiv);

                    }

                    // surface shear viscosity 
                    if (dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.SemiImplicit ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        var surfDeformRate = new BoussinesqScriven_SurfaceDeformationRate_GradU(d, D, muI * 0.5, penalty, false, dntParams.SurfStressTensor == SurfaceSressTensor.SemiImplicit);
                        XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(surfDeformRate);

                        if (dntParams.SurfStressTensor != SurfaceSressTensor.SemiImplicit) {
                            var surfDeformRateT = new BoussinesqScriven_SurfaceDeformationRate_GradUTranspose(d, D, muI * 0.5, penalty, false);
                            XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(surfDeformRateT);
                        }

                    }

                }


                // stabilization
                // =============

                switch(dntParams.STFstabilization) {
                    case DoNotTouchParameters.SurfaceTensionForceStabilization.surfaceDeformationRateLocal: {
                            XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(new SurfaceDeformationRate_LocalStabilization(d, D, false));
                            break;
                        }
                    case DoNotTouchParameters.SurfaceTensionForceStabilization.GradUxGradV: {
                            XOp.EquationComponents[CodName].Add(new LevelSetStabilization(d, D, 0.1));
                            break;
                        }
                    case DoNotTouchParameters.SurfaceTensionForceStabilization.surfaceDivergence: {
                            XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(new DynamicSurfaceTension_LB_SurfaceVelocityDivergence(d, D, 0.1));
                            break;
                        }
                    case DoNotTouchParameters.SurfaceTensionForceStabilization.EdgeDissipation: {
                            XOp.SurfaceElementOperator_Ls0.EquationComponents[CodName].Add(new DynamicSurfaceTension_LB_EdgeDissipation(d, D, sigma, 0.0));
                            break;
                        }
                    case DoNotTouchParameters.SurfaceTensionForceStabilization.None:
                        break;
                    default:
                        throw new NotImplementedException();
                        

                }

            }


            // artificial surface tension force 
            // ================================

            if (config.isPressureGradient && physParams.useArtificialSurfaceForce) {

                XOp.EquationComponents[CodName].Add(new SurfaceTension_ArfForceSrc(d, D));
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        /// <param name="degU"></param>
        /// <param name="NormalsRequired"></param>
        /// <param name="CurvatureRequired"></param>
        public static void AddSurfaceTensionForce(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, int degU,
            out bool NormalsRequired, out bool CurvatureRequired) {

            NormalsRequired = false;
            CurvatureRequired = false;

            for (int d = 0; d < D; d++) {
                AddSurfaceTensionForce_component(XOp, config, d, D, BcMap, LsTrk, degU, out NormalsRequired, out CurvatureRequired);
            }



        }

        #endregion


        //====================
        // Continuity equation
        //====================
        #region conti

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        public static void AddSpeciesContinuityEq(XSpatialOperatorMk2 XOp, INSE_Configuration config, int D,
            string spcName, SpeciesId spcId, IncompressibleMultiphaseBoundaryCondMap BcMap) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.ContinuityEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoSpc;
            switch (spcName) {
                case "A": { rhoSpc = physParams.rho_A; break; }
                case "B": { rhoSpc = physParams.rho_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            for (int d = 0; d < D; d++) {
                var src = new Operator.Continuity.DivergenceInSpeciesBulk_Volume(d, D, spcName, rhoSpc, -1, false);
                comps.Add(src);
                var flx = new Operator.Continuity.DivergenceInSpeciesBulk_Edge(d, BcMap, spcName, rhoSpc, -1, false);
                comps.Add(flx);
                //var stab = new PressureStabilizationInBulk(dntParams.PresPenalty2, physParams.reynolds_A, physParams.reynolds_B, spcName, spcId);
                //comps.Add(stab);
            }

        }

        /// <summary>
        /// 
        /// </summary>
        public static void AddInterfaceContinuityEq(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.ContinuityEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;

            // set components
            var comps = XOp.EquationComponents[CodName];

            var divPen = new Operator.Continuity.DivergenceAtLevelSet(D, rhoA, rhoB, config.isMatInt, -1, false);
            comps.Add(divPen);
            //var stabint = new PressureStabilizationAtLevelSet(LsTrk, dntParams.PresPenalty2, physParams.reynolds_A, physParams.reynolds_B);
            //comps.Add(stabint);
        }

        #endregion

    }

    /// <summary>
    /// base configuration options
    /// </summary>
    public interface ISolver_Configuration {

        /// <summary>
        /// physical parameters
        /// </summary>
        PhysicalParameters getPhysParams { get; }

        /// <summary>
        /// advanced operator configuration
        /// </summary>
        DoNotTouchParameters getDntParams { get; }


        /// <summary>
        /// Controls the domain variables that the operator should contain. <br/>
        /// This controls only the formal operator shape, not the actual components.
        /// </summary>
        bool[] getDomBlocks { get; }

        /// <summary>
        /// Controls the codomain variables that the operator should contain. <br/>
        ///This controls only the formal operator shape, not the actual components.
        /// </summary>
        bool[] getCodBlocks { get; }

    }


    /// <summary>
    /// configuration options for the bulk NSE including continuity equation
    /// </summary>
    public interface INSE_Configuration : ISolver_Configuration {

        /// <summary>
        /// include transport operator
        /// </summary>
        bool isTransport { get; }

        /// <summary>
        /// include viscous operator
        /// </summary>
        bool isViscous { get; set; }

        /// <summary>
        /// include pressure gradient
        /// </summary>
        bool isPressureGradient { get; }

        /// <summary>
        /// include continuity equation
        /// </summary>
        bool isContinuity { get; }

        /// <summary>
        /// Include a gravitation component
        /// </summary>
        bool isGravity { get; }

        /// <summary>
        /// Include a general volume force component
        /// </summary>
        bool isVolForce { get; }

        /// <summary>
        /// slip on the fluid-fluid interface
        /// </summary>
        bool isInterfaceSlip { get; }

    }


    /// <summary>
    /// extended configuration options for interface discretization
    /// </summary>
    public interface IXNSE_Configuration : INSE_Configuration {

        /// <summary>
        /// switch for moving mesh flux discretizations
        /// </summary>
        bool isMovingMesh { get; }

        /// <summary>
        /// true if the interface is a material interface
        /// </summary>
        bool isMatInt { get; }

        ///// <summary>
        ///// true if the interface pressure is prescribed i.e. for evaporation
        ///// </summary>
        //bool isPInterfaceSet { get; }

    }

}
