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

using BoSSS.Solution.XNSECommon.Operator.SurfaceTension;

namespace BoSSS.Solution.XNSECommon.newXSpatialOperator {
    
    /// <summary>
    /// 
    /// </summary>
    public static class XOperatorComponentsFactory {

        //=======================
        // Navier Stokes equation
        //=======================

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddSpeciesNSE_component(XSpatialOperatorMk2 XOp, string CodName, int d, int D, string spcName, SpeciesId spcId, 
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk, out bool U0meanrequired) {

            // check input
            if(XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.GetPhysParams;
            DoNotTouchParameters dntParams = config.GetDntParams;

            // set species arguments
            double rhoSpc, LFFSpc, muSpc;
            switch(spcName) {
                case "A": { rhoSpc = physParams.rho_A; LFFSpc = dntParams.LFFA; muSpc = physParams.mu_A; break; }
                case "B": { rhoSpc = physParams.rho_B; LFFSpc = dntParams.LFFB; muSpc = physParams.mu_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective operator
            // ===================
            U0meanrequired = false;
            if(physParams.IncludeConvection && config.isTransport) {
                var conv = new Operator.Convection.ConvectionInSpeciesBulk_LLF(D, BcMap, spcName, spcId, d, rhoSpc, LFFSpc, LsTrk);
                comps.Add(conv);
                U0meanrequired = true;
            }

            // pressure gradient
            // =================
            if(config.isPressureGradient) {
                var pres = new Operator.Pressure.PressureInSpeciesBulk(d, BcMap, spcName, spcId);
                comps.Add(pres);
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
                                dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);
                            comps.Add(Visc1);

                            if(dntParams.UseGhostPenalties) {
                                var Visc1Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                                    penalty, 0.0,
                                    BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);
                                XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc1Penalty);
                            }

                            break;
                        }
                    case ViscosityMode.FullySymmetric: {
                            // Bulk operator
                            var Visc1 = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                                dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);
                            comps.Add(Visc1);
                            
                            var Visc2 = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUtranspTerm(
                                dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                                BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);
                            comps.Add(Visc2);


                            if(dntParams.UseGhostPenalties) {
                                var Visc1Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUTerm(
                                    penalty, 0.0,
                                    BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);
                                var Visc2Penalty = new Operator.Viscosity.ViscosityInSpeciesBulk_GradUtranspTerm(
                                    penalty, 0.0,
                                    BcMap, spcName, spcId, d, D, physParams.mu_A, physParams.mu_B);

                                XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc1Penalty);
                                XOp.GhostEdgesOperator.EquationComponents[CodName].Add(Visc2Penalty);

                            }

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
        /// <param name="CodName"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddSpeciesNSE(XSpatialOperatorMk2 XOp, string[] CodName, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk, out bool U0meanrequired) {

            U0meanrequired = false;

            int D = CodName.Length;
            for(int d = 0; d < D; d++) {
                AddSpeciesNSE_component(XOp, CodName[d], d, D, spcName, spcId, BcMap, config, LsTrk, out U0meanrequired);
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_component(XSpatialOperatorMk2 XOp, string CodName, int d, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk) {

            // check input
            if(XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.GetPhysParams;
            DoNotTouchParameters dntParams = config.GetDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double LFFA = dntParams.LFFA;
            double LFFB = dntParams.LFFB;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;


            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective operator
            // ===================
            if(physParams.IncludeConvection && config.isTransport) {
                var conv = new Operator.Convection.ConvectionAtLevelSet_LLF(d, D, LsTrk, rhoA, rhoB, LFFA, LFFB, physParams.Material, BcMap, config.isMovingMesh);
                comps.Add(conv);
            }

            // pressure gradient
            // =================
            if(config.isPressureGradient) {
                var presLs = new Operator.Pressure.PressureFormAtLevelSet(d, D, LsTrk);
                comps.Add(presLs);
            }

            // viscous operator
            // ================
            if(config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {

                double penalty = dntParams.PenaltySafety;
                switch(dntParams.ViscosityMode) {
                    case ViscosityMode.Standard:
                        comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, true));
                        break;
                    case ViscosityMode.TransposeTermMissing:
                        comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, false));
                        break;                   
                    case ViscosityMode.FullySymmetric:
                        comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(LsTrk, muA, muB, penalty, d, dntParams.UseWeightedAverages));
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
        /// <param name="CodName"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE(XSpatialOperatorMk2 XOp, string[] CodName,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk) {

            int D = CodName.Length;
            for(int d = 0; d < D; d++) {
                AddInterfaceNSE_component(XOp, CodName[d], d, D, BcMap, config, LsTrk);
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="D"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddSurfaceTensionForce_component(XSpatialOperatorMk2 XOp, string CodName, int d, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk, int degU,
            out bool NormalsRequired, out bool CurvatureRequired) {

            // check input
            if(XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.GetPhysParams;
            DoNotTouchParameters dntParams = config.GetDntParams;

            // set arguments
            double sigma = physParams.Sigma;

            // set components
            //var comps = XOp.SurfaceElementOperator.EquationComponents[CodName];

            // surface stress tensor
            // =====================
            NormalsRequired = false;
            CurvatureRequired = false;
            if(config.isPressureGradient && physParams.Sigma != 0.0) {

                // isotropic part of the surface stress tensor
                // ===========================================

                if(dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux
                 || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local
                 || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {

                    if(dntParams.SST_isotropicMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                        IEquationComponent G = new SurfaceTension_LaplaceBeltrami_Surface(d, sigma * 0.5);
                        XOp.SurfaceElementOperator.EquationComponents[CodName].Add(G);
                        IEquationComponent H = new SurfaceTension_LaplaceBeltrami_BndLine(d, sigma * 0.5, dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
                        XOp.SurfaceElementOperator.EquationComponents[CodName].Add(H);
                    } else {
                        IEquationComponent isoSurfT = new IsotropicSurfaceTension_LaplaceBeltrami(d, D, sigma * 0.5, BcMap.EdgeTag2Type, physParams.theta_e, physParams.betaL);
                        XOp.SurfaceElementOperator.EquationComponents[CodName].Add(isoSurfT);
                    }
                    NormalsRequired = true;

                } else if(dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Projected
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_ClosestPoint
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_LaplaceBeltramiMean
                        || dntParams.SST_isotropicMode == SurfaceStressTensor_IsotropicMode.Curvature_Fourier) {

                    XOp.EquationComponents[CodName].Add(new CurvatureBasedSurfaceTension(d, D, LsTrk, sigma));

                    CurvatureRequired = true;

                } else {
                    throw new NotImplementedException("Not implemented.");
                }


                // dynamic part of the surface stress tensor
                // =========================================

                if(dntParams.SurfStressTensor != SurfaceSressTensor.Isotropic) {

                    double muI = physParams.mu_I;
                    double lamI = physParams.lambda_I;

                    double penalty_base = (degU + 1) * (degU + D) / D;
                    double penalty = penalty_base * dntParams.PenaltySafety;

                    // surface shear viscosity 
                    if(dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceRateOfDeformation ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.SemiImplicit ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        var surfDeformRate = new BoussinesqScriven_SurfaceDeformationRate_GradU(d, muI * 0.5, penalty);
                        XOp.SurfaceElementOperator.EquationComponents[CodName].Add(surfDeformRate);

                        if(dntParams.SurfStressTensor != SurfaceSressTensor.SemiImplicit) {
                            var surfDeformRateT = new BoussinesqScriven_SurfaceDeformationRate_GradUTranspose(d, muI * 0.5, penalty);
                            XOp.SurfaceElementOperator.EquationComponents[CodName].Add(surfDeformRateT);
                        }

                    }
                    // surface dilatational viscosity
                    if(dntParams.SurfStressTensor == SurfaceSressTensor.SurfaceVelocityDivergence ||
                        dntParams.SurfStressTensor == SurfaceSressTensor.FullBoussinesqScriven) {

                        var surfVelocDiv = new BoussinesqScriven_SurfaceVelocityDivergence(d, muI * 0.5, lamI * 0.5, penalty, BcMap.EdgeTag2Type);
                        XOp.SurfaceElementOperator.EquationComponents[CodName].Add(surfVelocDiv);

                    }
                }


                // stabilization
                // =============

                if(dntParams.UseLevelSetStabilization) {

                    XOp.EquationComponents[CodName].Add(new LevelSetStabilization(d, D, LsTrk));
                }

            }


            // artificial surface tension force 
            // ================================

            if(config.isPressureGradient && physParams.useArtificialSurfaceForce) {

                XOp.EquationComponents[CodName].Add(new SurfaceTension_ArfForceSrc(d, D, LsTrk));
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddSurfaceTensionForce(XSpatialOperatorMk2 XOp, string[] CodName,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk, int degU,
            out bool NormalsRequired, out bool CurvatureRequired) {

            NormalsRequired = false;
            CurvatureRequired = false;

            int D = CodName.Length;
            for(int d = 0; d < D; d++) {
                AddSurfaceTensionForce_component(XOp, CodName[d], d, D, BcMap, config, LsTrk, degU, out NormalsRequired, out CurvatureRequired);
            }



        }



        //====================
        // Continuity equation
        //====================

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="D"></param>
        /// <param name="spcName"></param>
        /// <param name="spcId"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddSpeciesContinuityEq(XSpatialOperatorMk2 XOp, string CodName, int D, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk) {

            // check input
            if(XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.GetPhysParams;
            DoNotTouchParameters dntParams = config.GetDntParams;

            // set species arguments
            double rhoSpc;
            switch(spcName) {
                case "A": { rhoSpc = physParams.rho_A; break; }
                case "B": { rhoSpc = physParams.rho_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            for(int d = 0; d < D; d++) {
                var src = new Operator.Continuity.DivergenceInSpeciesBulk_Volume(d, D, spcId, rhoSpc, dntParams.ContiSign, dntParams.RescaleConti);
                comps.Add(src);
                var flx = new Operator.Continuity.DivergenceInSpeciesBulk_Edge(d, BcMap, spcName, spcId, rhoSpc, dntParams.ContiSign, dntParams.RescaleConti);
                comps.Add(flx);
            }

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="CodName"></param>
        /// <param name="_D"></param>
        /// <param name="BcMap"></param>
        /// <param name="config"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceContinuityEq(XSpatialOperatorMk2 XOp, string CodName, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, IXNSE_Configuration config, LevelSetTracker LsTrk) {

            // check input
            if(XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");
            if(!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.GetPhysParams;
            DoNotTouchParameters dntParams = config.GetDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;

            // set components
            var comps = XOp.EquationComponents[CodName];

            var divPen = new Operator.Continuity.DivergenceAtLevelSet(D, LsTrk, rhoA, rhoB, config.isMatInt, dntParams.ContiSign, dntParams.RescaleConti);
            comps.Add(divPen);
        }

        //==============
        // Heat equation
        //==============

        public static void AddSpeciesHeatEq(XSpatialOperatorMk2 XOp) {

        }


        public static void AddInterfaceHeatEq(XSpatialOperatorMk2 XOp) {

        }


        //========================
        // Kinetic energy equation
        //========================

        public static void AddSpeciesKineticEnergyBalance(XSpatialOperatorMk2 XOp) {

        }

        public static void AddInterfaceKineticEnergyBalance(XSpatialOperatorMk2 XOp) {

        }

    }



    public interface ISolver_Configuration {

        PhysicalParameters GetPhysParams { get; }

        /// <summary>
        /// advanced operator configuration
        /// </summary>
        DoNotTouchParameters GetDntParams { get; }

    }


    /// <summary>
    /// 
    /// </summary>
    public interface IXNSE_Configuration : ISolver_Configuration {

        /// <summary>
        /// include transport operator
        /// </summary>
        bool isTransport { get; }

        /// <summary>
        /// include viscous operator
        /// </summary>
        bool isViscous { get; }

        /// <summary>
        /// include pressure gradient
        /// </summary>
        bool isPressureGradient { get; }

        /// <summary>
        /// include continuity equation
        /// </summary>
        bool isContinuity { get; }

        /// <summary>
        /// switch for moving mesh flux discretizations
        /// </summary>
        bool isMovingMesh { get; }

        /// <summary>
        /// true if the interface is a material interface
        /// </summary>
        bool isMatInt { get; }

    }
}
