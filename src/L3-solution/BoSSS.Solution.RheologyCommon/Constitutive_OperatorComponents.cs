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

using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.NSECommon;


namespace BoSSS.Solution.RheologyCommon {

    /// <summary>
    /// 
    /// </summary>
    public static partial class XConstitutiveOperatorComponentsFactory {

        //=======================
        // Navier Stokes equation
        //=======================

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
        public static void AddSpeciesConstitutive_component(XSpatialOperatorMk2 XOp, IRheology_Configuration config, int compRow, int compCol, int D, int stressDegree, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName;
            int d;

            if (compRow == 0 && compCol == 0) {
                CodName = EquationNames.ConstitutiveXX;
                d = 0;
            } else if (compRow == 0 && compCol == 1) {
                CodName = EquationNames.ConstitutiveXY;
                d = 1;
            } else if (compRow == 1 && compCol == 0) {
                CodName = EquationNames.ConstitutiveXY;
                d = 1;
            } else if (compRow == 1 && compCol == 1) {
                CodName = EquationNames.ConstitutiveYY;
                d = 2;
            } else {
                throw new NotSupportedException("Invalid column or row index (3D not supported). Row: " + compRow + " Col:" + compCol);
            }


            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            //double rhoSpc, LFFSpc, muSpc;
            //switch(spcName) {
            //    case "A": { rhoSpc = physParams.rho_A; LFFSpc = dntParams.LFFA; muSpc = physParams.mu_A; break; }
            //    case "B": { rhoSpc = physParams.rho_B; LFFSpc = dntParams.LFFB; muSpc = physParams.mu_B; break; }
            //     default: throw new ArgumentException("Unknown species.");
            // }

            // set components
            var comps = XOp.EquationComponents[CodName];

            // identity part
            // ===================
            var identity = new IdentityInBulk(d, spcName, spcId);
            comps.Add(identity);


            U0meanrequired = false;

            if (config.isOldroydB) {

                // convective operator
                // ===================
                var convective = new ConvectiveInBulk(compRow, compCol, BcMap, physParams.Weissenberg_a, dntParams.alpha, spcName, spcId);
                comps.Add(convective);
                U0meanrequired = true;

                // objective operator
                // ===================
                var objective = new ObjectiveInBulk(d, BcMap, physParams.Weissenberg_a, dntParams.ObjectiveParam, dntParams.StressPenalty, spcName, spcId);
                comps.Add(objective);

                // objective operator with Stress as param
                // ========================================
                var objective_Tparam = new ObjectiveInBulk_Tparam(d, BcMap, physParams.Weissenberg_a, dntParams.ObjectiveParam, spcName, spcId);
                comps.Add(objective_Tparam);
            }


            // viscous operator
            // ==================
            var viscosity = new ViscosityInBulk(d, BcMap, physParams.beta_a, dntParams.Penalty1, spcName, spcId);
            comps.Add(viscosity);

            // artificial diffusion
            // =====================
            if (config.isUseArtificialDiffusion == true) {
                //    if (d == 0)
                //        var diffusion = new DiffusionInBulk(stressDegree, D, ((GridData)GridData).Cells.cj, VariableNames.StressXX);
                //        comps.Add(diffusion);
                //    if (d == 1)
                //         var diffusion = new DiffusionInBulk(stressDegree, D, ((GridData)GridData).Cells.cj, VariableNames.StressXY);
                //        comps.Add(diffusion);
                //    if (d == 2)
                //         var diffusion = new DiffusionInBulk(stressDegree, D, ((GridData)GridData).Cells.cj, VariableNames.StressYY);
                //        comps.Add(diffusion);
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
        public static void AddSpeciesConstitutive(XSpatialOperatorMk2 XOp, IRheology_Configuration config, int D, int stressDegree, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            if (D > 2) {
                throw new NotSupportedException("Viscoelastic solver does not support 3D calculation. Only implemented for 2D cases!");
            }

            U0meanrequired = false;
            int compRow;
            int compCol;

            for (int d = 0; d < 3; d++) {

                if (d == 0) {
                    compRow = 0;
                    compCol = 0;
                    AddSpeciesConstitutive_component(XOp, config, compRow, compCol, D, stressDegree, spcName, spcId, BcMap, LsTrk, out U0meanrequired);
                } else if (d == 1) {
                    compRow = 0;
                    compCol = 1;
                    AddSpeciesConstitutive_component(XOp, config, compRow, compCol, D, stressDegree, spcName, spcId, BcMap, LsTrk, out U0meanrequired);
                } else if (d == 2) {
                    compRow = 1;
                    compCol = 1;
                    AddSpeciesConstitutive_component(XOp, config, compRow, compCol, D, stressDegree, spcName, spcId, BcMap, LsTrk, out U0meanrequired);
                }
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
        //public static void AddInterfaceNSE_component(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int d, int D,
        //    IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

        //    // check input
        //    if (XOp.IsCommited)
        //        throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

        //    string CodName = EquationNames.MomentumEquationComponent(d);
        //    if (!XOp.CodomainVar.Contains(CodName))
        //        throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

        //    PhysicalParameters physParams = config.getPhysParams;
        //    DoNotTouchParameters dntParams = config.getDntParams;

        //    // set species arguments
        //    double rhoA = physParams.rho_A;
        //    double rhoB = physParams.rho_B;
        //    double LFFA = dntParams.LFFA;
        //    double LFFB = dntParams.LFFB;
        //    double muA = physParams.mu_A;
        //    double muB = physParams.mu_B;


        //    // set components
        //    var comps = XOp.EquationComponents[CodName];

        //    // convective operator
        //    // ===================
        //    if(physParams.IncludeConvection && config.isTransport) {
        //        var conv = new Operator.Convection.ConvectionAtLevelSet_LLF(d, D, LsTrk, rhoA, rhoB, LFFA, LFFB, physParams.Material, BcMap, config.isMovingMesh);
        //        comps.Add(conv);
        //    }

        //    // pressure gradient
        //    // =================
        //    if(config.isPressureGradient) {
        //        var presLs = new Operator.Pressure.PressureFormAtLevelSet(d, D, LsTrk);
        //        comps.Add(presLs);
        //    }

        //    // viscous operator
        //    // ================
        //    if(config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {

        //        double penalty = dntParams.PenaltySafety;
        //        switch(dntParams.ViscosityMode) {
        //            case ViscosityMode.Standard:
        //                comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, true));
        //                break;
        //            case ViscosityMode.TransposeTermMissing:
        //                comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_Standard(LsTrk, muA, muB, penalty * 1.0, d, false));
        //                break;                   
        //            case ViscosityMode.FullySymmetric:
        //                comps.Add(new Operator.Viscosity.ViscosityAtLevelSet_FullySymmetric(LsTrk, muA, muB, penalty, d, dntParams.UseWeightedAverages));
        //                break;

        //            default:
        //                throw new NotImplementedException();
        //        }
        //    }

        //}


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="BcMap"></param>
        /// <param name="LsTrk"></param>
        //public static void AddInterfaceNSE(XSpatialOperatorMk2 XOp, IXNSE_Configuration config, int D,
        //    IncompressibleMultiphaseBoundaryCondMap BcMap,  LevelSetTracker LsTrk) {

        //    for(int d = 0; d < D; d++) {
        //        AddInterfaceNSE_component(XOp, config, d, D, BcMap, LsTrk);
        //    }
        //}


    }

    /// <summary>
    /// configuration options for the bulk NSE including continuity equation
    /// </summary>
    public interface IRheology_Configuration : ISolver_Configuration {

        /// <summary>
        /// include upper convected derivative
        /// </summary>
        bool isOldroydB { get; set; }

        /// <summary>
        /// include artificial diffusion term
        /// </summary>
        bool isUseArtificialDiffusion { get; set; }

        /// <summary>
        /// use Jacobian for linearization
        /// </summary>
        bool isUseJacobian { get; set; }

    }


    /// <summary>
    /// extended configuration options for interface discretizations
    /// </summary>
    public interface IXRheology_Configuration : IRheology_Configuration {

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
