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
        // Constitutive equations fo viscoelastic multiphase flow
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
        public static void AddSpeciesConstitutive_component(XDifferentialOperatorMk2 XOp, IRheology_Configuration config, int d, int D, int stressDegree, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName;

            if (d == 0) {
                CodName = EquationNames.ConstitutiveXX;
            } else if (d == 1) {
                CodName = EquationNames.ConstitutiveXY;
            } else if (d == 2) {
                CodName = EquationNames.ConstitutiveYY;
            } else {
                throw new NotSupportedException("Invalid index (3D not supported), d is: " + d);
            }


            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParametersRheology physParams = config.getPhysParams as PhysicalParametersRheology;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double WiSpc, betaSpc, giesekusfactorSpc;

            switch (spcName) {
                case "A": { WiSpc = physParams.Weissenberg_a; betaSpc = physParams.beta_a; giesekusfactorSpc = physParams.giesekusfactor_a; break; }
                case "B": { WiSpc = physParams.Weissenberg_b; betaSpc = physParams.beta_b; giesekusfactorSpc = physParams.giesekusfactor_b; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            // identity part
            // ===================
            var identity = new IdentityInBulk(d, giesekusfactorSpc, WiSpc, betaSpc, spcName, spcId);
            //var identity = new IdentityInBulk(d, spcName, spcId);
            comps.Add(identity);


            U0meanrequired = false;

            if (config.isOldroydB) {

                // convective operator
                // ===================
                var convective = new ConvectiveInBulk(d, BcMap, WiSpc, dntParams.alpha, spcName, spcId);
                comps.Add(convective);
                U0meanrequired = true;

                // objective operator
                // ===================
                var objective = new ObjectiveInBulk(d, BcMap, WiSpc, dntParams.StressPenalty, spcName, spcId);
                comps.Add(objective);

            }

            // viscous operator
            // ==================
            var viscosity = new ViscosityInBulk(d, BcMap, betaSpc, dntParams.Penalty1, spcName, spcId);
            comps.Add(viscosity);

            // artificial diffusion
            // =====================
            if (config.isUseArtificialDiffusion == true) {

                throw new NotImplementedException();
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
        public static void AddSpeciesConstitutive(XDifferentialOperatorMk2 XOp, IRheology_Configuration config, int D, int stressDegree, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            if (D > 2) {
                throw new NotSupportedException("Viscoelastic solver does not support 3D calculation. Only implemented for 2D cases!");
            }

            U0meanrequired = false;

            for (int d = 0; d < 3; d++) {
                    AddSpeciesConstitutive_component(XOp, config, d, D, stressDegree, spcName, spcId, BcMap, LsTrk, out U0meanrequired);
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
        public static void AddInterfaceConstitutive_component(XDifferentialOperatorMk2 XOp, IRheology_Configuration config, int d, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");


            string CodName;

            if (d == 0) {
                CodName = EquationNames.ConstitutiveXX;
            } else if (d == 1) {
                CodName = EquationNames.ConstitutiveXY;
            } else if (d == 2) {
                CodName = EquationNames.ConstitutiveYY;
            } else {
                throw new NotSupportedException("Invalid index (3D not supported), d is: " + d);
            }

            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParametersRheology physParams = config.getPhysParams as PhysicalParametersRheology;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set components
            var comps = XOp.EquationComponents[CodName];

            // identity part: local variable, therefore no contribution at interface!
            // ===================

            U0meanrequired = false;

            if (config.isOldroydB) {

                // convective operator
                // ===================
                var convective = new ConvectiveAtLevelSet(d, physParams.Weissenberg_a, physParams.Weissenberg_b, dntParams.alpha);
                comps.Add(convective);
                U0meanrequired = true;

                // objective operator
                // ===================
                var objective = new ObjectiveAtLevelSet( d, physParams.Weissenberg_a, physParams.Weissenberg_b);
                comps.Add(objective);
            }

            // viscous operator
            // ==================
            var viscosity = new ViscosityAtLevelSet(LsTrk, d, physParams.beta_a, physParams.beta_b, dntParams.Penalty1);
            comps.Add(viscosity);

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name = "XOp" ></param>
        /// <param name="config"></param>
        /// <param name = "BcMap" ></param>
        /// <param name="LsTrk"></param>
        /// <param name="U0meanrequired"></param>
        public static void AddInterfaceConstitutive(XDifferentialOperatorMk2 XOp, IRheology_Configuration config, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, out bool U0meanrequired) {

            U0meanrequired = false;

            for (int d = 0; d < 3; d++) {
                AddInterfaceConstitutive_component(XOp, config, d, D, BcMap, LsTrk, out U0meanrequired);
            }
        }


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
