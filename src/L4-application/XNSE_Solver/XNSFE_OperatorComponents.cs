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

using ilPSP;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XheatCommon;


namespace BoSSS.Application.XNSE_Solver {
    public static partial class XOperatorComponentsFactory {

        //==============================================================
        // additional interface components for NSFE interface components
        //==============================================================


        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="d"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_withEvaporation_component(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config,
            int d, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.MomentumEquationComponent(d);
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;
            double sigma = physParams.Sigma;


            // compute further arguments
            double hVapA = thermParams.hVap_A;
            double hVapB = thermParams.hVap_B;

            double Tsat = thermParams.T_sat;

            double f = config.thermParams.fc;
            double R = config.thermParams.Rc;
            double R_int = 0.0;
            if (config.thermParams.hVap_A > 0 && config.thermParams.hVap_B < 0) {
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoB * hVapA.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoA * hVapA.Pow2())));
            } else if (config.thermParams.hVap_A < 0 && config.thermParams.hVap_B > 0) {
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoA * hVapB.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoB * hVapB.Pow2())));
            }


            // set components
            var comps = XOp.EquationComponents[CodName];


            if (config.isTransport) {
                comps.Add(new ConvectionAtLevelSet_nonMaterialLLF(d, D, LsTrk, rhoA, rhoB, thermParams, R_int, sigma));
                comps.Add(new ConvectionAtLevelSet_Consistency(d, D, LsTrk, rhoA, rhoB, dntParams.ContiSign, dntParams.RescaleConti, thermParams, R_int, sigma));
            }

            if (config.isViscous) {
                comps.Add(new ViscosityAtLevelSet_FullySymmetric_withEvap(LsTrk, muA, muB, dntParams.PenaltySafety, d, rhoA, rhoB, thermParams, R_int, sigma));
            }

            comps.Add(new MassFluxAtInterface(d, D, LsTrk, rhoA, rhoB, thermParams, R_int, sigma));


        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceNSE_withEvaporation(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config, int D, LevelSetTracker LsTrk) {

            for (int d = 0; d < D; d++) {
                AddInterfaceNSE_withEvaporation_component(XOp, config, d, D, LsTrk);
            }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceContinuityEq_withEvaporation(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.ContinuityEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double sigma = physParams.Sigma;


            // compute further arguments
            double hVapA = thermParams.hVap_A;
            double hVapB = thermParams.hVap_B;

            double Tsat = thermParams.T_sat;

            double f = config.thermParams.fc;
            double R = config.thermParams.Rc;
            double R_int = 0.0;
            if (config.thermParams.hVap_A > 0 && config.thermParams.hVap_B < 0) {
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoB * hVapA.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoA * hVapA.Pow2())));
            } else if (config.thermParams.hVap_A < 0 && config.thermParams.hVap_B > 0) {
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoA * hVapB.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoB * hVapB.Pow2())));
            }


            // set components
            var comps = XOp.EquationComponents[CodName];

            var divEvap = new DivergenceAtLevelSet_withEvaporation(D, LsTrk, rhoA, rhoB, dntParams.ContiSign, dntParams.RescaleConti, thermParams, R_int, sigma);
            comps.Add(divEvap);

        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="XOp"></param>
        /// <param name="config"></param>
        /// <param name="D"></param>
        /// <param name="LsTrk"></param>
        public static void AddInterfaceHeatEq_withEvaporation(XSpatialOperatorMk2 XOp, XNSFE_OperatorConfiguration config, int D, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.HeatEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            ThermalParameters thermParams = config.getThermParams;

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double sigma = physParams.Sigma;


            // compute further arguments
            double hVapA = thermParams.hVap_A;
            double hVapB = thermParams.hVap_B;

            double Tsat = thermParams.T_sat;

            double f = config.thermParams.fc;
            double R = config.thermParams.Rc;
            double rho_l = 0.0;
            double R_int = 0.0;
            if (config.thermParams.hVap_A > 0 && config.thermParams.hVap_B < 0) {
                rho_l = rhoA;
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoB * hVapA.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoA * hVapA.Pow2())));
            } else if (config.thermParams.hVap_A < 0 && config.thermParams.hVap_B > 0) {
                rho_l = rhoB;
                R_int = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rhoA * hVapB.Pow2());
                //T_intMin = Tsat * (1 + (pc / (rhoB * hVapB.Pow2())));
            }


            // set components
            var comps = XOp.EquationComponents[CodName];


            // mass flux at interface
            // ======================
            comps.Add(new HeatFluxAtLevelSet(LsTrk, rho_l, thermParams, R_int, sigma));


            // convective part
            // ================
            if (thermParams.IncludeConvection) 
                comps.Add(new HeatConvectionAtLevelSet_Divergence(D, LsTrk, rhoA, rhoB, thermParams, R_int, sigma));


        }

        
    }
}
