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
using BoSSS.Solution.XNSECommon;

namespace BoSSS.Solution.XheatCommon {

    /// <summary>
    /// 
    /// </summary>
    public static partial class XOperatorComponentsFactory {

        //==============
        // Heat equation
        //==============

        public static void AddSpeciesHeatEq(XSpatialOperatorMk2 XOp, IHeat_Configuration config, int D, 
            string spcName, SpeciesId spcId, ThermalMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.HeatEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            if (config.isSeparated) {
                foreach (string cn in EquationNames.AuxHeatFlux(D)) {
                    if (!XOp.CodomainVar.Contains(cn))
                        throw new ArgumentException("CoDomain variable \"" + cn + "\" is not defined in Spatial Operator");
                }
            }

            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double capSpc, LFFSpc, kSpc;
            switch (spcName) {
                case "A": { capSpc = thermParams.rho_A * thermParams.c_A; LFFSpc = dntParams.LFFA; kSpc = thermParams.k_A; break; }
                case "B": { capSpc = thermParams.rho_B * thermParams.c_B; LFFSpc = dntParams.LFFB; kSpc = thermParams.k_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];


            // convective part
            // ================
            if (thermParams.IncludeConvection) {

                var conv = new HeatConvectionInSpeciesBulk(D, BcMap, spcName, spcId, capSpc, LFFSpc, LsTrk);
                comps.Add(conv);

            }


            // viscous operator (laplace)
            // ==========================
            if (!config.isSeparated) {
                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityInSpeciesBulk(
                    dntParams.UseGhostPenalties ? 0.0 : penalty, 1.0,
                    BcMap, D, spcName, spcId, thermParams.k_A, thermParams.k_B);

                comps.Add(Visc);

                if (dntParams.UseGhostPenalties) {
                    var ViscPenalty = new ConductivityInSpeciesBulk(penalty * 1.0, 0.0, BcMap, D,
                        spcName, spcId, thermParams.k_A, thermParams.k_B);
                    XOp.GhostEdgesOperator.EquationComponents[CodName].Add(ViscPenalty);
                }

            } else {

                comps.Add(new HeatFluxDivergenceInSpeciesBulk(D, BcMap, spcName, spcId));
                if (config.withStabilization)
                    comps.Add(new AuxiliaryStabilizationForm(D, BcMap, spcName, spcId));

                for (int d = 0; d < D; d++) {
                    comps = XOp.EquationComponents[EquationNames.AuxHeatFluxComponent(d)];

                    comps.Add(new AuxiliaryHeatFlux_Identity(d, spcId));   // cell local
                    comps.Add(new TemperatureGradientInSpeciesBulk(d, BcMap, spcName, spcId, kSpc));
                    if (config.withStabilization)
                        comps.Add(new TemperatureStabilizationForm(d, BcMap, spcName, spcId));
                }

            }

        }


        public static void AddInterfaceHeatEq(XSpatialOperatorMk2 XOp, IXHeat_Configuration config, int D, 
            ThermalMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommited)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.HeatEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            if (config.isSeparated) {
                foreach (string cn in EquationNames.AuxHeatFlux(D)) {
                    if (!XOp.CodomainVar.Contains(cn))
                        throw new ArgumentException("CoDomain variable \"" + cn + "\" is not defined in Spatial Operator");
                }
            }

            ThermalParameters thermParams = config.getThermParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            // set species arguments
            double capA = thermParams.rho_A * thermParams.c_A;
            double LFFA = dntParams.LFFA;
            double kA = thermParams.k_A;

            double capB = thermParams.rho_B * thermParams.c_B;
            double LFFB = dntParams.LFFB;
            double kB = thermParams.k_B;

            double Tsat = thermParams.T_sat;


            // set components
            var comps = XOp.EquationComponents[CodName];


            // convective part
            // ================
            if (thermParams.IncludeConvection) {

                var conv = new HeatConvectionAtLevelSet(D, LsTrk, capA, capB, LFFA, LFFB, BcMap, config.isMovingMesh, config.isEvaporation, Tsat);
                comps.Add(conv);

            }


            // viscous operator (laplace)
            // ==========================
            if (!config.isSeparated) {

                double penalty = dntParams.PenaltySafety;

                var Visc = new ConductivityAtLevelSet(LsTrk, kA, kB, penalty * 1.0, Tsat);
                comps.Add(Visc);

            } else {

                comps.Add(new HeatFluxDivergencetAtLevelSet(LsTrk, config.isEvaporation)); 
                //if(config.withStabilization)
                //    comps.Add(new AuxiliaryStabilizationFormAtLevelSet(LsTrk, config.isEvaporation));

                for (int d = 0; d < D; d++) {
                    comps = XOp.EquationComponents[EquationNames.AuxHeatFluxComponent(d)];

                    comps.Add(new TemperatureGradientAtLevelSet(d, LsTrk, kA, kB, config.isEvaporation, Tsat));
                    //if (config.withStabilization)
                    //    comps.Add(new TemperatureStabilizationFormAtLevelSet(d, LsTrk, config.isEvaporation, Tsat));
                }

            }

        }


    }


    /// <summary>
    /// configuration options for the bulk NSE including continuity equation
    /// </summary>
    public interface IHeat_Configuration : ISolver_Configuration {

        /// <summary>
        /// thermal parameters
        /// </summary>
        ThermalParameters getThermParams { get; }

        /// <summary>
        /// true if the heat equation is separated via the auxiliary heat flux formulation
        /// </summary>
        bool isSeparated { get; }

        /// <summary>
        /// additional stabilizations form via penalty terms
        /// </summary>
        bool withStabilization { get;  }

        /// <summary>
        /// include transport operator
        /// </summary>
        bool isHeatTransport { get; }

    }


    /// <summary>
    /// extended configuration options for interface discretizations
    /// </summary>
    public interface IXHeat_Configuration : IHeat_Configuration {

        /// <summary>
        /// switch for moving mesh flux discretizations
        /// </summary>
        bool isMovingMesh { get; }

        /// <summary>
        /// true if the interface allows for evaporative mass flux
        /// </summary>
        bool isEvaporation { get; }

    }


}
