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

namespace BoSSS.Solution.EnergyCommon {

    /// <summary>
    /// 
    /// </summary>
    public static partial class XOperatorComponentsFactory {


        //========================
        // Kinetic energy equation
        //========================
        #region energy

        public static void AddSpeciesKineticEnergyEquation(XDifferentialOperatorMk2 XOp, IEnergy_Configuration config, int D, string spcName, SpeciesId spcId,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.KineticEnergyEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            bool laplaceKinE = (config.getKinEviscousDiscretization == KineticEnergyViscousSourceTerms.laplaceKinE);
            bool divergenceP = (config.getKinEpressureDiscretization == KineticEnergyPressureSourceTerms.divergence);

            // set species arguments
            double rhoSpc, LFFSpc, muSpc;
            switch (spcName) {
                case "A": { rhoSpc = physParams.rho_A; LFFSpc = dntParams.LFFA; muSpc = physParams.mu_A; break; }
                case "B": { rhoSpc = physParams.rho_B; LFFSpc = dntParams.LFFB; muSpc = physParams.mu_B; break; }
                default: throw new ArgumentException("Unknown species.");
            }

            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective part
            // ================
            if (config.isTransport) {
                
                var convK = new KineticEnergyConvectionInSpeciesBulk(D, BcMap, spcName, spcId, rhoSpc, LFFSpc, LsTrk);
                //var convK = new KineticEnergyConvectionInSpeciesBulk_Upwind(D, BcMap, spcName, spcId, rhoSpc);
                comps.Add(convK);

            }

            // viscous terms
            // =============
            if (config.isViscous && !(muSpc == 0.0)) {

                double penalty = dntParams.PenaltySafety;

                // Laplace of kinetic energy
                // =========================
                if (laplaceKinE){

                    var kELap = new KineticEnergyLaplaceInSpeciesBulk(
                        penalty, //dntParams.UseGhostPenalties ? 0.0 : penalty, 
                        1.0,
                        BcMap, spcName, spcId, D, physParams.mu_A, physParams.mu_B);
                    comps.Add(kELap);

                    //if (dntParams.UseGhostPenalties) {
                    //    var kELapPenalty = new KineticEnergyLaplaceInSpeciesBulk(
                    //        penalty, 0.0,
                    //        BcMap, spcName, spcId, D, physParams.mu_A, physParams.mu_B);
                    //    XOp.GhostEdgesOperator.EquationComponents[CodName].Add(kELapPenalty);
                    //}
                }

                // Divergence of stress tensor
                // ===========================
                {
                    comps.Add(new StressDivergenceInSpeciesBulk(D, BcMap, spcName, spcId, muSpc, transposed: !laplaceKinE));

                    if (config.getKinEviscousDiscretization == KineticEnergyViscousSourceTerms.local) {
                        throw new ApplicationException("deprecated option");
                        //comps.Add(new StressDivergence_Local(D, muSpc, spcId, transposed: !laplaceKinE));
                    }

                }

                // Dissipation
                // ===========
                {
                    comps.Add(new Dissipation(D, muSpc, spcName, spcId, _withPressure: config.withPressureDissipation));
                }

            }

            // pressure term
            // =============
            if (config.isPressureGradient) {

                //if (divergenceP) {
                //    comps.Add(new DivergencePressureEnergyInSpeciesBulk(D, BcMap, spcName, spcId));
                //    //comps.Add(new ConvectivePressureTerm_LLF(D, BcMap, spcName, spcId, LFFSpc, LsTrk));
                //} else {
                //    comps.Add(new PressureGradientConvection(D, spcId));
                //}

            }

            // gravity (volume forces)
            // =======================
            {
                comps.Add(new PowerofGravity(D, spcName, spcId, rhoSpc));
            }


        }

        public static void AddInterfaceKineticEnergyEquation(XDifferentialOperatorMk2 XOp, IEnergy_Configuration config, int D,
            IncompressibleMultiphaseBoundaryCondMap BcMap, LevelSetTracker LsTrk, int degU) {

            // check input
            if (XOp.IsCommitted)
                throw new InvalidOperationException("Spatial Operator is already comitted. Adding of new components is not allowed");

            string CodName = EquationNames.KineticEnergyEquation;
            if (!XOp.CodomainVar.Contains(CodName))
                throw new ArgumentException("CoDomain variable \"" + CodName + "\" is not defined in Spatial Operator");

            PhysicalParameters physParams = config.getPhysParams;
            DoNotTouchParameters dntParams = config.getDntParams;

            bool laplaceKinE = (config.getKinEviscousDiscretization == KineticEnergyViscousSourceTerms.laplaceKinE);
            bool divergenceP = (config.getKinEpressureDiscretization == KineticEnergyPressureSourceTerms.divergence);

            // set species arguments
            double rhoA = physParams.rho_A;
            double rhoB = physParams.rho_B;
            double LFFA = dntParams.LFFA;
            double LFFB = dntParams.LFFB;
            double muA = physParams.mu_A;
            double muB = physParams.mu_B;
            double sigma = physParams.Sigma;

            double penalty_base = (degU + 1) * (degU + D) / D;
            double penalty = penalty_base * dntParams.PenaltySafety;

            // set components
            var comps = XOp.EquationComponents[CodName];

            // convective part
            // ================
            if (config.isTransport) {

                comps.Add(new KineticEnergyConvectionAtLevelSet(D, LsTrk, rhoA, rhoB, LFFA, LFFB, physParams.Material, BcMap, config.isMovingMesh));
            }

            // viscous terms
            // =============
            if (config.isViscous && (!(muA == 0.0) && !(muB == 0.0))) {

                if(laplaceKinE)
                    comps.Add(new KineticEnergyLaplaceAtInterface(LsTrk, muA, muB, penalty * 1.0));

                if(config.getKinEviscousDiscretization == KineticEnergyViscousSourceTerms.fluxFormulation)
                    comps.Add(new StressDivergenceAtLevelSet(LsTrk, muA, muB, transposed: !laplaceKinE));
            }


            // pressure term
            // =============
            if (config.isPressureGradient) {

                if (divergenceP) {
                    comps.Add(new DivergencePressureEnergyAtLevelSet(LsTrk.GridDat.SpatialDimension));
                    //comps.Add(new ConvectivePressureTermAtLevelSet_LLF(D, LsTrk, LFFA, LFFB, physParams.Material, BcMap, config.isMovingMesh));
                }
            }

            // surface energy
            // ==============
            {
                comps.Add(new SurfaceEnergy(D, sigma, rhoA, rhoB));
            }

        }

        #endregion       

    }

    public interface IEnergy_Configuration : IXNSE_Configuration {

        KineticEnergyViscousSourceTerms getKinEviscousDiscretization { get; }

        KineticEnergyPressureSourceTerms getKinEpressureDiscretization { get; }

        bool withPressureDissipation { get; }

    }


    public enum KineticEnergyViscousSourceTerms {

        /// <summary>
        /// all source terms are evaluated locally
        /// </summary>
        local,

        /// <summary>
        /// source terms in divergence form are discretized over the flux formulation
        /// </summary>
        fluxFormulation,

        /// <summary>
        /// the transposed stress divergncee term is rewritten with respect to the kinetic energy
        /// </summary>
        laplaceKinE

    }

    public enum KineticEnergyPressureSourceTerms {

        /// <summary>
        /// the divergence of the convective pressure term is discretized in flux formulation
        /// </summary>
        divergence,

        /// <summary>
        /// only the convective pressure gradient is discretized locally
        /// </summary>
        convectiveGradP

    }

}
