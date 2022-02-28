using BoSSS.Application.XNSFE_Solver;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// Solver using the mass fraction approach for calculation of combustion
    /// </summary>
    public class XNSEC_MixtureFraction : XNSEC {

        /// <summary>
        /// dirty hack...
        /// </summary>
        protected override IncompressibleBoundaryCondMap GetBcMap() {
            if(boundaryMap == null)
                base.boundaryMap = new LowMachMixtureFractionMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, new string[] { "A", "B" });
            return boundaryMap;
        }

        public override void DefineScalarEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            Console.WriteLine("==================");
            Console.WriteLine("Solving the mixture fraction equations!");
            Console.WriteLine("==================");

            //================================
            // Mixture fraction
            //================================
            opFactory.AddEquation(new LowMachMixtureFraction("A", D, boundaryMap, config, EoS_A));
            opFactory.AddEquation(new LowMachMixtureFraction("B", D, boundaryMap, config, EoS_B));
        }

        protected override void DefineMomentumEquations(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int d, int D, LevelSetUpdater lsUpdater) {
            opFactory.AddEquation(new LowMachNavierStokes_MixtureFractions("A", d, D, boundaryMap, config, EoS_A));
            opFactory.AddEquation(new LowMachNavierStokes_MixtureFractions("B", d, D, boundaryMap, config, EoS_B));
            opFactory.AddEquation(new NSEInterface_Newton("A", "B", d, D, boundaryMap, config, config.isMovingMesh));
            if (config.isEvaporation) {
                opFactory.AddEquation(new InterfaceNSE_Evaporation_Newton("A", "B", D, d, config));
            }
        }

        protected override void DefineContinuityEquation(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater) {
            opFactory.AddEquation(new LowMachContinuity_MixtureFractions(D, "A", config, boundaryMap, EoS_A));
            opFactory.AddEquation(new LowMachContinuity_MixtureFractions(D, "B", config, boundaryMap, EoS_B));
            opFactory.AddEquation(new InterfaceContinuityLowMach(config, D, LsTrk, config.isMatInt));
            //=== evaporation extension === //
            if (config.isEvaporation) {
                opFactory.AddEquation(new InterfaceContinuity_Evaporation_Newton_LowMach("A", "B", D, config));
            }

        }

        override protected void DefineAditionalParameters(OperatorFactory opFactory, XNSEC_OperatorConfiguration config, int D, LevelSetUpdater lsUpdater, int quadOrder) {
            // No parameters
        }

        protected override void AddMultigridConfigLevel(List<MultigridOperator.ChangeOfBasisConfig> configsLevel, int iLevel) {
            int D = this.GridData.SpatialDimension;
            int pVel = VelocityDegree();
            int pPrs = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.Pressure].Degree;
            int pMxtFrctn = this.Control.FieldOptions[BoSSS.Solution.NSECommon.VariableNames.MixtureFraction].Degree;

            // configurations for velocity
            for(int d = 0; d < D; d++) {
                var configVel_d = new MultigridOperator.ChangeOfBasisConfig() {
                    DegreeS = new int[] { pVel },
                    mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib_DropIndefinite,
                    //mode = MultigridOperator.Mode.Eye,
                    VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.VelocityVector(D)[d]) }
                };
                configsLevel.Add(configVel_d);
            }
            // configuration for pressure
            var configPres = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pPrs },
                mode = MultigridOperator.Mode.IdMass_DropIndefinite,
                //mode = MultigridOperator.Mode.Eye,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.Pressure) }
            };
            configsLevel.Add(configPres);

            // configuration for mixture fraction
            var confTemp = new MultigridOperator.ChangeOfBasisConfig() {
                DegreeS = new int[] { pMxtFrctn },
                mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                //mode = MultigridOperator.Mode.Eye,
                VarIndex = new int[] { this.XOperator.DomainVar.IndexOf(VariableNames.MixtureFraction) }
            };
            configsLevel.Add(confTemp);
        }
    }
}