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
using BoSSS.Solution;
using ilPSP.Tracing;
using NSE_SIMPLE.Multiphase;

namespace NSE_SIMPLE {

    /// <summary>
    /// SIMPLEStep for multiphase flows.
    /// </summary>
    public class SIMPLEStepMultiphase : BaseSIMPLEStepVariableDensity {

        OperatorFactoryLevelSet OperatorsLevelSet;
        MatrixFactoryLevelSet MatrixAssembliesLevelSet;
        
        MultiphaseSIMPLEControl MultiphaseControl;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="WorkingSet"></param>
        /// <param name="WorkingSetMatrices"></param>
        public SIMPLEStepMultiphase(SolverConfiguration SolverConf, VariableSet WorkingSet, VariableMatrices WorkingSetMatrices)
            : base(SolverConf, WorkingSet, WorkingSetMatrices) {

            this.SolverConf = SolverConf;
            this.MultiphaseControl = SolverConf.Control as MultiphaseSIMPLEControl;
            if (this.MultiphaseControl == null) {
                throw new ArgumentException("Invalid configuration", nameof(SolverConf)); 
            }

            // Construct SIMPLEOperators
            UnsetteledCoordinateMapping LevelSetMapping = new UnsetteledCoordinateMapping(WorkingSet.LevelSetBasis);

            OperatorsLevelSet = new OperatorFactoryLevelSet(LevelSetMapping,
                WorkingSet.Velocity.Current,
                WorkingSet.VelocityMean,
                SolverConf);

            // Construct matrix assemblies
            MatrixAssembliesLevelSet = new MatrixFactoryLevelSet(OperatorsLevelSet, LevelSetMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, SolverConf, base.BDF);
        }

        protected override OperatorFactoryFlowFieldVariableDensity GetOperatorsFlowField(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping, UnsetteledCoordinateMapping PressureMapping) {
            return (new OperatorFactoryFlowFieldMultiphase(VelocityMapping, VelocityVectorMapping, PressureMapping,
                base.SolverConf,
                base.WorkingSet.Velocity.Current, base.WorkingSet.VelocityMean, base.WorkingSet.Phi.Current, base.WorkingSet.PhiMean, base.WorkingSet.Eta));
        }

        public override void OverallIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            using (new FuncTrace()) {
                // Solve flow field and scalar field
                base.SolveFlowFieldEquations(ref SIMPLEStatus, dt, ResLogger);
                this.SolveScalarEquations(ref SIMPLEStatus, dt, ResLogger);

                // Check convergence
                if (ResLogger.Residuals["L2Norm p'"] < MultiphaseControl.L2NormPressureCorrection
                    && ResLogger.Residuals["L2Norm u_res"] < MultiphaseControl.L2NormVelocityResidual
                    && ResLogger.Residuals["L2Norm v_res"] < MultiphaseControl.L2NormVelocityResidual
                    && ResLogger.Residuals["L2Norm Phi_res"] < MultiphaseControl.L2NormLevelSetResidual) {

                    if (SolverConf.SpatialDimension == 2) {
                        SIMPLEStatus.IsConverged = true;
                    } else if (ResLogger.Residuals["L2Norm w_res"] < MultiphaseControl.L2NormVelocityResidual) {
                        SIMPLEStatus.IsConverged = true;
                    }
                }

                // Set SIMPLEStatus
                if (MultiphaseControl.MaxNoSIMPLEsteps == SIMPLEStatus.SIMPLEStepNo)
                    SIMPLEStatus.TerminateSIMPLE = true;

                if ((MultiphaseControl.Algorithm == SolutionAlgorithms.Steady_SIMPLE) && (SIMPLEStatus.SIMPLEStepNo % MultiphaseControl.SavePeriodSIMPLE == 0))
                    SIMPLEStatus.SaveStep = true;
            }
        }

        protected override ScalarFieldHistory<SinglePhaseField> GetScalarField() {
            return base.WorkingSet.Phi;
        }

        protected override SIMPLESolver GetCorrectorSolver() {
            SIMPLESolver Corrector;

            Corrector = new MultiphaseSolverCorrector(
                base.SolverConf,
                MultiphaseControl.CorrectorSolverFactory(),
                base.OperatorsFlowField.DivergenceConti,
                base.MatrixAssembliesFlowField.Corrector,
                base.WorkingSet.Velocity_Intrmed,
                base.WorkingSet.DivB4);

            return Corrector;
        }

        protected override void SolveScalarEquations(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            // Update level-set SIMPLE operators (first!) and matrix assemblies (second!)
            // ==========================================================================
            OperatorsLevelSet.LevelSetAdvection.Update();
            MatrixAssembliesLevelSet.LevelSet.Update();            
            if (base.UpdateApproximations && (MultiphaseControl.LevelSetRelaxationType == RelaxationTypes.Implicit)) {
                MatrixAssembliesLevelSet.LevelSetApprox.Update();                
            }

            // Level-Set solver
            // ================

            SIMPLESolver LevelSetSolver = new MultiphaseSolverLevelSet(
                SolverConf,
                MultiphaseControl.LevelSetSolverFactory(),
                MatrixAssembliesLevelSet.LevelSet,
                MatrixAssembliesLevelSet.LevelSetApprox,
                base.WorkingSet.Phi,
                base.BDF);

            base.WorkingSet.PhiRes.Clear();
            base.WorkingSet.PhiRes.Acc(1.0, base.WorkingSet.Phi.Current);

            LevelSetSolver.Solve(base.WorkingSet.Phi.Current.CoordinateVector, dt);
            LevelSetSolver.Dispose();

            // Update level-set variables
            // ==========================

            SIMPLEStepUtils.UpdateScalarFieldVariables(
                MultiphaseControl,
                MultiphaseControl.LevelSetRelaxationType,
                MultiphaseControl.RelaxationFactorLevelSet,
                base.WorkingSet.Phi,
                base.WorkingSet.PhiRes,
                base.WorkingSet.PhiMean,
                base.WorkingSet.Rho,
                base.WorkingSet.Eta,
                base.WorkingSetMatrices.Rho,
                MultiphaseControl.EoS,
                null);

            // Calculate residuals level-set
            // =============================

            ResLogger.ComputeL2Norm(base.WorkingSet.PhiRes);
        }        
    }
}
