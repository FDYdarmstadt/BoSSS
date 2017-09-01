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
using BoSSS.Solution.NSECommon;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// SIMPLEStep for Low-Mach flows.
    /// </summary>
    public class SIMPLEStepLowMach : BaseSIMPLEStepVariableDensity {

        OperatorFactoryTemperature OperatorsTemperature;
        MatrixFactoryTemperature MatrixAssembliesTemperature;
        LowMachSIMPLEControl LowMachControl;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="WorkingSet"></param>
        /// <param name="WorkingSetMatrices"></param>
        public SIMPLEStepLowMach(SolverConfiguration SolverConf, VariableSet WorkingSet, VariableMatrices WorkingSetMatrices)
            : base(SolverConf, WorkingSet, WorkingSetMatrices) {
            this.LowMachControl = SolverConf.Control as LowMachSIMPLEControl;
            if (this.LowMachControl == null) {
                throw new ArgumentException("Invalid config", nameof(SolverConf));
            }

            // Construct SIMPLEOperators
            UnsetteledCoordinateMapping TemperatureMapping = new UnsetteledCoordinateMapping(WorkingSet.TemperatureBasis);
            UnsetteledCoordinateMapping PressureMapping = new UnsetteledCoordinateMapping(WorkingSet.PressureBasis);

            OperatorsTemperature = new OperatorFactoryTemperature(
                TemperatureMapping,
                PressureMapping,
                base.WorkingSet.Velocity.Current,
                base.WorkingSet.VelocityMean,
                base.WorkingSet.Temperature.Current,
                base.WorkingSet.TemperatureMean,
                SolverConf);

            // Construct matrix assemblies
            MatrixAssembliesTemperature = new MatrixFactoryTemperature(
                OperatorsTemperature, TemperatureMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, WorkingSetMatrices.Rho.Matrix, SolverConf, base.BDF);
        }

        protected override OperatorFactoryFlowFieldVariableDensity GetOperatorsFlowField(UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping, UnsetteledCoordinateMapping PressureMapping) {
            return (new OperatorFactoryFlowFieldLowMach(VelocityMapping, VelocityVectorMapping, PressureMapping,
                SolverConf,
                base.WorkingSet.Velocity.Current, base.WorkingSet.VelocityMean, base.WorkingSet.Temperature.Current, base.WorkingSet.TemperatureMean, base.WorkingSet.Eta));
        }

        public override void OverallIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            if ((LowMachControl.MaxFlowSolverIterations == 1) && (LowMachControl.MaxTemperatureSolverIterations == 1))
                FullyCoupledIteration(ref SIMPLEStatus, dt, ResLogger);
            else
                SegregatedIteration(ref SIMPLEStatus, dt, ResLogger);
        }

        /// <summary>
        /// One fully coupled iteration,
        /// i.e. solving flow field and temperature field in one iteration.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt"></param>
        /// <param name="ResLogger"></param>
        private void FullyCoupledIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            // Solve flow field and temperature field
            base.SolveFlowFieldEquations(ref SIMPLEStatus, dt, ResLogger);
            this.SolveScalarEquations(ref SIMPLEStatus, dt, ResLogger);

            // Check convergence
            if (ResLogger.Residuals["L2Norm p'"] < LowMachControl.L2NormPressureCorrection
                && ResLogger.Residuals["L2Norm u_res"] < LowMachControl.L2NormVelocityResidual
                && ResLogger.Residuals["L2Norm v_res"] < LowMachControl.L2NormVelocityResidual
                && ResLogger.Residuals["L2Norm T_res"] < LowMachControl.L2NormTemperatureResidual) {

                if (base.SolverConf.SpatialDimension == 2) {
                    SIMPLEStatus.IsConverged = true;
                } else if (ResLogger.Residuals["L2Norm w_res"] < LowMachControl.L2NormVelocityResidual) {
                    SIMPLEStatus.IsConverged = true;
                }
            }

            // Set SIMPLEStatus
            if (LowMachControl.MaxNoSIMPLEsteps == SIMPLEStatus.SIMPLEStepNo)
                SIMPLEStatus.TerminateSIMPLE = true;

            if ((LowMachControl.Algorithm == SolutionAlgorithms.Steady_SIMPLE) && (SIMPLEStatus.SIMPLEStepNo % LowMachControl.SavePeriodSIMPLE == 0))
                SIMPLEStatus.SaveStep = true;

            /*
            // Solve only temperature equation (only used for debugging).              
            this.SolveScalarEquations(ref SIMPLEStatus, dt, ResLogger);
            SIMPLEStatus.SaveStep = true;

            // Check convergence temperature field
            if ((ResLogger.Residuals["L2Norm T_res"] < base.SolverConf.L2NormTemperatureRes)
                || (SolverConf.MaxNoSIMPLEsteps == SIMPLEStatus.SIMPLEStepNo)) {

                // Set SIMPLEStatus
                SIMPLEStatus.TerminateSIMPLE = true;
            }

            // Workaround for plotting residuals
            ResLogger.ComputeL2Norm(WorkingSet.Pressure_Correction, "p'");
            ResLogger.ComputeL2Norm(WorkingSet.VelRes);*/
        }

        /// <summary>
        /// Counter for inner iterations within one segregated step.
        /// </summary>
        int CntInnerIterSegregatedStep = 1;

        /// <summary>
        /// One segregated iteration,
        /// i.e. solving either flow field or temperature field in one iteration.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt"></param>
        /// <param name="ResLogger"></param>
        private void SegregatedIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            if (CntInnerIterSegregatedStep <= LowMachControl.MaxFlowSolverIterations) {
                // Solve equations flow field
                base.SolveFlowFieldEquations(ref SIMPLEStatus, dt, ResLogger);

                // Check convergence flow field
                if (ResLogger.Residuals["L2Norm p'"] < LowMachControl.L2NormPressureCorrection
                    && ResLogger.Residuals["L2Norm u_res"] < LowMachControl.L2NormVelocityResidual
                    && ResLogger.Residuals["L2Norm v_res"] < LowMachControl.L2NormVelocityResidual) {

                    if (base.SolverConf.SpatialDimension == 2) {
                        CntInnerIterSegregatedStep = LowMachControl.MaxFlowSolverIterations;
                    } else if (ResLogger.Residuals["L2Norm w_res"] < LowMachControl.L2NormVelocityResidual) {
                        CntInnerIterSegregatedStep = LowMachControl.MaxFlowSolverIterations;
                    }
                }

                // Set SIMPLEStatus
                if (CntInnerIterSegregatedStep == LowMachControl.MaxFlowSolverIterations)
                    SIMPLEStatus.SaveStep = true;

                // Workaround for plotting residuals
                ResLogger.ComputeL2Norm(base.WorkingSet.TemperatureRes);
            } else {
                // Solve temperature equation
                this.SolveScalarEquations(ref SIMPLEStatus, dt, ResLogger);

                // Check convergence temperature field
                if ((ResLogger.Residuals["L2Norm T_res"] < LowMachControl.L2NormTemperatureResidual)
                    || (CntInnerIterSegregatedStep == LowMachControl.MaxFlowSolverIterations + LowMachControl.MaxTemperatureSolverIterations)) {
                    CntInnerIterSegregatedStep = 0;
                    SIMPLEStatus.NextSegregatedStep();
                    // Set SIMPLEStatus
                    SIMPLEStatus.SaveStep = true;
                    if (SIMPLEStatus.Timestep[1] > LowMachControl.MaxNoSIMPLEsteps)
                        SIMPLEStatus.IsConverged = true;
                }

                // Workaround for plotting residuals
                ResLogger.ComputeL2Norm(WorkingSet.Pressure_Correction, "p'");
                ResLogger.ComputeL2Norm(WorkingSet.VelRes);
            }
            CntInnerIterSegregatedStep++;
        }

        protected override ScalarFieldHistory<SinglePhaseField> GetScalarField() {
            return base.WorkingSet.Temperature;
        }

        protected override SIMPLESolver GetCorrectorSolver() {
            SIMPLESolver Corrector;

            Corrector = new LowMachSolverCorrector(
                base.SolverConf,
                LowMachControl.CorrectorSolverFactory(),
                base.MatrixAssembliesFlowField.Corrector,
                base.OperatorsFlowField.DivergenceConti,
                base.WorkingSet.Velocity_Intrmed,
                base.WorkingSet.DivB4,
                base.BDF,
                base.WorkingSet.Temperature,
                LowMachControl.EoS);

            return Corrector;
        }

        protected override void SolveScalarEquations(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            // Update temperature field SIMPLE operators (first!) and matrix assemblies (second!)
            // ==================================================================================

            OperatorsTemperature.TemperatureConvection.Update();
            OperatorsTemperature.HeatConduction.Update();
            MatrixAssembliesTemperature.Temperature.Update();
            if ((LowMachControl.RelaxationModeTemperature == RelaxationTypes.Implicit)
                && !MatrixAssembliesTemperature.TemperatureApprox.IsConstant
                && ((SIMPLEStatus.SIMPLEStepNo - 1) % LowMachControl.PredictorApproximationUpdateCycle == 0))
                MatrixAssembliesTemperature.TemperatureApprox.Update();

            // Temperature solver
            // ==================

            SIMPLESolver TemperatureSolver = new LowMachSolverTemperature(
                this.SolverConf,
                LowMachControl.TemperatureSolverFactory(),
                base.WorkingSetMatrices.Rho.Matrix,
                this.MatrixAssembliesTemperature.Temperature,
                this.MatrixAssembliesTemperature.TemperatureApprox,
                base.WorkingSet.Temperature,
                base.BDF,
                LowMachControl.EoS,
                base.WorkingSet.ThermodynamicPressure);

            base.WorkingSet.TemperatureRes.Clear();
            base.WorkingSet.TemperatureRes.Acc(1.0, base.WorkingSet.Temperature.Current);

            TemperatureSolver.Solve(base.WorkingSet.Temperature.Current.CoordinateVector, dt);
            TemperatureSolver.Dispose();

            // Update temperature field variables
            // ==================================

            SIMPLEStepUtils.UpdateScalarFieldVariables(
                LowMachControl,
                LowMachControl.RelaxationModeTemperature,
                LowMachControl.RelexationFactorTemperature,
                base.WorkingSet.Temperature,
                base.WorkingSet.TemperatureRes,
                base.WorkingSet.TemperatureMean,
                base.WorkingSet.Rho,
                base.WorkingSet.Eta,
                base.WorkingSetMatrices.Rho,
                LowMachControl.EoS,
                base.WorkingSet.ThermodynamicPressure.Current);

            // Calculate residuals temperature
            // ===============================

            ResLogger.ComputeL2Norm(base.WorkingSet.TemperatureRes);
        }
    }
}
