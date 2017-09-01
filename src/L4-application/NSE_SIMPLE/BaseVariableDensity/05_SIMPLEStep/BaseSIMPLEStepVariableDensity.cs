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
using ilPSP.Tracing;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using NSE_SIMPLE.BaseVariableDensity;

namespace NSE_SIMPLE {

    /// <summary>
    /// Base class for SIMPLEStep for variable density flows,
    /// i.e. Low-Mach flows and multiphase.
    /// </summary>
    public abstract class BaseSIMPLEStepVariableDensity : ISIMPLEStep {

        protected SolverConfiguration SolverConf;

        // Variables
        // =========        
        protected VariableSet WorkingSet;
        protected VariableMatrices WorkingSetMatrices;

        protected BDFScheme BDF = null;

        // SIMPLEOperators
        // ===============
        protected OperatorFactoryFlowFieldVariableDensity OperatorsFlowField;

        // SIMPLEMatrixAssemblies
        // ======================
        protected MatrixFactoryVariableDensityFlowField MatrixAssembliesFlowField;

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="SolverConf"></param>
        /// <param name="WorkingSet"></param>        
        /// <param name="WorkingSetMatrices"></param>
        public BaseSIMPLEStepVariableDensity(SolverConfiguration SolverConf, VariableSet WorkingSet, VariableMatrices WorkingSetMatrices) {

            this.SolverConf = SolverConf;
            this.WorkingSet = WorkingSet;
            this.WorkingSetMatrices = WorkingSetMatrices;

            // Construct BDF scheme for unsteady solver
            if (SolverConf.Control.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE)
                BDF = new BDFScheme();

            // Construct SIMPLEOperators
            UnsetteledCoordinateMapping VelocityMapping = new UnsetteledCoordinateMapping(WorkingSet.VelBasis);
            UnsetteledCoordinateMapping PressureMapping = new UnsetteledCoordinateMapping(WorkingSet.PressureBasis);

            Basis[] VelBasisS = new Basis[SolverConf.SpatialDimension];
            for (int d = 0; d < SolverConf.SpatialDimension; d++) {
                VelBasisS[d] = WorkingSet.VelBasis;
            }
            UnsetteledCoordinateMapping VelocityVectorMapping = new UnsetteledCoordinateMapping(VelBasisS);

            OperatorsFlowField = GetOperatorsFlowField(VelocityMapping, VelocityVectorMapping, PressureMapping);

            // Construct matrix assemblies
            MatrixAssembliesFlowField = new MatrixFactoryVariableDensityFlowField(
                SolverConf,
                OperatorsFlowField,
                WorkingSetMatrices.Rho.Matrix,
                BDF,
                VelocityMapping,
                VelocityVectorMapping);
        }

        /// <summary>
        /// Get operators for flow field.
        /// </summary>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        /// <param name="PressureMapping"></param>  
        /// <returns></returns>
        protected abstract OperatorFactoryFlowFieldVariableDensity GetOperatorsFlowField(UnsetteledCoordinateMapping VelocityMapping,
            UnsetteledCoordinateMapping VelocityVectorMapping,
            UnsetteledCoordinateMapping PressureMapping);

        /// <summary>
        /// One SIMPLE iteration for low Mach number flows.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt"></param>
        /// <param name="ResLogger"></param>
        public abstract void OverallIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger);

        SIMPLESolver Corrector;
        protected bool UpdateApproximations;

        /// <summary>
        /// Solve flow field equations, i.e.
        /// Predictor for velocity and
        /// Corrector for pressure correction.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt"></param>
        /// <param name="ResLogger"></param>
        protected void SolveFlowFieldEquations(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            // Update flow field SIMPLE operators (first!) and matrix assemblies (second!)
            // ===========================================================================

            UpdateApproximations = (!MatrixAssembliesFlowField.PredictorApprox[0].IsConstant
                && ((SIMPLEStatus.SIMPLEStepNo - 1) % SolverConf.Control.PredictorApproximationUpdateCycle == 0)) ? true : false;

            for (int i = 0; i < SolverConf.SpatialDimension; i++) {
                OperatorsFlowField.Convection[i].Update();
                OperatorsFlowField.Visc[i].Update();
                OperatorsFlowField.Swip2[i].Update();
                if (SolverConf.Control.PhysicsMode == PhysicsMode.LowMach) {
                    OperatorsFlowField.Swip3[i].Update();
                    OperatorsFlowField.DivergenceConti[i].Update();
                }
                for (int j = 0; j < SolverConf.SpatialDimension; j++) {
                    MatrixAssembliesFlowField.ViscSplit[i, j].Update();
                }
                MatrixAssembliesFlowField.Predictor[i].Update();
            }

            if (UpdateApproximations) {
                for (int i = 0; i < SolverConf.SpatialDimension; i++) {
                    MatrixAssembliesFlowField.PredictorApprox[i].Update();
                    MatrixAssembliesFlowField.PredictorApproxInv[i].Update();
                }
                if (SolverConf.Control.PhysicsMode == PhysicsMode.Multiphase) {
                    MatrixAssembliesFlowField.Corrector.Update();
                }
            }

            if (SolverConf.Control.PhysicsMode == PhysicsMode.LowMach) {
                // For Low-Mach flows the corrector is never constant 
                // due to the density in the divergence operator.
                MatrixAssembliesFlowField.Corrector.Update();
            }

            // Predictor
            // =========

            for (int comp = 0; comp < SolverConf.SpatialDimension; comp++) {

                SIMPLESolver Predictor;

                VariableDensitySIMPLEControl varDensConf = SolverConf.Control as VariableDensitySIMPLEControl;
                Predictor = new VariableDensitySolverPredictor(
                    SolverConf,
                    SolverConf.Control.PredictorSolverFactory(),
                    WorkingSetMatrices.Rho.Matrix,
                    MatrixAssembliesFlowField.Predictor[comp],
                    MatrixAssembliesFlowField.PredictorApprox[comp],
                    OperatorsFlowField.PressureGradient,
                    WorkingSet.Pressure,
                    MatrixAssembliesFlowField.ViscSplit,
                    OperatorsFlowField.BuoyantForce,
                    WorkingSet.Velocity,
                    this.GetScalarField(),
                    varDensConf.EoS,
                    BDF);

                Predictor.Solve(WorkingSet.Velocity_Intrmed[comp].CoordinateVector, dt, comp);



                Predictor.Dispose();
            }

            // Corrector
            // =========

            if ((SIMPLEStatus.SIMPLEStepNoTotal == 1) || UpdateApproximations) {

                if (Corrector != null)
                    Corrector.Dispose();

                Corrector = this.GetCorrectorSolver();
            }
            Corrector.Solve(WorkingSet.Pressure_Correction.CoordinateVector, dt);

            // Update flow field variables
            // ===========================

            SIMPLEStepUtils.UpdateFlowFieldVariables(dt, SolverConf, OperatorsFlowField.PressureGradient, WorkingSet, BDF, MatrixAssembliesFlowField.PredictorApproxInv);
            if (WorkingSet.DivAfter != null) {
                WorkingSet.DivAfter.Clear();
                SolverUtils.CalculateMassDefect_Divergence(OperatorsFlowField.DivergenceConti, WorkingSet.Velocity.Current, WorkingSet.DivAfter);
            }

            // Calculate residuals flow field
            // ==============================

            ResLogger.ComputeL2Norm(WorkingSet.Pressure_Correction, "p'");
            ResLogger.ComputeL2Norm(WorkingSet.VelRes);
        }

        /// <summary>
        /// Get scalar field to calculate density and viscosity, i.e.
        /// level-set for multiphase and
        /// temperature for Low-Mach.
        /// </summary>
        /// <returns></returns>
        protected abstract ScalarFieldHistory<SinglePhaseField> GetScalarField();

        /// <summary>
        /// Get corrector solver for multiphase resp. Low-Mach.
        /// </summary>
        /// <returns></returns>
        protected abstract SIMPLESolver GetCorrectorSolver();

        /// <summary>
        /// Solve scalar equations, i.e.
        /// level-set equation for multiphase and
        /// temperature equation for Low-Mach.
        /// </summary>
        /// <param name="dt"></param>
        protected abstract void SolveScalarEquations(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger);

        ~BaseSIMPLEStepVariableDensity() {
            this.Dispose();
        }

        /// <summary>
        /// Dispose corrector solver.
        /// </summary>
        public void Dispose() {
            if (Corrector != null)
                Corrector.Dispose();
        }
    }
}
