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
using BoSSS.Solution;

namespace NSE_SIMPLE {

    /// <summary>
    /// SIMPLEStep for incompressible flows.
    /// </summary>
    public class SIMPLEStepIncompressible : ISIMPLEStep {

        SolverConfiguration m_SolverConf;

        VariableSet m_WorkingSet;

        OperatorFactoryFlowFieldIncompressible m_IncompressibleOperators;
        MatrixFactoryIncompressibleFlows m_IncompressibleMatrixAssemblies;

        BDFScheme m_BDF = null;

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="SolverConf"></param>
        /// <param name="WorkingSet"></param>
        public SIMPLEStepIncompressible(SolverConfiguration SolverConf, VariableSet WorkingSet) {

            m_SolverConf = SolverConf;
            m_WorkingSet = WorkingSet;

            // Construct BDF scheme for unsteady solver
            if (SolverConf.Control.Algorithm == SolutionAlgorithms.Unsteady_SIMPLE)
                m_BDF = new BDFScheme();

            // Construct all SIMPLEOperators, which are needed for incompressible flows
            UnsetteledCoordinateMapping VelocityMapping = new UnsetteledCoordinateMapping(WorkingSet.VelBasis);
            UnsetteledCoordinateMapping PressureMapping = new UnsetteledCoordinateMapping(WorkingSet.PressureBasis);

            m_IncompressibleOperators = new OperatorFactoryFlowFieldIncompressible(VelocityMapping,
                    PressureMapping,
                    SolverConf,
                    WorkingSet.Velocity.Current,
                    WorkingSet.VelocityMean);

            // Construct matrix assemblies
            m_IncompressibleMatrixAssemblies = new MatrixFactoryIncompressibleFlows(
                SolverConf, m_IncompressibleOperators, PressureMapping, WorkingSet.Pressure, m_BDF);
        }

        SIMPLESolver m_Corrector;

        /// <summary>
        /// One SIMPLE iteration for incompressible flows.
        /// </summary>
        /// <param name="SIMPLEStatus"></param>
        /// <param name="dt"></param>
        /// <param name="ResLogger"></param>
        public void OverallIteration(ref SIMPLEStepStatus SIMPLEStatus, double dt, ResidualLogger ResLogger) {
            using (new FuncTrace()) {

                //Update SIMPLE operators (first!) and matrix assemblies (second!)
                //================================================================

                for (int d = 0; d < m_SolverConf.SpatialDimension; d++) {
                    m_IncompressibleOperators.Convection[d].Update();
                    m_IncompressibleMatrixAssemblies.Predictor[d].Update();
                }

                //Determination whether approximation of predictor needs to be updated
                bool UpdatePredictorApprox;
                bool SetPredictorApproxAsConstant = false;
                switch (m_SolverConf.Control.PredictorApproximation) {
                    case PredictorApproximations.Identity:
                    case PredictorApproximations.Identity_IP1:
                        // For unsteady cases the approximation is 
                        // gamma * dt / (beta_0 + gamma * dt) * I,
                        // where gamma and beta_0 might not be constant for the first time steps
                        // depending on the current BDF order.
                        if (!m_IncompressibleMatrixAssemblies.PredictorApprox.IsConstant && (SIMPLEStatus.SIMPLEStepNo == 1)) {
                            if (m_SolverConf.BDFOrder < m_SolverConf.Control.TimeOrder) {
                                UpdatePredictorApprox = true;
                            } else if (m_SolverConf.BDFOrder == m_SolverConf.Control.TimeOrder) {
                                UpdatePredictorApprox = true;
                                SetPredictorApproxAsConstant = true;
                            } else {
                                throw new ApplicationException("Should not happen.");
                            }
                        } else {
                            UpdatePredictorApprox = false;
                        }
                        break;
                    case PredictorApproximations.Diagonal:
                    case PredictorApproximations.BlockDiagonal:
                        UpdatePredictorApprox = (!m_IncompressibleMatrixAssemblies.PredictorApprox.IsConstant
                            && ((SIMPLEStatus.SIMPLEStepNo - 1) % m_SolverConf.Control.PredictorApproximationUpdateCycle == 0)) ? true : false;
                        break;
                    default:
                        throw new NotImplementedException();
                }

                if (UpdatePredictorApprox) {
                    m_IncompressibleMatrixAssemblies.PredictorApprox.Update();
                    m_IncompressibleMatrixAssemblies.PredictorApproxInv.Update();
                    m_IncompressibleMatrixAssemblies.Corrector.Update();
                    if (SetPredictorApproxAsConstant) {
                        m_IncompressibleMatrixAssemblies.PredictorApprox.IsConstant = true;
                        m_IncompressibleMatrixAssemblies.PredictorApproxInv.IsConstant = true;
                        m_IncompressibleMatrixAssemblies.Corrector.IsConstant = true;
                    }
                }

                //Predictor
                //=========

                SIMPLESolver Predictor = new SolverPredictor(
                    m_SolverConf,
                    m_SolverConf.Control.PredictorSolverFactory(),
                    m_IncompressibleMatrixAssemblies.Predictor,
                    m_IncompressibleMatrixAssemblies.PredictorApprox,
                    m_IncompressibleOperators.PressureGradient,
                    m_BDF,
                    m_WorkingSet.Velocity,
                    m_WorkingSet.Pressure);

                for (int comp = 0; comp < m_SolverConf.SpatialDimension; comp++)
                    Predictor.Solve(m_WorkingSet.Velocity_Intrmed[comp].CoordinateVector, dt, comp);

                Predictor.Dispose();

                //Corrector
                //=========

                if ((SIMPLEStatus.SIMPLEStepNoTotal == 1) || UpdatePredictorApprox) {

                    if (m_Corrector != null)
                        m_Corrector.Dispose();

                    m_Corrector = new SolverCorrector(
                        m_SolverConf, m_SolverConf.Control.CorrectorSolverFactory(),
                        m_IncompressibleMatrixAssemblies.Corrector, m_IncompressibleOperators.VelocityDivergence,
                        m_BDF, m_WorkingSet.Velocity_Intrmed, m_WorkingSet.DivB4,
                        m_IncompressibleOperators.PressureStabilization, m_WorkingSet.Pressure);
                }

                m_Corrector.Solve(m_WorkingSet.Pressure_Correction.CoordinateVector, dt);

                //Update variables
                //================

                SIMPLEStepUtils.UpdateFlowFieldVariables(dt, m_SolverConf, m_IncompressibleOperators.PressureGradient, m_WorkingSet, m_BDF, m_IncompressibleMatrixAssemblies.PredictorApproxInv);
                if (m_WorkingSet.DivAfter != null) {
                    m_WorkingSet.DivAfter.Clear();
                    SolverUtils.CalculateMassDefect_Divergence(m_IncompressibleOperators.VelocityDivergence, m_WorkingSet.Velocity.Current, m_WorkingSet.DivAfter);
                }

                //Calculate residuals
                //===================

                ResLogger.ComputeL2Norm(m_WorkingSet.Pressure_Correction, "p'");
                ResLogger.ComputeL2Norm(m_WorkingSet.VelRes);

                //Check convergence
                //=================

                if (ResLogger.Residuals["L2Norm p'"] < m_SolverConf.Control.L2NormPressureCorrection
                    && ResLogger.Residuals["L2Norm u_res"] < m_SolverConf.Control.L2NormVelocityResidual
                    && ResLogger.Residuals["L2Norm v_res"] < m_SolverConf.Control.L2NormVelocityResidual) {

                    if (m_SolverConf.SpatialDimension == 2) {
                        SIMPLEStatus.IsConverged = true;
                    } else if (ResLogger.Residuals["L2Norm w_res"] < m_SolverConf.Control.L2NormVelocityResidual) {
                        SIMPLEStatus.IsConverged = true;
                    }
                }

                //Set SIMPLEStatus
                //================

                if (m_SolverConf.Control.MaxNoSIMPLEsteps == SIMPLEStatus.SIMPLEStepNo)
                    SIMPLEStatus.TerminateSIMPLE = true;

                if ((m_SolverConf.Control.Algorithm == SolutionAlgorithms.Steady_SIMPLE) && (SIMPLEStatus.SIMPLEStepNo % m_SolverConf.Control.SavePeriodSIMPLE == 0))
                    SIMPLEStatus.SaveStep = true;
            }
        }

        ~SIMPLEStepIncompressible() {
            this.Dispose();
        }

        /// <summary>
        /// Dispose corrector solver.
        /// </summary>
        public void Dispose() {
            if (m_Corrector != null)
                m_Corrector.Dispose();
        }
    }
}
