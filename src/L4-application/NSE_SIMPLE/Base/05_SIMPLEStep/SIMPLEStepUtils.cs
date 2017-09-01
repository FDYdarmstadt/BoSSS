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
using ilPSP.Tracing;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution;
using NSE_SIMPLE.LowMach;

namespace NSE_SIMPLE {

    /// <summary>
    /// Utility functions for SIMPLEStep, which are used for incompressible flows
    /// as well as low Mach number flows.
    /// </summary>
    public static class SIMPLEStepUtils {

        /// <summary>
        /// Update flow field variables after solving Predictor and Corrector.
        /// </summary>
        /// <param name="dt"></param>
        /// <param name="SolverConf"></param>
        /// <param name="PressureGradient"></param>        
        /// <param name="WorkingSet"></param>
        /// <param name="BDF"></param>
        /// <param name="PredictorApproxInv"></param>
        public static void UpdateFlowFieldVariables(double dt, SolverConfiguration SolverConf, SIMPLEOperator[] PressureGradient,
            VariableSet WorkingSet, BDFScheme BDF, params SIMPLEMatrixAssembly[] PredictorApproxInv) {

            using (var t = new FuncTrace()) {

                // Pressure - explicit under-relaxation
                // ====================================

                WorkingSet.Pressure.Acc(SolverConf.Control.RelaxationFactorPressure, WorkingSet.Pressure_Correction);

                if (!SolverConf.BcMap.DirichletPressureBoundary) {
                    if (!double.IsNaN(SolverConf.Control.PressureMeanValue)) {
                        // set specified mean value of pressure
                        double NegCalcMeanValuePressure = -1.0 * WorkingSet.Pressure.GetMeanValueTotal(null);
                        WorkingSet.Pressure.AccConstant(NegCalcMeanValuePressure);
                        WorkingSet.Pressure.AccConstant(SolverConf.Control.PressureMeanValue);
                    }
                    if ((SolverConf.Control.PressureReferencePoint != null) && double.IsNaN(SolverConf.Control.PressureMeanValue)) {
                        // set pressure to zero at the specified reference point
                        double ShiftPressure = -1.0 * WorkingSet.Pressure.ProbeAt(SolverConf.Control.PressureReferencePoint);
                        WorkingSet.Pressure.AccConstant(ShiftPressure);
                    }
                }

                // Velocity
                // ========

                for (int comp = 0; comp < SolverConf.SpatialDimension; comp++) {

                    // Velocity_Correction = - ApproximationPredictorMatrix^(-1) * PressureGradient * pri_Pressure_Correction
                    double[] GradientPressureCorrection = new double[PressureGradient[0].LocalLength];
                    PressureGradient[comp].OperatorMatrix.SpMVpara(1.0, WorkingSet.Pressure_Correction.CoordinateVector,
                        0.0, GradientPressureCorrection);
                    if (PredictorApproxInv.Length > 1) {
                        PredictorApproxInv[comp].AssemblyMatrix.SpMVpara(-1.0, GradientPressureCorrection,
                            0.0, WorkingSet.Velocity_Correction[comp].CoordinateVector);
                    } else {
                        PredictorApproxInv[0].AssemblyMatrix.SpMVpara(-1.0, GradientPressureCorrection,
                            0.0, WorkingSet.Velocity_Correction[comp].CoordinateVector);
                    }

                    for (int i = 0; i < WorkingSet.Velocity_Correction[comp].Mapping.LocalLength; i++) {
                        // VelRes = Vel(n) - Vel(n-1)
                        WorkingSet.VelRes[comp].CoordinateVector[i] =
                            WorkingSet.Velocity_Intrmed[comp].CoordinateVector[i] + WorkingSet.Velocity_Correction[comp].CoordinateVector[i]
                            - WorkingSet.Velocity.Current[comp].CoordinateVector[i];

                        // Vel = Vel^* + Vel'
                        WorkingSet.Velocity.Current[comp].CoordinateVector[i] =
                            WorkingSet.Velocity_Intrmed[comp].CoordinateVector[i] + WorkingSet.Velocity_Correction[comp].CoordinateVector[i];
                    }
                }

                // VelocityMean
                // ============

                WorkingSet.VelocityMean.Clear();
                WorkingSet.VelocityMean.AccLaidBack(1.0, WorkingSet.Velocity.Current);
            }
        }

        /// <summary>
        /// Update scalar field variables after solving scalar equation.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="ModeRelaxScalar"></param>
        /// <param name="relax_scalar"></param>
        /// <param name="Scalar"></param>
        /// <param name="ScalarRes"></param>
        /// <param name="ScalarMean"></param>
        /// <param name="Rho"></param>
        /// <param name="Eta"></param>
        /// <param name="RhoMatrix"></param>
        /// <param name="EoS"></param>
        /// <param name="ThermodynamicPressure">Null for multiphase flows.</param>
        public static void UpdateScalarFieldVariables(SIMPLEControl SolverConf, RelaxationTypes ModeRelaxScalar, double relax_scalar,
            ScalarFieldHistory<SinglePhaseField> Scalar, SinglePhaseField ScalarRes, SinglePhaseField ScalarMean,
            SinglePhaseField Rho, SinglePhaseField Eta, QuadratureMatrix RhoMatrix, MaterialLaw EoS,
            SinglePhaseField ThermodynamicPressure, bool UpdateRhoVisc = true) {
            using (new FuncTrace()) {

                // Explicit Under-Relaxation of scalar variable
                // ============================================

                if (ModeRelaxScalar == RelaxationTypes.Explicit) {
                    // phi = alpha * phi_new + (1-alpha) * phi_old
                    Scalar.Current.Scale(relax_scalar);
                    Scalar.Current.Acc((1.0 - relax_scalar), ScalarRes);
                }

                // Scalar residual
                // ===============

                ScalarRes.Scale(-1.0);
                ScalarRes.Acc(1.0, Scalar.Current);

                // ScalarMean
                // ==========

                ScalarMean.Clear();
                ScalarMean.AccLaidBack(1.0, Scalar.Current);

                // Thermodynamic pressure - only for Low-Mach number flows
                // =======================================================

                switch (SolverConf.PhysicsMode) {
                    case PhysicsMode.LowMach:
                        LowMachSIMPLEControl lowMachConf = SolverConf as LowMachSIMPLEControl;
                        if (lowMachConf.ThermodynamicPressureMode == ThermodynamicPressureMode.MassDetermined) {
                            ThermodynamicPressure.Clear();
                            ThermodynamicPressure.AccConstant(((MaterialLawLowMach)EoS).GetMassDeterminedThermodynamicPressure(lowMachConf.InitialMass.Value, Scalar.Current));
                        }
                        break;
                    case PhysicsMode.Multiphase:
                        break;
                    default:
                        throw new ApplicationException();
                }

                if (UpdateRhoVisc) {

                    // Density
                    // =======

                    Rho.Clear();
                    Rho.ProjectFunction(1.0, EoS.GetDensity, null, Scalar.Current);
                    RhoMatrix.Update();

                    // Viscosity
                    // =========

                    Eta.Clear();
                    Eta.ProjectFunction(1.0, EoS.GetViscosity, null, Scalar.Current);

                }
            }
        }
    }
}
