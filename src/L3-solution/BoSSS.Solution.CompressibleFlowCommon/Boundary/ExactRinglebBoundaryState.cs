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

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using System;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// A boundary condition that represents the exact solution of the Ringleb
    /// test case.
    /// </summary>
    public class ExactRinglebBoundaryState : BoundaryCondition {

        /// <summary>
        /// Constructs a new boundary state
        /// </summary>
        /// <param name="control"></param>
        public ExactRinglebBoundaryState(Material config)
            : base(config) {
        }

        /// <summary>
        /// Uses <see cref="RinglebExactSolution"/> in order to construct the
        /// boundary state corresponding to the exact solution of the Ringleb
        /// problem.
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="stateIn"></param>
        /// <returns></returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            RinglebMaterial ringlebMaterial = (RinglebMaterial)config;
            double kappa = ringlebMaterial.EquationOfState.HeatCapacityRatio;
            double pi = ringlebMaterial.RinglebReferencePressure;
            double a0 = ringlebMaterial.RinglebReferenceSpeedOfSound;
            double p0 = ringlebMaterial.RinglebReferenceTotalPressure;

            RinglebExactSolution.FlowState ringlebState = RinglebExactSolution.GetFlowState(
                x[0], x[1], kappa, pi, a0, p0);

            double rho = ringlebState.Density;
            Vector v = new Vector(ringlebState.Velocity[0], ringlebState.Velocity[1]);
            return new StateVector(
                stateIn.Material,
                rho,
                rho * v,
                (ringlebState.Pressure + kappa * pi) / (kappa - 1.0) + 0.5 * rho * v * v);
        }

        /// <summary>
        /// A hack to pass additional values required for the Ringleb solution
        /// </summary>
        public class RinglebMaterial : Material {

            /// <summary>
            /// Constructor
            /// </summary>
            public RinglebMaterial(IEquationOfState __EquationOfState, IViscosityLaw __ViscosityLaw, double __MachNumber, double __ReynoldsNumber, double __PrandtlNumber, double __FroudeNumber, double __ViscosityRatio)
                : base(__EquationOfState, __ViscosityLaw, __MachNumber, __ReynoldsNumber, __PrandtlNumber, __FroudeNumber, __ViscosityRatio) {
            }

            /// <summary>
            /// Reference pressure
            /// </summary>
            public double RinglebReferencePressure;

            /// <summary>
            /// Reference speed of sound
            /// </summary>
            public double RinglebReferenceSpeedOfSound;

            /// <summary>
            /// Reference total pressure
            /// </summary>
            public double RinglebReferenceTotalPressure;
        }

        /// <summary>
        /// Computes the exact solution of the Ringleb problem according to
        /// Emanuel2010 (with some extensions for stiffened gases that can be found
        /// in Mueller2014).
        /// </summary>
        public static class RinglebExactSolution {

            /// <summary>
            /// Constructs the complete flow state in the point
            /// (<paramref name="x"/>, <paramref name="y"/>) for the given set of
            /// parameters.
            /// </summary>
            /// <param name="x">x-coordinate</param>
            /// <param name="y">y-coordinate</param>
            /// <param name="kappa">Heat capacity ratio</param>
            /// <param name="pi">
            /// Pressure offset in case of a stiffened gas, cf.
            /// <see cref="MaterialProperty.StiffenedGas"/>
            /// </param>
            /// <param name="a0">
            /// The characteristic stagnation speed of sound of the flow.
            /// </param>
            /// <param name="p0">
            /// The characteristic stagnation pressure of the flow. Beware: In case
            /// of a stiffened gas, larger values are required in order to obtain a
            /// valid flow field
            /// </param>
            /// <returns>
            /// The exact solution to the Ringleb problem in
            /// (<paramref name="x"/>, <paramref name="y"/>).
            /// </returns>
            public static FlowState GetFlowState(double x, double y, double kappa, double pi, double a0, double p0) {
                double tau = ComputeTau(x, y, kappa);
                // Compute theta from y because using x may cause problems
                double theta = 0.5 * Math.Acos(2.0 * (y - EvaluateG(tau, kappa)) *
                    Math.Sqrt(2.0 / (kappa - 1.0)) * tau * Math.Pow(1.0 - tau, 1.0 / (kappa - 1.0)));

                double rho0 = kappa * (p0 + pi) / a0 / a0;
                var state = new FlowState(
                    V: Math.Sqrt(2.0 * a0 * a0 * tau / (kappa - 1.0)),
                    theta: theta,
                    tau: tau,
                    pressure: (p0 + pi) * Math.Pow(1.0 - tau, kappa / (kappa - 1.0)) - pi,
                    density: rho0 * Math.Pow(1.0 - tau, 1.0 / (kappa - 1.0)),
                    velocity_2d: a0 * Math.Sqrt(2.0 * tau / (kappa - 1.0)) *
                        new Vector(Math.Cos(theta), Math.Sin(theta)),
                    Mach: Math.Sqrt(2.0 / (kappa - 1.0) * tau / (1.0 - tau)),
                    kappa: kappa,
                    pi: pi);
                return state;
            }

            /// <summary>
            /// Numerically determines the root of the function
            /// \f$ 
            /// f = x^2 + (y-g(\tau))^2 - R(\tau)^2
            /// \f$ 
            /// where
            /// \f$ 
            /// R(tau) = \frac{1}{2} \frac{\sqrt{\frac{(\kappa - 1)}{2}}}{\tau (1-\tau)^{\frac{1}{1-\tau}}}
            /// \f$ 
            /// and \f$ g(\tau)\f$  is computed via
            /// <see cref="EvaluateG"/>. The method uses a simple bisection algorithm.
            /// </summary>
            /// <param name="x">x-coordinate</param>
            /// <param name="y">y-coordinate</param>
            /// <param name="kappa">Heat capacity ratio.</param>
            /// <returns>
            /// A root of \f$ f\f$ 
            /// </returns>
            private static double ComputeTau(double x, double y, double kappa) {
                double tauMin = 0.00000000001;
                double tauMax = 0.99999999999;

                int n = 0;
                while (true) {
                    if (n >= 100) {
                        throw new System.Exception("No convergence");
                    }

                    double tau = 0.5 * (tauMax + tauMin);

                    double g = EvaluateG(tau, kappa);
                    if (double.IsNaN(g)) {
                        throw new Exception(
                            "Encountered NaN in the function 'g' of the Ringleb solution");
                    }

                    double R = 0.5 * Math.Sqrt((kappa - 1.0) / 2.0) /
                        (tau * Math.Pow(1.0 - tau, 1.0 / (kappa - 1.0)));
                    double f = x * x + (y - g) * (y - g) - R * R;

                    if (Math.Abs(f) < 1e-13) {
                        break;
                    } else if (f > 0.0) {
                        tauMax = tau;
                    } else {
                        tauMin = tau;
                    }

                    n++;
                }

                return tauMin;
            }

            /// <summary>
            /// Evaluates the function \f$ g(\tau)\f$ 
            /// required by <see cref="ComputeTau"/>. The explicit formula for 
            /// \f$ g\f$  depends on <paramref name="kappa"/>.
            /// This method is currently implemented for <paramref name="kappa"/>
            /// = 1.4 and <paramref name="kappa"/> = 7.
            /// </summary>
            /// <param name="tau">Evaluation point</param>
            /// <param name="kappa">Heat capacity ratio</param>
            /// <returns>
            /// The function value \f$ g(\tau)\f$ 
            /// </returns>
            private static double EvaluateG(double tau, double kappa) {
                if (kappa == 1.4) {
                    double d = Math.Sqrt(1.0 - tau);
                    return -1.0 / Math.Sqrt(8.0 * (kappa - 1.0)) * (
                        (30.0 * tau * tau - 70.0 * tau + 46.0) / 15.0 / Math.Pow(1.0 - tau, 2.5)
                        - Math.Log((1.0 + d) / (1.0 - d)));
                } else if (kappa == 7.0) {
                    double d = Math.Pow(1.0 - tau, 1.0 / 6.0);
                    return -1.0 / Math.Sqrt(8.0 * (kappa - 1.0)) * (
                        6.0 / d
                        + Math.Sqrt(3.0) * (
                            Math.Atan((2.0 * d - 1.0) / Math.Sqrt(3.0)) +
                            Math.Atan((2.0 * d + 1.0) / Math.Sqrt(3.0)))
                        + Math.Log(1.0 - d)
                        - Math.Log(1.0 + d)
                        + 0.5 * Math.Log(1.0 - d + Math.Pow(1.0 - tau, 1.0 / 3.0))
                        - 0.5 * Math.Log(1.0 + d + Math.Pow(1.0 - tau, 1.0 / 3.0)));
                } else {
                    throw new NotImplementedException();
                }
            }

            /// <summary>
            /// Helper class to represent the flow state in a given point of the
            /// computational domain.
            /// </summary>
            public struct FlowState {

                /// <summary>
                /// Velocity magnitude
                /// </summary>
                public readonly double V;

                /// <summary>
                /// Flow angle
                /// </summary>
                public readonly double Theta;

                /// <summary>
                /// Characteristic Mach number
                /// </summary>
                public readonly double Tau;

                /// <summary>
                /// Pressure
                /// </summary>
                public readonly double Pressure;

                /// <summary>
                /// Density
                /// </summary>
                public readonly double Density;

                /// <summary>
                /// Velocity
                /// </summary>
                public readonly Vector Velocity;

                /// <summary>
                /// Momentum
                /// </summary>
                public readonly Vector Momentum;

                /// <summary>
                /// Local Mach number
                /// </summary>
                public readonly double Mach;

                /// <summary>
                /// Total energy per volume
                /// </summary>
                public readonly double Energy;

                /// <summary>
                /// Constructs a new flow state from the given data
                /// </summary>
                /// <param name="V"></param>
                /// <param name="theta"></param>
                /// <param name="tau"></param>
                /// <param name="pressure"></param>
                /// <param name="density"></param>
                /// <param name="velocity_2D"></param>
                /// <param name="Mach"></param>
                /// <param name="kappa"></param>
                /// <param name="pi"></param>
                public FlowState(double V, double theta, double tau, double pressure, double density, Vector velocity_2d, double Mach, double kappa, double pi) {
                    this.V = V;
                    this.Theta = theta;
                    this.Tau = tau;
                    this.Pressure = pressure;
                    this.Density = density;
                    this.Velocity = velocity_2d;
                    this.Mach = Mach;
                    this.Momentum = density * velocity_2d;
                    this.Energy = (pressure + kappa * pi) / (kappa - 1.0) + 0.5 * Momentum * velocity_2d;
                }
            }
        }
    }
}