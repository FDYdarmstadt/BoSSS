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
using BoSSS.Platform.LinAlg;
using CNS.Boundary;

namespace CNS.Tests.IsentropicVortex {

    /// <summary>
    /// Exact solution for an isentropic vortex in a uniform background flow
    /// </summary>
    public class IsentropicVortexExactSolution {

        private CNSControl control;

        private double vortexSpeed;

        private bool periodic;

        private double BoxSize;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="control"></param>
        /// <param name="vortexSpeed"></param>
        /// <param name="periodic"></param>
        /// <param name="BoxSize"></param>
        public IsentropicVortexExactSolution(CNSControl control, double vortexSpeed, bool periodic = false, double BoxSize=0.0) {
            this.control = control;
            this.vortexSpeed = vortexSpeed;
            this.periodic = periodic;
            this.BoxSize = BoxSize;
        }

        /// <summary>
        /// The x-coordinate w.r.t. origin in the initial configuration
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private double x(double[] X, double t) {      
            double x = X[0] - vortexSpeed * t;
            if (periodic) {
                while (x < -BoxSize/2.0) {
                    x += BoxSize;
                }
            }
            return x;
        }

        /// <summary>
        /// Radial coordinate w.r.t the initial configuration
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private double r(double[] X, double t) {
            return Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
        }

        /// <summary>
        /// Angular coordinate w.r.t the initial configuration
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private double phi(double[] X, double t) {
            return Math.Atan2(X[1], x(X, t));
        }

        /// <summary>
        /// Velocity magnitude
        /// </summary>
        /// <param name="X"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        private double uAbs(double[] X, double t) {
            return r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
        }

        /// <summary>
        /// Density
        /// </summary>
        /// <returns></returns>
        public Func<double[], double, double> rho() {
            double gamma = control.EquationOfState.HeatCapacityRatio;
            return (X,t) => Math.Pow(
                1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
        }

        /// <summary>
        /// x-velocity
        /// </summary>
        /// <returns></returns>
        public Func<double[], double, double> u() {
            return (X, t) => vortexSpeed - Math.Sin(phi(X, t)) * uAbs(X, t);
        }

        /// <summary>
        /// y-velocity
        /// </summary>
        /// <returns></returns>
        public Func<double[], double, double> v() {
            return (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);
        }

        /// <summary>
        /// Pressure
        /// </summary>
        /// <returns></returns>
        public Func<double[], double, double> p() {
            var rho = this.rho();
            return (X, t) => Math.Pow(rho(X, t), control.EquationOfState.HeatCapacityRatio);
        }
    }
}
