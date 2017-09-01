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

namespace CNS.Boundary {

    /// <summary>
    /// Implementation of a pressure boundary condition for subsonic inlets
    /// where the stagnation pressure and the stagnation temperature are
    /// prescribed.
    /// </summary>
    public class SubsonicPressureInlet : BoundaryCondition {

        /// <summary>
        /// The prescribed stagnation pressure
        /// </summary>
        public readonly Func<double[], double, double> TotalPressureFunction;

        /// <summary>
        /// The prescribed stagnation temperature
        /// </summary>
        public readonly Func<double[], double, double> TotalTemperatureFunction;

        /// <summary>
        /// <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="pressureFunction">
        /// The prescribed (stagnation) pressure
        /// </param>
        /// <param name="temperatureFunction">
        /// The prescribed (stagnation) temperature
        /// </param>
        public SubsonicPressureInlet(CNSControl config, Func<double[], double, double> pressureFunction, Func<double[], double, double> temperatureFunction)
            : base(config) {
            this.TotalPressureFunction = pressureFunction;
            this.TotalTemperatureFunction = temperatureFunction;
        }

        /// <summary>
        /// Calculates the state at the boundary from the prescribed stagnation
        /// pressure and the prescribed stagnation temperature. The calculation
        /// follows FerzigerPeric2001 (p. 315ff). The flow is prescribed to be
        /// normal to the edge.
        /// </summary>
        /// <param name="time">
        /// <see cref="BoundaryCondition.GetBoundaryState"/>
        /// </param>
        /// <param name="x">
        /// <see cref="BoundaryCondition.GetBoundaryState"/>
        /// </param>
        /// <param name="normal">
        /// <see cref="BoundaryCondition.GetBoundaryState"/>
        /// </param>
        /// <param name="stateIn">
        /// <see cref="BoundaryCondition.GetBoundaryState"/>
        /// </param>
        /// <returns>
        /// \f$ 
        /// T^+ = T_t(x) \frac{p_t(x)}{p^-}^{\frac{1.0 - \gamma}{\gamma}}
        /// \f$ 
        /// \f$ 
        /// |\vec{u^+}| = \sqrt{\frac{2 T^+}{(\gamma - 1) \mathrm{Ma}_\infty^2} \left(\frac{T_t(x)}{T^+} -1\right)}
        /// \f$ 
        /// \f$ \rho^+ = \frac{p^- }{T^+}\f$ 
        /// \f$ 
        /// (\rho^+, -\rho^+ |\vec{u^+}| \vec{n}, \frac{p}{\kappa - 1} + \frac{rho^+ |\vec{u^+}|^2}{2})^T
        /// \f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            double gamma = config.EquationOfState.HeatCapacityRatio;
            double Mach = config.MachNumber;
            Vector3D inwardNormal = new Vector3D();
            for (int i = 0; i < normal.Length; i++) {
                inwardNormal[i] = -normal[i];
            }

            double p0 = TotalPressureFunction(x, time);
            double T0 = TotalTemperatureFunction(x, time);

            double p = stateIn.Pressure;
            double T = TotalTemperatureFunction(x, time) * Math.Pow(p0 / p, (1.0 - gamma) / gamma);
            double rho = p / T;

            double VelocitySquare = 2.0 * T / ((gamma - 1.0) * (Mach* Mach)) * (T0/T - 1.0);
            Vector3D velocityOut = Math.Sqrt(VelocitySquare) * inwardNormal;

            return StateVector.FromPrimitiveQuantities(stateIn.Material, rho, velocityOut, p);
        }
    }
}
