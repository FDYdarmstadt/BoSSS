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
using CNS.MaterialProperty;

namespace CNS.Boundary {

    /// <summary>
    /// Calculation of boundary values for an outlet with a Mach number smaller
    /// than 1.0.
    /// </summary>
    public class SubsonicOutlet : BoundaryCondition {

        /// <summary>
        /// The function for the prescribed (dimensionless) pressure in the
        /// free stream.
        /// </summary>
        private Func<double[], double, double> pressureFunction;

        /// <summary>
        /// Sets <see cref="pressureFunction"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="pressureFunction">
        /// The function for the prescribed (dimensionless) pressure in the
        /// free stream.
        /// </param>
        public SubsonicOutlet(CNSControl config, Func<double[], double, double> pressureFunction)
            : base(config) {
            this.pressureFunction = pressureFunction;
        }

        /// <summary>
        /// Calculates the boundary values for a subsonic outlet (Mach number
        /// smaller than 1). We have to impose one condition (here, we choose
        /// the pressure <see cref="pressureFunction"/>) and extrapolate the
        /// other values following FerzigerPeric2001 (p. 315ff)
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
        /// (\rho^-, u^-, p^*)^T
        /// \f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            return StateVector.FromPrimitiveQuantities(
                stateIn.Material,
                stateIn.Density,
                stateIn.Velocity,
                pressureFunction(x, time));
        }
    }
}
