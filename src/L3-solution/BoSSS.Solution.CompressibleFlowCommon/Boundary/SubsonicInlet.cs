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
using BoSSS.Solution.CompressibleFlowCommon;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Implementation for boundary values for an inlet with a Mach number
    /// smaller than 1.
    /// </summary>
    public class SubsonicInlet : BoundaryCondition {

        /// <summary>
        /// A function specifying the density at the boundary.
        /// </summary>
        private readonly Func<double[], double, double> densityFunction;

        /// <summary>
        /// A function specifying the momentum at the boundary.
        /// </summary>
        private readonly Func<double[], double, double>[] velocityFunctions;

        /// <summary>
        /// Constructs a new subsonic inlet using the values defined by
        /// <paramref name="densityFunction"/> and
        /// <paramref name="velocityFunctions"/> as values at the boundary.
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="densityFunction">
        /// A function specifying the density at the boundary.
        /// </param>
        /// <param name="velocityFunctions">
        /// A function specifying the momentum at the boundary.
        /// </param>
        public SubsonicInlet(MaterialProperty.Material config, Func<double[], double, double> densityFunction, Func<double[], double, double>[] velocityFunctions)
            : base(config) {
            this.densityFunction = densityFunction;
            this.velocityFunctions = velocityFunctions;
        }

        /// <summary>
        /// Calculates the boundary values for a subsonic inlet (Mach number
        /// smaller than 1). Theoretically, we have to impose four conditions
        /// (in the inviscid as well as in the viscid case!) which is why we
        /// extrapolate the fifth value (in this case: the energy) to the
        /// boundary.
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
        /// \f$ (\rho^*, \vec{u}^*, p^-)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            double rhoOut = densityFunction(x, time);

            Vector uOut = new Vector();
            for (int i = 0; i < CNSEnvironment.NumberOfDimensions; i++) {
                uOut[i] = velocityFunctions[i](x, time);
            }

            return StateVector.FromPrimitiveQuantities(stateIn.Material, rhoOut, uOut, stateIn.Pressure);
        }
    }
}