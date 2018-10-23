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
    /// Implementation of the boundary condition for a supersonic inlet (i.e.
    /// Ma > 1.0). In this case, all values at the boundary have to be
    /// prescribed. For this matter, we choose the density, the momentum and
    /// the pressure as input.
    /// </summary>
    public class SupersonicInlet : BoundaryCondition {

        /// <summary>
        /// The prescribed density at the inlet
        /// </summary>
        public readonly Func<double[], double, double> DensityFunction;

        /// <summary>
        /// The prescribed velocity components at the inlet
        /// </summary>
        public readonly Func<double[], double, double>[] VelocityFunctions;

        /// <summary>
        /// The prescribed pressure at the inlet
        /// </summary>
        public Func<double[], double, double> PressureFunction;

        /// <summary>
        /// Constructs a new supersonic inlet using the values defined by
        /// <paramref name="densityFunction"/>,
        /// <paramref name="velocityFunctions"/> and
        /// <paramref name="pressureFunction"/> as values at the boundary.
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="densityFunction">
        /// The prescribed density at the inlet
        /// </param>
        /// <param name="velocityFunctions">
        /// The prescribed velocity components at the inlet
        /// </param>
        /// <param name="pressureFunction">
        /// The prescribed pressure at the inlet
        /// </param>
        public SupersonicInlet(MaterialProperty.Material config, Func<double[], double, double> densityFunction, Func<double[], double, double>[] velocityFunctions, Func<double[], double, double> pressureFunction)
            : base(config) {
            this.DensityFunction = densityFunction;
            this.VelocityFunctions = velocityFunctions;
            this.PressureFunction = pressureFunction;
        }

        /// <summary>
        /// Calculates the boundary values for a supersonic inlet (Mach number
        /// greater than 1). Theoretically, we have to impose five conditions
        /// (in the inviscid as well as in the viscid case!) which is why
        /// calculate the complete state from the given density, momentum and
        /// pressure
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
        /// \f$ (\rho^*, u_0^*[, u_1^*[, u_2^*]], p*)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            Vector uOut = new Vector();
            for (int i = 0; i < CNSEnvironment.NumberOfDimensions; i++) {
                uOut[i] = VelocityFunctions[i](x, time);
            }

            return StateVector.FromPrimitiveQuantities(
                stateIn.Material,
                DensityFunction(x, time),
                uOut,
                PressureFunction(x, time));
        }
    }
}
