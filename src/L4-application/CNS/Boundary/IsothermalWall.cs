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
using CNS.Exception;

namespace CNS.Boundary {

    /// <summary>
    /// Implementation of the boundary condition for a no-slip wall with a
    /// fixed wall temperature.
    /// </summary>
    public class IsothermalWall : BoundaryCondition {

        /// <summary>
        /// The predefined temperature
        /// </summary>
        public readonly Func<double[], double, double> TemperatureFunction;

        /// <summary>
        /// Optional velocities for a moving wall BC
        /// </summary>
        public readonly Func<double[], double, double>[] WallVelocities;

        /// <summary>
        /// <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="temperatureFunction">
        /// The predefined temperature
        /// </param>
        /// <param name="wallVelocities">
        /// Optional wall velocity (individual formula for each direction) 
        /// </param>
        public IsothermalWall(CNSControl config, Func<double[], double, double> temperatureFunction, Func<double[], double, double>[] wallVelocities = null)
            : base(config) {
            this.TemperatureFunction = temperatureFunction;
            this.WallVelocities = wallVelocities;

            if (wallVelocities != null) {
                if (wallVelocities.Length != CNSEnvironment.NumberOfDimensions) {
                    throw new ConfigurationException();
                }

                for (int d = 0; d < wallVelocities.Length; d++) {
                    if (wallVelocities[d] == null) {
                        throw new ConfigurationException();
                    }
                }
            }
        }

        /// <summary>
        /// On a no-slip wall, all velocity (and thus momentum) components
        /// vanish. The density can be extrapolated and we can impose a fixed
        /// temperature by setting \f$ T^* = T^-\f$  which is equivalent
        /// to \f$ e^* = e^-\f$ 
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
        /// \f$ (\rho^-, 0[, 0[, 0]], \rho^- e^*)^T\f$  where
        /// \f$ e* = (\gamma - 1.0) T^*\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            double gamma = config.EquationOfState.HeatCapacityRatio;
            double innerEnergy = TemperatureFunction(x, time) / (gamma - 1.0);
            double MachScaling = gamma * config.MachNumber * config.MachNumber;

            if (WallVelocities == null) {
                // Kinetic energy is zero at a no-slip boundary, we can omit it
                return new StateVector(
                    stateIn.Material,
                    stateIn.Density,
                    new Vector3D(),
                    stateIn.Density * innerEnergy);
            } else {
                Vector3D velocity = new Vector3D();
                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    velocity[d] = WallVelocities[d](x, time);
                }

#if DEBUG
                Vector3D n = new Vector3D(normal);
                if (Math.Abs(velocity * n) > 1e-10) {
                    throw new ConfigurationException(
                        "Wall velocity must be tangent to the wall");
                }
#endif

                // Kinetic energy is zero is solely determined by wall velocity
                return new StateVector(
                    stateIn.Material,
                    stateIn.Density,
                    stateIn.Density * velocity,
                    stateIn.Density * (innerEnergy + 0.5 * MachScaling * velocity * velocity));
            }
        }
    }
}
