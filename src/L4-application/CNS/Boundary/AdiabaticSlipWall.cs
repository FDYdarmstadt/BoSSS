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
    /// Implementation of the boundary value on an isolating wall with a slip
    /// condition for the velocity (e.g. for an inviscid fluid). Can also be
    /// used as a symmetry plane boundary condition in arbitrary flows
    /// </summary>
    public class AdiabaticSlipWall : BoundaryCondition {

        /// <summary>
        /// Optional velocities for a moving wall BC
        /// </summary>
        public readonly Func<double[], double, double>[] WallVelocities;

        /// <summary>
        /// <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="wallVelocities"></param>
        public AdiabaticSlipWall(CNSControl config, Func<double[], double, double>[] wallVelocities = null)
            : base(config) {
            this.WallVelocities = wallVelocities;
        }

        /// <summary>
        /// Calculates the boundary value for an isolating slip wall according
        /// to the first variant described in VegtVen2002. The basic idea is to
        /// impose \f$ \vec{u} \cdot \vec{n} = 0\f$  (slip wall) by
        /// setting the edge momentum to
        /// \f$ \vec{m}^+ = \vec{m}^- - 2((\rho \vec{u})^- \cdot \vec{n})\vec{n}\f$ 
        /// which mimics a mirrored flow at the other side of the wall.
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
        /// <see cref="BoundaryCondition.GetBoundaryState"/>
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            Vector3D normalVector = GetNormalVector(normal);

            StateVector stateOut;
            if (WallVelocities == null) {
                // VegtVen2002, page 14, second equation
                Vector3D mOut = stateIn.Momentum - 2.0 * (stateIn.Momentum * normalVector) * normalVector;
                stateOut = new StateVector(stateIn.Material, stateIn.Density, mOut, stateIn.Energy);
            } else {
                Vector3D uWall = new Vector3D();
                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    uWall[d] = WallVelocities[d](x, time);
                }

                Vector3D uOut = stateIn.Velocity
                    - 2.0 * ((stateIn.Velocity - uWall) * normalVector) * normalVector;
                stateOut = StateVector.FromPrimitiveQuantities(
                    stateIn.Material, stateIn.Density, uOut, stateIn.Pressure);
            }

            return stateOut;
        }
    }
}
