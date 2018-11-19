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
using BoSSS.Solution.CompressibleFlowCommon;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Implementation of the boundary condition for an isolating (i.e.
    /// adiabatic) wall with a no-slip condition for the velocity (i.e. the
    /// momentum at the wall is zero)
    /// </summary>
    public class AdiabaticWall : BoundaryCondition {

        /// <summary>
        /// <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        public AdiabaticWall(MaterialProperty.Material config)
            : base(config) {
        }

        /// <summary>
        /// On a no-slip wall, all velocity (and thus momentum) components
        /// vanish. The density can be extrapolated and we can impose a zero
        /// heat flux by setting \f$ T^* = T^-\f$  which is equivalent
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
        /// \f$ (\rho^-, 0[, 0[, 0]], \rho^- e^-)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            // Momentum is 0 at a no-slip boundary, thus kinetic energy is 0
            return new StateVector(
                stateIn.Material,
                stateIn.Density,
                new Vector(),
                stateIn.Density * stateIn.SpecificInnerEnergy);
        }
    }
}
