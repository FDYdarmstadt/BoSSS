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

namespace CNS.Boundary {

    /// <summary>
    /// Implementation of boundary condition for a supersonic outlet (i.e. an
    /// outlet with a local Mach number greater than 1.0)
    /// </summary>
    public class SupersonicOutlet : BoundaryCondition {

        /// <summary>
        /// <see cref="BoundaryCondition"/>
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        public SupersonicOutlet(CNSControl config)
            : base(config) {
        }

        /// <summary>
        /// At a supersonic outlet, we do not need to (and must not) impose
        /// any boundary values since all characteristics travel out of the
        /// domain and must thus be calculated
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
        /// \f$ (\rho^-, m_0^-[, m_1^-[, m_2^-]], (\rho E)^-)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, double[] x, double[] normal, StateVector stateIn) {
            return stateIn;
        }
    }
}
