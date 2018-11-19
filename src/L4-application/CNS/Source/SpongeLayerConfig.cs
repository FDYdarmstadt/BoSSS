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

using BoSSS.Solution.CompressibleFlowCommon;
using System;

namespace CNS.Source {

    /// <summary>
    /// Configuration options for a sponge layer source term which serves as a
    /// non-reflecting boundary condition (see <see cref="SpongeLayerSource"/>)
    /// </summary>
    [Serializable]
    public class SpongeLayerConfig {

        /// <summary>
        /// Reference state at the boundary of the domain, i.e. deviations from
        /// this state will be damped
        /// </summary>
        public StateVector referenceState = default(StateVector);

        /// <summary>
        /// Coordinate direction of the damping function (0 for x, 1 for y,
        /// 2 for z)
        /// </summary>
        public int layerDirection = -1;

        /// <summary>
        /// x-/y-/z-coordinate of the start of the sponge layer
        /// </summary>
        public double layerStart = double.NaN;

        /// <summary>
        /// x-/y-/z-coordinate of the end of the sponge layer
        /// </summary>
        public double layerEnd = double.NaN;

        /// <summary>
        /// Strength of the damping
        /// </summary>
        public double strength = double.NaN;
    }
}
