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

using System.Linq;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.EquationSystem;

namespace CNS.Source {

    /// <summary>
    /// Builder for the components of <see cref="Operators.Gravity"/>.
    /// </summary>
    public class GravityFluxBuilder : FluxBuilder {

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="control"></param>
        /// <param name="boundaryMap"></param>
        /// <param name="speciesMap"></param>
        public GravityFluxBuilder(CNSControl control, CompressibleBoundaryCondMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        /// <summary>
        /// Adds instances of <see cref="GravityMomentumSource"/> and
        /// <see cref="GravityEnergySource"/> to the corresponding components.
        /// Here, gravity is assumed to act into -y-direction in 2D and
        /// -z-direction in 3D
        /// </summary>
        /// <param name="mapping"></param>
        public override void BuildFluxes(Operator mapping) {
            mapping.MomentumComponents.Last().Add(new GravityMomentumSource(control));
            mapping.EnergyComponents.Add(new GravityEnergySource(control));
        }
    }
}
