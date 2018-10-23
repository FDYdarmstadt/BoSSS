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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.EquationSystem;

namespace CNS.Source {

    /// <summary>
    /// Provider for the optional source terms.
    /// </summary>
    public class CustomSourceBuilder : FluxBuilder {

        /// <summary>
        /// Builds a new source provider.
        /// </summary>
        /// <param name="control">
        /// Configuration options.
        /// </param>
        /// <param name="boundaryMap"></param>
        /// <param name="speciesMap"></param>
        public CustomSourceBuilder(CNSControl control, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        /// <summary>
        /// For each equation component: Adds custom sources (see
        /// <see cref="CNSControl"/>) to <paramref name="op"/> if they have
        /// been configured.
        /// </summary>
        /// <param name="op">
        /// <see cref="FluxBuilder.BuildFluxes"/>
        /// </param>
        public override void BuildFluxes(Operator op) {
            string[] argumentOrdering = CNSEnvironment.PrimalArgumentOrdering;

            foreach (var sourceFactory in control.CustomContinuitySources) {
                op.DensityComponents.Add(sourceFactory(speciesMap));
            }

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                foreach (var sourceFactory in control.CustomMomentumSources[d]) {
                    op.MomentumComponents[d].Add(sourceFactory(speciesMap));
                }
            }

            foreach (var sourceFactory in control.CustomEnergySources) {
                op.EnergyComponents.Add(sourceFactory(speciesMap));
            }
        }
    }
}