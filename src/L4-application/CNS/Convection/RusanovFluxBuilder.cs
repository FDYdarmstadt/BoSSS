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
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.EquationSystem;

namespace CNS.Convection {

    /// <summary>
    /// Flux builder for the <see cref="RusanovFlux"/>.
    /// </summary>
    public class RusanovFluxBuilder : FluxBuilder {

        /// <summary>
        /// <see cref="FluxBuilder"/>
        /// </summary>
        /// <param name="control"><see cref="FluxBuilder"/>
        /// <see cref="FluxBuilder"/>
        /// </param>
        /// <param name="boundaryMap"><see cref="FluxBuilder"/>
        /// <see cref="FluxBuilder"/>
        /// </param>
        /// <param name="speciesMap">
        /// <see cref="FluxBuilder"/>
        /// </param>
        public RusanovFluxBuilder(CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        /// <summary>
        /// Adds and instance of <see cref="RusanovFlux" /> to all components of
        /// the given <paramref name="op" />.
        /// </summary>
        /// <param name="op">
        /// The mapping to which the fluxes should be added.
        /// </param>
        public override void BuildFluxes(Operator op) {
            op.DensityComponents.Add(new RusanovFlux(
                control, boundaryMap, new EulerDensityComponent(), speciesMap));

            for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                op.MomentumComponents[d].Add(new RusanovFlux(
                    control, boundaryMap, new EulerMomentumComponent(d, control.EquationOfState.HeatCapacityRatio, control.MachNumber, CompressibleEnvironment.NumberOfDimensions), speciesMap));
            }

            op.EnergyComponents.Add(new RusanovFlux(
                control, boundaryMap, new EulerEnergyComponent(), speciesMap));
        }
    }
}
