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

using CNS.Boundary;
using CNS.EquationSystem;

namespace CNS.Convection {

    /// <summary>
    /// Builder for the fluxes based on the exact Riemann solver define in
    /// <see cref="GodunovFlux"/>.
    /// </summary>
    public class GodunovFluxBuilder : FluxBuilder {

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
        public GodunovFluxBuilder(CNSControl control, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
        }

        /// <summary>
        /// Adds appropriate instances of <see cref="GodunovFlux"/> to
        /// each equation component of <paramref name="op" />.
        /// </summary>
        /// <param name="op">
        /// The mapping to which the fluxes should be added.
        /// </param>
        public override void BuildFluxes(Operator op) {
            op.DensityComponents.Add(new GodunovFlux(
                control, boundaryMap, new EulerDensityComponent(), speciesMap));

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                op.MomentumComponents[d].Add(new GodunovFlux(
                    control,
                    boundaryMap,
                    new EulerMomentumComponent(d, control.EquationOfState.HeatCapacityRatio, control.MachNumber),
                    speciesMap));
            }

            op.EnergyComponents.Add(new GodunovFlux(
                control, boundaryMap, new EulerEnergyComponent(), speciesMap));
        }
    }
}
