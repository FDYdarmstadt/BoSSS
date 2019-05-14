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
    /// Specialized <see cref="FluxBuilder"/> for the HLL flux (see
    /// <see cref="HLLFlux"/>)
    /// </summary>
    public class HLLFluxBuilder : FluxBuilder {

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
        public HLLFluxBuilder(CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, double machNumber)
            : base(control, boundaryMap, speciesMap, machNumber) {
        }

        /// <summary>
        /// Adds instances of <see cref="HLLFlux" /> to all components of the
        /// <paramref name="op" />
        /// </summary>
        /// <param name="op">
        /// The mapping to which the fluxes should be added.
        /// </param>
        public override void BuildFluxes(Operator op) {
            op.DensityComponents.Add(new HLLFlux(
                control, boundaryMap, new EulerDensityComponent(), speciesMap));

            for (int i = 0; i < CNSEnvironment.NumberOfDimensions; i++) {
                op.MomentumComponents[i].Add(new HLLFlux(
                    control, boundaryMap, new EulerMomentumComponent(i, control.EquationOfState.HeatCapacityRatio, control.MachNumber, CNSEnvironment.NumberOfDimensions), speciesMap));
            }

            op.EnergyComponents.Add(new HLLFlux(
                control, boundaryMap, new EulerEnergyComponent(), speciesMap));
        }
    }
}
