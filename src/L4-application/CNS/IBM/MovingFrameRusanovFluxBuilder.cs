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
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.Convection;
using CNS.EquationSystem;

namespace CNS.IBM {

    /// <summary>
    /// Flux builder for the <see cref="MovingFrameRusanovFlux"/>.
    /// </summary>
    public class MovingFrameRusanovFluxBuilder : FluxBuilder {

        private ImmersedSpeciesMap ibmSpeciesMap;

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
        public MovingFrameRusanovFluxBuilder(CNSControl control, CompressibleBoundaryCondMap boundaryMap, ISpeciesMap speciesMap)
            : base(control, boundaryMap, speciesMap) {
            this.ibmSpeciesMap = speciesMap as ImmersedSpeciesMap;
            if (ibmSpeciesMap == null) {
                throw new System.Exception();
            }
        }

        /// <summary>
        /// Adds and instance of <see cref="MovingFrameRusanovFlux" /> to all components of
        /// the given <paramref name="op" />.
        /// </summary>
        /// <param name="op">
        /// The mapping to which the fluxes should be added.
        /// </param>
        public override void BuildFluxes(Operator op) {
            op.DensityComponents.Add(new MovingFrameRusanovFlux(
                control, boundaryMap, new EulerDensityComponent(), ibmSpeciesMap));

            for (int d = 0; d < CompressibleEnvironment.NumberOfDimensions; d++) {
                op.MomentumComponents[d].Add(new MovingFrameRusanovFlux(
                    control, boundaryMap, new EulerMomentumComponent(d, control.EquationOfState.HeatCapacityRatio, control.MachNumber, CompressibleEnvironment.NumberOfDimensions), ibmSpeciesMap));
            }

            op.EnergyComponents.Add(new MovingFrameRusanovFlux(
                control, boundaryMap, new EulerEnergyComponent(), ibmSpeciesMap));
        }
    }
}
