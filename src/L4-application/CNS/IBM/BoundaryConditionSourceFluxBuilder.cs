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
using System.Linq;

namespace CNS.IBM {

    /// <summary>
    /// Flux builder for the boundary conditions at immersed interfaces which
    /// are implemented via source terms. Note that currently, the convective
    /// fluxes are considered only.
    /// </summary>
    public class BoundaryConditionSourceFluxBuilder : FluxBuilder {

        private Operator standardOperator;

        private BoundaryCondition boundaryCondition;

        /// <summary>
        /// Constructs a new flux builder where the boundary conditions are
        /// evaluated using the convective fluxes defined by
        /// <paramref name="convectiveBuilder"/>.
        /// </summary>
        /// <param name="control"></param>
        /// <param name="boundaryMap"></param>
        /// <param name="speciesMap"></param>
        /// <param name="convectiveBuilder"></param>
        /// <param name="diffusiveBuilder"></param>
        public BoundaryConditionSourceFluxBuilder(
            IBMControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, FluxBuilder convectiveBuilder, FluxBuilder diffusiveBuilder)
            : base(control, boundaryMap, speciesMap) {
            standardOperator = new Operator(control);

            if (convectiveBuilder != null) {
                convectiveBuilder.BuildFluxes(standardOperator);
            }

            if (diffusiveBuilder != null) {
                diffusiveBuilder.BuildFluxes(standardOperator);
            }

            string levelSetBoundaryType = control.LevelSetBoundaryTag;
            boundaryCondition = boundaryMap.GetBoundaryCondition(levelSetBoundaryType);
        }

        /// <summary>
        /// Adds an instance of <see cref="BoundaryConditionSource"/>
        /// to each component of <paramref name="boundaryOperator"/>.
        /// </summary>
        /// <param name="boundaryOperator"></param>
        public override void BuildFluxes(Operator boundaryOperator) {
            foreach (var equationComponent in standardOperator.DensityComponents) {
                boundaryOperator.DensityComponents.Add(
                    equationComponent.CreateBoundaryConditionSource(control, speciesMap, boundaryCondition));
            }

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                foreach (var equationComponent in standardOperator.MomentumComponents[d]) {
                    boundaryOperator.MomentumComponents[d].Add(
                        equationComponent.CreateBoundaryConditionSource(control, speciesMap, boundaryCondition));
                }
            }

            foreach (var equationComponent in standardOperator.EnergyComponents) {
                boundaryOperator.EnergyComponents.Add(
                        equationComponent.CreateBoundaryConditionSource(control, speciesMap, boundaryCondition));
            }
        }
    }
}
