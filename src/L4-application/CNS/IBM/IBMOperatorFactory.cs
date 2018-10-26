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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon;
using CNS.Boundary;
using CNS.EquationSystem;
using System.Collections.Generic;

namespace CNS.IBM {

    /// <summary>
    /// Specialization of <see cref="OperatorFactory"/> for immersed boundary
    /// flows.
    /// </summary>
    public class IBMOperatorFactory : OperatorFactory {
        
        /// <summary>
        /// Builder for all terms only active at the immersed boundary
        /// </summary>
        protected readonly IList<FluxBuilder> immersedBoundaryFluxBuilders = new List<FluxBuilder>();

        /// <summary>
        /// Constructs a new operator factory which additionally implements
        /// the boundary conditions at immersed boundaries.
        /// </summary>
        /// <param name="control"></param>
        /// <param name="gridData"></param>
        /// <param name="workingSet"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryMap"></param>
        public IBMOperatorFactory(
            IBMControl control,
            IGridData gridData,
            CNSFieldSet workingSet,
            ISpeciesMap speciesMap,
            IBoundaryConditionMap boundaryMap)
            : base(control, gridData, workingSet, speciesMap, boundaryMap) {

            this.immersedBoundaryFluxBuilders.Add(new BoundaryConditionSourceFluxBuilder(
                control, boundaryMap, speciesMap, convectiveFluxBuilder, diffusiveFluxBuilder));
        }

        /// <summary>
        /// Returns an operator containing all terms that enforce boundary
        /// conditions at immersed boundaries
        /// </summary>
        /// <returns></returns>
        public Operator GetImmersedBoundaryOperator() {
            Operator op = new Operator(control);
            foreach (FluxBuilder builder in immersedBoundaryFluxBuilders) {
                builder.BuildFluxes(op);
            }

            return op;
        }
    }
}
