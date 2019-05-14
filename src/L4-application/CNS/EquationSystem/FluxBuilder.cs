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

using BoSSS.Foundation;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace CNS.EquationSystem {

    /// <summary>
    /// Base class for flux builders which consume a
    /// <see cref="SpatialOperator"/> and add all numerical fluxes required
    /// by a specific flux type.
    /// </summary>
    public abstract class FluxBuilder {

        /// <summary>
        /// Configuration options
        /// </summary>
        protected readonly CNSControl control;
        
        /// <summary>
        /// The boundary mapping which is required for the construction of the
        /// specific fluxes
        /// </summary>
        protected readonly BoundaryConditionMap boundaryMap;

        /// <summary>
        /// </summary>
        protected readonly ISpeciesMap speciesMap;

        protected readonly IEquationOfState equationOfState;

        protected readonly double machNumber;

        /// <summary>
        /// Constructs a new flux builder
        /// </summary>
        /// <param name="control">Configuration options</param>
        /// <param name="boundaryMap">
        /// The boundary mapping which is required for the construction of the
        /// specific fluxes.
        /// </param>
        /// <param name="speciesMap">
        /// The species mapping which is required to determine the active
        /// equation of state upon evaluation of the fluxes.
        /// </param>
        protected FluxBuilder(CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, double machNumber, IEquationOfState equationOfState = null) {
            this.control = control;
            this.boundaryMap = boundaryMap;
            this.speciesMap = speciesMap;
            this.machNumber = machNumber;
            this.equationOfState = equationOfState;
        }

        /// <summary>
        /// Implement this method such that it creates flux functions for
        /// affected component of <paramref name="op"/>.
        /// </summary>
        /// <param name="op">
        /// The mapping to which the fluxes should be <b>added</b>.
        /// </param>
        abstract public void BuildFluxes(Operator op);
    }
}
