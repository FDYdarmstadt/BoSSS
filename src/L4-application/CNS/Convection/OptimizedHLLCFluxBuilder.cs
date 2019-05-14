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
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using CNS.EquationSystem;

namespace CNS.Convection {

    /// <summary>
    /// Builder for the optimized variant of the HLLC flux
    /// </summary>
    public class OptimizedHLLCFluxBuilder : FluxBuilder {

        /// <summary>
        /// <see cref="FluxBuilder.FluxBuilder"/>
        /// </summary>
        /// <param name="control"></param>
        /// <param name="boundaryMap"></param>
        /// <param name="speciesMap"></param>
        public OptimizedHLLCFluxBuilder(CNSControl control, BoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IEquationOfState equationOfState, double machNumber)
            : base(control, boundaryMap, speciesMap, machNumber, equationOfState) {
        }

        /// <summary>
        /// Maps appropriate instances of
        /// <see cref="OptimizedHLLCDensityFlux"/>,
        /// <see cref="OptimizedHLLCMomentumFlux"/> and
        /// <see cref="OptimizedHLLCEnergyFlux "/> to each component of the
        /// given <paramref name="op"/>.
        /// </summary>
        /// <param name="op">
        /// <see cref="FluxBuilder.BuildFluxes"/>
        /// </param>
        public override void BuildFluxes(Operator op) {
            op.DensityComponents.Add(new OptimizedHLLCDensityFlux(speciesMap, boundaryMap, equationOfState, machNumber));

            for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                op.MomentumComponents[d].Add(new OptimizedHLLCMomentumFlux(speciesMap, boundaryMap, d, equationOfState, machNumber));
            }

            op.EnergyComponents.Add(new OptimizedHLLCEnergyFlux(speciesMap, boundaryMap, equationOfState, machNumber));
        }
    }
}
