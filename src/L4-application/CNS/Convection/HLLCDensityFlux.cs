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
using BoSSS.Solution.CompressibleFlowCommon.Boundary;

namespace CNS.Convection {

    /// <summary>
    /// Implements the parts of the HLLCFlux specific to the continuity
    /// equation
    /// </summary>
    public class HLLCDensityFlux : HLLCFlux {

        /// <summary>
        /// <see cref="HLLCFlux"/>
        /// </summary>
        /// <param name="config"><see cref="HLLCFlux"/></param>
        /// <param name="boundaryMap"><see cref="HLLCFlux"/></param>
        /// <param name="equationComponent"><see cref="HLLCFlux"/></param>
        /// <param name="speciesMap"><see cref="HLLCFlux"/></param>
        public HLLCDensityFlux(CNSControl config, IBoundaryConditionMap boundaryMap, EulerDensityComponent equationComponent, ISpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
        }

        /// <summary>
        /// See Toro2009, equation 10.73
        /// </summary>
        /// <param name="state">
        /// <see cref="HLLCFlux.GetModifiedVariableValue"/>
        /// </param>
        /// <param name="cellWaveSpeed">
        /// <see cref="HLLCFlux.GetModifiedVariableValue"/>
        /// </param>
        /// <param name="cellNormalVelocity">
        /// <see cref="HLLCFlux.GetModifiedVariableValue"/>
        /// </param>
        /// <param name="intermediateWaveSpeed">
        /// <see cref="HLLCFlux.GetModifiedVariableValue"/>
        /// </param>
        /// <param name="normal">
        /// <see cref="HLLCFlux.GetModifiedVariableValue"/>
        /// </param>
        /// <returns>See Toro2009, equation 10.73</returns>
        protected override double GetModifiedVariableValue(StateVector state, double cellWaveSpeed, double cellNormalVelocity, double intermediateWaveSpeed, ref Vector normal) {
            return state.Density
                * (cellWaveSpeed - cellNormalVelocity)
                / (cellWaveSpeed - intermediateWaveSpeed);
        }
    }
}
