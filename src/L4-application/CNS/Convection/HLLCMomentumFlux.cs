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
using CNS.Boundary;

namespace CNS.Convection {

    /// <summary>
    /// Implements the parts of the HLLCFlux specific to the momentum
    /// equation
    /// </summary>
    public class HLLCMomentumFlux : HLLCFlux {

        /// <summary>
        /// The index of the momentum component represented by a specific
        /// instance of this class
        /// </summary>
        private int momentumComponent;

        /// <summary>
        /// <see cref="HLLCFlux"/>
        /// </summary>
        /// <param name="config"><see cref="HLLCFlux"/></param>
        /// <param name="boundaryMap"><see cref="HLLCFlux"/></param>
        /// <param name="equationComponent"><see cref="HLLCFlux"/></param>
        /// <param name="speciesMap"><see cref="HLLCFlux"/></param>
        public HLLCMomentumFlux(CNSControl config, IBoundaryConditionMap boundaryMap, EulerMomentumComponent equationComponent, ISpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
            this.momentumComponent = equationComponent.MomentumComponent;
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
        /// <remarks>
        /// Always remember: DON'T multiply with intermediateWaveSpeed here
        /// like Toro does! He is working in a rotated coordinate system! In
        /// our setup, this makes result dependent on the direction of the
        /// normal vector. Correct version follows from BattenEtAl1997, 
        /// equation 37
        /// </remarks>
        protected override double GetModifiedVariableValue(StateVector state, double cellWaveSpeed, double cellNormalVelocity, double intermediateWaveSpeed, ref Vector3D normal) {
            return state.Density
                * (cellWaveSpeed - cellNormalVelocity)
                / (cellWaveSpeed - intermediateWaveSpeed)
                * (state.Velocity[momentumComponent]
                    + (intermediateWaveSpeed - cellNormalVelocity) * normal[momentumComponent]);
        }
    }
}
