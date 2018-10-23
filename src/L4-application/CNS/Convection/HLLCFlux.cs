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
using CNS.Boundary;
using CNS.MaterialProperty;
using System;

namespace CNS.Convection {

    /// <summary>
    /// Abstract base class for the components of the HLLC flux as described in
    /// Toro2009 (p. 331f)
    /// </summary>
    abstract public class HLLCFlux : EulerFlux {

        /// <summary>
        /// <see cref="EulerFlux"/>
        /// </summary>
        /// <param name="config"><see cref="EulerFlux"/></param>
        /// <param name="boundaryMap"><see cref="EulerFlux"/></param>
        /// <param name="equationComponent"><see cref="EulerFlux"/></param>
        /// <param name="speciesMap"><see cref="EulerFlux"/></param>
        protected HLLCFlux(CNSControl config, IBoundaryConditionMap boundaryMap, IEulerEquationComponent equationComponent, ISpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
            if (config.EquationOfState is IdealGas == false) {
                throw new Exception("HLLC flux currently only works for ideal gases");
            }
        }

        /// <summary>
        /// Evaluates the HLLC flux according to Toro2009
        /// </summary>
        /// <param name="x">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <param name="time">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <param name="stateIn">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <param name="stateOut">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <param name="normal">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <param name="edgeIndex">
        /// <see cref="InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>
        /// </param>
        /// <returns>see Toro2009, equation 10.71</returns>
        protected internal override double InnerEdgeFlux(double[] x, double time, StateVector stateIn, StateVector stateOut, ref Vector normal, int edgeIndex) {
            double waveSpeedIn;
            double waveSpeedOut;
            EstimateWaveSpeeds(stateIn, stateOut, ref normal, out waveSpeedIn, out waveSpeedOut);

            double MachScaling = config.EquationOfState.HeatCapacityRatio * config.MachNumber * config.MachNumber;

            double normalVelocityIn = stateIn.Velocity * normal;
            double normalVelocityOut = stateOut.Velocity * normal;

            double cIn = stateIn.Density * (waveSpeedIn - normalVelocityIn);
            double cOut = stateOut.Density * (waveSpeedOut - normalVelocityOut);

            // cf. Toro2009, equation 10.70
            // corrected according to dimensionless equations 
            double intermediateWaveSpeed = (cOut * normalVelocityOut - cIn * normalVelocityIn
                + (stateIn.Pressure - stateOut.Pressure) / MachScaling) / (cOut - cIn);

            double edgeFlux = 0.0;

            // cf. Toro2009, equation 10.71
            if (waveSpeedIn > 0.0) {
                edgeFlux = equationComponent.Flux(stateIn) * normal;
            } else if (waveSpeedIn <= 0.0 && 0.0 < intermediateWaveSpeed) {
                edgeFlux = equationComponent.Flux(stateIn) * normal + waveSpeedIn * (
                    GetModifiedVariableValue(stateIn, waveSpeedIn, normalVelocityIn, intermediateWaveSpeed, ref normal)
                    - equationComponent.VariableValue(stateIn));
            } else if (intermediateWaveSpeed <= 0.0 && 0.0 <= waveSpeedOut) {
                edgeFlux = equationComponent.Flux(stateOut) * normal + waveSpeedOut * (
                    GetModifiedVariableValue(stateOut, waveSpeedOut, normalVelocityOut, intermediateWaveSpeed, ref normal)
                    - equationComponent.VariableValue(stateOut));
            } else if (waveSpeedOut < 0.0) {
                edgeFlux = equationComponent.Flux(stateOut) * normal;
            } else {
                throw new Exception("Inconsistency in HLLC flux detected. Might the flow state be invalid?");
            }

            return edgeFlux;
        }

        /// <summary>
        /// Implement this method for the different components of the Euler
        /// system by calculating the specific correction factor associated
        /// with the so called star region in the HLLC flux (see Toro2009)
        /// </summary>
        /// <param name="state">The flow state inside a cell</param>
        /// <param name="cellWaveSpeed">
        /// The fastest signal velocity in a cell
        /// </param>
        /// <param name="cellNormalVelocity">
        /// The velocity normal to the edge
        /// </param>
        /// <param name="intermediateWaveSpeed">
        /// An estimation for the signal velocity in the so called star region
        /// </param>
        /// <param name="normal">A unit vector normal to the edge</param>
        /// <returns>
        /// An estimate for the variable value in the star region
        /// </returns>
        abstract protected double GetModifiedVariableValue(StateVector state, double cellWaveSpeed, double cellNormalVelocity, double intermediateWaveSpeed, ref Vector normal);
    }
}
