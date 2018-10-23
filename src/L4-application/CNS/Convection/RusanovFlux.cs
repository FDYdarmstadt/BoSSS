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
using System;

namespace CNS.Convection {

    /// <summary>
    /// Implementation of the Rusanov flux (also known as the local
    /// Lax-Friedrichs flux) according to Toro2009 (p. 329).
    /// </summary>
    public class RusanovFlux : EulerFlux {

        /// <summary>
        /// <see cref="EulerFlux"/>
        /// </summary>
        /// <param name="config"><see cref="EulerFlux"/></param>
        /// <param name="boundaryMap"><see cref="EulerFlux"/></param>
        /// <param name="equationComponent"><see cref="EulerFlux"/></param>
        /// <param name="speciesMap"><see cref="EulerFlux"/></param>
        public RusanovFlux(CNSControl config, IBoundaryConditionMap boundaryMap, IEulerEquationComponent equationComponent, ISpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
        }

        /// <summary>
        /// Evaluates the Rusanov flux (also known as the local Lax-Friedrichs
        /// flux as stated in Toro2009, equations 10.55 and 10.56.
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
        /// <returns>
        /// \f$ 
        /// \frac{1}{2} (F_L + F_R - S^+ (U_R - U_L))
        /// \f$ 
        /// where
        /// \f$ 
        /// S^+ = \max \{|u_L| + a_L, |u_r| + a_R\}
        /// \f$ 
        /// </returns>
        protected internal override double InnerEdgeFlux(double[] x, double time, StateVector stateIn, StateVector stateOut, ref Vector normal, int edgeIndex) {
            double waveSpeedIn = Math.Abs(stateIn.Velocity * normal) + stateIn.SpeedOfSound;
            double waveSpeedOut = Math.Abs(stateOut.Velocity * normal) + stateOut.SpeedOfSound;
            double penalty = Math.Max(waveSpeedIn, waveSpeedOut);

            double valueIn = equationComponent.VariableValue(stateIn);
            double fluxIn = equationComponent.Flux(stateIn) * normal;
            double valueOut = equationComponent.VariableValue(stateOut);
            double fluxOut = equationComponent.Flux(stateOut) * normal;
            return 0.5 * (fluxIn + fluxOut - penalty * (valueOut - valueIn));
        }
    }
}
