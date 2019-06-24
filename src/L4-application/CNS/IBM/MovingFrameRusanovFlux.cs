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

using System;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;

namespace CNS.IBM {

    public class MovingFrameRusanovFlux : EulerFlux {

        private Func<double[], double, Vector> levelSetVelocity;

        public MovingFrameRusanovFlux(CompressibleControl config, IBoundaryConditionMap boundaryMap, IEulerEquationComponent equationComponent, ImmersedSpeciesMap speciesMap)
            : base(config, boundaryMap, equationComponent, speciesMap) {
            this.levelSetVelocity = speciesMap.Control.LevelSetVelocity;
        }

        protected override double InnerEdgeFlux(double[] x, double time, StateVector stateIn, StateVector stateOut, ref Vector normal, int edgeIndex) {
            // Version: Subtract -s * flux from upwind
            double uIn = stateIn.Velocity * normal;
            double uOut = stateOut.Velocity * normal;
            double correction = 0.0;
            if (edgeIndex < 0) {
                double s = levelSetVelocity(x, time) * normal;
                double relativeSpeed = 0.5 * (uIn + uOut) - s;
                if (relativeSpeed > 0.0) {
                    correction = s * equationComponent.VariableValue(stateIn);
                } else {
                    correction = s * equationComponent.VariableValue(stateOut);
                }
            }

            double waveSpeedIn = uIn + stateIn.SpeedOfSound;
            double waveSpeedOut = uOut + stateOut.SpeedOfSound;
            double penalty = Math.Max(waveSpeedIn, waveSpeedOut);

            double valueIn = equationComponent.VariableValue(stateIn);
            double fluxIn = equationComponent.Flux(stateIn) * normal;
            double valueOut = equationComponent.VariableValue(stateOut);
            double fluxOut = equationComponent.Flux(stateOut) * normal;

            return 0.5 * (fluxIn + fluxOut) - 0.5 * penalty * (valueOut - valueIn) - correction;
        }
    }
}
