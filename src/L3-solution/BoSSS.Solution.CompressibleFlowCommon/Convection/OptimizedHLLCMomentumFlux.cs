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
using ilPSP;
using System;

namespace BoSSS.Solution.CompressibleFlowCommon.Convection {

    /// <summary>
    /// Highly optimized version of the HLLC flux for ideal gases.
    /// </summary>
    public class OptimizedHLLCMomentumFlux : OptimizedHLLCFlux {

        /// <summary>
        /// <see cref="OptimizedHLLCMomentumFlux.OptimizedHLLCMomentumFlux"/>
        /// </summary>
        private int component;

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        /// <param name="component">
        /// The component of the momentum vector.
        /// </param>
        public OptimizedHLLCMomentumFlux(IBoundaryConditionMap boundaryMap, int component, Material material)
            : base(boundaryMap, material) {
            this.component = component;
        }

        /// <summary>
        /// <see cref="INonlinearFlux.InnerEdgeFlux"/>
        /// </summary>
        /// <param name="time"></param>
        /// <param name="jEdge"></param>
        /// <param name="x"></param>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <param name="Offset"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        public override void InnerEdgeFlux(
            double time,
            int jEdge,
            MultidimensionalArray x,
            MultidimensionalArray normal,
            MultidimensionalArray[] Uin,
            MultidimensionalArray[] Uout,
            int Offset,
            int Lenght,
            MultidimensionalArray Output) {

            int NoOfNodes = Uin[0].GetLength(1);
            int D = CompressibleEnvironment.NumberOfDimensions;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double MachScaling = gamma * Mach * Mach;

            for (int e = 0; e < Lenght; e++) {
                for (int n = 0; n < NoOfNodes; n++) {
                    double densityIn = Uin[0][e + Offset, n];
                    double densityOut = Uout[0][e + Offset, n];

                    double momentumSquareIn = 0.0;
                    double momentumSquareOut = 0.0;
                    double normalVelocityIn = 0.0;
                    double normalVelocityOut = 0.0;
                    for (int d = 0; d < D; d++) {
                        momentumSquareIn += Uin[d + 1][e + Offset, n] * Uin[d + 1][e + Offset, n];
                        momentumSquareOut += Uout[d + 1][e + Offset, n] * Uout[d + 1][e + Offset, n];
                        normalVelocityIn += Uin[d + 1][e + Offset, n] * normal[e + Offset, n, d];
                        normalVelocityOut += Uout[d + 1][e + Offset, n] * normal[e + Offset, n, d];
                    }
                    normalVelocityIn /= densityIn;
                    normalVelocityOut /= densityOut;

                    double energyIn = Uin[D + 1][e + Offset, n];
                    double energyOut = Uout[D + 1][e + Offset, n];

                    double pIn = (gamma - 1.0) * (energyIn - MachScaling * 0.5 * momentumSquareIn / densityIn);
                    double pOut = (gamma - 1.0) * (energyOut - MachScaling * 0.5 * momentumSquareOut / densityOut);

                    double speedOfSoundIn = Math.Sqrt(pIn / densityIn) / Mach;
                    double speedOfSoundOut = Math.Sqrt(pOut / densityOut) / Mach;

                    double densityMean = 0.5 * (densityIn + densityOut);
                    double pressureMean = 0.5 * (pIn + pOut);
                    double speedOfSoundMean = 0.5 * (speedOfSoundIn + speedOfSoundOut);
                    double velocityJump = normalVelocityOut - normalVelocityIn;

                    // Calculate the pressure estimate at the edge and use it to
                    // calculate the components of the correction factor qIn and qOut
                    double pStar = pressureMean - 0.5 * MachScaling * velocityJump * speedOfSoundMean * densityMean;
                    pStar = Math.Max(0.0, pStar);

                    double qIn = 1.0;
                    if (pStar > pIn) {
                        qIn = Math.Sqrt(1.0 + 0.5 * (gamma + 1.0) * (pStar / pIn - 1.0) / gamma);
                    }

                    double qOut = 1.0;
                    if (pStar > pOut) {
                        qOut = Math.Sqrt(1.0 + 0.5 * (gamma + 1.0) * (pStar / pOut - 1.0) / gamma);
                    }

                    // Determine the wave speeds
                    double waveSpeedIn = normalVelocityIn - speedOfSoundIn * qIn;
                    double waveSpeedOut = normalVelocityOut + speedOfSoundOut * qOut;

                    double cIn = densityIn * (waveSpeedIn - normalVelocityIn);
                    double cOut = densityOut * (waveSpeedOut - normalVelocityOut);

                    // cf. Toro2009, equation 10.70
                    double intermediateWaveSpeed =
                        (cOut * normalVelocityOut - cIn * normalVelocityIn + (pIn - pOut) / MachScaling) /
                        (cOut - cIn);

                    double edgeFlux = 0.0;
                    // cf. Toro2009, equation 10.71
                    // flux = state.Momentum[MomentumComponent] * state.Velocity
                    //   + state.Pressure * ComponentVector;
                    if (intermediateWaveSpeed > 0.0) {
                        edgeFlux = normalVelocityIn * Uin[component + 1][e + Offset, n]
                            + pIn/MachScaling * normal[e + Offset, n, component];
                        if (waveSpeedIn <= 0.0) {
                            double modifiedMomentum = densityIn *
                                (waveSpeedIn - normalVelocityIn) /
                                (waveSpeedIn - intermediateWaveSpeed) *
                                (Uin[component + 1][e + Offset, n] / densityIn +
                                    (intermediateWaveSpeed - normalVelocityIn) * normal[e + Offset, n, component]);
                            edgeFlux += waveSpeedIn * (modifiedMomentum - Uin[component + 1][e + Offset, n]);
                        }
                    } else {
                        edgeFlux = normalVelocityOut * Uout[component + 1][e + Offset, n]
                            + pOut/MachScaling * normal[e + Offset, n, component];
                        if (waveSpeedOut >= 0.0) {
                            double modifiedMomentum = densityOut *
                                (waveSpeedOut - normalVelocityOut) /
                                (waveSpeedOut - intermediateWaveSpeed) *
                                (Uout[component + 1][e + Offset, n] / densityOut +
                                    (intermediateWaveSpeed - normalVelocityOut) * normal[e + Offset, n, component]);
                            edgeFlux += waveSpeedOut * (modifiedMomentum - Uout[component + 1][e + Offset, n]);
                        }
                    }

                    Output[e + Offset, n] += edgeFlux;
                }
            }
        }

        /// <summary>
        /// <see cref="INonlinearFlux.Flux"/>
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="Offset"></param>
        /// <param name="Length"></param>
        /// <param name="Output"></param>
        public override void Flux(
            double time,
            MultidimensionalArray x,
            MultidimensionalArray[] U,
            int Offset,
            int Length,
            MultidimensionalArray Output) {
            int D = CompressibleEnvironment.NumberOfDimensions;
            int NoOfNodes = Output.GetLength(1);
            int L = U.Length;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;

            for (int e = 0; e < Length; e++) {
                for (int n = 0; n < NoOfNodes; n++) {
                    double density = U[0][e + Offset, n];
                    double energy = U[D + 1][e + Offset, n];

                    //state.Momentum[MomentumComponent] * state.Velocity
                    //    + 1/(gamma * Mach^2) * state.Pressure * ComponentVector;
                    double momentumSquared = 0.0;
                    for (int d = 0; d < D; d++) {
                        Output[e + Offset, n, d] += U[component + 1][e + Offset, n] *
                            U[d + 1][e + Offset, n] / density;
                        momentumSquared += U[d + 1][e + Offset, n] * U[d + 1][e + Offset, n];
                    }
                    Output[e + Offset, n, component] +=
                        ((gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density)) / gammaMachSquared;
                }
            }
        }
    }
}
