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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using System;
using System.Diagnostics;

namespace XESF.Fluxes {

    /// <summary>
    /// Highly optimized version of the HLLC flux for ideal gases.
    /// </summary>
    public class HLLCMomentumFlux : HLLCFlux {

        /// <summary>
        /// <see cref="HLLCMomentumFlux.HLLCMomentumFlux"/>
        /// </summary>
        protected int component;

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        /// <param name="component">
        /// The component of the momentum vector.
        /// </param>
        public HLLCMomentumFlux(IBoundaryConditionMap boundaryMap, int component, Material material)
            : base(boundaryMap, material) {
            this.component = component;
        }


        /// <summary>
        /// <see cref="INonlinearFlux.InnerEdgeFlux"/>
        /// </summary>
        public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double MachScaling = gamma * Mach * Mach;

            double densityIn = Uin[0];
            double densityOut = Uout[0];

            double momentumSquareIn = 0.0;
            double momentumSquareOut = 0.0;
            double normalVelocityIn = 0.0;
            double normalVelocityOut = 0.0;
            for(int d = 0; d < D; d++) {
                momentumSquareIn += Uin[d + 1] * Uin[d + 1];
                momentumSquareOut += Uout[d + 1] * Uout[d + 1];
                normalVelocityIn += Uin[d + 1] * normal[d];
                normalVelocityOut += Uout[d + 1] * normal[d];
            }
            normalVelocityIn /= densityIn;
            normalVelocityOut /= densityOut;

            double energyIn = Uin[D + 1];
            double energyOut = Uout[D + 1];

            double pIn = (gamma - 1.0) * (energyIn - MachScaling * 0.5 * momentumSquareIn / densityIn);
            double pOut = (gamma - 1.0) * (energyOut - MachScaling * 0.5 * momentumSquareOut / densityOut);

            double speedOfSoundIn = Math.Sqrt(pIn / densityIn) / Mach;
            double speedOfSoundOut = Math.Sqrt(pOut / densityOut) / Mach;
            //double speedOfSoundIn = (pIn / densityIn).Pow2() / Mach;
            //double speedOfSoundOut = (pOut / densityOut).Pow2() / Mach;


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
            double speedDiff = cOut * normalVelocityOut - cIn * normalVelocityIn;
            //if (Math.Abs(speedDiff) < 1e-13) {
            //    speedDiff = 0.0;
            //}

            double pIn_minus_pOut = pIn - pOut;
            //if (Math.Abs(pIn_minus_pOut) < 1e-13) {
            //    pIn_minus_pOut = 0.0;
            //}

            double intermediateWaveSpeed =
                (speedDiff + pIn_minus_pOut / MachScaling) /
                (cOut - cIn);

            //double intermediateWaveSpeed =
            //    (cOut * normalVelocityOut - cIn * normalVelocityIn + (pIn - pOut) / MachScaling) /
            //    (cOut - cIn);

            double edgeFlux = 0.0;
            // cf. Toro2009, equation 10.71
            // flux = state.Momentum[MomentumComponent] * state.Velocity
            //   + state.Pressure * ComponentVector;
            if (intermediateWaveSpeed > 0.0) {
                edgeFlux = normalVelocityIn * Uin[component + 1]
                    + pIn / MachScaling * normal[component];
                if (waveSpeedIn <= 0.0) {
                    double modifiedMomentum = densityIn *
                        (waveSpeedIn - normalVelocityIn) /
                        (waveSpeedIn - intermediateWaveSpeed) *
                        (Uin[component + 1] / densityIn +
                            (intermediateWaveSpeed - normalVelocityIn) * normal[ component]);
                    edgeFlux += waveSpeedIn * (modifiedMomentum - Uin[component + 1]);
                }
            } else {
                edgeFlux = normalVelocityOut * Uout[component + 1]
                    + pOut / MachScaling * normal[component];
                if (waveSpeedOut >= 0.0) {
                    double modifiedMomentum = densityOut *
                        (waveSpeedOut - normalVelocityOut) /
                        (waveSpeedOut - intermediateWaveSpeed) *
                        (Uout[component + 1] / densityOut +
                            (intermediateWaveSpeed - normalVelocityOut) * normal[component]);
                    edgeFlux += waveSpeedOut * (modifiedMomentum - Uout[component + 1]);
                }
            }

            //Debug.Assert(!double.IsNaN(edgeFlux) || double.IsInfinity(edgeFlux));
            return edgeFlux;
        
        }

        /// <summary>
        /// <see cref="INonlinearFlux.Flux"/>
        /// </summary>
        public Vector Flux(double[] U, int D) {
            Vector Output = new Vector(D);
            int L = U.Length;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;


            double density = U[0];
            double energy = U[D + 1];

            //state.Momentum[MomentumComponent] * state.Velocity
            //    + 1/(gamma * Mach^2) * state.Pressure * ComponentVector;
            double momentumSquared = 0.0;
            for (int d = 0; d < D; d++) {
                Output[d] += U[component + 1] * U[d + 1] / density;
                momentumSquared += U[d + 1] * U[d + 1];
            }

            double flux = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density) / gammaMachSquared;
            //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
            Output[component] += flux;
            return Output;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D) * (_vIN - _vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector fluxvector = Flux(U, cpv.D);
            return -1.0 * fluxvector * GradV;

        }
    }

}
