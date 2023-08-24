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
using System.Diagnostics;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;

namespace XESF.Fluxes {

    /// <summary>
    /// Optimized version of the HLLC density flux for ideal gases.
    /// </summary>
    public class HLLCDensityFlux : HLLCFlux {

        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public HLLCDensityFlux(IBoundaryConditionMap boundaryMap, Material material)
            : base(boundaryMap, material) {
        }

        public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {
            double ret = 0;
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
            //pIn = Math.Max(pIn, ilPSP.Utils.BLAS.MachineEps);
            //pOut = Math.Max(pOut, ilPSP.Utils.BLAS.MachineEps);

            double speedOfSoundIn = Math.Sqrt(pIn / densityIn) / Mach;
            double speedOfSoundOut = Math.Sqrt(pOut / densityOut) / Mach;
            //Debug.Assert(!double.IsNaN(speedOfSoundIn) || double.IsInfinity(speedOfSoundIn));
            //Debug.Assert(!double.IsNaN(speedOfSoundOut) || double.IsInfinity(speedOfSoundOut));



            double densityMean = 0.5 * (densityIn + densityOut);
            double pressureMean = 0.5 * (pIn + pOut);
            double speedOfSoundMean = 0.5 * (speedOfSoundIn + speedOfSoundOut);
            double velocityJump = normalVelocityOut - normalVelocityIn;

            // Calculate the pressure estimate at the edge and use it to
            // calculate the components of the correction factor qIn and qOut
            double pStar = pressureMean - 0.5 * MachScaling * velocityJump * speedOfSoundMean * densityMean;
            pStar = Math.Max(0.0, pStar);

            double qIn = 1.0;
            if(pStar > pIn) {
                qIn = Math.Sqrt(1.0 + 0.5 * (gamma + 1.0) * (pStar / pIn - 1.0) / gamma);
            }

            double qOut = 1.0;
            if(pStar > pOut) {
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
            if(intermediateWaveSpeed > 0.0) {
                edgeFlux = normalVelocityIn * densityIn;
                if(waveSpeedIn <= 0.0) {
                    double modifiedDensity = densityIn * (waveSpeedIn - normalVelocityIn)
                        / (waveSpeedIn - intermediateWaveSpeed);
                    edgeFlux += waveSpeedIn * (modifiedDensity - densityIn);
                }
            } else {
                edgeFlux = normalVelocityOut * densityOut;
                if(waveSpeedOut >= 0.0) {
                    double modifiedDensity = densityOut * (waveSpeedOut - normalVelocityOut)
                        / (waveSpeedOut - intermediateWaveSpeed);
                    edgeFlux += waveSpeedOut * (modifiedDensity - densityOut);
                }
            }

            //Debug.Assert(!double.IsNaN(edgeFlux) || double.IsInfinity(edgeFlux));
            ret += edgeFlux;
            return ret;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            return InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D)*(_vIN-_vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {


            Vector fluxvector = new Vector(cpv.D);
            for(int d = 0; d < cpv.D; d++) {
                // state.Momentum
                double flux = U[d + 1];
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                fluxvector[d] = flux;
            }
            return -1.0*fluxvector* GradV;

        }
    }
}
