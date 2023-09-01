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
using XESF.Fluxes;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESF.Fluxes {

    public class HLLCEnergyFlux_XDG : HLLCEnergyFlux, ISpeciesFilter {

        private SpeciesId speciesId;

        public HLLCEnergyFlux_XDG(IBoundaryConditionMap boundaryMap, Material material, string spcName, SpeciesId speciesId) : base(boundaryMap, material) {
            this.speciesId = speciesId;
            this.speciesName = spcName;
            this.ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }

        /// <summary>
        /// Duplicated code in
        /// <see cref="RoeDensityFlux.BorderEdgeFlux(double, int, MultidimensionalArray, MultidimensionalArray, bool, byte[], int, MultidimensionalArray[], int, int, MultidimensionalArray)"/>,
        /// <see cref="HLLCMomentumFlux_XDG.BorderEdgeFlux(double, int, MultidimensionalArray, MultidimensionalArray, bool, byte[], int, MultidimensionalArray[], int, int, MultidimensionalArray)"/>,
        /// <see cref="HLLCEnergyFlux_XDG.BorderEdgeFlux(double, int, MultidimensionalArray, MultidimensionalArray, bool, byte[], int, MultidimensionalArray[], int, int, MultidimensionalArray)"/>.
        /// In the XDG context, the <see cref="AdiabaticSlipWall"/> boundary condition is sensitive to round-off errors --> The HLLC flux computation is not robust in this case (intermediate wave speed has to be
        /// exactly zero), so the call of the inner edge fluxes <see cref="HLLCDensityFlux.InnerEdgeFlux(double, int, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray[], MultidimensionalArray[], int, int, MultidimensionalArray)"/> etc. does not give the correct results.
        /// In order not to make things more complex and confusing with the usage of inheritance etc., this solution was considered to be best (and easy to implement).
        /// </summary>
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            double[] Uout = GetExternalState(ref inp, Uin);

            XDGCompressibleBoundaryCondMap boundaryMap = (XDGCompressibleBoundaryCondMap)this.boundaryMap;
            // Get boundary condition on this edge
            BoundaryCondition boundaryCondition = boundaryMap.GetBoundaryConditionForSpecies(inp.EdgeTag, this.speciesName);

            if(boundaryCondition is AdiabaticSlipWall) {
                return XDGBorderEdgeFlux(inp.Normal, Uin, Uout, inp.D) * Vin;
            } else {
                double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
                return ret * (Vin);
            }
        }

        /// <summary>
        /// Helper function to compensate round-off problems (computation of the intermediate wave speed in all HLLC fluxes)
        /// This has been only experienced in the XDG context so far.
        /// </summary>
        private double XDGBorderEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) {

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

            // cf. Toro2009, equation 10.70
            double intermediateWaveSpeed = 0.0;

            double edgeFlux = 0.0;
            // cf. Toro2009, equation 10.71
            // flux = state.Velocity * (state.Energy + state.Pressure);
            if(intermediateWaveSpeed > 0.0) {
                edgeFlux = normalVelocityIn * (energyIn + pIn);
                if(waveSpeedIn <= 0.0) {
                    double factor = densityIn * intermediateWaveSpeed * MachScaling
                        + pIn / (waveSpeedIn - normalVelocityIn);
                    double modifiedEnergy = (waveSpeedIn - normalVelocityIn) /
                        (waveSpeedIn - intermediateWaveSpeed) *
                        (energyIn + factor * (intermediateWaveSpeed - normalVelocityIn));
                    edgeFlux += waveSpeedIn * (modifiedEnergy - energyIn);
                }
            } else {
                edgeFlux = normalVelocityOut * (energyOut + pOut);
                if(waveSpeedOut >= 0.0) {
                    double factor = densityOut * intermediateWaveSpeed * MachScaling
                        + pOut / (waveSpeedOut - normalVelocityOut);
                    double modifiedEnergy = (waveSpeedOut - normalVelocityOut) /
                        (waveSpeedOut - intermediateWaveSpeed) *
                        (energyOut + factor * (intermediateWaveSpeed - normalVelocityOut));
                    edgeFlux += waveSpeedOut * (modifiedEnergy - energyOut);
                }
            }

            //Debug.Assert(!double.IsNaN(edgeFlux) || double.IsInfinity(edgeFlux));
            return edgeFlux;

        }
    }
}
