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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XESF.Fluxes {
    public class HLLCDensityFlux_XDG_Wall : ILevelSetForm, ISupportsJacobianComponent {

        #region ILevelSetForm members
        public int LevelSetIndex {
            get;
            private set;
        }

        public string PositiveSpecies {
            get;
            private set;
        }

        public string NegativeSpecies {
            get;
            private set;
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        //private static StreamWriter writer;

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            Vector normalVec_B = new Vector(-inp.Normal[0], -inp.Normal[1]);
            double[] uWall = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_B, new StateVector(uB, this.material)).ToArray();

#if DEBUG
            Vector normalVec_A = new Vector(inp.Normal[0], inp.Normal[1]);
            double[] test = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_A, new StateVector(uB, this.material)).ToArray();
            for (int i = 0; i < uWall.Length; i++) {
                if (uWall[i] != test[i]) {
                    throw new NotSupportedException("Normal vector should not matter for adiabatic slip wall flux computation");
                }
            }
#endif
            // Define MultidimensionalArrays for INonLinearFlux usage
            int D = inp.D;

            double[] Uin = uB;
            double[] Uout = uWall;
            Vector normal = normalVec_B;

            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double MachScaling = gamma * Mach * Mach;

            double densityIn = Uin[0];
            double densityOut = Uout[0];

            double momentumSquareIn = 0.0;
            double momentumSquareOut = 0.0;
            double normalVelocityIn = 0.0;
            double normalVelocityOut = 0.0;
            for (int d = 0; d < D; d++) {
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

            //double cIn = densityIn * (waveSpeedIn - normalVelocityIn);
            //double cOut = densityOut * (waveSpeedOut - normalVelocityOut);

            // cf. Toro2009, equation 10.70
            //double speedDiff = cOut * normalVelocityOut - cIn * normalVelocityIn;
            //speedDiff = 0.0;

            //double pIn_minus_pOut = pIn - pOut;
            //pIn_minus_pOut = 0.0;

            //double intermediateWaveSpeed =
            //    (speedDiff + pIn_minus_pOut / MachScaling) /
            //    (cOut - cIn);

            double intermediateWaveSpeed = 0.0;

            //double intermediateWaveSpeed =
            //    (cOut * normalVelocityOut - cIn * normalVelocityIn + (pIn - pOut) / MachScaling) /
            //    (cOut - cIn);

            double edgeFlux = 0.0;
            // cf. Toro2009, equation 10.71
            if (intermediateWaveSpeed > 0.0) {
                edgeFlux = normalVelocityIn * densityIn;
                if (waveSpeedIn <= 0.0) {
                    double modifiedDensity = densityIn * (waveSpeedIn - normalVelocityIn)
                        / (waveSpeedIn - intermediateWaveSpeed);
                    edgeFlux += waveSpeedIn * (modifiedDensity - densityIn);
                }
            } else {
                edgeFlux = normalVelocityOut * densityOut;
                if (waveSpeedOut >= 0.0) {
                    double modifiedDensity = densityOut * (waveSpeedOut - normalVelocityOut)
                        / (waveSpeedOut - intermediateWaveSpeed);
                    edgeFlux += waveSpeedOut * (modifiedDensity - densityOut);
                }
            }

            Debug.Assert(!double.IsNaN(edgeFlux) || double.IsInfinity(edgeFlux));

            //double fluxA;
            //double fluxB;

            //if (vA != 0.0) {
            //    fluxA = edgeFlux;
            //    fluxB = 0.0;

            //    normalVec_B[0] = 0.0;
            //    normalVec_B[1] = 0.0;
            //} else {
            //    fluxA = 0.0;
            //    fluxB = edgeFlux;

            //    normalVec_A[0] = 0.0;
            //    normalVec_A[1] = 0.0;
            //}

            // StreamWriter
            //if (writer == null) {
            //    writer = new StreamWriter("XDG_Flux_rho.txt");
            //    writer.WriteLine("bulkFlux \t\t\t x \t y \t \t \t \t  ### \t n_x_Neg \t n_y_Neg \t flux_Neg \t ### \t n_x_Pos \t n_y_Pos \t flux_Pos \t ### \t Uin[0] \t Uin[1] \t Uin[2] \t Uin[3] \t ### \t Uout[0] \t Uout[1] \t Uout[2] \t Uout[3]");
            //}

            //string bulkFluxName = "rho";

            //string resultLine;
            //if (vB != 0.0) {
            //    resultLine = String.Format(bulkFluxName + "\t\t ### \t {0:0.000000} \t {1:0.000000} \t  ### \t {2:0.000000} \t {3:0.000000} \t {4:0.000000} \t ### \t {5:0.000000} \t {6:0.000000} \t {7:0.000000} \t ### \t {8:0.000000} \t {9:0.000000} \t {10:0.000000} \t {11:0.000000} \t ### \t {12:0.000000} \t {13:0.000000} \t {14:0.000000} \t {15:0.000000}", inp.X[0], inp.X[1], normalVec_A[0], normalVec_A[1], fluxA, normalVec_B[0], normalVec_B[1], fluxB, uB[0], uB[1], uB[2], uB[3], uWall[0], uWall[1], uWall[2], uWall[3]);
            //    writer.WriteLine(resultLine);
            //    writer.Flush();
            //}

            return edgeFlux * vB;
        }
        #endregion

        private readonly Material material;

        private readonly AdiabaticSlipWall adiabaticSlipWall;

        public HLLCDensityFlux_XDG_Wall(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;
            this.adiabaticSlipWall = new AdiabaticSlipWall(material);

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }
        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");
            return new IEquationComponent[] {
                new LevelSetFormDifferentiator(this,SpatialDimension)
            };
        }
    }
}
