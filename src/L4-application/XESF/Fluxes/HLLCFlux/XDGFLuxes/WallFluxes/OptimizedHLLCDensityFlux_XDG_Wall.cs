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
using BoSSS.Solution.CompressibleFlowCommon.Convection;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace XESF.Fluxes {
    public class OptimizedHLLCDensityFlux_XDG_Wall : ILevelSetForm, INonlinLevelSetForm_V, IEquationComponentChecking {

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

        #region INonlinLevelSetForm_V members
        public void NonlinInternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
            Debug.Assert(inp.Len == 1);
            Debug.Assert(Koeff_Vin.GetLength(0) == inp.Len);
            Debug.Assert(Koeff_Vot.GetLength(0) == inp.Len);
            MultidimensionalArray fA = Koeff_Vin;
            MultidimensionalArray fB = Koeff_Vot;
            Debug.Assert(fA.Lengths.ListEquals(fB.Lengths, (a, b) => a == b));

            MultidimensionalArray Output = MultidimensionalArray.Create(fA.Lengths);

            #region Flux computation
            int Length = inp.Len;
            int Offset = 0;
            int NoOfNodes = uA[0].GetLength(1);

            int D = CompressibleEnvironment.NumberOfDimensions;
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double MachScaling = gamma * Mach * Mach;

            // It is only called for one cell
            int e = 0;
            //for (int e = 0; e < Length; e++) {

            MultidimensionalArray[] uWall = new MultidimensionalArray[uB.Length];
            for (int i = 0; i < uB.Length; i++) {
                uWall[i] = MultidimensionalArray.Create(uB[i].GetLength(0), uB[i].GetLength(1));
            }
            this.adiabaticSlipWall.GetBoundaryState(uWall, inp.time, inp.Nodes, inp.Normals, uB, e + Offset, Length, normalFlipped: true, material: material);

            // Loop over nodes
            for (int n = 0; n < NoOfNodes; n++) {
                double densityIn = uB[0][e + Offset, n];
                double densityOut = uWall[0][e + Offset, n];

                double momentumSquareIn = 0.0;
                double momentumSquareOut = 0.0;
                double normalVelocityIn = 0.0;
                double normalVelocityOut = 0.0;
                for (int d = 0; d < D; d++) {
                    momentumSquareIn += uB[d + 1][e + Offset, n] * uB[d + 1][e + Offset, n];
                    momentumSquareOut += uWall[d + 1][e + Offset, n] * uWall[d + 1][e + Offset, n];
                    normalVelocityIn += uB[d + 1][e + Offset, n] * -inp.Normals[e + Offset, n, d];
                    normalVelocityOut += uWall[d + 1][e + Offset, n] * -inp.Normals[e + Offset, n, d];
                }
                normalVelocityIn /= densityIn;
                normalVelocityOut /= densityOut;

                double energyIn = uB[D + 1][e + Offset, n];
                double energyOut = uWall[D + 1][e + Offset, n];

                double pIn = (gamma - 1.0) * (energyIn - MachScaling * 0.5 * momentumSquareIn / densityIn);
                double pOut = (gamma - 1.0) * (energyOut - MachScaling * 0.5 * momentumSquareOut / densityOut);

                double speedOfSoundIn = Math.Sqrt(pIn / densityIn) / Mach;
                double speedOfSoundOut = Math.Sqrt(pOut / densityOut) / Mach;
                Debug.Assert(!double.IsNaN(speedOfSoundIn) || double.IsInfinity(speedOfSoundIn));
                Debug.Assert(!double.IsNaN(speedOfSoundOut) || double.IsInfinity(speedOfSoundOut));

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
                Output[e + Offset, n] += edgeFlux;
            }
            //}
            #endregion  

            //fA.Acc(+1.0, Output);
            fB.Acc(1.0, Output);
        }
        #endregion

        #region IEquationComponentChecking members
        public bool IgnoreVectorizedImplementation => false;
        #endregion

        //private readonly LevelSetTracker levelSetTracker;

        private readonly Material material;

        private readonly AdiabaticSlipWall adiabaticSlipWall;

        public OptimizedHLLCDensityFlux_XDG_Wall(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;
            this.adiabaticSlipWall = new AdiabaticSlipWall(material);

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }
    }
}
