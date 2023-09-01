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
using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace XESF.Fluxes {
    public class RoeMomentumFlux_XDG_Wall : ILevelSetForm, ISupportsJacobianComponent, IEquationComponentChecking {

        #region ILevelSetForm members

        #region presets
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
        #endregion presets

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            
            #region Main

            // IMPORT SETUP
            ///
            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals,
                MultidimensionalArray V0_inv, double[] uWall) = Setup(inp, uA, uB);

            int D = inp.D;
            double[] Uin = uB;
            double[] Uout = uWall;
            Vector normal = new Vector(inp.Normal); //xxx schauen ob das die richtige Referenz ist

            double[] Udiff = Uin.CloneAs();
            Udiff.AccV(-1, Uout);
            var result = V0_inv.MatVecMul(1.0, Udiff);

            Vector n = new Vector(normal);
            double FL = Flux(Uin, D) * n;
            double FR = Flux(Uout, D) * n;

            // +++++++++++++++ V0 +++++++++++++++
            double[] V0 = new double[D + 2];
            V0[0] = vAverage[component] - speedOfSoundAverage * normal[component];
            //xxx wenn Error, dann u.a. hier prüfen
            for(int i = 1; i < D + 1; i++) {
                if(component == i - 1) {
                    V0[i] = (vAverage[component] - normal[component]) * normal[i - 1] + 1;
                } else {
                    V0[i] = (vAverage[component] - normal[component]) * normal[i - 1];
                }
            }
            V0[D + 1] = vAverage[component] + speedOfSoundAverage * normal[component];
            // +++++++++++++++ end +++++++++++++++

            double k_j = 0;
            for(int i = 0; i < D + 2; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }

            double momentumFlux = 0.5 * (FL + FR) + 0.5 * k_j;

            if(momentumFlux.IsNaN()) {
                Console.WriteLine("*************** Fehler im Code: Momentum Flux ***************");
            }
            return momentumFlux;




            #endregion Main


        }

        private (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv, double[] uWall) Setup(CommonParams inp, double[] uA, double[] uB) {
            throw new NotImplementedException();
        }
        #endregion

        #region Flux-Vector
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
            for(int d = 0; d < D; d++) {
                Output[d] += U[component + 1] * U[d + 1] / density;
                momentumSquared += U[d + 1] * U[d + 1];
            }

            double flux = (gamma - 1.0) * (energy - gammaMachSquared * 0.5 * momentumSquared / density) / gammaMachSquared;
            //Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
            Output[component] += flux;
            return Output;
        }
        #endregion Flux-Vector


        #region //right now this method is not used
        //#region INonlinLevelSetForm_V members
        //public void NonlinInternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
        //    Debug.Assert(inp.Len == 1);
        //    Debug.Assert(Koeff_Vin.GetLength(0) == inp.Len);
        //    Debug.Assert(Koeff_Vot.GetLength(0) == inp.Len);
        //    MultidimensionalArray fA = Koeff_Vin;
        //    MultidimensionalArray fB = Koeff_Vot;
        //    Debug.Assert(fA.Lengths.ListEquals(fB.Lengths, (a, b) => a == b));

        //    MultidimensionalArray Output = MultidimensionalArray.Create(fA.Lengths);

        //    #region Flux computation
        //    int Length = inp.Len;
        //    int Offset = 0;
        //    int NoOfNodes = uA[0].GetLength(1);

        //    int D = CompressibleEnvironment.NumberOfDimensions;
        //    double gamma = this.material.EquationOfState.HeatCapacityRatio;
        //    double Mach = this.material.MachNumber;
        //    double MachScaling = gamma * Mach * Mach;

        //    // It is only called for one cell
        //    int e = 0;
        //    //for (int e = 0; e < Length; e++) {

        //    MultidimensionalArray[] uWall = new MultidimensionalArray[uB.Length];
        //    for(int i = 0; i < uB.Length; i++) {
        //        uWall[i] = MultidimensionalArray.Create(uB[i].GetLength(0), uB[i].GetLength(1));
        //    }
        //    this.adiabaticSlipWall.GetBoundaryState(uWall, inp.time, inp.Nodes, inp.Normals, uB, e + Offset, Length, normalFlipped: true, material: material);

        //    // Loop over nodes
        //    for(int n = 0; n < NoOfNodes; n++) {
        //        double densityIn = uB[0][e + Offset, n];
        //        double densityOut = uWall[0][e + Offset, n];

        //        double momentumSquareIn = 0.0;
        //        double momentumSquareOut = 0.0;
        //        double normalVelocityIn = 0.0;
        //        double normalVelocityOut = 0.0;
        //        for(int d = 0; d < D; d++) {
        //            momentumSquareIn += uB[d + 1][e + Offset, n] * uB[d + 1][e + Offset, n];
        //            momentumSquareOut += uWall[d + 1][e + Offset, n] * uWall[d + 1][e + Offset, n];
        //            normalVelocityIn += uB[d + 1][e + Offset, n] * -inp.Normals[e + Offset, n, d];
        //            normalVelocityOut += uWall[d + 1][e + Offset, n] * -inp.Normals[e + Offset, n, d];
        //        }
        //        normalVelocityIn /= densityIn;
        //        normalVelocityOut /= densityOut;

        //        double energyIn = uB[D + 1][e + Offset, n];
        //        double energyOut = uWall[D + 1][e + Offset, n];

        //        double pIn = (gamma - 1.0) * (energyIn - MachScaling * 0.5 * momentumSquareIn / densityIn);
        //        double pOut = (gamma - 1.0) * (energyOut - MachScaling * 0.5 * momentumSquareOut / densityOut);

        //        double speedOfSoundIn = Math.Sqrt(pIn / densityIn) / Mach;
        //        double speedOfSoundOut = Math.Sqrt(pOut / densityOut) / Mach;
        //        Debug.Assert(!double.IsNaN(speedOfSoundIn) || double.IsInfinity(speedOfSoundIn));
        //        Debug.Assert(!double.IsNaN(speedOfSoundOut) || double.IsInfinity(speedOfSoundOut));

        //        double densityMean = 0.5 * (densityIn + densityOut);
        //        double pressureMean = 0.5 * (pIn + pOut);
        //        double speedOfSoundMean = 0.5 * (speedOfSoundIn + speedOfSoundOut);
        //        double velocityJump = normalVelocityOut - normalVelocityIn;

        //        // Calculate the pressure estimate at the edge and use it to
        //        // calculate the components of the correction factor qIn and qOut
        //        double pStar = pressureMean - 0.5 * MachScaling * velocityJump * speedOfSoundMean * densityMean;
        //        pStar = Math.Max(0.0, pStar);

        //        double qIn = 1.0;
        //        if(pStar > pIn) {
        //            qIn = Math.Sqrt(1.0 + 0.5 * (gamma + 1.0) * (pStar / pIn - 1.0) / gamma);
        //        }

        //        double qOut = 1.0;
        //        if(pStar > pOut) {
        //            qOut = Math.Sqrt(1.0 + 0.5 * (gamma + 1.0) * (pStar / pOut - 1.0) / gamma);
        //        }

        //        // Determine the wave speeds
        //        double waveSpeedIn = normalVelocityIn - speedOfSoundIn * qIn;
        //        double waveSpeedOut = normalVelocityOut + speedOfSoundOut * qOut;

        //        //double cIn = densityIn * (waveSpeedIn - normalVelocityIn);
        //        //double cOut = densityOut * (waveSpeedOut - normalVelocityOut);

        //        // cf. Toro2009, equation 10.70
        //        //double speedDiff = cOut * normalVelocityOut - cIn * normalVelocityIn;
        //        //speedDiff = 0.0;

        //        //double pIn_minus_pOut = pIn - pOut;
        //        //pIn_minus_pOut = 0.0;

        //        //double intermediateWaveSpeed =
        //        //    (speedDiff + pIn_minus_pOut / MachScaling) /
        //        //    (cOut - cIn);

        //        double intermediateWaveSpeed = 0.0;

        //        //double intermediateWaveSpeed =
        //        //    (cOut * normalVelocityOut - cIn * normalVelocityIn + (pIn - pOut) / MachScaling) /
        //        //    (cOut - cIn);

        //        double edgeFlux = 0.0;
        //        // cf. Toro2009, equation 10.71
        //        // flux = state.Momentum[MomentumComponent] * state.Velocity
        //        //   + state.Pressure * ComponentVector;
        //        if(intermediateWaveSpeed > 0.0) {
        //            edgeFlux = normalVelocityIn * uB[component + 1][e + Offset, n]
        //                + pIn / MachScaling * -inp.Normals[e + Offset, n, component];
        //            if(waveSpeedIn <= 0.0) {
        //                double modifiedMomentum = densityIn *
        //                    (waveSpeedIn - normalVelocityIn) /
        //                    (waveSpeedIn - intermediateWaveSpeed) *
        //                    (uB[component + 1][e + Offset, n] / densityIn +
        //                        (intermediateWaveSpeed - normalVelocityIn) * -inp.Normals[e + Offset, n, component]);
        //                edgeFlux += waveSpeedIn * (modifiedMomentum - uB[component + 1][e + Offset, n]);
        //            }
        //        } else {
        //            edgeFlux = normalVelocityOut * uWall[component + 1][e + Offset, n]
        //                + pOut / MachScaling * -inp.Normals[e + Offset, n, component];
        //            if(waveSpeedOut >= 0.0) {
        //                double modifiedMomentum = densityOut *
        //                    (waveSpeedOut - normalVelocityOut) /
        //                    (waveSpeedOut - intermediateWaveSpeed) *
        //                    (uWall[component + 1][e + Offset, n] / densityOut +
        //                        (intermediateWaveSpeed - normalVelocityOut) * -inp.Normals[e + Offset, n, component]);
        //                edgeFlux += waveSpeedOut * (modifiedMomentum - uWall[component + 1][e + Offset, n]);
        //            }
        //        }

        //        Debug.Assert(!double.IsNaN(edgeFlux) || double.IsInfinity(edgeFlux));
        //        Output[e + Offset, n] += edgeFlux;
        //    }
        //    //}
        //    #endregion

        //    //fA.Acc(+1.0, Output);
        //    fB.Acc(1.0, Output);
        //}
        #endregion

        #region IEquationComponentChecking members
        public bool IgnoreVectorizedImplementation => false;
        #endregion

        //private readonly LevelSetTracker levelSetTracker;

        private readonly Material material;

        private readonly AdiabaticSlipWall adiabaticSlipWall;

        private readonly int component;

        public RoeMomentumFlux_XDG_Wall(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component = int.MinValue, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;
            this.adiabaticSlipWall = new AdiabaticSlipWall(material);
            this.component = component;

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
