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
    public class RoeDensityFlux_XDG_Wall : ILevelSetForm, ISupportsJacobianComponent, IEquationComponentChecking {

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


       /* #region Setup 
                // ----------------------------------------------------------------------
                Vector normalVec_B = new Vector(-inp.Normal[0], -inp.Normal[1]);
                double[] uWall = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_B, new StateVector(uB, this.material)).ToArray();

#if DEBUG
                Vector normalVec_A = new Vector(inp.Normal[0], inp.Normal[1]);
                double[] test = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_A, new StateVector(uB, this.material)).ToArray();
                for(int i = 0; i < uWall.Length; i++) {
                    if(uWall[i] != test[i]) {
                        throw new NotSupportedException("Normal vector should not matter for adiabatic slip wall flux computation");
                    }
                }
#endif

                int D = inp.D;
                double[] Uin = uB;
                double[] Uout = uWall;
                Vector normal = normalVec_B;

                // ----------------------------------------------------------------------

                MultidimensionalArray V0_inv = MultidimensionalArray.Create(D + 2, D + 2);
                double densityIn = Uin[0];
                double densityOut = Uout[0];
                double[] vIn = new double[D];
                double[] vOut = new double[D];

                double momentumSquaredIn = 0;
                double momentumSquaredOut = 0;
                for(int d = 0; d < D; d++) {
                    vIn[d] = Uin[d + 1] / densityIn;
                    vOut[d] = Uout[d + 1] / densityOut;
                    momentumSquaredIn += Math.Pow(Uin[d + 1], 2);
                    momentumSquaredOut += Math.Pow(Uout[d + 1], 2);
                }
                double energyIn = Uin[D + 1];
                double energyOut = Uout[D + 1];

                // copied from HLLCEnergyFlux.cs, to calculate the pressure
                double gamma = this.material.EquationOfState.HeatCapacityRatio;
                double Mach = this.material.MachNumber;
                double gammaMachSquared = gamma * Mach * Mach;

                double pressureIn = (gamma - 1.0) * (energyIn - gammaMachSquared * 0.5 * momentumSquaredIn / densityIn);
                double pressureOut = (gamma - 1.0) * (energyOut - gammaMachSquared * 0.5 * momentumSquaredOut / densityOut);

                // Enthalpy: H = (E + pressure) / roh ; Toro p.347
                double enthalpyIn = (energyIn + pressureIn) / densityIn;
                double enthalpyOut = (energyOut + pressureOut) / densityOut;

                Vector vAverage = new Vector(D);

                for(int i = 0; i < D; i++) {
                    // 0:uAverage, 1:vAverage, 2:wAverage
                    vAverage[i] = (Math.Sqrt(densityIn) * vIn[i] + Math.Sqrt(densityOut) * vOut[i])
                        / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
                }
                double velocityAverage = vAverage.L2Norm();

                double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                    / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));

                double speedOfSoundAverage = Math.Sqrt((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)));

                // Toro p.374, eq. (11.58), EigenValues[0...4]
                double vN = vAverage * (new Vector(normal));

                //for Eigenvals
                double s_alpha = 10;
                double SmoothedAbs(double x) {
                    return x * Math.Tanh(s_alpha * x);
                }

                //Eigenvals 
                double[] eigenVals = new double[D + 2];
                eigenVals[0] = SmoothedAbs(vN - speedOfSoundAverage);
                for(int i = 1; i < D + 1; i++) {
                    // if D=2 or D=3 set the values; else it's = 0
                    eigenVals[i] = SmoothedAbs(vN);
                }
                eigenVals[D + 1] = SmoothedAbs(vN + speedOfSoundAverage);

            #region V0_inv 
                //first row
                double cdivgam = speedOfSoundAverage / (gamma - 1);
                V0_inv[0, 0] = vAverage.L2NormPow2() / 2 + vN * cdivgam;
                for(int j = 1; j < D + 1; j++) {
                    V0_inv[0, j] = -vAverage[j - 1] - cdivgam * normal[j - 1];
                }
                V0_inv[0, D + 1] = 1;

                //midle rows
                for(int i = 1; i < D + 1; i++) {
                    //first col
                    V0_inv[i, 0] = (2 * speedOfSoundAverage * cdivgam) * ((vN + 1) * normal[i - 1] - vAverage[i - 1]) - vAverage.L2NormPow2() * normal[i - 1];

                    //middle part
                    for(int j = 1; j < D + 1; j++) {
                        if(i == j) {
                            V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +                //xxx i oder j jetzt I'm confused
                                (2 * speedOfSoundAverage * cdivgam) * (1 - normal[i - 1] * normal[j - 1]);
                        } else {
                            V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +
                                (2 * speedOfSoundAverage * cdivgam) * (-normal[i - 1] * normal[j - 1]);
                        }

                    }
                    //last col
                    V0_inv[i, D + 1] = -2 * normal[i - 1];
                }

                //last row
                V0_inv[D + 1, 0] = vAverage.L2NormPow2() / 2 - vN * cdivgam;
                for(int j = 1; j < D + 1; j++) {
                    V0_inv[D + 1, j] = -vAverage[j - 1] + cdivgam * normal[j - 1];
                }
                V0_inv[D + 1, D + 1] = 1;

                V0_inv.Scale((gamma - 1) / (2 * speedOfSoundAverage * speedOfSoundAverage));
            #endregion V0_inv
        #endregion Setup
       */
            
            #region Main

            // IMPORT SETUP
            ///
            (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, 
                MultidimensionalArray V0_inv, double[] uWall) = Setup(ref inp, uA, uB);
                        
            int D = inp.D;
            double[] Uin = uB;
            double[] Uout = uWall; 
            
            double vN = vAverage * (new Vector(inp.Normal)); //xxx Error wegen ref Schlüsselwort ?? 
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            ///

            Vector n = new Vector(inp.Normal);
            double FL = Flux(uA, D) * n;
            double FR = Flux(uB, D) * n;
            double[] Udiff = Uin.CloneAs();
            Udiff.AccV(-1, Uout);
            var result = V0_inv.MatVecMul(1.0, Udiff);

            // +++++++++++++++ V0 +++++++++++++++
            double[] V0 = new double[D + 2];
            V0[0] = 1;
            for(int i = 0; i < D; i++) {
                V0[i + 1] = inp.Normal[i];
            }
            V0[D + 1] = 1;
            // +++++++++++++++ end +++++++++++++++

            double k_j = 0;
            for(int i = 0; i < D + 2; i++) {
                k_j += V0[i] * eigenVals[i] * result[i];
            }

            double densityFlux = 0.5 * (FL + FR) + 0.5 * k_j;

            if(densityFlux.IsNaN()) {
                Console.WriteLine("*************** Fehler im Code: Density Flux ***************");
            }
            return densityFlux;
        #endregion Main
        }
    #endregion ILevelSetForm members


         //SETUP FUNCTION ---------------------------------------------------------------
             #region Setup 
             //public (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) Setup(double[] normal, double[] Uin, double[] Uout, int D) {
             public (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv, double[] uWall) Setup(ref CommonParams inp, double[] uA, double[] uB) {

                 // ----------------------------------------------------------------------
                 Vector normalVec_B = new Vector(-inp.Normal[0], -inp.Normal[1]);
                 double[] uWall = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_B, new StateVector(uB, this.material)).ToArray();

        #if DEBUG
                 Vector normalVec_A = new Vector(inp.Normal[0], inp.Normal[1]);
                 double[] test = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_A, new StateVector(uB, this.material)).ToArray();
                 for(int i = 0; i < uWall.Length; i++) {
                     if(uWall[i] != test[i]) {
                         throw new NotSupportedException("Normal vector should not matter for adiabatic slip wall flux computation");
                     }
                 }
        #endif

                 int D = inp.D;
                 double[] Uin = uB;
                 double[] Uout = uWall;
                 Vector normal = normalVec_B;

                 // ----------------------------------------------------------------------

                 MultidimensionalArray V0_inv = MultidimensionalArray.Create(D + 2, D + 2);
                 double densityIn = Uin[0];
                 double densityOut = Uout[0];
                 double[] vIn = new double[D];
                 double[] vOut = new double[D];

                 double momentumSquaredIn = 0;
                 double momentumSquaredOut = 0;
                 for(int d = 0; d < D; d++) {
                     vIn[d] = Uin[d + 1] / densityIn;
                     vOut[d] = Uout[d + 1] / densityOut;
                     momentumSquaredIn += Math.Pow(Uin[d + 1], 2);
                     momentumSquaredOut += Math.Pow(Uout[d + 1], 2);
                 }
                 double energyIn = Uin[D + 1];
                 double energyOut = Uout[D + 1];

                 // copied from HLLCEnergyFlux.cs, to calculate the pressure
                 double gamma = this.material.EquationOfState.HeatCapacityRatio;
                 double Mach = this.material.MachNumber;
                 double gammaMachSquared = gamma * Mach * Mach;

                 double pressureIn = (gamma - 1.0) * (energyIn - gammaMachSquared * 0.5 * momentumSquaredIn / densityIn);
                 double pressureOut = (gamma - 1.0) * (energyOut - gammaMachSquared * 0.5 * momentumSquaredOut / densityOut);

                 // Enthalpy: H = (E + pressure) / roh ; Toro p.347
                 double enthalpyIn = (energyIn + pressureIn) / densityIn;
                 double enthalpyOut = (energyOut + pressureOut) / densityOut;

                 Vector vAverage = new Vector(D);

                 for(int i = 0; i < D; i++) {
                     // 0:uAverage, 1:vAverage, 2:wAverage
                     vAverage[i] = (Math.Sqrt(densityIn) * vIn[i] + Math.Sqrt(densityOut) * vOut[i])
                         / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
                 }
                 double velocityAverage = vAverage.L2Norm();

                 double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                     / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));

                 double speedOfSoundAverage = Math.Sqrt((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)));

                 // Toro p.374, eq. (11.58), EigenValues[0...4]
                 double vN = vAverage * (new Vector(normal));

                 //for Eigenvals
                 double s_alpha = 10;
                 double SmoothedAbs(double x) {
                     return x * Math.Tanh(s_alpha * x);
                 }

                 //Eigenvals 
                 double[] eigenVals = new double[D + 2];
                 eigenVals[0] = SmoothedAbs(vN - speedOfSoundAverage);
                 for(int i = 1; i < D + 1; i++) {
                     // if D=2 or D=3 set the values; else it's = 0
                     eigenVals[i] = SmoothedAbs(vN);
                 }
                 eigenVals[D + 1] = SmoothedAbs(vN + speedOfSoundAverage);

                 #region V0_inv 
                 //first row
                 double cdivgam = speedOfSoundAverage / (gamma - 1);
                 V0_inv[0, 0] = vAverage.L2NormPow2() / 2 + vN * cdivgam;
                 for(int j = 1; j < D + 1; j++) {
                     V0_inv[0, j] = -vAverage[j - 1] - cdivgam * normal[j - 1];
                 }
                 V0_inv[0, D + 1] = 1;

                 //midle rows
                 for(int i = 1; i < D + 1; i++) {
                     //first col
                     V0_inv[i, 0] = (2 * speedOfSoundAverage * cdivgam) * ((vN + 1) * normal[i - 1] - vAverage[i - 1]) - vAverage.L2NormPow2() * normal[i - 1];

                     //middle part
                     for(int j = 1; j < D + 1; j++) {
                         if(i == j) {
                             V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +                //xxx i oder j jetzt I'm confused
                                 (2 * speedOfSoundAverage * cdivgam) * (1 - normal[i - 1] * normal[j - 1]);
                         } else {
                             V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +
                                 (2 * speedOfSoundAverage * cdivgam) * (-normal[i - 1] * normal[j - 1]);
                         }

                     }
                     //last col
                     V0_inv[i, D + 1] = -2 * normal[i - 1];
                 }

                 //last row
                 V0_inv[D + 1, 0] = vAverage.L2NormPow2() / 2 - vN * cdivgam;
                 for(int j = 1; j < D + 1; j++) {
                     V0_inv[D + 1, j] = -vAverage[j - 1] + cdivgam * normal[j - 1];
                 }
                 V0_inv[D + 1, D + 1] = 1;

                 V0_inv.Scale((gamma - 1) / (2 * speedOfSoundAverage * speedOfSoundAverage));
                 #endregion V0_inv


                 return (vAverage, speedOfSoundAverage, enthalpyAverage, eigenVals, V0_inv, uWall);
             }
             #endregion Setup
        

     #region Flux-Vector
        public Vector Flux(double[] U, int D) {
            Vector fluxvector = new Vector(D);
            for(int d = 0; d < D; d++) {
                // state.Momentum
                double flux = U[d + 1];
                Debug.Assert(!double.IsNaN(flux) || double.IsInfinity(flux));
                fluxvector[d] = flux;
            }
            return fluxvector;
        }
     #endregion Flux-Vector

     #region IEquationComponentChecking members
     public bool IgnoreVectorizedImplementation => false;
     #endregion

     //private readonly LevelSetTracker levelSetTracker;

     private readonly Material material;

     private readonly AdiabaticSlipWall adiabaticSlipWall;

     public RoeDensityFlux_XDG_Wall(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
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
