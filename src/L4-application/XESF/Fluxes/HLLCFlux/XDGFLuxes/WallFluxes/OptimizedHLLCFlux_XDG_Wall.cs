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
    public class OptimizedHLLCFlux_XDG_Wall : ILevelSetForm,
        //INonlinLevelSetForm_V,
        IEquationComponentChecking {

        #region ILevelSetForm members
        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
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

        public bool IgnoreVectorizedImplementation => false;

        //private static StreamWriter writer;

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            // Adiabatic slip wall
            //Vector normalVec = inp.Normal;
            Vector normalVec_A = new Vector(inp.Normal[0], inp.Normal[1]);
            Vector normalVec_B = new Vector(-inp.Normal[0], -inp.Normal[1]);

            for (int i = 0; i < uB.Length; i++) {
                if (i != 2) {
                    if (uB[i] == 0.0) {
                        throw new NotSupportedException("Some component of uB is zero (which is not m1)");
                    }
                }
            }

            double[] uWall = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_A, new StateVector(uB, this.material)).ToArray();
            double[] test = this.adiabaticSlipWall.GetBoundaryState(inp.time, inp.X, normalVec_B, new StateVector(uB, this.material)).ToArray();

            for (int i = 0; i < uWall.Length; i++) {
                if (uWall[i] != test[i]) {
                    throw new NotSupportedException("Normal vector does matter for adiabatic slip wall calculation");
                }
            }

            // Define MultidimensionalArrays for INonLinearFlux usage
            int D = inp.D;

            // normal
            MultidimensionalArray normal_A = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray normal_tmp = MultidimensionalArray.CreateWrapper(normalVec_A, normalVec_A.Dim);
            normal_A.SetSubArray(normal_tmp, 0, 0, -1);

            // normal flipped
            MultidimensionalArray normal_B = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray normalFlipped_tmp = MultidimensionalArray.CreateWrapper(normalVec_B, normalVec_B.Dim);
            normal_B.SetSubArray(normalFlipped_tmp, 0, 0, -1);

            MultidimensionalArray x = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray x_tmp = MultidimensionalArray.CreateWrapper(inp.X, inp.X.Dim);
            x.SetSubArray(x_tmp, 0, 0, -1);

            MultidimensionalArray[] uWall_mda = new MultidimensionalArray[D + 2];
            MultidimensionalArray[] uB_mda = new MultidimensionalArray[D + 2];
            for (int d = 0; d < D + 2; d++) {
                uWall_mda[d] = MultidimensionalArray.Create(1, 1);
                uWall_mda[d][0, 0] = uWall[d];
                uB_mda[d] = MultidimensionalArray.Create(1, 1);
                uB_mda[d][0, 0] = uB[d];
            }

            double xCoord = inp.X[0];
            double yCoord = inp.X[1];

            // Flux across the interface
            MultidimensionalArray outputNeg = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray outputPos = MultidimensionalArray.Create(1, 1);
            this.bulkFlux.InnerEdgeFlux(inp.time, int.MinValue, x, normal_A, uWall_mda, uB_mda, 0, 1, outputNeg);
            this.bulkFlux.InnerEdgeFlux(inp.time, int.MinValue, x, normal_B, uB_mda, uWall_mda, 0, 1, outputPos);

            //if (Math.Abs(outputNeg[0, 0] - outputPos[0, 0]) > 1e-14) {
            //    throw new NotSupportedException(String.Format("Flux inconsistency detected: Diff is {0}", outputNeg[0, 0] - outputPos[0, 0]));
            //}

            //double flux = outputNeg[0, 0] * (0.0 - vB);

            double flux;
            if (vA != 0.0) {
                flux = outputNeg[0, 0];
            } else {
                flux = outputPos[0, 0];
            }

            double fluxA;
            double fluxB;

            if (vA != 0.0) {
                fluxA = flux;
                fluxB = 0.0;

                normalVec_B[0] = 0.0;
                normalVec_B[1] = 0.0;
            } else {
                fluxA = 0.0;
                fluxB = flux;

                normalVec_A[0] = 0.0;
                normalVec_A[1] = 0.0;
            }

            //Console.WriteLine(String.Format("Flux at ({0:0.00000000}, {1:0.00000000}) = {2:0.00000000} \t normal = ({3:0.00000000}, {4:0.00000000})", inp.X[0], inp.X[1], flux, normalVec[0], normalVec[1]));

            // StreamWriter
            //if (writer == null) {
            //    writer = new StreamWriter("XDG_Flux.txt");
            //    writer.WriteLine("bulkFlux \t\t\t x \t y \t \t \t \t  ### \t n_x_Neg \t n_y_Neg \t flux_Neg \t ### \t n_x_Pos \t n_y_Pos \t flux_Pos \t ### \t Uin[0] \t Uin[1] \t Uin[2] \t Uin[3] \t ### \t Uout[0] \t Uout[1] \t Uout[2] \t Uout[3]");
            //}

            //string bulkFluxName = null;
            //if (this.bulkFlux is OptimizedHLLCDensityFlux) {
            //    bulkFluxName = "rho";
            //} else if (this.bulkFlux is OptimizedHLLCMomentumFlux tmp) {
            //    bulkFluxName = "m";
            //} else if (this.bulkFlux is OptimizedHLLCEnergyFlux) {
            //    bulkFluxName = "rhoE";
            //}

            //string resultLine;
            //if (vB != 0.0) {
            //    resultLine = String.Format(bulkFluxName + "\t\t ### \t {0:0.000000} \t {1:0.000000} \t  ### \t {2:0.000000} \t {3:0.000000} \t {4:0.000000} \t ### \t {5:0.000000} \t {6:0.000000} \t {7:0.000000} \t ### \t {8:0.000000} \t {9:0.000000} \t {10:0.000000} \t {11:0.000000} \t ### \t {12:0.000000} \t {13:0.000000} \t {14:0.000000} \t {15:0.000000}", inp.X[0], inp.X[1], normalVec_A[0], normalVec_A[1], fluxA, normalVec_B[0], normalVec_B[1], fluxB, uB[0], uB[1], uB[2], uB[3], uWall[0], uWall[1], uWall[2], uWall[3]);
            //    writer.WriteLine(resultLine);
            //    writer.Flush();
            //}

            //return flux;
            return flux * vB;
        }
        #endregion

        #region INonlinLevelSetForm_V members
        //public void LevelSetForm_V(LevSetIntParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_V) {
        //    int NumOfCells = inp.Len;
        //    Debug.Assert(NumOfCells == 1);
        //    Debug.Assert(Koeff_V.GetLength(0) == NumOfCells);

        //    // Split coefficient array into species A and B
        //    MultidimensionalArray fA = Koeff_V.ExtractSubArrayShallow(-1, -1, 0);
        //    MultidimensionalArray fB = Koeff_V.ExtractSubArrayShallow(-1, -1, 1);
        //    Debug.Assert(fA.Lengths.ListEquals(fB.Lengths, (a, b) => a == b));

        //    var tmp = MultidimensionalArray.Create(fB.Lengths);

        //    // Pass only those things that are really needed
        //    bulkFlux.InnerEdgeFlux(
        //        time: double.NaN,
        //        jEdge: -1,
        //        x: inp.Nodes,
        //        normal: inp.Normals,
        //        Uin: uA,
        //        Uout: uB,
        //        Offset: 0,
        //        Length: inp.Len,
        //        Output: tmp);

        //    fA.Acc(+1.0, tmp);
        //    fB.Acc(-1.0, tmp);
        //}
        #endregion

        //private readonly LevelSetTracker levelSetTracker;

        protected readonly OptimizedHLLCFlux bulkFlux;

        private readonly Material material;

        private readonly AdiabaticSlipWall adiabaticSlipWall;

        public OptimizedHLLCFlux_XDG_Wall(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component = int.MinValue) {
            //this.levelSetTracker = levelSetTracker;
            this.material = material;
            this.adiabaticSlipWall = new AdiabaticSlipWall(material);

            switch (flux) {
                case FluxVariables.Density:
                    this.bulkFlux = new OptimizedHLLCDensityFlux(boundaryConditionMap, material);
                    break;
                case FluxVariables.Momentum:
                    this.bulkFlux = new OptimizedHLLCMomentumFlux(boundaryConditionMap, component, material);
                    break;
                case FluxVariables.Energy:
                    this.bulkFlux = new OptimizedHLLCEnergyFlux(boundaryConditionMap, material);
                    break;
                default:
                    throw new NotImplementedException();
            }
        }
    }
}
