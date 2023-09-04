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

namespace XESF.Fluxes {
    public class OptimizedHLLCFlux_XDG_Interface : ILevelSetForm, INonlinLevelSetForm_V, IEquationComponentChecking {

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

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            // Define MultidimensionalArrays for INonLinearFlux usage
            int D = inp.D;

            MultidimensionalArray normal = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray normal_tmp = MultidimensionalArray.CreateWrapper(inp.Normal, inp.Normal.Dim);
            normal.SetSubArray(normal_tmp, 0, 0, -1);

            MultidimensionalArray x = MultidimensionalArray.Create(1, 1, D);
            MultidimensionalArray x_tmp = MultidimensionalArray.CreateWrapper(inp.X, inp.X.Dim);
            x.SetSubArray(x_tmp, 0, 0, -1);

            MultidimensionalArray[] Uin = new MultidimensionalArray[D + 2];
            MultidimensionalArray[] Uout = new MultidimensionalArray[D + 2];
            for (int d = 0; d < D + 2; d++) {
                Uin[d] = MultidimensionalArray.Create(1, 1);
                Uin[d][0, 0] = uA[d];
                Uout[d] = MultidimensionalArray.Create(1, 1);
                Uout[d][0, 0] = uB[d];
            }

            // Flux across the interface
            MultidimensionalArray outputNeg = MultidimensionalArray.Create(1, 1);
            this.bulkFlux.InnerEdgeFlux(inp.time, int.MinValue, x, normal, Uin, Uout, 0, 1, outputNeg);
            double flux = outputNeg[0, 0];

            return flux * (vA - vB);

            //throw new NotSupportedException("Has to be checked again.");
        }
        #endregion

        #region INonlinLevelSetForm_V members
        public void NonlinInternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
            Debug.Assert(inp.Len == 1);
            Debug.Assert(Koeff_Vin.GetLength(0) == inp.Len);
            Debug.Assert(Koeff_Vot.GetLength(0) == inp.Len);

            // Split coefficient array into species A and B
            MultidimensionalArray fA = Koeff_Vin;//.ExtractSubArrayShallow(-1, -1, 0);
            MultidimensionalArray fB = Koeff_Vot;//.ExtractSubArrayShallow(-1, -1, 1);
            Debug.Assert(fA.Lengths.ListEquals(fB.Lengths, (a, b) => a == b));

            MultidimensionalArray tmp = MultidimensionalArray.Create(fA.Lengths);

            // Pass only those things that are really needed
            bulkFlux.InnerEdgeFlux(
                time: double.NaN,
                jEdge: -1,
                x: inp.Nodes,
                normal: inp.Normals,
                Uin: uA,
                Uout: uB,
                Offset: 0,
                Length: inp.Len,
                Output: tmp);

            fA.Acc(+1.0, tmp);
            fB.Acc(-1.0, tmp);
        }
        #endregion

        #region IEquationComponentChecking
        public bool IgnoreVectorizedImplementation => false;
        #endregion

        //private readonly LevelSetTracker levelSetTracker;

        protected readonly OptimizedHLLCFlux bulkFlux;

        public OptimizedHLLCFlux_XDG_Interface(LevelSetTracker levelSetTracker, IBoundaryConditionMap boundaryConditionMap, Material material, FluxVariables flux, int component = int.MinValue, int levelSetIndex = 0, string posSpecies = "B", string negSpecies = "A") {
            //this.levelSetTracker = levelSetTracker;

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

            LevelSetIndex = levelSetIndex;
            PositiveSpecies = posSpecies;
            NegativeSpecies = negSpecies;
        }
    }
}
