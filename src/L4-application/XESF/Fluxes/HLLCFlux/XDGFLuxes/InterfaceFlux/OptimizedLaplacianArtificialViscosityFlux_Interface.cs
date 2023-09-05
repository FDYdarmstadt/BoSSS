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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.ShockCapturing;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace XDGShock.Fluxes {

    public class OptimizedLaplacianArtificialViscosityFlux_Interface : INonlinLevelSetForm_V, INonlinLevelSetForm_GradV, IEquationComponentChecking {

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
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ArgumentOrdering {
            get;
            private set;
        }

        public IList<string> ParameterOrdering {
            get {
                return new string[] { "artificialViscosity" };
            }
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            // Convert CommonParams to CommonParams
            CommonParams commonParams = new CommonParams {
                Normal = inp.Normal,
                X = inp.X,
                Parameters_IN = inp.Parameters_IN,
                Parameters_OUT = inp.Parameters_OUT,
                iEdge = int.MaxValue, // Alternative: use cell index as edge index --> Might this cause problems?
                GridDat = levelSetTracker.GridDat,
                time = inp.time
            };

            // Call bulk flux
            return ((IEdgeForm)bulkFlux).InnerEdgeForm(ref commonParams, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
            //throw new NotImplementedException();
        }
        #endregion

        #region INonlinLevelSetForm_V members
        public void NonlinInternalEdge_V(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_Vin, MultidimensionalArray Koeff_Vot) {
            // Convert LevSetIntParams to EdgeFormParams
            EdgeFormParams edgeParams = inp;
            MultidimensionalArray fin = Koeff_Vin;
            MultidimensionalArray fout = Koeff_Vot;

            // Call bulk flux
            //((INonlinEdgeForm_V)bulkFlux).InternalEdge(ref edgeParams, uA, uB, Grad_uA, Grad_uB, fin, fout);

            int NumOfCells = edgeParams.Len;
            Debug.Assert(NumOfCells == 1);
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fout.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            double Penalty = bulkFlux.Penalties[edgeParams.e0];
            int cell = 0;

            for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                                                            // SIPG Flux Loops
                double viscosityIn = edgeParams.ParameterVars_IN[0][cell, node];
                double viscosityOut = edgeParams.ParameterVars_OUT[0][cell, node];

                double flux = 0.0;
                for (int d = 0; d < levelSetTracker.GridDat.SpatialDimension; d++) {
                    double n = edgeParams.Normals[cell, node, d];
                    flux -= 0.5 * (viscosityIn * Grad_uA[0][cell, node, d] + viscosityOut * Grad_uB[0][cell, node, d]) * n;
                }
                flux += Math.Max(viscosityIn, viscosityOut) * (uA[0][cell, node] - uB[0][cell, node]) * Penalty;

                fin[cell, node] += flux;
                fout[cell, node] -= flux;
            }
        }
        #endregion

        #region INonlinLevelSetForm_GradV members
        public void NonlinInternalEdge_GradV(ref EdgeFormParams inp, MultidimensionalArray[] uA, MultidimensionalArray[] uB, MultidimensionalArray[] Grad_uA, MultidimensionalArray[] Grad_uB, MultidimensionalArray Koeff_GradVin, MultidimensionalArray Koeff_GradVot) {
            // Convert LevSetIntParams to EdgeFormParams
            EdgeFormParams edgeParams = inp;
            MultidimensionalArray fin = Koeff_GradVin;
            MultidimensionalArray fout = Koeff_GradVot;

            // Call bulk flux
            //((INonlinEdgeForm_GradV)bulkFlux).InternalEdge(ref edgeParams, uA, uB, Grad_uA, Grad_uB, fin, fout);

            int NumOfCells = edgeParams.Len;
            Debug.Assert(NumOfCells == 1);
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fout.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            int cell = 0;

            for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                double uJump = 0.5 * (uA[0][cell, node] - uB[0][cell, node]);
                double fluxIn = edgeParams.ParameterVars_IN[0][cell, node] * uJump;
                double fluxOut = edgeParams.ParameterVars_OUT[0][cell, node] * uJump;

                for (int d = 0; d < levelSetTracker.GridDat.SpatialDimension; d++) {
                    double n = edgeParams.Normals[cell, node, d];
                    fin[cell, node, d] -= fluxIn * n;
                    fout[cell, node, d] -= fluxOut * n;
                }
            }
        }
        #endregion
               
        #region IEquationComponentChecking members
        public bool IgnoreVectorizedImplementation => false;
        #endregion

        private readonly LevelSetTracker levelSetTracker;

        private readonly OptimizedLaplacianArtificialViscosityFlux bulkFlux;

        public OptimizedLaplacianArtificialViscosityFlux_Interface(LevelSetTracker levelSetTracker, string ArgumentVarName, double penaltySafetyFactor, double penaltyFactor, Dictionary<SpeciesId, MultidimensionalArray> cellLengthScales) {
            this.levelSetTracker = levelSetTracker;
            this.bulkFlux = new OptimizedLaplacianArtificialViscosityFlux(levelSetTracker, ArgumentVarName, penaltySafetyFactor, penaltyFactor, cellLengthScales);
            this.ArgumentOrdering = new string[] { ArgumentVarName };
        }
    }
}
