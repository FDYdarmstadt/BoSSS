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
using BoSSS.Foundation.Grid.Classic;
using CNS.Diffusion;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace CNS.ShockCapturing {

    class OptimizedLaplacianArtificialViscosityFlux : INonlinear2ndOrderForm {

        private GridData gridData;

        private int dimension;

        private string ArgumentName;

        public bool AdiabaticWall { get; set; }

        private MultidimensionalArray cellLengthScale;

        private double penaltyFactor;

        public OptimizedLaplacianArtificialViscosityFlux(GridData gridData, double order, int dimension, string ArgumentVarName) {
            this.dimension = dimension;
            cellLengthScale = gridData.Cells.cj;
            this.gridData = gridData;
            penaltyFactor = (order + 1) * (order + (double)dimension) / (double)dimension;
            ArgumentName = ArgumentVarName;
        }

        #region IEquationComponent Members
        IList<string> IEquationComponent.ArgumentOrdering {
            get {
                return new string[] { ArgumentName };
            }
        }

        IList<string> IEquationComponent.ParameterOrdering {
            get {
                return new string[] { Variables.ArtificialViscosity };
            }
        }
        #endregion

        #region IEdgeComponent Members
        TermActivationFlags IEdgeForm.BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        TermActivationFlags IEdgeForm.InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV);
            }
        }

        double IEdgeForm.InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            throw new NotImplementedException();
        }

        double IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotImplementedException();
        }
        #endregion

        #region INonlineEdgeform_GradV Members
        void INonlinEdgeForm_GradV.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    double uJump = 0.5 * (Uin[0][cell, node] - Uout[0][cell, node]);
                    double fluxIn = efp.ParameterVars_IN[0][cell, node] * uJump;
                    double fluxOut = efp.ParameterVars_OUT[0][cell, node] * uJump;

                    for (int d = 0; d < dimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        fin[cell, node, d] -= fluxIn * n;
                        fot[cell, node, d] -= fluxOut * n;
                    }
                }
            }
        }

        void INonlinEdgeForm_GradV.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {

            //Boundary value is zero neumann boundary!
            // i.e. do nothing!

            //int NumOfCells = efp.Len;
            //Debug.Assert(fin.GetLength(0) == NumOfCells);
            //int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            //int NumOfArguments = this.ArgumentOrdering.Count;
            //Debug.Assert(NumOfArguments == Uin.Length);
            //Debug.Assert(NumOfArguments == GradUin.Length);

            //double[] U_in = new double[NumOfArguments];
            //double[] U_ot = new double[NumOfArguments];

            //for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
            //    int iEdge = efp.e0 + cell;

            //    for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

            //        for (int na = 0; na < NumOfArguments; na++) {
            //            U_in[na] = Uin[na][cell, node];
            //        }

            //        //// Reference implementation
            //        //for (int d = 0; d < dimension; d++) {
            //        //        double n = efp.Normals[cell, node, d];
            //        //         fin[cell, node, d] -= 0.5 * efp.ParameterVars_IN[0][cell, node] * n;

            //        //}


            //    }
            //}
        }
        #endregion

        #region INonlineEdgeform_V Members
        void INonlinEdgeForm_V.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);

                int jCellIn = gridData.Edges.CellIndices[iEdge, 0];
                int jCellOut = gridData.Edges.CellIndices[iEdge, 1];
                double Penalty = penaltyFactor * Math.Max(cellLengthScale[jCellIn], cellLengthScale[jCellOut]);

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    // SIPG Flux Loops
                    double viscosityIn = efp.ParameterVars_IN[0][cell, node];
                    double viscosityOut = efp.ParameterVars_OUT[0][cell, node];

                    double flux = 0.0;
                    for (int d = 0; d < dimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        flux -= 0.5 * (viscosityIn * GradUin[0][cell, node, d] + viscosityOut * GradUout[0][cell, node, d]) * n;
                    }
                    flux += Math.Max(viscosityIn, viscosityOut) * (Uin[0][cell, node] - Uout[0][cell, node]) * Penalty;
                    
                    fin[cell, node] += flux;
                    fot[cell, node] -= flux;
                }
            }
        }

        void INonlinEdgeForm_V.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {
            //int NumOfCells = efp.Len;
            //Debug.Assert(fin.GetLength(0) == NumOfCells);
            //int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            //int NumOfArguments = this.ArgumentOrdering.Count;
            //Debug.Assert(NumOfArguments == Uin.Length);
            //Debug.Assert(NumOfArguments == GradUin.Length);

            //double[] U_in = new double[NumOfArguments];
            //double[] U_ot = new double[NumOfArguments];
            //double[,] GradU_in = new double[dimension, NumOfArguments];

            //for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
            //    int iEdge = efp.e0 + cell;

            //    for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

            //        for (int na = 0; na < NumOfArguments; na++) {
            //            U_in[na] = Uin[na][cell, node];

            //            for (int d = 0; d < dimension; d++) {
            //                GradU_in[d, na] = GradUin[na][cell, node, d];
            //            }
            //        }

            //        // SIPG Flux Loops
            //        //double flux = 0.0;
            //        //for (int k = 0; k < dimension; k++) {
            //        //    double nk = efp.Normals[cell, node, k];
            //        //    for (int l = 0; l < dimension; l++) {
            //        //        double nl = efp.Normals[cell, node, l];
            //        //        for (int j = 0; j < NumOfArguments; j++) {
            //        //            // consistency
            //        //            flux -= (GTensorOut[k, l, j] * GradU_in[l,j]) * nk;
            //        //            // penalty
            //        //            flux += (GTensorOut[k, l, j]) * (U_in[j] - U_ot[j]) * nl * Penalty * nk;
            //        //        }
            //        //    }
            //        //}
            //        //fin[cell, node] += flux;
            //    }
            //}
        }
        #endregion

        #region IVolumeForm Members
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV);
            }
        }

        double IVolumeForm.VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            throw new NotImplementedException();
        }
        #endregion

        #region INonLinearVolumeForm_GradV Members
        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm,
            MultidimensionalArray[] U, MultidimensionalArray[] GradU,
            MultidimensionalArray f) {
            int NumofCells = prm.Len;
            Debug.Assert(f.GetLength(0) == NumofCells);
            int NumOfNodes = f.GetLength(1); // no of nodes per cell

            for (int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    double viscosity = prm.ParameterVars[0][cell, node];
                    for (int d = 0; d < dimension; d++) {
                        //acc -= GradU[0, d] * GradV[d] * this.Nu(cpv.Xglobal, cpv.Parameters, cpv.jCell) * this.m_alpha;
                        f[cell, node, d] += viscosity * GradU[0][cell, node, d];
                    }
                }
            }
        }
        #endregion
    }
}
