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
using BoSSS.Solution.CompressibleFlowCommon.Diffusion;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    public class OptimizedLaplacianArtificialViscosityFlux : INonlinear2ndOrderForm {

        public GridData GridData {
            private set;
            get;
        }

        public bool AdiabaticWall {
            get; set;
        }

        public MultidimensionalArray CellLengthScale {
            private set;
            get;
        }

        public double PenaltyFactor {
            private set;
            get;
        }

        private readonly int dimension;

        private readonly string ArgumentName;

        private readonly double? dirichlet;

        public OptimizedLaplacianArtificialViscosityFlux(GridData gridData, double order, int dimension, string ArgumentVarName, double? dirichlet = null) {
            this.dimension = dimension;
            this.CellLengthScale = gridData.Cells.cj;
            this.GridData = gridData;
            this.PenaltyFactor = (order + 1) * (order + (double)dimension) / (double)dimension;
            this.ArgumentName = ArgumentVarName;
            this.dirichlet = dirichlet;
        }

        #region IEquationComponent Members
        IList<string> IEquationComponent.ArgumentOrdering {
            get {
                return new string[] { ArgumentName };
            }
        }

        IList<string> IEquationComponent.ParameterOrdering {
            get {
                return new string[] { "artificialViscosity" };
            }
        }
        #endregion

        #region IEdgeForm Members
        TermActivationFlags IEdgeForm.BoundaryEdgeTerms {
            get {
                //return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
                return TermActivationFlags.AllOn;
            }
        }

        TermActivationFlags IEdgeForm.InnerEdgeTerms {
            get {
                //return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
                return TermActivationFlags.AllOn;
            }
        }

        double IEdgeForm.InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double penalty = this.PenaltyFactor * Math.Max(this.CellLengthScale[inp.jCellIn], this.CellLengthScale[inp.jCellOut]);
            return InternalEdgeForXDG_nonOptimized(ref inp, _uA, _uB, _Grad_uA, _Grad_uB, _vA, _vB, _Grad_vA, _Grad_vB, penalty);
        }

        double IEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double penalty = this.PenaltyFactor * this.CellLengthScale[inp.jCellIn];
            return BoundaryEdgeForXDG_nonOptimized(ref inp, _uA, _Grad_uA, _vA, _Grad_vA, penalty);
        }

        public double InternalEdgeForXDG_nonOptimized(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB, double penalty) {
            double Acc = 0.0;

            double nuA = inp.Parameters_IN[0];
            double nuB = inp.Parameters_OUT[0];

            for (int d = 0; d < inp.D; d++) {
                Acc += 0.5 * (nuA * _Grad_uA[0, d] + nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }

            double nuMax = (Math.Abs(nuA) > Math.Abs(nuB)) ? nuA : nuB;

            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * penalty * nuMax; // penalty term

            return Acc;
        }

        public double BoundaryEdgeForXDG_nonOptimized(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA, double penalty) {
            if (dirichlet != null) {
                double Acc = 0.0;
                if (Math.Abs(inp.X[0]) < 1e-14) {
                    double nuA = inp.Parameters_IN[0];

                    for (int d = 0; d < inp.D; d++) {
                        Acc += nuA * _Grad_uA[0, d] * _vA * inp.Normale[d];  // consistency term
                        Acc += nuA * _Grad_vA[d] * (_uA[0] - (double)dirichlet) * inp.Normale[d];  // symmetry term
                    }

                    Acc -= (_uA[0] - (double)dirichlet) * _vA * penalty * nuA; // penalty term
                }
                return Acc;
            } else {
                //Boundary value is zero neumann boundary, i.e. do nothing
                return 0.0;
            }
        }
        #endregion

        #region IVolumeForm Members
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        double IVolumeForm.VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double ret = 0;
            double viscosity = cpv.Parameters[0];
            for (int d = 0; d < dimension; d++) {
                ret += viscosity * GradU[0, d] * GradV[d];
            }
            return ret;
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
            if (dirichlet != null) {
                int NumOfCells = efp.Len;
                Debug.Assert(fin.GetLength(0) == NumOfCells);
                int NumOfNodes = fin.GetLength(1); // no of nodes per cell

                for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                    int iEdge = efp.e0 + cell;

                    for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                        double[] X = efp.NodesGlobal.ExtractSubArrayShallow(cell, node, -1).To1DArray();
                        if (Math.Abs(X[0]) <= 1e-14) {
                            double uJump = 1.0 * (Uin[0][cell, node] - (double)dirichlet);
                            double fluxIn = efp.ParameterVars_IN[0][cell, node] * uJump;

                            for (int d = 0; d < dimension; d++) {
                                double n = efp.Normals[cell, node, d];
                                fin[cell, node, d] -= fluxIn * n;
                            }
                        }
                    }
                }
            } else {
                // Boundary value is zero neumann boundary, i.e. do nothing
            }
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

            InternalEdgeForXDG(efp, Uin, Uout, GradUin, GradUout, fin, fot, null);
        }

        void INonlinEdgeForm_V.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {

            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            BoundaryEdgeForXDG(efp, Uin, GradUin, fin, null);

        }

        /// <summary>
        /// Hack in order to use the bulk AV flux as AV interface flux in XDG simulations.
        /// The question is how to pass the penalty factors in a suitable way. According to Florian,
        /// this is a general issue and has to be tackle in a global manner.
        /// Currently, we do not have any better idea.
        /// </summary>
        public void InternalEdgeForXDG(EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout, MultidimensionalArray fin, MultidimensionalArray fot, double[] penalties) {
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int NumOfCells = efp.Len;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);
                int iEdge = efp.e0 + cell;

                double penalty;
                if (penalties != null) {
                    // Use given penalties at the interface
                    Debug.Assert(iEdge == efp.e0, "Cell index should be the edge index at the interface, int cell should always be zero");
                    penalty = penalties[iEdge];
                } else {
                    // Penalty is calculated according to the formula in background cells
                    int jCellIn = this.GridData.Edges.CellIndices[iEdge, 0];
                    int jCellOut = this.GridData.Edges.CellIndices[iEdge, 1];
                    penalty = this.PenaltyFactor * Math.Max(this.CellLengthScale[jCellIn], this.CellLengthScale[jCellOut]);
                }

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                    // SIPG Flux Loops
                    double viscosityIn = efp.ParameterVars_IN[0][cell, node];
                    double viscosityOut = efp.ParameterVars_OUT[0][cell, node];

                    double flux = 0.0;
                    for (int d = 0; d < dimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        flux -= 0.5 * (viscosityIn * GradUin[0][cell, node, d] + viscosityOut * GradUout[0][cell, node, d]) * n;    // Consistency term
                    }
                    flux += Math.Max(viscosityIn, viscosityOut) * (Uin[0][cell, node] - Uout[0][cell, node]) * penalty; // Penalty term

                    fin[cell, node] += flux;
                    fot[cell, node] -= flux;
                }
            }
        }

        public void BoundaryEdgeForXDG(EdgeFormParams efp, MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin, MultidimensionalArray fin, double[] penalties) {
            if (dirichlet != null) {
                int NumOfNodes = fin.GetLength(1); // no of nodes per cell
                int NumOfCells = efp.Len;

                for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                                                                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);
                    int iEdge = efp.e0 + cell;

                    double penalty;
                    if (penalties != null) {
                        // Use given penalties at the interface
                        Debug.Assert(iEdge == efp.e0, "Cell index should be the edge index at the interface, int cell should always be zero");
                        penalty = penalties[iEdge];
                    } else {
                        // Penalty is calculated according to the formula in background cells
                        int jCellIn = this.GridData.Edges.CellIndices[iEdge, 0];
                        penalty = this.PenaltyFactor * this.CellLengthScale[jCellIn];
                    }

                    for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...
                                                                    // SIPG Flux Loops

                        double[] X = efp.NodesGlobal.ExtractSubArrayShallow(cell, node, -1).To1DArray();
                        if (Math.Abs(X[0]) <= 1e-14) {
                            double viscosityIn = efp.ParameterVars_IN[0][cell, node];

                            double flux = 0.0;
                            for (int d = 0; d < dimension; d++) {
                                double n = efp.Normals[cell, node, d];
                                flux -= viscosityIn * GradUin[0][cell, node, d] * n;    // Consistency term
                            }
                            flux += viscosityIn * (Uin[0][cell, node] - (double)dirichlet) * penalty; // Penalty term

                            fin[cell, node] += flux;
                        }
                    }
                }
            } else {
                // Boundary value is zero neumann boundary, i.e. do nothing
            }
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
