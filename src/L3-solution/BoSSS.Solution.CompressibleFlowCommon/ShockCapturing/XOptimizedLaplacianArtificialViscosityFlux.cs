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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.Diffusion;
using BoSSS.Solution.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Solution.CompressibleFlowCommon.ShockCapturing {

    /// <summary>
    /// Implements the negative Laplace operator
    /// </summary>
    public class XOptimizedLaplacianArtificialViscosityFlux : INonlinear2ndOrderForm {

        public GridData GridData {
            get;
            private set;
        }

        public bool AdiabaticWall {
            get;
            set;
        }

        private readonly BoundaryCondMap<XDGHeatBcType> boundaryCondMap;

        private readonly double[] penalties;

        private readonly string ArgumentName;

        public XOptimizedLaplacianArtificialViscosityFlux(BoundaryCondMap<XDGHeatBcType> boundaryCondMap, LevelSetTracker levelSetTracker, string ArgumentVarName, double penaltySafetyFactor, double penaltyFactor, Dictionary<SpeciesId, MultidimensionalArray> inverseLengthScales) {
            this.GridData = levelSetTracker.GridDat;
            this.ArgumentName = ArgumentVarName;
            this.boundaryCondMap = boundaryCondMap;

            // Calculate penalties
            CellMask cutCells = levelSetTracker.Regions.GetCutCellMask();
            CellMask speciesAWithOutCutCells = levelSetTracker.Regions.GetSpeciesMask("A").Except(cutCells);
            CellMask speciesBWithOutCutCells = levelSetTracker.Regions.GetSpeciesMask("B").Except(cutCells);

            double[] lengthScales_A = inverseLengthScales[levelSetTracker.GetSpeciesId("A")].To1DArray();
            double[] lengthScales_B = inverseLengthScales[levelSetTracker.GetSpeciesId("B")].To1DArray();

            this.penalties = new double[lengthScales_A.Length];

            foreach (int cell in speciesAWithOutCutCells.ItemEnum) {
                this.penalties[cell] = penaltySafetyFactor * penaltyFactor * lengthScales_A[cell];
            }

            foreach (int cell in speciesBWithOutCutCells.ItemEnum) {
                this.penalties[cell] = penaltySafetyFactor * penaltyFactor * lengthScales_B[cell];
            }

            foreach (int cell in cutCells.ItemEnum) {
                this.penalties[cell] = double.NaN;
            }
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
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        TermActivationFlags IEdgeForm.InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        double IInnerEdgeForm.InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            return this.InnerEdgeFormHack(ref inp, _uA, _uB, _Grad_uA, _Grad_uB, _vA, _vB, _Grad_vA, _Grad_vB);
        }

        public double InnerEdgeFormHack(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;

            double pnlty = Math.Max(this.penalties[inp.jCellIn], this.penalties[inp.jCellOut]);
            double nuA = inp.Parameters_IN[0];
            double nuB = inp.Parameters_OUT[0];

            for (int d = 0; d < inp.D; d++) {
                Acc -= 0.5 * (nuA * _Grad_uA[0, d] + nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
                Acc -= 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            }

            //Acc *= this.m_alpha;

            double nuMax = (Math.Abs(nuA) > Math.Abs(nuB)) ? nuA : nuB;

            Acc += (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * nuMax; // penalty term

            return Acc;
        }

        double IBoundaryEdgeForm.BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            double Acc = 0.0;

            double pnlty = 2 * this.penalties[inp.jCellIn];
            double nuA = inp.Parameters_IN[0];

            XDGHeatBcType edgeType = this.boundaryCondMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case XDGHeatBcType.Dirichlet:
                    Func<double[], double, double> dirichletFunction = this.boundaryCondMap.bndFunction["u"][inp.EdgeTag];
                    double g_D = dirichletFunction(inp.X, inp.time);

                    for (int d = 0; d < inp.D; d++) {
                        double nd = inp.Normal[d];
                        Acc -= (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc -= (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }

                    Acc += nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                    break;

                case XDGHeatBcType.ZeroNeumann:
                    double g_N = 0.0;

                    Acc -= nuA * g_N * _vA;     // consistency
                    break;

                default:
                    throw new NotSupportedException();
            }

            return Acc;
        }
        #endregion

        #region IVolumeForm Members
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return TermActivationFlags.GradUxGradV;
            }
        }

        double IVolumeForm.VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                acc += GradU[0, d] * GradV[d] * cpv.Parameters[0];
            return acc;
        }
        #endregion

        #region INonlineEdgeform_GradV Members
        void INonlinInnerEdgeForm_GradV.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            InternalEdge_GradV(ref efp, Uin, Uout, GradUin, GradUout, fin, fot);
        }

        public void InternalEdge_GradV(ref EdgeFormParams efp,
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

                    for (int d = 0; d < GridData.SpatialDimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        fin[cell, node, d] -= fluxIn * n;
                        fot[cell, node, d] -= fluxOut * n;
                    }
                }
            }
        }

        /// <summary>
        /// Symmetry term
        /// </summary>
        void INonlinBoundaryEdgeForm_GradV.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int dimension = efp.GridDat.SpatialDimension;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;

                byte edgeTag = this.GridData.iGeomEdges.EdgeTags[iEdge];
                XDGHeatBcType edgeType = this.boundaryCondMap.EdgeTag2Type[edgeTag];
                Func<double[], double, double> dirichletFunction = this.boundaryCondMap.bndFunction["u"][edgeTag];

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    // Global node coordinates
                    double[] X = new double[dimension];
                    for (int i = 0; i < dimension; i++) {
                        X[i] = efp.Nodes[cell, node, i];
                    }

                    switch (edgeType) {
                        case XDGHeatBcType.Dirichlet:
                            double g_D = dirichletFunction(X, efp.time);
                            //double uJump = 0.5 * (Uin[0][cell, node] - g_D);
                            double uJump = 1.0 * (Uin[0][cell, node] - g_D);
                            double fluxIn = efp.ParameterVars_IN[0][cell, node] * uJump;

                            for (int d = 0; d < GridData.SpatialDimension; d++) {
                                double n = efp.Normals[cell, node, d];
                                fin[cell, node, d] -= fluxIn * n;
                            }
                            break;

                        default:
                            // Boundary value is zero neumann boundary, i.e. do nothing
                            break;
                    }
                }
            }
        }
        #endregion

        #region INonlineEdgeform_V Members
        void INonlinInnerEdgeForm_V.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            InternalEdge_V(ref efp, Uin, Uout, GradUin, GradUout, fin, fot);
        }

        public void InternalEdge_V(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {

            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;

                // Penalty is calculated according to the formula in background cells
                int jCellIn = this.GridData.Edges.CellIndices[iEdge, 0];
                int jCellOut = this.GridData.Edges.CellIndices[iEdge, 1];
                double penalty = Math.Max(this.penalties[jCellIn], this.penalties[jCellOut]);

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    // SIPG Flux Loops
                    double viscosityIn = efp.ParameterVars_IN[0][cell, node];
                    double viscosityOut = efp.ParameterVars_OUT[0][cell, node];

                    double flux = 0.0;
                    for (int d = 0; d < GridData.SpatialDimension; d++) {
                        double n = efp.Normals[cell, node, d];
                        flux -= 0.5 * (viscosityIn * GradUin[0][cell, node, d] + viscosityOut * GradUout[0][cell, node, d]) * n;    // Consistency term
                    }
                    flux += Math.Max(viscosityIn, viscosityOut) * (Uin[0][cell, node] - Uout[0][cell, node]) * penalty; // Penalty term

                    fin[cell, node] += flux;
                    fot[cell, node] -= flux;
                }
            }
        }

        /// <summary>
        /// Consistency and penalty term
        /// </summary>
        void INonlinBoundaryEdgeForm_V.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {

            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell
            int dimension = efp.GridDat.SpatialDimension;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;

                // Penalty is calculated according to the formula in background cells
                int jCellIn = this.GridData.Edges.CellIndices[iEdge, 0];
                double penalty = 2 * this.penalties[jCellIn];

                byte edgeTag = this.GridData.iGeomEdges.EdgeTags[iEdge];
                XDGHeatBcType edgeType = this.boundaryCondMap.EdgeTag2Type[edgeTag];
                Func<double[], double, double> dirichletFunction = this.boundaryCondMap.bndFunction["u"][edgeTag];

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    // Global node coordinates
                    double[] X = new double[dimension];
                    for (int i = 0; i < dimension; i++) {
                        X[i] = efp.Nodes[cell, node, i];
                    }

                    // SIPG Flux Loops
                    double viscosityIn = efp.ParameterVars_IN[0][cell, node];

                    switch (edgeType) {
                        case XDGHeatBcType.Dirichlet:
                            double g_D = dirichletFunction(X, efp.time);

                            double flux = 0.0;
                            for (int d = 0; d < GridData.SpatialDimension; d++) {
                                double n = efp.Normals[cell, node, d];
                                flux -= viscosityIn * GradUin[0][cell, node, d] * n;    // Consistency term
                            }
                            flux += viscosityIn * (Uin[0][cell, node] - g_D) * penalty; // Penalty term
                            fin[cell, node] += flux;
                            break;
                        default:
                            // Boundary value is zero neumann boundary, i.e. do nothing
                            break;
                    }
                }
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
                    for (int d = 0; d < GridData.SpatialDimension; d++) {
                        f[cell, node, d] += viscosity * GradU[0][cell, node, d];
                    }
                }
            }
        }
        #endregion
    }
}
