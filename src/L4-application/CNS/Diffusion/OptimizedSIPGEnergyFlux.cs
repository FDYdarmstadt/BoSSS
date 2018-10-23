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
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CNS.Diffusion {

    /// <summary>
    /// Optimized implementation of the SIPG Energy Flux for ideal gases
    /// </summary>
    class OptimizedSIPGEnergyFlux : INonlinear2ndOrderForm {

        private CNSControl config;

        private ISpeciesMap speciesMap;

        private IBoundaryConditionMap boundaryMap;

        private IGridData gridDat;

        public bool AdiabaticWall { get; set; }

        double penaltyFactor;
        
        /// <summary>
        /// Dictionary, especially needed adiabatic wall
        /// - index: edge tag
        /// - value: some switch
        /// </summary>
        protected bool[] edgeTagBool = new bool[byte.MaxValue];

        // [Dimension, Dimension, NumOfArguments]
        // [ k , l , j] --> indices according to Hartmann2008 or AnnualReport2014_SKE (i doesn't exist)
        double[,,] GTensorIn;
        double[,,] GTensorOut;
        int dimension;
        Material material;

        private MultidimensionalArray cellMetricBack;

        private Func<MultidimensionalArray> cellMetricFunc;


        /// <summary>
        /// Metric for each cell, calculation depends on IBM or non IBM case
        /// necessary for the calculation of the penalty factor
        /// </summary>
        protected MultidimensionalArray cellMetric {
            get {
                if (cellMetricBack == null) {
                    cellMetricBack = cellMetricFunc();
                }
                return cellMetricBack;
            }
        }

        public OptimizedSIPGEnergyFlux(CNSControl config, IBoundaryConditionMap boundaryMap, ISpeciesMap speciesMap, IGridData gridData, Func<MultidimensionalArray> cellMetricFunc) {
            this.config = config;
            this.speciesMap = speciesMap;
            this.boundaryMap = boundaryMap;
            this.gridDat = gridData;
            this.dimension = gridDat.SpatialDimension;
            this.material = speciesMap.GetMaterial(double.NaN);
            this.cellMetricFunc = cellMetricFunc;

            double p = new int[] { config.DensityDegree, config.MomentumDegree, config.EnergyDegree }.Max();
            penaltyFactor = config.SIPGPenaltyScaling * p * p;

            foreach (byte edgeTag in gridData.iGeomEdges.EdgeTags) {
                if (boundaryMap.EdgeTagNames[edgeTag].StartsWith("adiabaticWall", StringComparison.InvariantCultureIgnoreCase)) {
                    edgeTagBool[edgeTag] = true;
                } else {
                    edgeTagBool[edgeTag] = false;
                }
            }

            // [NumOfArguments, dimension, dimension]
            // [ k , l , j] --> indices according to Hartmann2008 or AnnualReport2014_SKE (i doesn't exist)
            GTensorIn = new double[dimension, dimension, dimension + 2];
            GTensorOut = new double[dimension, dimension, dimension + 2];
        }

        #region IEquationComponent Members
        /// <summary>
        /// <see cref="CNSEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return CNSEnvironment.PrimalArgumentOrdering;
            }
        }

        /// <summary>
        /// Empty (i.e., no parameters are used)
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return null;
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
            bool adiaWall = this.AdiabaticWall;
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            int NumOfArguments = this.ArgumentOrdering.Count;
            Debug.Assert(NumOfArguments == Uin.Length);
            Debug.Assert(NumOfArguments == GradUin.Length);
            Debug.Assert(NumOfArguments == Uout.Length);
            Debug.Assert(NumOfArguments == GradUout.Length);

            double[] U_in = new double[NumOfArguments];
            double[] U_ot = new double[NumOfArguments];

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    for (int na = 0; na < NumOfArguments; na++) {
                        U_in[na] = Uin[na][cell, node];
                        U_ot[na] = Uout[na][cell, node];
                    }

                    unsafe
                    {
                        fixed (double* pGin = GTensorIn, pGOut = GTensorOut)
                        {
                            UpdateBoundaryTensorComponent(U_in, adiaWall, dimension, pGin, material, cell);
                            UpdateBoundaryTensorComponent(U_ot, adiaWall, dimension, pGOut, material, cell);
                        }
                    }

                    //UpdateTensorComponent(U_in, dimension, GTensorIn, material);
                    //UpdateTensorComponent(U_ot, dimension, GTensorOut, material);

                    //for (int k = 0; k < dimension; k++) {    
                    //    for (int l = 0; l < dimension; l++) {
                    //        double n = efp.Normals[cell, node, l];
                    //        for (int j = 0; j < NumOfArguments; j++) {
                    //            double uJump = U_in[j] - U_ot[j];
                    //            fin[cell, node, k] -= 0.5 * GTensorIn[k, l, j] * uJump * n;
                    //            fot[cell, node, k] -= 0.5 * GTensorOut[k, l, j] * uJump * n;
                    //        }
                    //    }
                    //}
                    unsafe
                    {
                        fixed (double* pGin = GTensorIn, pGot = GTensorOut, pUin = U_in, pUot = U_ot)
                        {
                            double* pGinVar = pGin;
                            double* pGotVar = pGot;
                            for (int k = 0; k < dimension; k++) {
                                double accIn = 0;
                                double accOut = 0;
                                for (int l = 0; l < dimension; l++) {
                                    double* pUinVar = pUin;
                                    double* pUotVar = pUot;
                                    double n = efp.Normals[cell, node, l];
                                    for (int j = 0; j < NumOfArguments; j++) {
                                        double uJump = *pUinVar - *pUotVar;
                                        accIn -= 0.5 * *pGinVar * uJump * n;
                                        accOut -= 0.5 * *pGotVar * uJump * n;

                                        //Increment Pointer
                                        pGinVar++;
                                        pGotVar++;
                                        pUinVar++;
                                        pUotVar++;
                                    }
                                }
                                fin[cell, node, k] += accIn;
                                fot[cell, node, k] += accOut;
                            }
                        }
                    }
                }
            }
        }

        void INonlinEdgeForm_GradV.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            int NumOfArguments = this.ArgumentOrdering.Count;
            Debug.Assert(NumOfArguments == Uin.Length);
            Debug.Assert(NumOfArguments == GradUin.Length);

            double[] U_in = new double[NumOfArguments];
            double[] U_ot = new double[NumOfArguments];

            StateVector stateIn;
            StateVector stateOut;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    for (int na = 0; na < NumOfArguments; na++) {
                        U_in[na] = Uin[na][cell, node];
                    }
                    //Preparing for BoundaryState
                    double[] X = new double[dimension];
                    double[] Normale = new double[dimension];
                    for (int i = 0; i < dimension; i++) {
                        X[i] = efp.NodesGlobal[cell, node, i];
                        Normale[i] = efp.Normals[cell, node, i];
                    }

                    stateIn = new StateVector(U_in, speciesMap.GetMaterial(double.NaN));
                    stateOut = boundaryMap.GetBoundaryState(efp.GridDat.iGeomEdges.EdgeTags[iEdge], 0.0, X, Normale, stateIn);
                    U_ot = stateOut.ToArray();

                    bool adiabaticWall = edgeTagBool[efp.GridDat.iGeomEdges.EdgeTags[iEdge]];

                    unsafe
                    {
                        fixed (double* pGOut = GTensorOut)
                        {
                            UpdateBoundaryTensorComponent(U_ot, adiabaticWall, dimension, pGOut, material, cell);
                        }
                    }


                    //// Reference implementation
                    //for (int k = 0; k < dimension; k++) {
                    //    for (int l = 0; l < dimension; l++) {
                    //        double n = efp.Normals[cell, node, l];
                    //        for (int j = 0; j < NumOfArguments; j++) {
                    //            fin[cell, node, k] -= GTensorOut[k, l, j] * (U_in[j] - U_ot[j]) * n;
                    //        }
                    //    }
                    //}

                    unsafe
                    {
                        fixed (double* pGot = GTensorOut, gUin = U_in, gUot = U_ot)
                        {
                            double* pGotVar = pGot;
                            for (int k = 0; k < dimension; k++) {
                                for (int l = 0; l < dimension; l++) {
                                    double n = efp.Normals[cell, node, l];
                                    double* pUinVar = gUin;
                                    double* pUotVar = gUot;
                                    for (int j = 0; j < NumOfArguments; j++) {
                                        fin[cell, node, k] -= *pGotVar * (*pUinVar - *pUotVar) * n;
                                        //Increment Pointer
                                        pGotVar++;
                                        pUinVar++;
                                        pUotVar++;
                                    }
                                }
                            }
                        }
                    }

                }
            }
        }
        #endregion

        public static int EVIL_HACK_CELL_INDEX = -1;

        #region INonlineEdgeform_V Members
        void INonlinEdgeForm_V.InternalEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] Uout, MultidimensionalArray[] GradUin, MultidimensionalArray[] GradUout,
            MultidimensionalArray fin, MultidimensionalArray fot) {
            bool adiaWall = this.AdiabaticWall;
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            Debug.Assert(fot.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell

            int NumOfArguments = this.ArgumentOrdering.Count;
            Debug.Assert(NumOfArguments == Uin.Length);
            Debug.Assert(NumOfArguments == GradUin.Length);
            Debug.Assert(NumOfArguments == Uout.Length);
            Debug.Assert(NumOfArguments == GradUout.Length);


            double[] U_in = new double[NumOfArguments];
            double[] U_ot = new double[NumOfArguments];
            double[,] GradU_in = new double[dimension, NumOfArguments];
            double[,] GradU_ot = new double[dimension, NumOfArguments];

            int[,] E2Cl = gridDat.iGeomEdges.LogicalCellIndices;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);
                int jCellIn = E2Cl[iEdge, 0];
                int jCellOut = E2Cl[iEdge, 1];
                double Penalty = penaltyFactor * Math.Max(cellMetric[jCellIn], cellMetric[jCellOut]);
                if (EVIL_HACK_CELL_INDEX >= 0) {
                    Penalty = penaltyFactor * cellMetric[EVIL_HACK_CELL_INDEX];
                }

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    for (int na = 0; na < NumOfArguments; na++) {
                        U_in[na] = Uin[na][cell, node];
                        U_ot[na] = Uout[na][cell, node];

                        for (int d = 0; d < dimension; d++) {
                            GradU_in[d, na] = GradUin[na][cell, node, d];
                            GradU_ot[d, na] = GradUout[na][cell, node, d];
                        }
                    }
                    unsafe
                    {
                        fixed (double* pGin = GTensorIn, pGOut = GTensorOut)
                        {
                            UpdateBoundaryTensorComponent(U_in, adiaWall, dimension, pGin, material, cell);
                            UpdateBoundaryTensorComponent(U_ot, adiaWall, dimension, pGOut, material, cell);
                        }
                    }
                    //UpdateTensorComponent(U_in, dimension, GTensorIn, material);
                    //UpdateTensorComponent(U_ot, dimension, GTensorOut, material);

                    // SIPG Flux Loops
                    //double flux = 0.0;
                    //for (int k = 0; k < dimension; k++) {
                    //    double nk = efp.Normals[cell, node, k];
                    //    for (int l = 0; l < dimension; l++) {
                    //        double nl = efp.Normals[cell, node, l];
                    //        for (int j = 0; j < NumOfArguments; j++) {
                    //            // consistency
                    //            flux -= 0.5 * (GTensorIn[k, l, j] * GradU_in[l,j] + GTensorOut[k, l, j] * GradU_ot[l,j]) * nk;
                    //            // penalty
                    //            flux += 0.5 * (GTensorIn[k, l, j] + GTensorOut[k, l, j]) * (U_in[j] - U_ot[j]) * nl * Penalty * nk;
                    //        }
                    //    }
                    //}
                    //fin[cell, node] += flux;
                    //fot[cell, node] -= flux;
                    double flux = 0.0;
                    unsafe
                    {
                        fixed (double* pGin = GTensorIn, pGot = GTensorOut, pGradUin = GradU_in, pGradUot = GradU_ot, pUin = U_in, pUot = U_ot)
                        {
                            double* pGinVar = pGin;
                            double* pGotVar = pGot;
                            for (int k = 0; k < dimension; k++) {
                                double* pGradUinVar = pGradUin;
                                double* pGradUotVar = pGradUot;
                                double nk = efp.Normals[cell, node, k];
                                for (int l = 0; l < dimension; l++) {
                                    double factor = efp.Normals[cell, node, l] * Penalty * nk;
                                    double* pUinVar = pUin;
                                    double* pUotVar = pUot;
                                    for (int j = 0; j < NumOfArguments; j++) {
                                        // consistency
                                        flux -= 0.5 * (*pGinVar * *pGradUinVar + *pGotVar * *pGradUotVar) * nk;
                                        // penalty
                                        flux += 0.5 * (*pGinVar + *pGotVar) * (*pUinVar - *pUotVar) * factor;

                                        //Increment Pointer
                                        pGinVar++;
                                        pGotVar++;
                                        pGradUinVar++;
                                        pGradUotVar++;
                                        pUinVar++;
                                        pUotVar++;
                                    }
                                }
                            }
                        }
                    }
                    fin[cell, node] += flux;
                    fot[cell, node] -= flux;
                }
            }
        }

        void INonlinEdgeForm_V.BoundaryEdge(ref EdgeFormParams efp,
            MultidimensionalArray[] Uin, MultidimensionalArray[] GradUin,
            MultidimensionalArray fin) {
            int NumOfCells = efp.Len;
            Debug.Assert(fin.GetLength(0) == NumOfCells);
            int NumOfNodes = fin.GetLength(1); // no of nodes per cell


            int NumOfArguments = this.ArgumentOrdering.Count;
            Debug.Assert(NumOfArguments == Uin.Length);
            Debug.Assert(NumOfArguments == GradUin.Length);


            double[] U_in = new double[NumOfArguments];
            double[] U_ot = new double[NumOfArguments];
            double[,] GradU_in = new double[dimension, NumOfArguments];

            StateVector stateIn;
            StateVector stateOut;

            int[,] E2Cl = gridDat.iGeomEdges.LogicalCellIndices;

            for (int cell = 0; cell < NumOfCells; cell++) { // loop over cells...
                int iEdge = efp.e0 + cell;
                //double Penalty = penalty(gridDat.Edges.CellIndices[iEdge, 0], gridDat.Edges.CellIndices[iEdge, 1], gridDat.Cells.cj);
                int jCellIn = E2Cl[iEdge, 0];
                double Penalty = penaltyFactor * cellMetric[jCellIn];

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    for (int na = 0; na < NumOfArguments; na++) {
                        U_in[na] = Uin[na][cell, node];

                        for (int d = 0; d < dimension; d++) {
                            GradU_in[d, na] = GradUin[na][cell, node, d];
                        }
                    }
                    //Preparing for BoundaryState
                    double[] X = new double[dimension];
                    double[] Normale = new double[dimension];
                    for (int i = 0; i < dimension; i++) {
                        X[i] = efp.NodesGlobal[cell, node, i];
                        Normale[i] = efp.Normals[cell, node, i];
                    }

                    stateIn = new StateVector(U_in, speciesMap.GetMaterial(double.NaN));
                    stateOut = boundaryMap.GetBoundaryState(efp.GridDat.iGeomEdges.EdgeTags[iEdge], 0.0, X, Normale, stateIn);
                    U_ot = stateOut.ToArray();

                    bool adiabaticWall = edgeTagBool[efp.GridDat.iGeomEdges.EdgeTags[iEdge]];

                    unsafe
                    {
                        fixed (double* pGOut = GTensorOut)
                        {
                            UpdateBoundaryTensorComponent(U_ot, adiabaticWall, dimension, pGOut, material, cell);
                        }
                    }

                    // SIPG Flux Loops
                    //double flux = 0.0;
                    //for (int k = 0; k < dimension; k++) {
                    //    double nk = efp.Normals[cell, node, k];
                    //    for (int l = 0; l < dimension; l++) {
                    //        double nl = efp.Normals[cell, node, l];
                    //        for (int j = 0; j < NumOfArguments; j++) {
                    //            // consistency
                    //            flux -= (GTensorOut[k, l, j] * GradU_in[l,j]) * nk;
                    //            // penalty
                    //            flux += (GTensorOut[k, l, j]) * (U_in[j] - U_ot[j]) * nl * Penalty * nk;
                    //        }
                    //    }
                    //}
                    //fin[cell, node] += flux;
                    double flux = 0.0;
                    unsafe
                    {
                        fixed (double* pGout = GTensorOut, pGradUin = GradU_in, pUin = U_in, pUot = U_ot)
                        {
                            double* pGoutVar = pGout;
                            for (int k = 0; k < dimension; k++) {
                                double* pGradUinVar = pGradUin;
                                double nk = efp.Normals[cell, node, k];
                                for (int l = 0; l < dimension; l++) {
                                    double factor = efp.Normals[cell, node, l] * Penalty * nk;
                                    double* pUinVar = pUin;
                                    double* pUotVar = pUot;
                                    for (int j = 0; j < NumOfArguments; j++) {
                                        // consistency
                                        flux -= *pGoutVar * *pGradUinVar * nk;
                                        // penalty
                                        flux += *pGoutVar * (*pUinVar - *pUotVar) * factor;
                                        //Increment Pointer
                                        pGradUinVar++;
                                        pGoutVar++;
                                        pUinVar++;
                                        pUotVar++;
                                    }
                                }
                            }
                        }
                    }
                    fin[cell, node] += flux;
                }
            }
        }
        #endregion


        #region IVolumeForm Members
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV | TermActivationFlags.UxGradV);
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


            int NumOfArguments = this.ArgumentOrdering.Count;
            Debug.Assert(NumOfArguments == U.Length);
            Debug.Assert(NumOfArguments == GradU.Length);

            double[,] GradU_in = new double[dimension, NumOfArguments];
            double[] U_in = new double[NumOfArguments];

            for (int cell = 0; cell < NumofCells; cell++) { // loop over cells...
                int jCell = prm.j0 + cell;

                for (int node = 0; node < NumOfNodes; node++) { // loop over nodes...

                    for (int na = 0; na < NumOfArguments; na++) {
                        U_in[na] = U[na][cell, node];

                        for (int d = 0; d < dimension; d++) {
                            GradU_in[d, na] = GradU[na][cell, node, d];
                        }
                    }

                    unsafe
                    {
                        fixed (double* pGin = GTensorIn)
                        {
                            UpdateTensorComponent(U_in, dimension, pGin, material, cell);
                        }
                    }
                    //UpdateTensorComponent(U_in, dimension, GTensorIn, material);


                    //for (int k = 0; k < dimension; k++) {
                    //    for (int l = 0; l < dimension; l++) {
                    //        for (int j = 0; j < NumOfArguments; j++) {
                    //            f[cell, node, k] += GTensorIn[k, l, j] * GradU_in[l,j];
                    //        }
                    //    }
                    //}
                    unsafe
                    {
                        fixed (double* pGin = GTensorIn, pGradUin = GradU_in)
                        {
                            double* pGinVar = pGin;
                            for (int k = 0; k < dimension; k++) {
                                double acc = 0;
                                double* pGradUinVar = pGradUin;
                                for (int l = 0; l < dimension; l++) {
                                    for (int j = 0; j < NumOfArguments; j++) {
                                        acc += *pGinVar * *pGradUinVar;

                                        //Increment Pointer
                                        pGinVar++;
                                        pGradUinVar++;
                                    }
                                }
                                f[cell, node, k] += acc;
                            }
                        }
                    }

                }
            }
            #endregion
        }



        #region SIPG Flux Helper

        unsafe protected void UpdateTensorComponent(double[] state, int dimension, double* Tensor, Material material, int cellIndex) {

            double gamma = config.EquationOfState.HeatCapacityRatio;
            double alpha = config.ViscosityRatio;
            double alphaPlus43 = alpha + (4.0 / 3.0);
            double alphaPlus13 = alpha + (1.0 / 3.0);
            double alphaMinus23 = alpha - (2.0 / 3.0);

            double gamma_Pr = gamma / config.PrandtlNumber;
            double MachScaling = gamma * config.MachNumber * config.MachNumber;

            switch (dimension) {
                case 1:
                    double v1 = state[1] / state[0];
                    double E = state[2] / state[0];
                    double VelocitySquared = v1 * v1;

                    double Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    double mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];


                    //Tensor[0, 0, 0] = -((alpha43 - 1.0) * u * u + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * mu_rhoRe;
                    //Tensor[0, 0, 1] = (alpha43 - lambdaKappa_etaPr) * u * mu_rhoRe;
                    //Tensor[0, 0, 2] = lambdaKappa_etaPr * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 2] = gamma_Pr * mu_rhoRe;
                    break;
                case 2:
                    v1 = state[1] / state[0];
                    double v2 = state[2] / state[0];
                    E = state[3] / state[0];
                    VelocitySquared = v1 * v1 + v2 * v2;

                    Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];

                    // G_xx
                    //Tensor[0, 0, 0] = -((alpha43 - 1.0) * u * u + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * eta_rhoRe;
                    //Tensor[0, 0, 1] = (alpha43 - lambdaKappa_etaPr) * u * eta_rhoRe;
                    //Tensor[0, 0, 2] = (1 - lambdaKappa_etaPr) * v * eta_rhoRe;
                    //Tensor[0, 0, 3] = lambdaKappa_etaPr * eta_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 3] = gamma_Pr * mu_rhoRe;


                    // G_xy
                    //Tensor[0, 1, 0] = -(alpha43 - 1.0) * u * v * eta_rhoRe;
                    //Tensor[0, 1, 1] = v * eta_rhoRe;
                    //Tensor[0, 1, 2] = alpha23 * u * eta_rhoRe;
                    ////Tensor[0, 1, 3] = 0 * eta_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 0] = -MachScaling *  alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 1] = MachScaling * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 2] = MachScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0 * eta_rhoRe;


                    // G_yx
                    //Tensor[1, 0, 0] = -(alpha43 - 1.0) * u * v * eta_rhoRe;
                    //Tensor[1, 0, 1] = alpha23 * v * eta_rhoRe;
                    //Tensor[1, 0, 2] = u * eta_rhoRe;
                    ////Tensor[1, 0, 3] = 0 * eta_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 1] = MachScaling * alphaMinus23 * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 2] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0 * eta_rhoRe;


                    // G_yy
                    //Tensor[1, 1, 0] = -((alpha43 - 1.0) * v * v + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * eta_rhoRe;
                    //Tensor[1, 1, 1] = (1 - lambdaKappa_etaPr) * u * eta_rhoRe;
                    //Tensor[1, 1, 2] = (alpha43 - lambdaKappa_etaPr) * v * eta_rhoRe;
                    //Tensor[1, 1, 3] = lambdaKappa_etaPr * eta_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 0] = -(MachScaling * (alphaPlus13 * v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 2] = MachScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 3] = gamma_Pr * mu_rhoRe;
                    break;
                case 3:
                    v1 = state[1] / state[0];
                    v2 = state[2] / state[0];
                    double v3 = state[3] / state[0];
                    E = state[4] / state[0];
                    VelocitySquared = v1 * v1 + v2 * v2 + v3 * v3;

                    Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];

                    // G_xx
                    Tensor[0 * 3 * 5 + 0 * 5 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 3] = MachScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    // G_xy
                    Tensor[0 * 3 * 5 + 1 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 1 * 5 + 1] = MachScaling * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 1 * 5 + 2] = MachScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                    //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                    // G_xz
                    Tensor[0 * 3 * 5 + 2 * 5 + 0] = - MachScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 2 * 5 + 1] = MachScaling * v3 * mu_rhoRe;
                    //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                    Tensor[0 * 3 * 5 + 2 * 5 + 3] = MachScaling * alphaMinus23 * v1 * mu_rhoRe;
                    //Tensor[0, 2, 4] = 0 * eta_rhoRe;


                    // G_yx
                    Tensor[1 * 3 * 5 + 0 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 0 * 5 + 1] = MachScaling * alphaMinus23 * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 0 * 5 + 2] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                    //Tensor[1, 0, 4] = 0 * eta_rhoRe;


                    // G_yy
                    Tensor[1 * 3 * 5 + 1 * 5 + 0] = -(MachScaling * (alphaPlus13* v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 2] = MachScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 3] = MachScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    // G_yz
                    Tensor[1 * 3 * 5 + 2 * 5 + 0] = -MachScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                    Tensor[1 * 3 * 5 + 2 * 5 + 2] = MachScaling * v3 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 2 * 5 + 3] = MachScaling * alphaMinus23 * v2 * mu_rhoRe;
                    //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                    // G_zx
                    Tensor[2 * 3 * 5 + 0 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 0 * 5 + 1] = MachScaling * alphaMinus23 * v3 * mu_rhoRe;
                    //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                    Tensor[2 * 3 * 5 + 0 * 5 + 3] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                    // G_zy
                    Tensor[2 * 3 * 5 + 1 * 5 + 0] = -MachScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                    Tensor[2 * 3 * 5 + 1 * 5 + 2] = MachScaling * alphaMinus23 * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 1 * 5 + 3] = MachScaling * v2 * mu_rhoRe;
                    //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                    // G_zz
                    Tensor[2 * 3 * 5 + 2 * 5 + 0] = -(MachScaling * (alphaPlus13 * v3 * v3 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 3] = MachScaling * (alphaPlus43 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    break;
                default:
                    throw new ArgumentOutOfRangeException("dimension");
            }
        }


        unsafe protected void UpdateBoundaryTensorComponent(double[] state, bool adiabaticWall, int dimension, double* Tensor, Material material, int cellIndex) {

            double gamma = config.EquationOfState.HeatCapacityRatio;
            double alpha = config.ViscosityRatio;
            double alphaPlus43 = alpha + (4.0 / 3.0);
            double alphaPlus13 = alpha + (1.0 / 3.0);
            double alphaMius23 = alpha - (2.0 / 3.0);

            double gamma_Pr = adiabaticWall ? 0 : gamma / config.PrandtlNumber;
            double MachScaling = gamma * config.MachNumber * config.MachNumber;

            switch (dimension) {
                case 1:
                    double v1 = state[1] / state[0];
                    double E = state[2] / state[0];
                    double VelocitySquared = v1 * v1;

                    double Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    double mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];                

                    //Tensor[0, 0, 0] = -((alpha43 - 1.0) * u * u + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * eta_rhoRe;
                    //Tensor[0, 0, 1] = (alpha43 - lambdaKappa_etaPr) * u * eta_rhoRe;
                    //Tensor[0, 0, 2] = lambdaKappa_etaPr * eta_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 2] = gamma_Pr * mu_rhoRe;
                    break;
                case 2:
                    v1 = state[1] / state[0];
                    double v2 = state[2] / state[0];
                    E = state[3] / state[0];
                    VelocitySquared = v1 * v1 + v2 * v2;

                    Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];

                    // G_xx
                    //Tensor[0, 0, 0] = -((alpha43 - 1.0) * u * u + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * eta_rhoRe;
                    //Tensor[0, 0, 1] = (alpha43 - lambdaKappa_etaPr) * u * eta_rhoRe;
                    //Tensor[0, 0, 2] = (1 - lambdaKappa_etaPr) * v * eta_rhoRe;
                    //Tensor[0, 0, 3] = lambdaKappa_etaPr * eta_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 0 * 4 + 3] = gamma_Pr * mu_rhoRe;


                    // G_xy
                    //Tensor[0, 1, 0] = -(alpha43 - 1.0) * u * v * eta_rhoRe;
                    //Tensor[0, 1, 1] = v * eta_rhoRe;
                    //Tensor[0, 1, 2] = alpha23 * u * eta_rhoRe;
                    ////Tensor[0, 1, 3] = 0 * eta_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 1] = MachScaling * v2 * mu_rhoRe;
                    Tensor[0 * 2 * 4 + 1 * 4 + 2] = MachScaling * alphaMius23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0 * eta_rhoRe;


                    // G_yx
                    //Tensor[1, 0, 0] = -(alpha43 - 1.0) * u * v * eta_rhoRe;
                    //Tensor[1, 0, 1] = alpha23 * v * eta_rhoRe;
                    //Tensor[1, 0, 2] = u * eta_rhoRe;
                    ////Tensor[1, 0, 3] = 0 * eta_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 0] = -MachScaling *  alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 1] = MachScaling * alphaMius23 * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 0 * 4 + 2] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0 * eta_rhoRe;


                    // G_yy
                    //Tensor[1, 1, 0] = -((alpha43 - 1.0) * v * v + TwoKineticEnergy + lambdaKappa_etaPr * (E - TwoKineticEnergy)) * eta_rhoRe;
                    //Tensor[1, 1, 1] = (1 - lambdaKappa_etaPr) * u * eta_rhoRe;
                    //Tensor[1, 1, 2] = (alpha43 - lambdaKappa_etaPr) * v * eta_rhoRe;
                    //Tensor[1, 1, 3] = lambdaKappa_etaPr * eta_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 0] = -(MachScaling * (alphaPlus13 * v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 2] = MachScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1 * 2 * 4 + 1 * 4 + 3] = gamma_Pr * mu_rhoRe;
                    break;
                case 3:
                    v1 = state[1] / state[0];
                    v2 = state[2] / state[0];
                    double v3 = state[3] / state[0];
                    E = state[4] / state[0];
                    VelocitySquared = v1 * v1 + v2 * v2 + v3 * v3;

                    Viscosity = material.ViscosityLaw.GetViscosity((gamma - 1.0) * (E - 0.5 * MachScaling * VelocitySquared), cellIndex);
                    mu_rhoRe = Viscosity / config.ReynoldsNumber / state[0];

                    // G_xx
                    Tensor[0 * 3 * 5 + 0 * 5 + 0] = -(MachScaling * (alphaPlus13 * v1 * v1 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 1] = MachScaling * (alphaPlus43 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 3] = MachScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 0 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    // G_xy
                    Tensor[0 * 3 * 5 + 1 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 1 * 5 + 1] = MachScaling * v2 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 1 * 5 + 2] = MachScaling * alphaMius23 * v1 * mu_rhoRe;
                    //Tensor[0, 1, 3] = 0 * eta_rhoRe;
                    //Tensor[0, 1, 4] = 0 * eta_rhoRe;

                    // G_xz
                    Tensor[0 * 3 * 5 + 2 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[0 * 3 * 5 + 2 * 5 + 1] = MachScaling * v3 * mu_rhoRe;
                    //Tensor[0, 2, 2] = 0 * eta_rhoRe;
                    Tensor[0 * 3 * 5 + 2 * 5 + 3] = MachScaling * alphaMius23 * v1 * mu_rhoRe;
                    //Tensor[0, 2, 4] = 0 * eta_rhoRe;


                    // G_yx
                    Tensor[1 * 3 * 5 + 0 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 0 * 5 + 1] = MachScaling * alphaMius23 * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 0 * 5 + 2] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[1, 0, 3] = 0 * eta_rhoRe;
                    //Tensor[1, 0, 4] = 0 * eta_rhoRe;


                    // G_yy
                    Tensor[1 * 3 * 5 + 1 * 5 + 0] = -(MachScaling * (alphaPlus13 * v2 * v2 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 2] = MachScaling * (alphaPlus43 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 3] = MachScaling * (1 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 1 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    // G_yz
                    Tensor[1 * 3 * 5 + 2 * 5 + 0] = -MachScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[1, 2, 1] = 0 * eta_rhoRe;
                    Tensor[1 * 3 * 5 + 2 * 5 + 2] = MachScaling * v3 * mu_rhoRe;
                    Tensor[1 * 3 * 5 + 2 * 5 + 3] = MachScaling * alphaMius23 * v2 * mu_rhoRe;
                    //Tensor[1, 2, 4] = 0 * eta_rhoRe;

                    // G_zx
                    Tensor[2 * 3 * 5 + 0 * 5 + 0] = -MachScaling * alphaPlus13 * v1 * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 0 * 5 + 1] = MachScaling * alphaMius23 * v3 * mu_rhoRe;
                    //Tensor[2, 0, 2] = 0 * eta_rhoRe;
                    Tensor[2 * 3 * 5 + 0 * 5 + 3] = MachScaling * v1 * mu_rhoRe;
                    //Tensor[2, 0, 4] = 0 * eta_rhoRe;

                    // G_zy
                    Tensor[2 * 3 * 5 + 1 * 5 + 0] = -MachScaling * alphaPlus13 * v2 * v3 * mu_rhoRe;
                    //Tensor[2, 1, 1] = 0 * eta_rhoRe;
                    Tensor[2 * 3 * 5 + 1 * 5 + 2] = MachScaling * alphaMius23 * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 1 * 5 + 3] = MachScaling * v2 * mu_rhoRe;
                    //Tensor[2, 1, 4] = 0 * eta_rhoRe;

                    // G_zz
                    Tensor[2 * 3 * 5 + 2 * 5 + 0] = -(MachScaling * (alphaPlus13 * v3 * v3 + VelocitySquared) + gamma_Pr * (E - MachScaling * VelocitySquared)) * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 1] = MachScaling * (1 - gamma_Pr) * v1 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 2] = MachScaling * (1 - gamma_Pr) * v2 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 3] = MachScaling * (alphaPlus43 - gamma_Pr) * v3 * mu_rhoRe;
                    Tensor[2 * 3 * 5 + 2 * 5 + 4] = gamma_Pr * mu_rhoRe;

                    break;
                default:
                    throw new ArgumentOutOfRangeException("dimension");
            }
        }

        #endregion

    }
}
