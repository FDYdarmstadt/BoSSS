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
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Implementation of a LocalTimeStepping Method with Adams-Bashforth time
    /// integration, according to: 
    /// Winters, A. R., &amp; Kopriva, D. A. (2013). High-Order Local Time Stepping
    /// on Moving DG Spectral Element Meshes. Journal of Scientific Computing.
    /// DOI: 10.1007/s10915-013-9730-z
    /// </summary>
    public class AdamsBashforthLTS : AdamsBashforth {

        /// <summary>
        /// Number of local time steps for each sub-grid
        /// [index=0] --> largest sub-grid == 1 time step
        /// </summary>
        public List<int> NumOfLocalTimeSteps {
            get;
            protected set;
        }

        /// <summary>
        /// List containing all sub-grids
        /// </summary>
        protected List<SubGrid> subgridList;

        /// <summary>
        /// Local evolvers
        /// </summary>
        protected ABevolve[] localABevolve;

        /// <summary>
        /// Number of maximum local time-steps, 
        /// i.e. ratio between first sub-grid (largest time step) and last sub-grid (smallest time step)
        /// </summary>
        protected int MaxLocalTS;

        /// <summary>
        /// Number of sub-grids which the whole grid is subdivided into 
        /// </summary>
        protected int numOfSubgrids;

        /// <summary>
        /// Saves Boundary Topology of sub-grids
        /// </summary>
        protected int[,] BoundaryTopology;

        /// <summary>
        /// Helper array of sub-grids, which only stores the boundary elements for each sub-grid
        /// </summary>
        protected SubGrid[] BoundarySgrds;

        /// <summary>
        /// Stores the local cell indices of each boundary sub-grid.
        /// Needed to avoid MPI communication errors.
        /// </summary>
        protected int[][] jSub2jCell;

        /// <summary>
        /// Helper Field, just needed for visualization of the individual sub-grids
        /// </summary>
        public DGField SgrdField;

        /// <summary>
        /// Information about the grid
        /// </summary>
        protected IGridData gridData;

        /// <summary>
        /// Bool for triggering the fluss correction
        /// </summary>
        protected bool conservative;

        private bool IBM = false;

        Queue<double>[] historyTime_Q;

        /// <summary>
        /// another ctor
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fields"></param>
        /// <param name="order"></param>
        /// <param name="NumOfSgrd"></param>
        /// <remarks>
        /// This constructor does not support parameter fields and no restriction to a sub-grid
        /// </remarks>
        public AdamsBashforthLTS(SpatialOperator spatialOp, int order, int NumOfSgrd, params DGField[] Fields)
            : this(spatialOp, new CoordinateMapping(Fields), null, order, NumOfSgrd) {
        }

        /// <summary>
        /// LTS Constructor
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if
        /// <paramref name="spatialOp"/> contains no parameters; must match the
        /// parameter field list of <paramref name="spatialOp"/>, see
        /// <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="order">LTS order</param>
        /// <param name="NumOfSgrd">Number of elements groups</param>
        /// <param name="timeStepConstraints"></param>
        /// <param name="sgrd">
        /// optional restriction to computational domain
        /// </param>
        /// <param name="conservative">Bool for triggering the fluss correction</param>
        /// <remarks>
        /// Uses the k-Mean clustering, see
        /// <see cref="BoSSS.Solution.Utils.ClusteringKmean"/>, to generate
        /// the element groups
        /// </remarks>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int NumOfSgrd, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null, bool conservative = true)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, sgrd) {

            this.numOfSubgrids = NumOfSgrd;
            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.conservative = conservative;
            this.SgrdField = new SinglePhaseField(new Basis(gridData, 0));

            // NumOfSgrd can be changed by CreateSubGrids, if less significant different element sizes than NumOfSgrd exist
            NumOfLocalTimeSteps = new List<int>(this.numOfSubgrids);

            CreateSubGrids();
            CalculateNumberOfLocalTS();
            //if (this.NumOfSgrd == 1)
            //    throw new ArgumentException("Clustering yields only to one sub-grid, LTS is not possible! Element sizes of your grid are too similar");
            localABevolve = new ABevolve[this.numOfSubgrids];

            // i == "Grid Id"
            for (int i = 0; i < subgridList.Count; i++) {
                localABevolve[i] = new ABevolve(spatialOp, Fieldsmap, Parameters, order, sgrd: subgridList[i]);
            }

            GetBoundaryTopology();

            for (int i = 0; i < this.numOfSubgrids; i++) {
                Console.WriteLine("LTS: id=" + i + " -> sub-steps=" + NumOfLocalTimeSteps[i] + " and elements=" + subgridList[i].GlobalNoOfCells);
            }

            //if (this.NumOfSgrd < 2)
            //    throw new ArgumentException("Local Time Stepping is only possible with 2 or more subgrids, but is used with " + this.NumOfSgrd);


            //RungeKuttaScheme = new RungeKutta(
            //    RungeKutta.RungeKuttaSchemes.ExplicitEuler,
            //    spatialOp,
            //    Fieldsmap,
            //    Parameters,
            //    timeStepConstraints,
            //    sgrd);
        }

        /// <summary>
        /// another ctor.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="NumOfSgrd"></param>
        /// <param name="spatialOp"></param>
        /// <param name="Fields"></param>
        public AdamsBashforthLTS(int order, int NumOfSgrd, SpatialOperator spatialOp, params DGField[] Fields)
            : this(spatialOp, new CoordinateMapping(Fields), null, order, NumOfSgrd) {
        }

        /// <summary>
        /// LTS Constructor
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if
        /// <paramref name="spatialOp"/> contains no parameters; must match the
        /// parameter field list of <paramref name="spatialOp"/>, see
        /// <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="order">LTS order</param>
        /// <param name="sgrdArray">
        /// Array with sub-grids on which different time steps are performed
        /// </param>
        /// <param name="timeStepConstraints"></param>
        /// <param name="sgrd">
        /// optional restriction to computational domain
        /// </param>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, SubGrid[] sgrdArray, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, sgrd) {

            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.subgridList = sgrdArray.ToList();
            this.numOfSubgrids = sgrdArray.Length;
            localABevolve = new ABevolve[numOfSubgrids];
            NumOfLocalTimeSteps = new List<int>(numOfSubgrids);
            CalculateNumberOfLocalTS();

            // i == "Grid Id"
            for (int i = 0; i < sgrdArray.Length; i++) {
                localABevolve[i] = new ABevolve(spatialOp, Fieldsmap, Parameters, order, sgrd: sgrdArray[i]);
            }
            GetBoundaryTopology();
        }

        /// <summary>
        /// ctor for LTS with IBM, currently under development!
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters"></param>
        /// <param name="order"></param>
        /// <param name="NumOfSgrd"></param>
        /// <param name="IBM"></param>
        /// <param name="timeStepConstraints"></param>
        /// <param name="sgrd"></param>
        /// <param name="convervative">Bool for triggering the fluss correction</param>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int NumOfSgrd, bool IBM, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null, bool convervative = true)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, sgrd) {

            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.numOfSubgrids = NumOfSgrd;
            this.IBM = IBM;
            this.conservative = convervative;

            SgrdField = new SinglePhaseField(new Basis(gridData, 0));
        }

        /// <summary>
        /// Gives the minimal euclidean distance between two vertices in each cell
        /// </summary>
        /// <returns></returns>
        protected virtual MultidimensionalArray GetSmallestDistanceInCells() {
            return gridData.iGeomCells.h_min;
        }

        /// <summary>
        /// Return a fixed metric for testing
        /// </summary>
        /// <returns></returns>
        protected virtual MultidimensionalArray GetTestMetric() {
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(gridData.iLogicalCells.NoOfLocalUpdatedCells);
            // 2 sub-grids

            //Debugger.Break();
            //if (gridData.CellPartitioning.MpiRank == 0)
            //    cellMetric[0] = 1;
            //else if (gridData.CellPartitioning.MpiRank == 1)
            //    cellMetric[0] = 0.5;
            //else if (gridData.CellPartitioning.MpiRank == 2) {
            //    cellMetric[0] = 1;
            //    cellMetric[1] = 0.5;
            //}


            //cellMetric[0] = 1;
            //cellMetric[1] = 0.1;
            //cellMetric[2] = 1;
            //cellMetric[3] = 0.05;

            //if (gridData.CellPartitioning.MpiRank == 0) {
            //    cellMetric[0] = 1;
            //    cellMetric[1] = 0.1;
            //}

            //if (gridData.CellPartitioning.MpiRank == 1) {
            //    cellMetric[0] = 1;
            //    cellMetric[1] = 0.05;
            //}

            if (gridData.CellPartitioning.MpiRank == 0) {
                cellMetric[0] = 1;
            }

            if (gridData.CellPartitioning.MpiRank == 1) {
                cellMetric[0] = 0.1;
            }

            if (gridData.CellPartitioning.MpiRank == 2) {
                cellMetric[0] = 1;
                cellMetric[1] = 0.05;
            }

            //cellMetric[0] = 1;
            //for (int i = 1; i < gridData.iGeomCells.NoOfCells; i++) {
            //    cellMetric[i] = 1.0 - gridData.CellPartitioning.MpiRank * 0.5;
            //}

            //cellMetric[0] = 1;
            //cellMetric[1] = 1;
            //cellMetric[2] = 1;
            //cellMetric[3] = 0.5;


            //cellMetric[0] = 1;
            //cellMetric[1] = 0.5;
            //cellMetric[2] = 2;

            //cellMetric[0] = 0.1;
            //cellMetric[1] = 0.1;
            //cellMetric[2] = 0.1;
            //cellMetric[3] = 0.1;
            //cellMetric[4] = 0.1;
            //cellMetric[5] = 0.1;
            //cellMetric[6] = 0.1;
            //cellMetric[7] = 1;
            //cellMetric[8] = 1;
            //cellMetric[9] = 1;

            // 3 sub-grids
            //cellMetric[0] = 1;
            //cellMetric[1] = 1;
            //cellMetric[2] = 1;
            //cellMetric[3] = 1;
            //cellMetric[4] = 0.5;
            //cellMetric[5] = 0.5;
            //cellMetric[6] = 0.5;
            //cellMetric[7] = 0.1;
            //cellMetric[8] = 0.1;
            //cellMetric[9] = 0.1;

            return cellMetric;
        }

        /// <summary>
        /// Returns a metric for the Kmeans algorithm (cell clustering)
        /// </summary>
        /// <returns></returns>
        protected virtual MultidimensionalArray GetCellMetric() {
            return GetSmallestDistanceInCells();
            //return GetTestMetric();
            //return GetTimeStepsConstraintsInCells();
        }

        /// <summary>
        /// Creates the sub-grids for the LTS algorithm
        /// </summary>
        protected void CreateSubGrids() {
            int NumOfCells = gridData.iLogicalCells.NoOfLocalUpdatedCells;

            MultidimensionalArray cellMetric = GetCellMetric();
            MultidimensionalArray means = CreateMeans(cellMetric);

            ClusteringKmean Kmean = new ClusteringKmean(cellMetric.To1DArray(), numOfSubgrids, means.To1DArray());
            // The corresponding sub-grid IDs
            int[] clustered = Kmean.Cluster();
            int[] ClusterCount = Kmean.ClusterCount;

            unsafe
            {
                int[] globalCC = new int[numOfSubgrids];
                // send = means[]
                // receive = globalMeans[]
                fixed (int* pSend = &ClusterCount[0], pRcv = &globalCC[0]) {
                    csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), numOfSubgrids, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                }
                ClusterCount = globalCC;
            }

            int counter = numOfSubgrids;
            for (int i = 0; i < numOfSubgrids; i++) {
                if (ClusterCount[i] == 0) {
                    System.Console.WriteLine("Sub-grid/Cluster " + (i + 1) + ", with mean value " + Kmean.means[i] + ", is empty and not used anymore!");
                    counter--;
                }
            }

            subgridList = new List<SubGrid>(counter);

            // Generating BitArray for all Subgrids, even for those which are empty, i.e ClusterCount == 0
            BitArray[] baMatrix = new BitArray[numOfSubgrids];
            for (int i = 0; i < numOfSubgrids; i++) {
                baMatrix[i] = new BitArray(NumOfCells);
            }

            // Filling the BitArrays
            SgrdField.Clear();
            for (int i = 0; i < NumOfCells; i++) {
                if (clustered[i] != -1) { // Happens only in the IBM case for void cells
                    baMatrix[clustered[i]][i] = true;
                    // For Debugging: Visualizes the clusters in a field
                    SgrdField.SetMeanValue(i, clustered[i] + 0 * gridData.CellPartitioning.MpiRank);
                }
            }

            // Generating the sub-grids
            int j = 0;
            for (int i = 0; i < numOfSubgrids; i++) {
                // Generating only the sub-grids which are not empty
                if (ClusterCount[i] != 0) {
                    BitArray ba = baMatrix[i];
                    subgridList.Add(new SubGrid(new CellMask(gridData, ba)));
                    j++;
                }
            }
            numOfSubgrids = counter;
        }

        /// <summary>
        /// Creates an array with an tanh spaced distributions of the mean
        /// values between maximum and minimum value of a given cell metric, 
        /// e.g minimal distance between two nodes in a cell <see cref="GridData.CellData.h_min"/>
        /// </summary>
        /// <param name="cellMetric">the given cell metric</param>
        /// <returns>Double array of size NumOfSgrd</returns>
        protected MultidimensionalArray CreateMeans(MultidimensionalArray cellMetric) {
            //MultidimensionalArray means = MultidimensionalArray.Create(NumOfSgrd);
            double h_min = cellMetric.Min(d => double.IsNaN(d) ? double.MaxValue : d); // .Where(d => !double.IsNaN(d)).ToArray().Min();
            double h_max = cellMetric.Max();
            Console.WriteLine("LTS: Create Means");
            // Getting global h_min and h_max
            ilPSP.MPICollectiveWatchDog.Watch();
            h_min = h_min.MPIMin();
            h_max = h_max.MPIMax();


            if (h_min == h_max)
                h_max += 0.1 * h_max; // Dirty Hack for IBM cases with equidistant grids


            // Tanh Spacing, which yields to more cell cluster for smaller cells
            var means = Grid1D.TanhSpacing(h_min, h_max, numOfSubgrids, 4.0, true).Reverse().ToArray();

            // Equidistant spacing, in general not the best choice
            //means = GenericBlas.Linspace(h_min, h_max, NumOfSgrd).Reverse().ToArray();

            return MultidimensionalArray.CreateWrapper(means, numOfSubgrids);
        }

        /// <summary>
        /// Checks if the mean values of the k-Mean clustering are in descending order
        /// </summary>
        /// <param name="means">Mean-values after k-Mean Clustering</param>
        protected void CheckMeans(double[] means) {
            for (int i = 0; i < means.Length - 1; i++) {
                if (means[i] < means[i + 1])
                    throw new ArgumentException("k-Mean clustering: Mean-Values are not in descending order");
            }
        }

        /// <summary>
        /// Performs one time step
        /// </summary>
        /// <param name="dt">size of time step</param>
        public override double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {
                if (localABevolve[0].HistoryChangeRate.Count >= order - 1) {

                    if (timeStepConstraints != null) {
                        dt = CalculateTimeStep();
                    }

                    double[,] CorrectionMatrix = new double[this.numOfSubgrids, this.numOfSubgrids];

                    // Saves the results at t_n
                    double[] y0 = new double[Mapping.LocalLength];
                    DGCoordinates.CopyTo(y0, 0);

                    double time0 = m_Time;
                    double time1 = m_Time + dt;

                    // evolve function
                    // evolves each sub-grid with its own time step: only one step
                    // the result is not written to m_DGCoordinates!!!
                    for (int i = 0; i < numOfSubgrids; i++) {
                        //localABevolve[i].completeBndFluxes.Clear();
                        //if (localABevolve[i].completeBndFluxes.Any(x => x != 0.0)) Console.WriteLine("Not all Bnd fluxes were used in correction step!!!");
                        localABevolve[i].Perform(dt / (double)NumOfLocalTimeSteps[i]);
                    }

                    // After evolving each cell update the time with dt_min
                    // Update AB_LTS.Time
                    m_Time = m_Time + dt / (double)NumOfLocalTimeSteps[numOfSubgrids - 1];

                    // Saves the History of DG_Coordinates for each sub-grid
                    Queue<double[]>[] historyDGC_Q = new Queue<double[]>[numOfSubgrids];
                    for (int i = 0; i < numOfSubgrids; i++) {
                        historyDGC_Q[i] = localABevolve[i].HistoryDGCoordinate;
                    }

                    // Saves DtHistory for each sub-grid
                    historyTime_Q = new Queue<double>[numOfSubgrids];
                    for (int i = 0; i < numOfSubgrids; i++) {
                        historyTime_Q[i] = localABevolve[i].HistoryTime;
                    }

                    // Perform the local time steps
                    for (int localTS = 1; localTS < MaxLocalTS; localTS++) {
                        for (int id = 1; id < numOfSubgrids; id++) {
                            //Evolve Condition: Is "ABevolve.Time" at "AB_LTS.Time"?
                            if ((localABevolve[id].Time - m_Time) < 1e-10) {
                                double localDt = dt / NumOfLocalTimeSteps[id];
                                DGCoordinates.Clear();
                                DGCoordinates.CopyFrom(historyDGC_Q[id].Last(), 0);

                                double[] interpolatedCells = InterpolateBoundaryValues(historyDGC_Q, id, localABevolve[id].Time);
                                DGCoordinates.axpy<double[]>(interpolatedCells, 1);

                                localABevolve[id].Perform(localDt);

                                m_Time = localABevolve.Min(s => s.Time);
                            }

                            // Are we at an (intermediate -) syncronization levels ?
                            // For conservatvity, we have to correct the values of the larger cell cluster
                            for (int idCoarse = 0; idCoarse < id; idCoarse++) {
                                if (Math.Abs(localABevolve[id].Time - localABevolve[idCoarse].Time) < 1e-10 &&
                                     !(Math.Abs(localABevolve[idCoarse].Time - CorrectionMatrix[idCoarse, id]) < 1e-10)) {
                                    if (conservative) {
                                        CorrectFluxes(idCoarse, id, historyDGC_Q);
                                    }
                                    CorrectionMatrix[idCoarse, id] = localABevolve[idCoarse].Time;
                                }
                            }
                        }
                    }

                    // Finalize step
                    DGCoordinates.Clear();
                    for (int id = 0; id < numOfSubgrids; id++) {
                        DGCoordinates.axpy<double[]>(historyDGC_Q[id].Last(), 1);
                    }

                    // Update time
                    m_Time = time0 + dt;

                } else {

                    double[] currentChangeRate = new double[Mapping.LocalLength];
                    double[] upDGC = new double[Mapping.LocalLength];

                    if (localABevolve[0].HistoryTime.Count == 0) {
                        for (int i = 0; i < numOfSubgrids; i++) {
                            localABevolve[i].HistoryTime.Enqueue(m_Time);
                        }
                    }

                    // Needed for the history
                    for (int i = 0; i < subgridList.Count; i++) {
                        double[] localCurrentChangeRate = new double[currentChangeRate.Length];
                        double[] edgeFlux = new double[gridData.iGeomEdges.Count * Mapping.Fields.Count];
                        localABevolve[i].ComputeChangeRate(localCurrentChangeRate, m_Time, 0, edgeFlux);
                        localABevolve[i].HistoryChangeRate.Enqueue(localCurrentChangeRate);
                        localABevolve[i].HistoryBndFluxes.Enqueue(edgeFlux);
                    }

                    dt = RungeKuttaScheme.Perform(dt);

                    DGCoordinates.CopyTo(upDGC, 0);

                    // Saves ChangeRateHistory for AB LTS
                    // Only entries for the specific sub-grid
                    for (int i = 0; i < subgridList.Count; i++) {
                        localABevolve[i].HistoryDGCoordinate.Enqueue(OrderValuesBySubgrid(subgridList[i], upDGC));
                        localABevolve[i].HistoryTime.Enqueue(RungeKuttaScheme.Time);
                    }

                    // RK is a global timeStep
                    // -> time update for all other timeStepper with rk.Time
                    m_Time = RungeKuttaScheme.Time;
                    foreach (ABevolve ab in localABevolve) {
                        ab.ResetTime(m_Time);
                    }
                }
            }
            return dt;
        }

        /// <summary>
        /// To achieve a conservative time stepping scheme, we have to correct the DG coordinates of 
        /// large interface cells
        /// </summary>
        /// <param name="coarseID">cluster ID of the large cells</param>
        /// <param name="fineID">cluster ID of the small cells</param>
        /// <param name="historyDGC_Q"></param>
        protected void CorrectFluxes(int coarseID, int fineID, Queue<double[]>[] historyDGC_Q) {
            // Gather edgeFlux data
            double[] edgeBndFluxCoarse = localABevolve[coarseID].completeBndFluxes;
            double[] edgeBndFluxFine = localABevolve[fineID].completeBndFluxes;

            int[] LocalCellIdx2SubgridIdx = subgridList[coarseID].LocalCellIndex2SubgridIndex;

            CellMask CellMaskCoarse = subgridList[coarseID].VolumeMask;

            //Only the edges of coarseID and fineID are needed
            EdgeMask unionEdgeMask = subgridList[coarseID].BoundaryEdgesMask.Intersect(subgridList[fineID].BoundaryEdgesMask);

            int cellCoarse;
            int cellFine;

            MultidimensionalArray basisScale = gridData.ChefBasis.Scaling;
            int noOfFields = Mapping.Fields.Count;

            //loop over all BoundaryEdges of the coarse sgrd
            foreach (Chunk chunk in unionEdgeMask) {
                foreach (int edge in chunk.Elements) {

                    int cell1 = gridData.iLogicalEdges.CellIndices[edge, 0];
                    int cell2 = gridData.iLogicalEdges.CellIndices[edge, 1];

                    if (LocalCellIdx2SubgridIdx[cell1] >= 0) { // <- check for MPI save operation
                        //cell1 is coarse
                        cellCoarse = cell1;
                        cellFine = cell2;
                    } else {
                        // cell2 is coarse
                        cellCoarse = cell2;
                        cellFine = cell1;
                    }

                    // Just do correction on real coarseCells, no ghost cells
                    if (CellMaskCoarse.Contains(cellCoarse)) {
                        // cell at boundary
                        // f== each field
                        // n== basis polynomial
                        for (int f = 0; f < noOfFields; f++) {
                            int n = 0; //only 0-mode is accumulated, the rest is not needed

                            int indexCoarse = Mapping.LocalUniqueCoordinateIndex(f, cellCoarse, n);
                            int indexEdge = noOfFields * edge + f;

                            double basisScaling = basisScale[cellCoarse] / basisScale[cellFine];

                            // edgefluxCorrection    = Flux from coarse -> fine                    + Flux from fine -> coarse  
                            double edgeDG_Correction = edgeBndFluxCoarse[indexEdge] * basisScaling + edgeBndFluxFine[indexEdge];


                            // Update Fluxes
                            double fluxScaling = 1.0 / localABevolve[coarseID].ABCoefficients[0];
                            localABevolve[coarseID].CurrentChangeRate[indexCoarse] += fluxScaling * edgeDG_Correction;
                            localABevolve[coarseID].HistoryBndFluxes.Last()[indexEdge] -= fluxScaling * edgeDG_Correction / basisScaling;
                            //Update local DGCoordinates
                            historyDGC_Q[coarseID].Last()[indexCoarse] -= edgeDG_Correction;

                            //We used this edge, now clear data in completeBndFluxes
                            //otherwise it will be used again in a next intermediate synchronization
                            edgeBndFluxCoarse[indexEdge] = 0.0;
                            edgeBndFluxFine[indexEdge] = 0.0;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Resets the time for all individual time stepper within the LTS algorithm, 
        /// i.e all ABevolve helper and Runge-Kutta time stepper. 
        /// It is needed after a restart, such that all time stepper restart from the 
        /// same common simulation time. 
        /// </summary>
        /// <param name="NewTime"></param>
        public override void ResetTime(double NewTime) {
            base.ResetTime(NewTime);
            RungeKuttaScheme.ResetTime(NewTime);

            foreach (var ABevolve in localABevolve) {
                ABevolve.ResetTime(NewTime);
            }

            RungeKuttaScheme.ResetTime(NewTime);
        }

        /// <summary>
        /// Calculated topology of the grid, i.e creating the boundarySubgrids
        /// for each sub-grid
        /// </summary>
        protected void GetBoundaryTopology() {
            // NumOfSgrd - 1, because largest grid (id=0) don't need a boundary cells
            BoundaryTopology = new int[numOfSubgrids - 1, gridData.iLogicalCells.NoOfLocalUpdatedCells];
            ArrayTools.SetAll(BoundaryTopology, -1);
            BoundarySgrds = new SubGrid[numOfSubgrids - 1];
            jSub2jCell = new int[numOfSubgrids - 1][];
            int[][] LocalCells2SubgridIndex = new int[numOfSubgrids][];
            BitArray[] SgrdWithGhostCells = new BitArray[numOfSubgrids];
            // prepare the calculation and  save temporarily all array which involve MPI communication
            for (int id = 0; id < numOfSubgrids; id++) {
                LocalCells2SubgridIndex[id] = subgridList[id].LocalCellIndex2SubgridIndex;
                SgrdWithGhostCells[id] = subgridList[id].VolumeMask.GetBitMaskWithExternal();
            }

            for (int id = 1; id < numOfSubgrids; id++) {
                SubGrid sgrd = subgridList[id];
                BitArray BoBA = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);

                //BitArray SgrdWithGhostCell = sgrd.VolumeMask.GetBitMaskWithExternal();
                //int[] LocalCellIndex2SubgridIndex = sgrd.LocalCellIndex2SubgridIndex;

                foreach (BoSSS.Foundation.Grid.Chunk chunk in sgrd.BoundaryEdgesMask) {
                    foreach (int edge in chunk.Elements) {

                        int cell1 = gridData.iLogicalEdges.CellIndices[edge, 0];
                        int cell2 = gridData.iLogicalEdges.CellIndices[edge, 1];

                        if (cell2 >= gridData.iLogicalCells.NoOfLocalUpdatedCells) { //special case: cell2 is "ghost-cell" at MPI border
                            if (SgrdWithGhostCells[id][cell2]) {
                                int gridId = GetSubgridIdOf(cell1, LocalCells2SubgridIndex);
                                if (gridId != -1) { // cell is not in void area of IBM
                                    BoundaryTopology[id - 1, cell1] = gridId;
                                    BoBA[cell1] = true;
                                }
                            }
                        } else if (cell1 >= 0 && cell2 >= 0 && LocalCells2SubgridIndex[id][cell1] >= 0 && LocalCells2SubgridIndex[id][cell2] < 0) {
                            //BoT[id - 1, cell2] = getSgrdIdOf(cell2, LocalCells2SubgridIndex);
                            //BoBA[cell2] = true;
                            int gridId = GetSubgridIdOf(cell2, LocalCells2SubgridIndex);
                            if (gridId != -1) { // cell is not in void area of IBM
                                BoundaryTopology[id - 1, cell2] = gridId;
                                BoBA[cell2] = true;
                            }

                        } else if (cell1 >= 0 && cell2 >= 0 && LocalCells2SubgridIndex[id][cell2] >= 0 && LocalCells2SubgridIndex[id][cell1] < 0) {
                            //BoT[id - 1, cell1] = getSgrdIdOf(cell1, LocalCells2SubgridIndex);
                            //BoBA[cell1] = true;
                            int gridId = GetSubgridIdOf(cell1, LocalCells2SubgridIndex);
                            if (gridId != -1) { // cell is not in void area of IBM
                                BoundaryTopology[id - 1, cell1] = gridId;
                                BoBA[cell1] = true;
                            }
                        }
                    }
                }
                //Creating the Boundary sub-grid
                BoundarySgrds[id - 1] = new SubGrid(new CellMask(gridData, BoBA));
                jSub2jCell[id - 1] = BoundarySgrds[id - 1].SubgridIndex2LocalCellIndex;
            }

            // Debugging the boundary topology with MPI
            //int ii = 0;
            //SgrdField.Clear();
            //foreach (SubGrid Sgrd in BoSgrd) {
            //    foreach (BoSSS.Foundation.Grid.Chunk chunk in Sgrd.VolumeMask) {
            //        foreach (int cell in chunk.Elements) {
            //            SgrdField.SetMeanValue(cell, ii + 1);
            //        }
            //    }
            //    ii++;
            //}
        }

        /// <summary>
        /// Calculates to which sub-grid the cell belongs. Each cell belongs only to one sub-grid!
        /// </summary>
        /// <param name="cell">
        /// Cell-ID</param>
        /// <param name="LocalCells2SubgridIndex">
        /// storage of all <see cref="SubGrid.LocalCellIndex2SubgridIndex"/> arrays for each LTS sub-grid</param>
        /// <returns>LTS sub-grid ID to which the cell belong</returns>
        protected int GetSubgridIdOf(int cell, int[][] LocalCells2SubgridIndex) {
            int id = -1;
            for (int i = 0; i < numOfSubgrids; i++) {
                if (LocalCells2SubgridIndex[i][cell] >= 0)
                    id = i;
            }
            return id;
        }

        /// <summary>
        /// Caluclates the particular local timesteps
        /// </summary>
        /// <returns>the largest stable timestep</returns>
        protected override double CalculateTimeStep() {
            if (timeStepConstraints.First().dtMin != timeStepConstraints.First().dtMax) {
                double[] localDts = new double[numOfSubgrids];
                for (int i = 0; i < numOfSubgrids; i++) {
                    // Use "harmonic sum" of step - sizes, see
                    // WatkinsAsthanaJameson2016 for the reasoning
                    double dt = 1.0 / timeStepConstraints.Sum(
                            c => 1.0 / c.GetGloballyAdmissibleStepSize(subgridList[i]));
                    if (dt == 0.0) {
                        throw new ArgumentException(
                            "Time-step size is exactly zero.");
                    } else if (double.IsNaN(dt)) {
                        throw new ArgumentException(
                            "Could not determine stable time-step size in Subgrid " + i + ". This indicates illegal values in some cells.");
                    }

                    // restrict timesteps
                    dt = Math.Min(dt, timeStepConstraints.First().Endtime - Time);
                    dt = Math.Min(Math.Max(dt, timeStepConstraints.First().dtMin), timeStepConstraints.First().dtMax);

                    localDts[i] = dt;
                }

                if (IBM) {
                    //localDts[NumOfSgrd - 1] *= 1.0;
                    //localDts[NumOfSgrd - 1] *= 1.0/ timeStepConstraints.First().dtFraction;
                }

                NumOfLocalTimeSteps.Clear();
                for (int i = 0; i < numOfSubgrids; i++) {
                    double fraction = localDts[0] / localDts[i];
                    //Accounting for roundoff errors
                    double eps = 1.0e-1;
                    int subSteps;
                    if (fraction > Math.Floor(fraction) + eps) {
                        subSteps = (int)Math.Ceiling(fraction);
                    } else {
                        subSteps = (int)Math.Floor(fraction);
                    }

                    NumOfLocalTimeSteps.Add(subSteps);
                }
                MaxLocalTS = NumOfLocalTimeSteps.Last();

                // Prints for each timestep the substeps
                //int ii = 0;
                //foreach (int i in NumOfLocalTimeSteps) {
                //    if (ii != 0) {
                //        Console.Write("[{0}] with {1} steps({2:0.####E-00}) and dt={3:0.####E-00} \n", ii, NumOfLocalTimeSteps[ii], (localDts[0] / (double)NumOfLocalTimeSteps[ii]), localDts[ii]);
                //    } else {
                //        Console.Write("\n[{0}] with {1} steps({2:0.####E-00}) and dt={3:0.####E-00} \n", ii, NumOfLocalTimeSteps[ii], (localDts[0] / (double)NumOfLocalTimeSteps[ii]), localDts[ii]);
                //    }
                //    ii++;
                //}

                return localDts[0];

            } else {
                MaxLocalTS = NumOfLocalTimeSteps.Last();
                double dt = timeStepConstraints.First().dtMin;
                dt = Math.Min(dt, timeStepConstraints.First().Endtime - Time);
                return dt;
            }
        }


        /// <summary>
        /// Calculates the number of sub-steps for each sub-grid
        /// </summary>
        protected void CalculateNumberOfLocalTS() {
            NumOfLocalTimeSteps.Clear();

            double[] sendHmin = new double[numOfSubgrids];
            double[] rcvHmin = new double[numOfSubgrids];

            MultidimensionalArray cellMetric = GetCellMetric();
            for (int i = 0; i < numOfSubgrids; i++) {
                double h_min = double.MaxValue;
                CellMask volumeMask = subgridList[i].VolumeMask;
                foreach (Chunk c in volumeMask) {
                    int JE = c.JE;
                    for (int j = c.i0; j < JE; j++) {
                        h_min = Math.Min(cellMetric[j], h_min);
                    }
                }
                sendHmin[i] = h_min;
            }

            // MPI to ensure that each processor has the same NumLocalTS
            unsafe
            {
                fixed (double* pSend = sendHmin, pRcv = rcvHmin) {
                    csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), numOfSubgrids, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }
            }

            int jj = 0; // Counter for the current position in the Lists
            for (int i = 0; i < numOfSubgrids; i++) {
                double fraction = rcvHmin[0] / rcvHmin[i];
                // Accounting for roundoff errors
                double eps = 1.0e-2;
                int subSteps;
                if (fraction > Math.Floor(fraction) + eps) {
                    subSteps = (int)Math.Ceiling(fraction);
                } else {
                    subSteps = (int)Math.Floor(fraction);
                }
                if (i > 0 && subSteps == NumOfLocalTimeSteps[jj - 1]) {
                    // Combine both subgrids and remove the old ones
                    SubGrid combinedSubgrid = new SubGrid(subgridList[jj].VolumeMask.Union(subgridList[jj - 1].VolumeMask));
                    subgridList.RemoveRange(jj - 1, 2);
                    subgridList.Insert(jj - 1, combinedSubgrid);
                    Console.WriteLine("Clustering leads to sub-grids which are too similar, i.e. they have the same local time step size. They are combined.");
                } else {
                    NumOfLocalTimeSteps.Add(subSteps);
                    jj++;
                }

            }
            Debug.Assert(NumOfLocalTimeSteps.Count == subgridList.Count);
            numOfSubgrids = NumOfLocalTimeSteps.Count;
            MaxLocalTS = NumOfLocalTimeSteps[numOfSubgrids - 1];
        }

        /// <summary>
        /// Interpolates the boundary elements for sub-grid of "id"
        /// </summary>
        /// <param name="historyDG">History of DG coordinates</param>
        /// <param name="id">ID of the sub-grid</param>
        /// <param name="interpolTime">Interpolation time</param>
        /// <returns>
        /// Array of the complete grid, which has only non-zero entries for
        /// the interpolated cells
        /// </returns>
        protected double[] InterpolateBoundaryValues(Queue<double[]>[] historyDG, int id, double interpolTime) {
            double[] result = new double[Mapping.LocalLength];
            SubGrid sgrd = BoundarySgrds[id - 1];

            //calculation
            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                //int cell = sgrd.SubgridIndex2LocalCellIndex[j]; --> changed to local-only operation
                int cell = jSub2jCell[id - 1][j];
                int BoundaryGridId = BoundaryTopology[id - 1, cell];
                // cell at boundary
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        double[] valueHist = new double[order];
                        int k = 0;
                        foreach (double[] histArray in historyDG[BoundaryGridId]) {
                            valueHist[k] = histArray[index];
                            k++;
                        }
                        double[] timeHistory = GetBoundaryCellTimeHistory(BoundaryGridId, cell);
                        result[index] = Interpolate(timeHistory, valueHist, interpolTime, order);
                    }
                }

            }
            return result;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="clusterId"></param>
        /// <param name="cell"></param>
        /// <returns></returns>
        virtual protected double[] GetBoundaryCellTimeHistory(int clusterId, int cell) {
            return historyTime_Q[clusterId].ToArray();
        }

        /// <summary>
        /// Interpolates a y-values for the X values with the given (x,y) pairs.
        /// Only second an third order is needed
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="X"></param>
        /// <param name="order"></param>
        /// <returns>Interpolated y value at X</returns>
        private double Interpolate(double[] x, double[] y, double X, int order) {
            switch (order) {
                case 1:
                    return y[0];

                case 2:
                    double a = (y[1] - y[0]) / (x[1] - x[0]);
                    double b = y[0] - a * x[0];
                    return a * X + b;
                case 3:
                    return y[0] * ((X - x[1]) * (X - x[2])) / ((x[0] - x[1]) * (x[0] - x[2])) +
                           y[1] * ((X - x[0]) * (X - x[2])) / ((x[1] - x[0]) * (x[1] - x[2])) +
                           y[2] * ((X - x[0]) * (X - x[1])) / ((x[2] - x[0]) * (x[2] - x[1]));
                default:
                    throw new ArgumentException("LTS is only supported for order 2 and 3, but was " + order);
            }
        }

        /// <summary>
        /// Takes a double[] with results for the global grid and gives an
        /// array with only entries for the specific sub-grid
        /// </summary>
        /// <param name="sgrd"></param>
        /// <param name="results">Result for the complete grid</param>
        /// <returns>Array with entries only for the sgrd-cells</returns>
        protected double[] OrderValuesBySubgrid(SubGrid sgrd, double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                int cell = sgrd.SubgridIndex2LocalCellIndex[j];
                // cell in sgrd
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        ordered[index] = results[index];
                    }
                }
            }
            return ordered;
        }
    }
}
