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
using BoSSS.Foundation.IO;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
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

        //protected Clusterer.Clustering CurrentClustering;

        /// <summary>
        /// Local evolvers
        /// </summary>
        protected ABevolve[] localABevolve;

        /// <summary>
        /// Number of maximum local time-steps, 
        /// i.e. ratio between first sub-grid (largest time step) and last sub-grid (smallest time step)
        /// </summary>
        protected int maxNumOfLocalTimeSteps;

        /// <summary>
        /// Number of sub-grids which the whole grid is subdivided into 
        /// </summary>
        protected int numOfClusters;

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
        /// Information about the grid
        /// </summary>
        protected IGridData gridData;

        /// <summary>
        /// Bool for triggering the flux correction
        /// </summary>
        protected bool fluxCorrection;

        /// <summary>
        /// Bool for triggering the adaptive LTS version
        /// </summary>
        protected bool adaptive;

        /// <summary>
        /// Constant number of sub-grids specified by the user
        /// </summary>
        static int numOfSubgridsInit;

        /// <summary>
        /// localABevolve objects from the previous time step needed for copying the histories
        /// </summary>
        protected ABevolve[] localABevolvePrevious;

        /// <summary>
        /// Interval for reclustering when LTS is used in adpative mode
        /// 0: standard LTS, e.g., 10: <see cref="Clusterer.CheckForNewClustering(Clusterer.Clustering, Clusterer.Clustering)"/> is called in every tenth time step
        /// </summary>
        protected int reclusteringInterval;

        /// <summary>
        /// <see cref="Clusterer"/> based on the time step constraint in each cell. Devides the grids into sub-grids.
        /// </summary>
        protected Clusterer clusterer;

        Queue<double>[] historyTime_Q;

        private int timeStepCount;

        public Clusterer.Clustering CurrentClustering {
            get;
            protected set;
        }

        //################# Hack for saving to database in every (A)LTS sub-step
        private Action<TimestepNumber, double> saveToDBCallback;
        //################# Hack for saving to database in every (A)LTS sub-step

        /// <summary>
        /// Standard constructor for the (adaptive) local time stepping algorithm
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap">Coordinate mapping for the variable fields</param>
        /// <param name="Parameters">optional parameter fields, can be null if <paramref name="spatialOp"/> contains no parameters; must match the parameter field list of <paramref name="spatialOp"/>, see <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/></param>
        /// <param name="order">LTS/AB order</param>
        /// <param name="numOfSubgrids">Amount of sub-grids/clusters to be used for LTS</param>
        /// <param name="timeStepConstraints">Time step constraints for later usage as metric</param>
        /// <param name="subGrid">Sub-grids, e.g., from previous time steps</param>
        /// <param name="fluxCorrection">Bool for triggering the fluss correction</param>
        /// <param name="reclusteringInterval">Interval for potential reclustering</param>
        /// <remarks>Uses the k-Mean clustering, see <see cref="BoSSS.Solution.Utils.Kmeans"/>, to generate the element groups</remarks>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int numOfSubgrids, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid subGrid = null, bool fluxCorrection = true, int reclusteringInterval = 0, Action<TimestepNumber, double> saveToDBCallback = null)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, subGrid) {

            if (reclusteringInterval != 0) {
                numOfSubgridsInit = numOfSubgrids;
                this.timeStepCount = 1;
                this.adaptive = true;
                RungeKuttaScheme.OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);
            }

            this.reclusteringInterval = reclusteringInterval;
            this.numOfClusters = numOfSubgrids;
            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.fluxCorrection = fluxCorrection;

            NumOfLocalTimeSteps = new List<int>(this.numOfClusters);

            clusterer = new Clusterer(this.gridData, this.timeStepConstraints);
            CurrentClustering = clusterer.CreateClustering(this.numOfClusters, subGrid);    // Might remove clusters when their centres are too close

            for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                Console.WriteLine("(A)LTS: id=" + i + " -> elements=" + CurrentClustering.Clusters[i].GlobalNoOfCells);
            }

            CurrentClustering = CalculateNumberOfLocalTS(CurrentClustering); // Might remove clusters when time step sizes are too similar

            localABevolve = new ABevolve[this.numOfClusters];

            // i == "Grid Id"
            for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                localABevolve[i] = new ABevolve(spatialOp, Fieldsmap, Parameters, order, adaptive: this.adaptive, sgrd: CurrentClustering.Clusters[i]);
                localABevolve[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);
            }

            GetBoundaryTopology();

            for (int i = 0; i < this.numOfClusters; i++) {
                Console.WriteLine("(A)LTS: id=" + i + " -> sub-steps=" + NumOfLocalTimeSteps[i] + " and elements=" + CurrentClustering.Clusters[i].GlobalNoOfCells);
            }

            // Saving time steps in subgrids
            //this.saveToDBCallback = saveToDBCallback;
        }

        /// <summary>
        /// Constructor for LTS with IBM, currently under development!
        /// </summary>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int NumOfSgrd, bool IBM, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null, bool fluxCorrection = true)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, sgrd) {

            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.numOfClusters = NumOfSgrd;
            this.fluxCorrection = fluxCorrection;
        }

        /// <summary>
        /// Performs one time step
        /// </summary>
        /// <param name="dt">Time step size that equals -1, if no fixed time step is prescribed</param>
        public override double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {

                if (localABevolve[0].HistoryChangeRate.Count >= order - 1) {
                    bool reclustered = false;
                    if (adaptive) {
                        if (timeStepCount % reclusteringInterval == 0) {
                            // Necessary in order to use the number of sub-grids specified by the user for the reclustering in each time step
                            // Otherwise the value could be changed by the constructor of the parent class (AdamsBashforthLTS.cs) --> CreateSubGrids()
                            numOfClusters = numOfSubgridsInit;

                            // Fix for update problem of artificial viscosity
                            RaiseOnBeforComputechangeRate(Time, dt);

                            var oldClustering = CurrentClustering;
                            CurrentClustering = clusterer.CreateClustering(numOfClusters, null);
                            CurrentClustering = CalculateNumberOfLocalTS(CurrentClustering); // Might remove sub-grids when time step sizes are too similar
                            reclustered = clusterer.CheckForNewClustering(oldClustering, CurrentClustering);

                            // After the intitial phase, activate adaptive mode for all ABevolve objects
                            foreach (ABevolve abE in localABevolve)
                                abE.adaptive = true;

                            if (reclustered) {
                                // Store all localAbevolve objects from the last time step for copying the histories
                                ShortenHistories(localABevolve);
                                localABevolvePrevious = localABevolve;

                                // Create array of Abevolve objects based on the new clustering
                                localABevolve = new ABevolve[this.numOfClusters];

                                for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                                    localABevolve[i] = new ABevolve(Operator, Mapping, ParameterMapping, order, adaptive: true, sgrd: CurrentClustering.Clusters[i]);
                                    localABevolve[i].ResetTime(m_Time);
                                    localABevolve[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);
                                }

                                CopyHistoriesOfABevolver();

                            } else {
                                //Console.WriteLine("#####Clustering has NOT changed in timestep{0}#####", timeStepCount);
                            }


                            GetBoundaryTopology();
                        }
                    }

                    if (timeStepConstraints != null) {
                        dt = CalculateTimeStep();
                    }

                    // Number of substeps could have changed
                    if (reclustered) {
                        for (int i = 0; i < this.numOfClusters; i++) {
                            Console.WriteLine("LTS: id=" + i + " -> sub-steps=" + NumOfLocalTimeSteps[i] + " and elements=" + CurrentClustering.Clusters[i].GlobalNoOfCells);
                        }
                    }

                    double[,] CorrectionMatrix = new double[this.numOfClusters, this.numOfClusters];

                    // Saves the results at t_n
                    double[] y0 = new double[Mapping.LocalLength];
                    DGCoordinates.CopyTo(y0, 0);

                    double time0 = m_Time;
                    double time1 = m_Time + dt;

                    TimestepNumber subTimestep = new TimestepNumber(timeStepCount - 1);

                    // Evolves each sub-grid with its own time step: only one step
                    // the result is not written to m_DGCoordinates!!!
                    for (int i = 0; i < numOfClusters; i++) {
                        //localABevolve[i].completeBndFluxes.Clear();
                        //if (localABevolve[i].completeBndFluxes.Any(x => x != 0.0)) Console.WriteLine("Not all Bnd fluxes were used in correction step!!!");
                        localABevolve[i].Perform(dt / (double)NumOfLocalTimeSteps[i]);
                    }

                    // After evolving each cell update the time with dt_min
                    // Update AB_LTS.Time
                    m_Time = m_Time + dt / (double)NumOfLocalTimeSteps[numOfClusters - 1];

                    if (saveToDBCallback != null) {
                        subTimestep = subTimestep.NextIteration();
                        saveToDBCallback(subTimestep, m_Time);
                    }

                    // Saves the History of DG_Coordinates for each sub-grid
                    Queue<double[]>[] historyDGC_Q = new Queue<double[]>[numOfClusters];
                    for (int i = 0; i < numOfClusters; i++) {
                        historyDGC_Q[i] = localABevolve[i].HistoryDGCoordinate;
                    }

                    if (!adaptive) {
                        // Saves DtHistory for each sub-grid
                        historyTime_Q = new Queue<double>[numOfClusters];
                        for (int i = 0; i < numOfClusters; i++) {
                            historyTime_Q[i] = localABevolve[i].HistoryTime;
                        }
                    }

                    // Perform the local time steps
                    for (int localTS = 1; localTS < maxNumOfLocalTimeSteps; localTS++) {
                        for (int id = 1; id < numOfClusters; id++) {
                            //Evolve Condition: Is "ABevolve.Time" at "AB_LTS.Time"?
                            if ((localABevolve[id].Time - m_Time) < 1e-10) {
                                double localDt = dt / NumOfLocalTimeSteps[id];

                                foreach (Chunk chunk in CurrentClustering.Clusters[id].VolumeMask) {
                                    foreach (int cell in chunk.Elements) {
                                        // f == each field
                                        // n == basis polynomial
                                        foreach (DGField f in Mapping.Fields) {
                                            for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                                                int coordinateIndex = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                                                DGCoordinates[coordinateIndex] = historyDGC_Q[id].Last()[coordinateIndex];
                                            }
                                        }
                                    }
                                }

                                Dictionary<int, double> myDic = InterpolateBoundaryValues(historyDGC_Q, id, localABevolve[id].Time);

                                foreach (KeyValuePair<int, double> kvp in myDic) {
                                    DGCoordinates[kvp.Key] = kvp.Value;
                                }

                                localABevolve[id].Perform(localDt);

                                m_Time = localABevolve.Min(s => s.Time);

                                if (saveToDBCallback != null) {
                                    subTimestep = subTimestep.NextIteration();
                                    saveToDBCallback(subTimestep, m_Time);
                                }
                            }

                            // Are we at an (intermediate-) syncronization levels?
                            // For conservatvity, we have to correct the values of the larger cell cluster
                            if (fluxCorrection) {
                                for (int idCoarse = 0; idCoarse < id; idCoarse++) {
                                    if (Math.Abs(localABevolve[id].Time - localABevolve[idCoarse].Time) < 1e-10 &&
                                         !(Math.Abs(localABevolve[idCoarse].Time - CorrectionMatrix[idCoarse, id]) < 1e-10)) {
                                        if (fluxCorrection) {
                                            CorrectFluxes(idCoarse, id, historyDGC_Q);
                                        }
                                        CorrectionMatrix[idCoarse, id] = localABevolve[idCoarse].Time;
                                    }
                                }
                            }
                        }
                    }

                    // Finalize step
                    // Use unmodified values in history of DGCoordinates (DGCoordinates could have been modified by
                    // InterpolateBoundaryValues, should be resetted afterwards) 
                    DGCoordinates.Clear();
                    for (int id = 0; id < numOfClusters; id++) {
                        DGCoordinates.axpy<double[]>(historyDGC_Q[id].Last(), 1);
                    }

                    // Update time
                    m_Time = time0 + dt;

                } else {

                    double[] currentChangeRate = new double[Mapping.LocalLength];
                    double[] upDGC = new double[Mapping.LocalLength];

                    // Save history: Time
                    if (adaptive) {
                        for (int i = 0; i < numOfClusters; i++) {
                            double[] currentTime = new double[localABevolve[i].sgrd.LocalNoOfCells];
                            for (int j = 0; j < currentTime.Length; j++) {
                                currentTime[j] = m_Time;
                            }
                            localABevolve[i].historyTimePerCell.Enqueue(currentTime);
                        }
                    } else {
                        if (localABevolve[0].HistoryTime.Count == 0) {
                            for (int i = 0; i < numOfClusters; i++) {
                                localABevolve[i].HistoryTime.Enqueue(m_Time);
                            }
                        }
                    }

                    // Needed for the history
                    for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
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
                    for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                        localABevolve[i].HistoryDGCoordinate.Enqueue(OrderValuesBySubgrid(CurrentClustering.Clusters[i], upDGC));
                        if (!adaptive)
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

            if (adaptive)
                timeStepCount++;

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

            int[] LocalCellIdx2SubgridIdx = CurrentClustering.Clusters[coarseID].LocalCellIndex2SubgridIndex;

            CellMask CellMaskCoarse = CurrentClustering.Clusters[coarseID].VolumeMask;

            //Only the edges of coarseID and fineID are needed
            EdgeMask unionEdgeMask = CurrentClustering.Clusters[coarseID].BoundaryEdgesMask.Intersect(CurrentClustering.Clusters[fineID].BoundaryEdgesMask);

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
        /// <param name="NewTime">Time to be set</param>
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
            BoundaryTopology = new int[numOfClusters - 1, gridData.iLogicalCells.NoOfLocalUpdatedCells];
            ArrayTools.SetAll(BoundaryTopology, -1);
            BoundarySgrds = new SubGrid[numOfClusters - 1];
            jSub2jCell = new int[numOfClusters - 1][];
            int[][] LocalCells2SubgridIndex = new int[numOfClusters][];
            BitArray[] SgrdWithGhostCells = new BitArray[numOfClusters];
            // prepare the calculation and  save temporarily all array which involve MPI communication
            for (int id = 0; id < numOfClusters; id++) {
                LocalCells2SubgridIndex[id] = CurrentClustering.Clusters[id].LocalCellIndex2SubgridIndex;
                SgrdWithGhostCells[id] = CurrentClustering.Clusters[id].VolumeMask.GetBitMaskWithExternal();
            }

            for (int id = 1; id < numOfClusters; id++) {
                SubGrid sgrd = CurrentClustering.Clusters[id];
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
        /// <param name="cell">Cell-ID</param>
        /// <param name="LocalCells2SubgridIndex">Storage of all <see cref="SubGrid.LocalCellIndex2SubgridIndex"/>
        /// arrays for each LTS sub-grid</param>
        /// <returns>LTS sub-grid ID which the cell belong to</returns>
        protected int GetSubgridIdOf(int cell, int[][] LocalCells2SubgridIndex) {
            int id = -1;
            for (int i = 0; i < numOfClusters; i++) {
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
                double[] localDts = new double[numOfClusters];
                for (int i = 0; i < numOfClusters; i++) {
                    // Use "harmonic sum" of step - sizes, see
                    // WatkinsAsthanaJameson2016 for the reasoning
                    double dt = 1.0 / timeStepConstraints.Sum(
                            c => 1.0 / c.GetGloballyAdmissibleStepSize(CurrentClustering.Clusters[i]));
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

                NumOfLocalTimeSteps.Clear();
                for (int i = 0; i < numOfClusters; i++) {
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
                maxNumOfLocalTimeSteps = NumOfLocalTimeSteps.Last();

                // Prints the substeps for each timestep 
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
                maxNumOfLocalTimeSteps = NumOfLocalTimeSteps.Last();
                double dt = timeStepConstraints.First().dtMin;
                dt = Math.Min(dt, timeStepConstraints.First().Endtime - Time);
                return dt;
            }
        }


        /// <summary>
        /// Calculates the number of sub-steps for each sub-grid
        /// </summary>
        protected Clusterer.Clustering CalculateNumberOfLocalTS(Clusterer.Clustering clustering) {
            NumOfLocalTimeSteps.Clear();

            double[] sendHmin = new double[clustering.NumberOfClusters];
            double[] rcvHmin = new double[clustering.NumberOfClusters];

            MultidimensionalArray cellMetric = clusterer.GetCellMetric(clustering.SubGrid);
            for (int i = 0; i < clustering.NumberOfClusters; i++) {
                double h_min = double.MaxValue;
                CellMask volumeMask = clustering.Clusters[i].VolumeMask;
                foreach (Chunk c in volumeMask) {
                    int JE = c.JE;
                    for (int j = c.i0; j < JE; j++) {
                        h_min = Math.Min(cellMetric[clustering.SubGrid.LocalCellIndex2SubgridIndex[j]], h_min);
                    }
                }
                sendHmin[i] = h_min;
            }

            // MPI to ensure that each processor has the local time step sizes
            unsafe {
                fixed (double* pSend = sendHmin, pRcv = rcvHmin) {
                    csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), clustering.NumberOfClusters, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }
            }

            int jj = 0; // Counter for the current position in the Lists
            List<SubGrid> newClusters = new List<SubGrid>();
            for (int i = 0; i < clustering.NumberOfClusters; i++) {
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
                    SubGrid combinedSubgrid = new SubGrid(clustering.Clusters[jj].VolumeMask.Union(clustering.Clusters[jj - 1].VolumeMask));
                    newClusters.RemoveAt(jj - 1);   // (can be done better without removing subgrids)
                    newClusters.Add(combinedSubgrid);
                    Console.WriteLine("CalculateNumberOfLocalTS: Clustering leads to sub-grids which are too similar, i.e. they have the same local time step size. They are combined.");
                } else {
                    newClusters.Add(clustering.Clusters[jj]);
                    NumOfLocalTimeSteps.Add(subSteps);
                    jj++;
                }
            }

            this.numOfClusters = newClusters.Count;
            this.maxNumOfLocalTimeSteps = NumOfLocalTimeSteps.Last();

            return new Clusterer.Clustering(newClusters, clustering.SubGrid);
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
        protected Dictionary<int, double> InterpolateBoundaryValues(Queue<double[]>[] historyDG, int id, double interpolTime) {
            SubGrid sgrd = BoundarySgrds[id - 1];
            Dictionary<int, double> myDic = new Dictionary<int, double>();

            for (int j = 0; j < sgrd.LocalNoOfCells; j++) {
                //int cell = sgrd.SubgridIndex2LocalCellIndex[j]; --> changed to local-only operation
                int cell = jSub2jCell[id - 1][j];
                int BoundaryGridId = BoundaryTopology[id - 1, cell];
                // cell at boundary
                // f== each field
                // n== basis polynomial
                foreach (DGField f in Mapping.Fields) {
                    for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                        int cellIndex = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                        double[] valueHist = new double[order];
                        int k = 0;
                        foreach (double[] histArray in historyDG[BoundaryGridId]) {
                            valueHist[k] = histArray[cellIndex];
                            k++;
                        }
                        double[] timeHistory = GetBoundaryCellTimeHistory(BoundaryGridId, cell);
                        double interpolatedValue = Interpolate(timeHistory, valueHist, interpolTime, order);
                        myDic.Add(cellIndex, interpolatedValue);
                    }
                }
            }
            return myDic;
        }

        /// <summary>
        /// Returns the update times of the boundary cells of a particular cluster.
        /// Needed for the flux interpolation when updating the cells of smaller clusters.
        /// </summary>
        /// <param name="clusterId">Id of the cluster</param>
        /// <param name="cell">Boundary cell</param>
        /// <returns>Array of update times of the boundary cells of a cluster</returns>
        virtual protected double[] GetBoundaryCellTimeHistory(int clusterId, int cell) {
            if (adaptive) {
                Queue<double[]> historyTimePerCell = localABevolve[clusterId].historyTimePerCell;

                // Mapping from local cell index to subgrid index
                int subgridIndex = localABevolve[clusterId].sgrd.LocalCellIndex2SubgridIndex[cell]; // Könnte im Parallelen zu Problemen führen (Stephan)

                // Add times from history
                double[] result = new double[order];
                int i = 0;
                foreach (double[] timePerCell in historyTimePerCell) {
                    result[i] = timePerCell[subgridIndex];
                    i++;
                }

                // Add current time
                result[i] = localABevolve[clusterId].Time;

                return result;
            } else
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

        /// <summary>
        /// Copies the histories from all ABevolve objects from the last time step
        /// to the new ABevolve objects from the current time step.
        /// The information in the ABevolve objects is "cluster-based".
        /// Therefore, all information is first copied in a "cell-based" intermediate state
        /// and then redistributed to the ABevolve objects based on the new clustering.
        /// </summary>
        private void CopyHistoriesOfABevolver() {

            // Previous ABevolve objects: Link "LocalCells --> SubgridIndex"
            int[][] LocalCells2SubgridIndexPrevious = new int[localABevolvePrevious.Length][];
            for (int id = 0; id < localABevolvePrevious.Length; id++) {
                LocalCells2SubgridIndexPrevious[id] = localABevolvePrevious[id].sgrd.LocalCellIndex2SubgridIndex;
            }

            // New ABevolve objects: Link "SubgridIndex --> LocalCells"
            //int[][] LocalCells2SubgridIndex = new int[localABevolve.Length][];
            int[][] Subgrid2LocalCellsIndex = new int[localABevolve.Length][];
            for (int id = 0; id < localABevolve.Length; id++) {
                //LocalCells2SubgridIndex[id] = localABevolve[id].sgrd.LocalCellIndex2SubgridIndex;
                Subgrid2LocalCellsIndex[id] = localABevolve[id].sgrd.SubgridIndex2LocalCellIndex;
            }

            // Helper arrays
            double[] changeRatesPrevious = new double[Mapping.LocalLength];
            double[] DGCoordinatesPrevious = new double[Mapping.LocalLength];
            double[] timesPerCellPrevious = new double[Mapping.LocalLength];

            // Very likely, this is unnecessary (already done in Perform(dt))
            ShortenHistories(localABevolvePrevious);

            // Copy histories from previous to new ABevolve objects (loop over all cells) depending on the LTS order
            for (int ord = 0; ord < order - 1; ord++) {
                for (int cell = 0; cell < gridData.iLogicalCells.NoOfLocalUpdatedCells; cell++) {
                    // Previous subgrid of the cell
                    int oldClusterID = GetSubgridIdOfPrevious(cell, LocalCells2SubgridIndexPrevious);

                    // Previous subgrid index of the cell
                    int subgridIndex = LocalCells2SubgridIndexPrevious[oldClusterID][cell];

                    // Store time-history of the cell
                    timesPerCellPrevious[cell] = localABevolvePrevious[oldClusterID].historyTimePerCell.ElementAt(ord)[subgridIndex];

                    foreach (DGField f in Mapping.Fields) {
                        for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                            // f == field, n == basis polynomial
                            int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                            changeRatesPrevious[index] = localABevolvePrevious[oldClusterID].HistoryChangeRate.ElementAt(ord)[index];
                            DGCoordinatesPrevious[index] = localABevolvePrevious[oldClusterID].HistoryDGCoordinate.ElementAt(ord)[index];
                        }
                    }
                }

                // Fill histories of the new ABevolve objects
                for (int id = 0; id < localABevolve.Length; id++) {
                    localABevolve[id].historyTimePerCell.Enqueue(OrderValuesBySubgridLength(localABevolve[id].sgrd, timesPerCellPrevious, Subgrid2LocalCellsIndex[id]));
                    localABevolve[id].HistoryChangeRate.Enqueue(OrderValuesBySubgrid(localABevolve[id].sgrd, changeRatesPrevious));
                    localABevolve[id].HistoryDGCoordinate.Enqueue(OrderValuesBySubgrid(localABevolve[id].sgrd, DGCoordinatesPrevious));
                }
            }
        }

        /// <summary>
        /// Links values in an array to their specific position in an subgrids
        /// where the array has exactly the length of the subgrid 
        /// </summary>
        /// <param name="subgrid">The particular subgrids</param>
        /// <param name="values">Values to order</param>
        /// <param name="SubgridIndex2LocalCellIndex">Indices: Subgrid cells --> global cells</param>
        /// <returns></returns>
        private double[] OrderValuesBySubgridLength(SubGrid subgrid, double[] values, int[] SubgridIndex2LocalCellIndex) {
            double[] result = new double[subgrid.LocalNoOfCells];

            for (int i = 0; i < subgrid.LocalNoOfCells; i++) {
                result[i] = values[SubgridIndex2LocalCellIndex[i]];
            }

            return result;
        }

        /// <summary>
        /// Deletes unnecessary entries in the histories
        /// </summary>
        private void ShortenHistories(ABevolve[] abEArray) {
            foreach (ABevolve abE in abEArray) {
                if (abE.historyTimePerCell.Count > order - 1)
                    abE.historyTimePerCell.Dequeue();

                if (abE.HistoryChangeRate.Count > order - 1)
                    abE.HistoryChangeRate.Dequeue();

                if (abE.HistoryDGCoordinate.Count > order - 1)
                    abE.HistoryDGCoordinate.Dequeue();
            }
        }

        /// <summary>
        /// Calculates to which sub-grid the cell belongs. Each cell belongs only to one sub-grid!
        /// </summary>
        /// <param name="cell">
        /// Cell-ID</param>
        /// <param name="LocalCells2SubgridIndex">
        /// storage of all <see cref="SubGrid.LocalCellIndex2SubgridIndex"/> arrays for each LTS sub-grid</param>
        /// <returns>LTS sub-grid ID to which the cell belong</returns>
        private int GetSubgridIdOfPrevious(int cell, int[][] LocalCells2SubgridIndex) {
            int id = -1;
            for (int i = 0; i < localABevolvePrevious.Length; i++) {            // Might be adapted in parent class, looping over ABevolve objects instead of NumOfSubgrids
                if (LocalCells2SubgridIndex[i][cell] >= 0)
                    id = i;
            }
            return id;
        }
    }
}