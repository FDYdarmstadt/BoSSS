﻿/* =======================================================================
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
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Implementation of a local time stepping algorithm using with Adams-Bashforth time integration, according to: 
    /// Winters, A. R., &amp; Kopriva, D. A. (2013). High-Order Local Time Stepping
    /// on Moving DG Spectral Element Meshes. Journal of Scientific Computing.
    /// DOI: 10.1007/s10915-013-9730-z
    /// </summary>
    public class AdamsBashforthLTS : AdamsBashforth {

        /// <summary>
        /// Local evolvers for each cluster
        /// </summary>
        protected ABevolve[] ABevolver;

        /// <summary>
        /// Boundary topology of the clusters
        /// </summary>
        protected int[,] BoundaryTopology;

        /// <summary>
        /// Helper array of sub-grids which only stores the boundary elements for each sub-grid
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
        /// Bool for triggering the adaptive LTS version with reclustering
        /// </summary>
        protected bool adaptive;

        /// <summary>
        /// Constant number of sub-grids specified by the user
        /// </summary>
        public int NumberOfClustersInitial {
            get;
            protected set;
        }

        /// <summary>
        /// Interval for reclustering when LTS is used in adpative mode
        /// 0: standard LTS, e.g., 10: <see cref="Clusterer.CheckForNewClustering(Clusterer.Clustering, Clusterer.Clustering)"/> is called in every tenth time step
        /// </summary>
        protected int reclusteringInterval;

        /// <summary>
        /// A <see cref="Clusterer"/> that divides the grids into sub-grids.
        /// </summary>
        protected Clusterer clusterer;

        /// <summary>
        /// Time history queue for each cluster (only used in LTS runs)
        /// </summary>
        Queue<double>[] historyTime_Q;

        /// <summary>
        /// The current <see cref="Clusterer.Clustering"/>
        /// </summary>
        public Clusterer.Clustering CurrentClustering {
            get;
            protected set;
        }

        private bool reclusteredByGridRedist = false;

        public void SetReclusteredByGridRedist(bool result) {
            this.reclusteredByGridRedist = result;
        }

        private bool forceReclustering;

        private bool Logging;

        public double[] log_clusterDts {
            get;
            private set;
        }

        public int[] log_clusterElements {
            get;
            private set;
        }

        public int[] log_clusterSubSteps {
            get;
            private set;
        }

        public bool ConsoleOutput {
            get;
            private set;
        }

        ////################# Hack for saving to database in every (A)LTS sub-step
        //private Action<TimestepNumber, double> saveToDBCallback;
        ////################# Hack for saving to database in every (A)LTS sub-step

        /// <summary>
        /// Standard constructor for the (adaptive) local time stepping algorithm
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap">Coordinate mapping for the variable fields</param>
        /// <param name="Parameters">optional parameter fields, can be null if <paramref name="spatialOp"/> contains no parameters; must match the parameter field list of <paramref name="spatialOp"/>, see <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/></param>
        /// <param name="order">LTS/AB order</param>
        /// <param name="numOfClusters">Amount of sub-grids/clusters to be used for LTS</param>
        /// <param name="timeStepConstraints">Time step constraints for later usage as metric</param>
        /// <param name="subGrid">Sub-grids, e.g., from previous time steps</param>
        /// <param name="fluxCorrection">Bool for triggering the fluss correction</param>
        /// <param name="reclusteringInterval">Interval for potential reclustering</param>
        /// <param name="saveToDBCallback">Hack for plotting all sub-steps</param>
        /// <param name="consoleOutput"></param>
        /// <param name="forceReclustering"></param>
        /// <param name="logging"></param>
        /// <param name="maxNumOfSubSteps"></param>
        /// <remarks>Uses the k-Mean clustering, see <see cref="BoSSS.Solution.Utils.Kmeans"/>, to generate the element groups</remarks>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int numOfClusters, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid subGrid = null, bool fluxCorrection = true, int reclusteringInterval = 0, Action<TimestepNumber, double> saveToDBCallback = null, int maxNumOfSubSteps = 0, bool forceReclustering = false, bool logging = false, bool consoleOutput = false)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, subGrid) {

            this.forceReclustering = forceReclustering;
            this.Logging = logging;
            this.ConsoleOutput = consoleOutput;

            if (reclusteringInterval != 0) {
                NumberOfClustersInitial = numOfClusters;
                this.adaptive = true;
            }

            // Add OnBeforeComputeChangeRate (AV) to start-up phase time stepper
            RungeKuttaScheme.OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);

            this.reclusteringInterval = reclusteringInterval;
            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.fluxCorrection = fluxCorrection;

            if (ConsoleOutput) {
                Console.WriteLine("### This is ABLTS ctor ###");
            }

            clusterer = new Clusterer(this.gridData, maxNumOfSubSteps);
            CurrentClustering = clusterer.CreateClustering(numOfClusters, this.TimeStepConstraints, this.SubGrid);    // Might remove clusters when their centres are too close
            CurrentClustering = clusterer.TuneClustering(CurrentClustering, Time, this.TimeStepConstraints); // Might remove clusters when their time step sizes are too similar

            ABevolver = new ABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new ABevolve(spatialOp, Fieldsmap, Parameters, order, adaptive: this.adaptive, sgrd: CurrentClustering.Clusters[i]);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }

            GetBoundaryTopology();

            // Saving time steps in subgrids
            //this.saveToDBCallback = saveToDBCallback;
        }

        /// <summary>
        /// Constructor for (A)LTS with IBM
        /// </summary>
        public AdamsBashforthLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int numOfClusters, bool IBM, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid subGrid = null, bool fluxCorrection = true, int reclusteringInterval = 0, bool forceReclustering = false, bool logging = false, bool consoleOutput = false)
            : base(spatialOp, Fieldsmap, Parameters, order, timeStepConstraints, subGrid) {

            this.forceReclustering = forceReclustering;
            this.Logging = logging;
            this.ConsoleOutput = consoleOutput;
            this.gridData = Fieldsmap.Fields.First().GridDat;
            this.fluxCorrection = fluxCorrection;

            if (reclusteringInterval != 0) {
                NumberOfClustersInitial = numOfClusters;
                this.adaptive = true;
            }
            this.reclusteringInterval = reclusteringInterval;

            // Add OnBeforeComputeChangeRate (AV) to start-up phase time stepper
            RungeKuttaScheme.OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
        }

        /// <summary>
        /// Performs a time step
        /// </summary>
        /// <param name="dt">Time step size that equals -1, if no fixed time step is prescribed</param>
        public override double Perform(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {

                if (ABevolver[0].HistoryChangeRate.Count >= order - 1) {
                    // +++++++++++++++++++++++++++++++++++++++++++
                    // Standard case -- sufficiently large history
                    // +++++++++++++++++++++++++++++++++++++++++++
                    using (var bt = new BlockTrace("AB_LTS_standard", tr)) {

                        if (!reclusteredByGridRedist) {
                            TryNewClustering(dt);
                        }

                        List<int> numberOfLocalTimeSteps = new List<int>();
                        double[] clusterDts = new double[CurrentClustering.NumberOfClusters];

                        // Set the number of sub steps (is calculated in every time step, regardless of whether a reclustering has been performed or not)
                        if (TimeStepConstraints != null) {
                            //dt = CalculateTimeStep();
                            // If no dtFixed is set
                            if (TimeStepConstraints.First().dtMin != TimeStepConstraints.First().dtMax) {
                                (clusterDts, numberOfLocalTimeSteps) = clusterer.GetPerCluster_dtHarmonicSum_SubSteps(CurrentClustering, Time, TimeStepConstraints, eps: 1.0e-1);
                                dt = clusterDts[0];
                            } else {    // dtFixed is set
                                        // Not nice, but working
                                dt = CalculateTimeStep();
                                numberOfLocalTimeSteps = CurrentClustering.SubStepsInitial;
                                for (int i = 0; i < numberOfLocalTimeSteps.Count; i++) {
                                    clusterDts[i] = dt / numberOfLocalTimeSteps[i];
                                }
                            }
                        }

                        // Log time info
                        if (Logging) {
                            this.log_clusterDts = clusterDts;
                            this.log_clusterSubSteps = numberOfLocalTimeSteps.ToArray();
                            this.log_clusterElements = CurrentClustering.Clusters.Select(s => s.GlobalNoOfCells).ToArray();
                        }

                        if (ConsoleOutput) {
                            for (int i = 0; i < numberOfLocalTimeSteps.Count; i++) {
                                Console.WriteLine("Perform(dt):\t\t id={0} -> sub-steps={1}\telements={2}\tdt={3:0.#######E-00}", i, numberOfLocalTimeSteps[i], CurrentClustering.Clusters[i].GlobalNoOfCells, clusterDts[i]);
                                //Console.WriteLine("Perform(dt):\t\t id={0} -> sub-steps={1}\telements={2}\tdt={3}", i, numberOfLocalTimeSteps[i], CurrentClustering.Clusters[i].GlobalNoOfCells, clusterDts[i]);
                            }

                            if (numberOfLocalTimeSteps.Last() > (clusterer.MaxSubSteps + 1) && clusterer.Restrict) {
                                throw new Exception(String.Format("Number of local time steps is larger than {0}! Restriction failed!", clusterer.MaxSubSteps));
                            }
                        }

                        double[,] CorrectionMatrix = new double[CurrentClustering.NumberOfClusters, CurrentClustering.NumberOfClusters];

                        // Test code
                        //double[] bla = new double[numberOfLocalTimeSteps.Count];
                        //for (int i = 0; i < bla.Length; i++) {
                        //    bla[i] = dt / numberOfLocalTimeSteps[i];
                        //    if (bla[i] != clusterDts[i]) {
                        //        throw new Exception("clusterDts wrong");
                        //    }
                        //}

                        // Saves the results at t_n
                        double[] y0 = new double[Mapping.LocalLength];
                        CurrentState.CopyTo(y0, 0);

                        double time0 = m_Time;
                        double time1 = m_Time + clusterDts[0];

                        // Evolves each sub-grid with its own time step (only one step)
                        // (The result is not written to m_DGCoordinates!)
                        for (int i = 0; i < ABevolver.Length; i++) {
                            //localABevolve[i].completeBndFluxes.Clear();
                            //if (localABevolve[i].completeBndFluxes.Any(x => x != 0.0)) Console.WriteLine("Not all Bnd fluxes were used in correction step!!!");
                            ABevolver[i].Perform(clusterDts[i]);
                        }

                        // After evolving each cell update the time with dt_min
                        m_Time = m_Time + clusterDts.Last();

                        TimestepNumber subTimestep = new TimestepNumber(TimeInfo.TimeStepNumber - 1);

                        //if (saveToDBCallback != null) {
                        //    subTimestep = subTimestep.NextIteration();
                        //    saveToDBCallback(subTimestep, m_Time);
                        //}

                        // Saves the history of DG_Coordinates for each cluster
                        Queue<double[]>[] historyDGC_Q = new Queue<double[]>[CurrentClustering.NumberOfClusters];
                        for (int i = 0; i < historyDGC_Q.Length; i++) {
                            historyDGC_Q[i] = ABevolver[i].HistoryDGCoordinate;
                        }

                        if (!adaptive) {
                            // Saves DtHistory for each cluster
                            historyTime_Q = new Queue<double>[CurrentClustering.NumberOfClusters];
                            for (int i = 0; i < historyTime_Q.Length; i++) {
                                historyTime_Q[i] = ABevolver[i].HistoryTime;
                            }
                        }

                        // Perform the local time steps
                        for (int localTS = 1; localTS < numberOfLocalTimeSteps.Last(); localTS++) {
                            for (int id = 1; id < CurrentClustering.NumberOfClusters; id++) {
                                //Evolve Condition: Is "ABevolve.Time" at "AB_LTS.Time"?
                                if ((ABevolver[id].Time - m_Time) < 1e-10) {
                                    foreach (Chunk chunk in CurrentClustering.Clusters[id].VolumeMask) {
                                        foreach (int cell in chunk.Elements) {
                                            // f == each field
                                            // n == basis polynomial
                                            foreach (DGField f in Mapping.Fields) {
                                                for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                                                    int coordinateIndex = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                                                    CurrentState[coordinateIndex] = historyDGC_Q[id].Last()[coordinateIndex];
                                                }
                                            }
                                        }
                                    }

                                    Dictionary<int, double> myDic = InterpolateBoundaryValues(historyDGC_Q, id, ABevolver[id].Time);

                                    foreach (KeyValuePair<int, double> kvp in myDic) {
                                        CurrentState[kvp.Key] = kvp.Value;
                                    }

                                    ABevolver[id].Perform(clusterDts[id]);

                                    m_Time = ABevolver.Min(s => s.Time);

                                    //if (saveToDBCallback != null) {
                                    //    subTimestep = subTimestep.NextIteration();
                                    //    saveToDBCallback(subTimestep, m_Time);
                                    //}
                                }

                                // Are we at an (intermediate-) syncronization levels?
                                // For conservatvity, we have to correct the values of the larger cell cluster
                                if (fluxCorrection) {
                                    for (int idCoarse = 0; idCoarse < id; idCoarse++) {
                                        if (Math.Abs(ABevolver[id].Time - ABevolver[idCoarse].Time) < 1e-10 &&
                                             !(Math.Abs(ABevolver[idCoarse].Time - CorrectionMatrix[idCoarse, id]) < 1e-10)) {
                                            if (fluxCorrection) {
                                                CorrectFluxes(idCoarse, id, historyDGC_Q);
                                            }
                                            CorrectionMatrix[idCoarse, id] = ABevolver[idCoarse].Time;
                                        }
                                    }
                                }
                            }
                        }

                        // Finalize step
                        // Use unmodified values in history of DGCoordinates (DGCoordinates could have been modified by
                        // InterpolateBoundaryValues, should be resetted afterwards) 
                        CurrentState.Clear();
                        for (int id = 0; id < historyDGC_Q.Length; id++) {
                            CurrentState.axpy<double[]>(historyDGC_Q[id].Last(), 1);
                        }

                        // Update time
                        m_Time = time0 + clusterDts[0];
                    }
                } else {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Startup - use Runge Rutta until history is sufficient
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                    using (var rkPhase = new BlockTrace("AB_LTS_Rkstartup", tr)) {

                        double[] currentChangeRate = new double[Mapping.LocalLength];
                        double[] upDGC = new double[Mapping.LocalLength];

                        // Save time history
                        if (adaptive) {
                            for (int i = 0; i < ABevolver.Length; i++) {
                                double[] currentTime = new double[ABevolver[i].ABSubGrid.LocalNoOfCells];
                                for (int j = 0; j < currentTime.Length; j++) {
                                    currentTime[j] = m_Time;
                                }
                                ABevolver[i].historyTimePerCell.Enqueue(currentTime);
                            }
                        } else {
                            if (ABevolver[0].HistoryTime.Count == 0) {
                                for (int i = 0; i < ABevolver.Length; i++) {
                                    ABevolver[i].HistoryTime.Enqueue(m_Time);
                                }
                            }
                        }

                        // Needed for the history
                        for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                            double[] localCurrentChangeRate = new double[currentChangeRate.Length];
                            double[] edgeFlux = new double[gridData.iGeomEdges.Count * Mapping.Fields.Count];
                            ABevolver[i].ComputeChangeRate(localCurrentChangeRate, m_Time, 0, edgeFlux);
                            ABevolver[i].HistoryChangeRate.Enqueue(localCurrentChangeRate);
                            ABevolver[i].HistoryBoundaryFluxes.Enqueue(edgeFlux);
                        }

                        dt = RungeKuttaScheme.Perform(dt);

                        CurrentState.CopyTo(upDGC, 0);

                        // Saves ChangeRateHistory for AB LTS
                        // Only entries for the specific cluster
                        for (int i = 0; i < CurrentClustering.NumberOfClusters; i++) {
                            ABevolver[i].HistoryDGCoordinate.Enqueue(OrderValuesByCluster(CurrentClustering.Clusters[i], upDGC));
                            if (!adaptive)
                                ABevolver[i].HistoryTime.Enqueue(RungeKuttaScheme.Time);
                        }

                        // RK is a global timeStep
                        // -> time update for all other timeStepper with rk.Time
                        m_Time = RungeKuttaScheme.Time;
                        foreach (ABevolve ab in ABevolver) {
                            ab.ResetTime(m_Time, TimeInfo.TimeStepNumber);
                        }
                    }
                }
            }

            return dt;
        }

        /// <summary>
        /// To achieve a conservative time stepping scheme, we have to correct the DG coordinates of 
        /// large interface cells
        /// </summary>
        /// <param name="coarseID">Cluster ID of the large cells</param>
        /// <param name="fineID">Cluster ID of the small cells</param>
        /// <param name="historyDGC_Q">History of the DG coordinates</param>
        protected void CorrectFluxes(int coarseID, int fineID, Queue<double[]>[] historyDGC_Q) {
            // Gather edgeFlux data
            double[] edgeBndFluxCoarse = ABevolver[coarseID].CompleteBoundaryFluxes;
            double[] edgeBndFluxFine = ABevolver[fineID].CompleteBoundaryFluxes;

            int[] LocalCellIdx2SubgridIdx = CurrentClustering.Clusters[coarseID].LocalCellIndex2SubgridIndex;

            CellMask CellMaskCoarse = CurrentClustering.Clusters[coarseID].VolumeMask;

            //Only the edges of coarseID and fineID are needed
            EdgeMask unionEdgeMask = CurrentClustering.Clusters[coarseID].BoundaryEdgesMask.Intersect(CurrentClustering.Clusters[fineID].BoundaryEdgesMask);

            int cellCoarse;
            int cellFine;

            MultidimensionalArray basisScale = gridData.ChefBasis.Scaling;
            int noOfFields = Mapping.Fields.Count;

            //loop over all BoundaryEdges of the coarse sgrd
            {
                foreach (int edge in unionEdgeMask.ItemEnum) {

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
                            double fluxScaling = 1.0 / ABevolver[coarseID].ABCoefficients[0];
                            ABevolver[coarseID].CurrentChangeRate[indexCoarse] += fluxScaling * edgeDG_Correction;
                            ABevolver[coarseID].HistoryBoundaryFluxes.Last()[indexEdge] -= fluxScaling * edgeDG_Correction / basisScaling;
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
        public override void ResetTime(double NewTime, int timestepNumber) {
            base.ResetTime(NewTime, timestepNumber);
            RungeKuttaScheme.ResetTime(NewTime, timestepNumber);

            foreach (var ABevolve in ABevolver) {
                ABevolve.ResetTime(NewTime, timestepNumber);
            }

            RungeKuttaScheme.ResetTime(NewTime, timestepNumber);
        }

        /// <summary>
        /// Calculated topology of the grid, i.e creating the boundarySubgrids
        /// for each sub-grid
        /// </summary>
        protected void GetBoundaryTopology() {
            // NumOfSgrd - 1, because largest grid (id=0) don't need a boundary cells
            BoundaryTopology = new int[CurrentClustering.NumberOfClusters - 1, gridData.iLogicalCells.NoOfLocalUpdatedCells];
            ArrayTools.SetAll(BoundaryTopology, -1);
            BoundarySgrds = new SubGrid[CurrentClustering.NumberOfClusters - 1];
            jSub2jCell = new int[CurrentClustering.NumberOfClusters - 1][];
            int[][] LocalCells2SubgridIndex = new int[CurrentClustering.NumberOfClusters][];
            BitArray[] SgrdWithGhostCells = new BitArray[CurrentClustering.NumberOfClusters];
            // prepare the calculation and  save temporarily all array which involve MPI communication
            for (int id = 0; id < CurrentClustering.NumberOfClusters; id++) {
                LocalCells2SubgridIndex[id] = CurrentClustering.Clusters[id].LocalCellIndex2SubgridIndex;
                SgrdWithGhostCells[id] = CurrentClustering.Clusters[id].VolumeMask.GetBitMaskWithExternal();
            }

            for (int id = 1; id < CurrentClustering.NumberOfClusters; id++) {
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
                                int gridId = GetClusterIDOf(cell1, ABevolver);
                                if (gridId != -1) { // cell is not in void area of IBM
                                    BoundaryTopology[id - 1, cell1] = gridId;
                                    BoBA[cell1] = true;
                                }
                            }
                        } else if (cell1 >= 0 && cell2 >= 0 && LocalCells2SubgridIndex[id][cell1] >= 0 && LocalCells2SubgridIndex[id][cell2] < 0) {
                            //BoT[id - 1, cell2] = getSgrdIdOf(cell2, LocalCells2SubgridIndex);
                            //BoBA[cell2] = true;
                            int gridId = GetClusterIDOf(cell2, ABevolver);
                            if (gridId != -1) { // cell is not in void area of IBM
                                BoundaryTopology[id - 1, cell2] = gridId;
                                BoBA[cell2] = true;
                            }

                        } else if (cell1 >= 0 && cell2 >= 0 && LocalCells2SubgridIndex[id][cell2] >= 0 && LocalCells2SubgridIndex[id][cell1] < 0) {
                            //BoT[id - 1, cell1] = getSgrdIdOf(cell1, LocalCells2SubgridIndex);
                            //BoBA[cell1] = true;
                            int gridId = GetClusterIDOf(cell1, ABevolver);
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
        /// Returns the cluster ID a cell belongs to. Each cell belongs only to one sub-grid!
        /// </summary>
        /// <param name="cell">The cell of interest</param>
        /// <param name="ABevolver">The clusters that have to be checked</param>
        /// <returns>Cluster ID of the cell </returns>
        protected int GetClusterIDOf(int cell, ABevolve[] ABevolver) {
            int id = -1;
            for (int i = 0; i < ABevolver.Length; i++) {
                if (ABevolver[i].ABSubGrid.LocalCellIndex2SubgridIndex[cell] >= 0)
                    id = i;
            }
            return id;
        }

        /// <summary>
        /// Caluclates the local time step sizes
        /// </summary>
        /// <returns>the largest stable timestep</returns>
        protected override double CalculateTimeStep() {
            if (TimeStepConstraints.First().dtMin != TimeStepConstraints.First().dtMax) {
                var result = clusterer.GetPerCluster_dtHarmonicSum_SubSteps(CurrentClustering, Time, TimeStepConstraints, eps: 1.0e-1);
                double[] localDts = result.Item1;
                return localDts[0];
            } else {
                double dt = TimeStepConstraints.First().dtMin;
                dt = Math.Min(dt, TimeStepConstraints.First().Endtime - Time);
                return dt;
            }
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
            SubGrid subGrid = BoundarySgrds[id - 1];
            Dictionary<int, double> myDic = new Dictionary<int, double>();

            for (int j = 0; j < subGrid.LocalNoOfCells; j++) {
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
                Queue<double[]> historyTimePerCell = ABevolver[clusterId].historyTimePerCell;

                // Mapping from local cell index to sub-grid index
                int subGridIndex = ABevolver[clusterId].ABSubGrid.LocalCellIndex2SubgridIndex[cell]; // Könnte im Parallelen zu Problemen führen (Stephan)

                // Add times from history
                double[] result = new double[order];
                int i = 0;
                foreach (double[] timePerCell in historyTimePerCell) {
                    result[i] = timePerCell[subGridIndex];
                    i++;
                }

                // Add current time
                result[i] = ABevolver[clusterId].Time;

                return result;
            } else {
                return historyTime_Q[clusterId].ToArray();
            }
        }

        /// <summary>
        /// Interpolates y-values for the x-values with the given (x,y)-pairs.
        /// Only second and third order is needed.
        /// </summary>
        /// <param name="x">x-values </param>
        /// <param name="y">y-values</param>
        /// <param name="X">(x,y)-pair</param>
        /// <param name="order">interpolation order</param>
        /// <returns>Interpolated y-value at (x,y)</returns>
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
        /// <param name="cluster">The cluster</param>
        /// <param name="results">Result for the complete grid</param>
        /// <returns>Array with entries only for the sgrd-cells</returns>
        protected double[] OrderValuesByCluster(SubGrid cluster, double[] results) {
            double[] ordered = new double[Mapping.LocalLength];

            for (int j = 0; j < cluster.LocalNoOfCells; j++) {
                int cell = cluster.SubgridIndex2LocalCellIndex[j];
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
        protected void CopyHistoriesOfABevolver(ABevolve[] oldABevolver) {

            // Previous ABevolve objects: Link "LocalCells --> SubgridIndex"
            int[][] oldLocalCells2SubgridIndex = new int[oldABevolver.Length][];
            for (int id = 0; id < oldABevolver.Length; id++) {
                oldLocalCells2SubgridIndex[id] = oldABevolver[id].ABSubGrid.LocalCellIndex2SubgridIndex;
            }

            // New ABevolve objects: Link "SubgridIndex --> LocalCells"
            int[][] Subgrid2LocalCellsIndex = new int[ABevolver.Length][];
            for (int id = 0; id < ABevolver.Length; id++) {
                Subgrid2LocalCellsIndex[id] = ABevolver[id].ABSubGrid.SubgridIndex2LocalCellIndex;
            }

            // More helper arrays
            double[] oldChangeRates = new double[Mapping.LocalLength];
            double[] oldDGCoordinates = new double[Mapping.LocalLength];
            double[] oldHistoryTimePerCell = new double[Mapping.LocalLength];

            // Copy histories from previous to new ABevolve objects (loop over all cells) depending on the LTS order
            for (int ord = 0; ord < order - 1; ord++) {
                foreach (Chunk chunk in this.SubGrid.VolumeMask) {
                    foreach (int cell in chunk.Elements) {
                        // Previous cluster of the cell
                        int oldClusterID = GetClusterIDOf(cell, oldABevolver);

                        // Previous cluster index of the cell
                        int oldClusterIndex = oldLocalCells2SubgridIndex[oldClusterID][cell];

                        // Previous time history of the cell
                        oldHistoryTimePerCell[cell] = oldABevolver[oldClusterID].historyTimePerCell.ElementAt(ord)[oldClusterIndex];

                        // Previous change rates and DG coordinates of the cell
                        foreach (DGField f in Mapping.Fields) {
                            for (int n = 0; n < f.Basis.GetLength(cell); n++) {
                                // f == field, n == basis polynomial
                                int index = Mapping.LocalUniqueCoordinateIndex(f, cell, n);
                                oldChangeRates[index] = oldABevolver[oldClusterID].HistoryChangeRate.ElementAt(ord)[index];
                                oldDGCoordinates[index] = oldABevolver[oldClusterID].HistoryDGCoordinate.ElementAt(ord)[index];
                            }
                        }
                    }
                }

                // Fill histories of the new ABevolve objects
                for (int id = 0; id < ABevolver.Length; id++) {
                    ABevolver[id].historyTimePerCell.Enqueue(OrderValuesByClusterLength(ABevolver[id].ABSubGrid, oldHistoryTimePerCell));
                    ABevolver[id].HistoryChangeRate.Enqueue(OrderValuesByCluster(ABevolver[id].ABSubGrid, oldChangeRates));
                    ABevolver[id].HistoryDGCoordinate.Enqueue(OrderValuesByCluster(ABevolver[id].ABSubGrid, oldDGCoordinates));
                }
            }
        }

        /// <summary>
        /// Links values in an array to their specific position in a subgrid
        /// where the array has exactly the length of the subgrid 
        /// </summary>
        /// <param name="subGrid"><see cref="SubGrid"/> of interest</param>
        /// <param name="values">Values to order</param>
        private double[] OrderValuesByClusterLength(SubGrid subGrid, double[] values) {
            double[] result = new double[subGrid.LocalNoOfCells];

            for (int i = 0; i < subGrid.LocalNoOfCells; i++) {
                result[i] = values[subGrid.SubgridIndex2LocalCellIndex[i]];
            }

            return result;
        }

        /// <summary>
        /// Deletes unnecessary entries in the histories
        /// </summary>
        protected void ShortenHistories(ABevolve[] ABevolver) {
            foreach (ABevolve abE in ABevolver) {
                if (abE.historyTimePerCell.Count > order - 1)
                    abE.historyTimePerCell.Dequeue();

                if (abE.HistoryChangeRate.Count > order - 1)
                    abE.HistoryChangeRate.Dequeue();

                if (abE.HistoryDGCoordinate.Count > order - 1)
                    abE.HistoryDGCoordinate.Dequeue();
            }
        }

        /// <summary>
        /// Creates new <see cref="ABevolve"/> objects when a reclustering is performed
        /// </summary>
        protected virtual void CreateNewABevolver() {
            ABevolver = new ABevolve[CurrentClustering.NumberOfClusters];

            for (int i = 0; i < ABevolver.Length; i++) {
                ABevolver[i] = new ABevolve(Operator, Mapping, ParameterMapping, order, adaptive: true, sgrd: CurrentClustering.Clusters[i]);
                ABevolver[i].ResetTime(m_Time, TimeInfo.TimeStepNumber);
                ABevolver[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforeComputechangeRate(t1, t2);
            }
        }

        public bool TryNewClustering(double dt, bool calledByMPIRedist = false) {
            //using (var tr = new ilPSP.Tracing.FuncTrace()) {

                bool reclustered = false;

                if (ABevolver[0].HistoryChangeRate.Count >= order - 1) {
                    if (adaptive) {
                        if (TimeInfo.TimeStepNumber % reclusteringInterval == 0) {
                            // Fix for update problem of artificial viscosity
                            RaiseOnBeforeComputechangeRate(Time, dt);

                            if (ConsoleOutput) {
                                Console.WriteLine("\n-------------------------------------------------------------------------------------------");
                                Console.WriteLine("### BUILD NEW CLUSTERING FOR TESTING ###");
                            }

                            // Necessary in order to use the number of sub-grids specified by the user for the reclustering in each time step
                            // Otherwise the value could be changed by the constructor of the parent class (AdamsBashforthLTS.cs) --> CreateSubGrids()
                            Clusterer.Clustering newClustering = clusterer.CreateClustering(NumberOfClustersInitial, this.TimeStepConstraints, this.SubGrid);
                            newClustering = clusterer.TuneClustering(newClustering, Time, this.TimeStepConstraints); // Might remove sub-grids when their time step sizes are too similar

                            if (calledByMPIRedist || forceReclustering) {
                                reclustered = true;
                            } else {
                                reclustered = clusterer.CheckForNewClustering(CurrentClustering, newClustering);
                            }

                            if (ConsoleOutput) {
                                Console.WriteLine("-------------------------------------------------------------------------------------------");
                            }

                            // After the intitial phase, activate adaptive mode for all ABevolve objects
                            foreach (ABevolve abE in ABevolver) {
                                abE.adaptive = true;
                            }

                            if (reclustered) {

                                if (ConsoleOutput) {
                                    Console.WriteLine("### RECLUSTERING ###");
                                }

                                CurrentClustering = newClustering;
                                ShortenHistories(ABevolver);
                                ABevolve[] oldABevolver = ABevolver;
                                CreateNewABevolver();
                                CopyHistoriesOfABevolver(oldABevolver);
                            }

                            GetBoundaryTopology();
                        }
                    }
                }

                return reclustered;
            }
        }
    }
//}