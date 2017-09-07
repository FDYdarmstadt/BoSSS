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
using BoSSS.Solution.Timestepping;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Solution {

    /// <summary>
    /// Implementation of an adaptive local time stepping algorithm (A-LTS)
    /// that features a dynamic (re-)clustering in all cells in each time step
    /// depending on a temporarly changing metric, e.g., a time step constraint
    /// <para>This class is based on <see cref="AdamsBashforthLTS"/>.</para>
    /// </summary>
    public class AdamsBashforthAdaptiveLTS : AdamsBashforthLTS {

        /// <summary>
        /// Hack for testing A-LTS with the scalar transport equation
        /// </summary>
        public ChangeRateCallback UpdateSensorAndAV {
            get;
            private set;
        }

        /// <summary>
        /// Constant number of sub-grids specified by the user
        /// </summary>
        static int numOfSubgridsInit;

        /// <summary>
        /// localABevolve objects from the previous time step needed for copying the histories
        /// </summary>
        protected ABevolve[] localABevolvePrevious;

        public MultidimensionalArray MetricOne {
            get;
            set;
        }

        public MultidimensionalArray MetricTwo {
            get;
            set;
        }

        public double ChangeMetricTime {
            get;
            set;
        }

        public bool IsNUnitTest {
            get;
            set;
        }

        /// <summary>
        /// Standard constructor for the adaptive local time stepping algorithm
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap">Coordinate mapping for the variable fields</param>
        /// <param name="Parameters">Coordinate mapping for the parameter fields</param>
        /// <param name="order">LTS/AB order</param>
        /// <param name="numOfSubgrids">Amount of sub-grids/clusters to be used for LTS</param>
        /// <param name="timeStepConstraints">Time step constraints for later usage as metric</param>
        /// <param name="sgrd">Sub-grids, e.g., from previous time steps</param>
        public AdamsBashforthAdaptiveLTS(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, int numOfSubgrids, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null, ChangeRateCallback test = null) : base(spatialOp, Fieldsmap, Parameters, order, numOfSubgrids, timeStepConstraints, sgrd, conservative: false) {

            numOfSubgridsInit = numOfSubgrids;
            UpdateSensorAndAV = test;

            // Choose Explicit Euler als time stepping scheme for start-up phase
            //RungeKuttaScheme = new RungeKutta(
            //    RungeKutta.RungeKuttaSchemes.ExplicitEuler,
            //    spatialOp,
            //    Fieldsmap,
            //    Parameters,
            //    timeStepConstraints,
            //    sgrd);


            // ########################################### CNS
            //RungeKuttaScheme.OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);

            //for (int i = 0; i < subgridList.Count; i++) {
            //    localABevolve[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);
            //}
            // ########################################### CNS


            // Hack for scalar transport
            if (test != null) {
                for (int i = 0; i < subgridList.Count; i++) {
                    localABevolve[i].OnBeforeComputeChangeRate += UpdateSensorAndAV;
                }
                RungeKuttaScheme.OnBeforeComputeChangeRate += UpdateSensorAndAV;
            }
        }

        /// <summary>
        /// Creates a cell metric based on the active time step constraints
        /// </summary>
        /// <returns>Returns a cell metric as multi-dimensional array</returns>
        protected virtual MultidimensionalArray GetTimeStepsConstraintsInCells() {
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(gridData.iLogicalCells.NoOfLocalUpdatedCells);

            // Adapted from Variables.cs --> DerivedVariable CFL
            for (int i = 0; i < gridData.iLogicalCells.NoOfLocalUpdatedCells; i++) {
                cellMetric[i] = timeStepConstraints.Min(c => c.GetLocalStepSize(i, 1));
            }

            return cellMetric;
        }

        /// <summary>
        /// Creates a cell metric based on 1/artificial viscosity
        /// </summary>
        /// <returns>Returns a cell metric as multi-dimensional array</returns>
        protected virtual MultidimensionalArray GetOneOverAVInCells() {
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(gridData.iGeomCells.NoOfCells);
            DGField avField = ParameterMapping.Fields.Where(c => c.Identification.Equals("viscosity")).Single();

            foreach (Chunk chunk in CellMask.GetFullMask(avField.GridDat)) {
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;
                    cellMetric[cell] = 1 / (avField.GetMeanValue(cell) + 1);
                }
            }

            return cellMetric;
        }

        /// <summary>
        /// Metric used by the Kmeans algorithm for cell clustering
        /// </summary>
        /// <returns> Returns the used cell metric as multi-dimensional array</returns>
        protected override MultidimensionalArray GetCellMetric() {
            if (IsNUnitTest)
                return GetNUnitMetric();
            else
                //return GetTimeStepsConstraintsInCells();
                return GetSmallestDistanceInCells();
            //return GetOneOverAVInCells();
            //return GetTestMetric();
        }

        /// <summary>
        /// Return a fixed metric for testing
        /// </summary>
        /// <returns></returns>
        protected override MultidimensionalArray GetTestMetric() {
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(gridData.iLogicalCells.NoOfLocalUpdatedCells);

            //    if (m_Time < 2e-5) {
            //        cellMetric[0] = 1;
            //        cellMetric[1] = 1;
            //        cellMetric[2] = 1;
            //        cellMetric[3] = 0.1;
            //        cellMetric[4] = 0.1;
            //        cellMetric[5] = 0.1;
            //        cellMetric[6] = 0.1;
            //        cellMetric[7] = 0.1;
            //        cellMetric[8] = 0.1;
            //        cellMetric[9] = 0.1;
            //} else {
            //        cellMetric[0] = 0.1;
            //        cellMetric[1] = 0.1;
            //        cellMetric[2] = 0.1;
            //        cellMetric[3] = 0.1;
            //        cellMetric[4] = 0.1;
            //        cellMetric[5] = 0.1;
            //        cellMetric[6] = 0.1;
            //        cellMetric[7] = 1;
            //        cellMetric[8] = 1;
            //        cellMetric[9] = 1;
            //    }

            if (m_Time < 4e-5) {

                //for (int i = 0; i < gridData.iGeomCells.NoOfCells; i++) {
                //    cellMetric[i] = 1.0 - gridData.CellPartitioning.MpiRank * 0.5;
                //}

                //cellMetric[0] = 1;
                //cellMetric[1] = 0.5;
                //cellMetric[2] = 1;
                //cellMetric[3] = 0.5;


                //cellMetric[0] = 1;
                //for (int i = 1; i < gridData.iGeomCells.NoOfCells; i++) {
                //    cellMetric[i] = 1.0 - gridData.CellPartitioning.MpiRank * 0.5;
                //}

                //cellMetric[0] = 2;
                //cellMetric[1] = 0.5;
                //cellMetric[2] = 0.5;
                //cellMetric[3] = 1.0;

                // 2 sub-grids
                //cellMetric[0] = 1;
                //cellMetric[1] = 1;
                //cellMetric[2] = 1;
                //cellMetric[3] = 0.5;

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
                //cellMetric[1] = 0.1;
                //cellMetric[2] = 1;
                //cellMetric[3] = 0.05;


            } else {
                // 3 sub-grids
                //cellMetric[0] = 0.5;
                //cellMetric[1] = 0.5;
                //cellMetric[2] = 1;
                //cellMetric[3] = 2;

                //cellMetric[0] = 1;
                //cellMetric[1] = 0.5;
                //cellMetric[2] = 2;

                // 2 sub-grids
                //cellMetric[0] = 0.5;
                //cellMetric[1] = 0.5;
                //cellMetric[2] = 0.5;
                //cellMetric[3] = 1;
            }

            return cellMetric;
        }

        /// <summary>
        /// Performs one time step
        /// </summary>
        /// <param name="dt">Time step size</param>
        public override double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {

                if (localABevolve[0].HistoryChangeRate.Count >= order - 1) {

                    #region AV+LTS
                    // Necessary in order to use the number of sub-grids specified by the user for the reclustering in each time step
                    // Otherwise the value could be changed by the constructor of the parent class (AdamsBashforthLTS.cs) --> CreateSubGrids()
                    numOfSubgrids = numOfSubgridsInit;

                    CreateSubGrids();
                    CalculateNumberOfLocalTS(); // Might remove sub-grids when time step sizes are too similar

                    bool reclustered = CheckForNewClustering();
                    //bool reclustered = true;

                    // After intitial phase, activate adaptive mode for all ABevolve objects
                    foreach (ABevolve abE in localABevolve)
                        abE.adaptive = true;

                    //if (order != 1 && reclustered) {
                    if (reclustered) {
                        //CalculateNumberOfLocalTS();   // Change: is called after CreateSubGrids(), might be simplified

                        // Store all localAbevolve objects from the last time step for copying the histories
                        ShortenHistories(localABevolve);
                        localABevolvePrevious = localABevolve;

                        // Create array of Abevolve objects based on the new clustering
                        localABevolve = new ABevolve[this.numOfSubgrids];

                        for (int i = 0; i < subgridList.Count; i++) {
                            localABevolve[i] = new ABevolve(Operator, Mapping, ParameterMapping, order, adaptive: true, sgrd: subgridList[i]);
                            localABevolve[i].ResetTime(m_Time);
                            //localABevolve[i].OnBeforeComputeChangeRate += (t1, t2) => this.RaiseOnBeforComputechangeRate(t1, t2);
                            if (UpdateSensorAndAV != null)
                                localABevolve[i].OnBeforeComputeChangeRate += UpdateSensorAndAV;     // Scalar transport
                        }

                        CopyHistoriesOfABevolver();
                    } else
                        Console.WriteLine("#####Clustering has NOT changed#####");

                    GetBoundaryTopology();
                    #endregion

                    if (timeStepConstraints != null) {
                        dt = CalculateTimeStep();
                    }

                    for (int i = 0; i < this.numOfSubgrids; i++) {
                        Console.WriteLine("LTS: id=" + i + " -> sub-steps=" + NumOfLocalTimeSteps[i] + " and elements=" + subgridList[i].GlobalNoOfCells);
                    }

                    // Saves the results at t_n
                    double[] y0 = new double[Mapping.LocalLength];
                    DGCoordinates.CopyTo(y0, 0);

                    double time0 = m_Time;
                    double time1 = m_Time + dt;

                    // Evolves each sub-grid with its own time step: only one step
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
                    //historyTime_Q = new Queue<double>[numOfSubgrids];
                    //for (int i = 0; i < numOfSubgrids; i++) {
                    //    historyTime_Q[i] = localABevolve[i].HistoryTime;
                    //}

                    // Perform the local time steps
                    for (int localTS = 1; localTS < MaxLocalTS; localTS++) {
                        for (int id = 1; id < numOfSubgrids; id++) {

                            //Evolve condition: Is "ABevolve.Time" at "AB_LTS.Time"?
                            //if (Math.Abs(localABevolve[id].Time - m_Time) < 1e-10) {
                            if ((localABevolve[id].Time - m_Time) < 1e-10) {
                                double localDt = dt / NumOfLocalTimeSteps[id];
                                DGCoordinates.Clear();
                                DGCoordinates.CopyFrom(historyDGC_Q[id].Last(), 0);

                                double[] interpolatedCells = InterpolateBoundaryValues(historyDGC_Q, id, localABevolve[id].Time);
                                DGCoordinates.axpy<double[]>(interpolatedCells, 1);

                                localABevolve[id].Perform(localDt);

                                m_Time = localABevolve.Min(s => s.Time);
                            }

                            // Are we at an (intermediate-) syncronization levels?
                            // For conservatvity, we have to correct the values of the larger cell cluster
                            //double[,] CorrectionMatrix = new double[this.NumOfSgrd, this.NumOfSgrd];

                            //for (int idCoarse = 0; idCoarse < id; idCoarse++) {
                            //    if (Math.Abs(localABevolve[id].Time - localABevolve[idCoarse].Time) < 1e-10 &&
                            //         !(Math.Abs(localABevolve[idCoarse].Time - CorrectionMatrix[idCoarse, id]) < 1e-10)) {
                            //        CorrectFluxes(idCoarse, id, historyDGC_Q);
                            //        CorrectionMatrix[idCoarse, id] = localABevolve[idCoarse].Time;
                            //    }
                            //}
                        }
                    }

                    // Finalize Step
                    DGCoordinates.Clear();
                    for (int id = 0; id < numOfSubgrids; id++) {
                        DGCoordinates.axpy<double[]>(historyDGC_Q[id].Last(), 1);
                    }

                    m_Time = time0 + dt;

                } else {

                    // Start-Up phase
                    double[] currentChangeRate = new double[Mapping.LocalLength];
                    double[] upDGC = new double[Mapping.LocalLength];


                    // Save history: Time
                    for (int i = 0; i < numOfSubgrids; i++) {
                        double[] currentTime = new double[localABevolve[i].sgrd.LocalNoOfCells];
                        for (int j = 0; j < currentTime.Length; j++) {
                            currentTime[j] = m_Time;
                        }
                        localABevolve[i].historyTimePerCell.Enqueue(currentTime);
                    }

                    // Save histories: Change rate, boundary fluxes
                    for (int i = 0; i < subgridList.Count; i++) {
                        double[] localCurrentChangeRate = new double[currentChangeRate.Length];
                        double[] edgeFlux = new double[gridData.iGeomEdges.Count * Mapping.Fields.Count];
                        localABevolve[i].ComputeChangeRate(localCurrentChangeRate, m_Time, 0, edgeFlux);
                        localABevolve[i].HistoryChangeRate.Enqueue(localCurrentChangeRate);
                        localABevolve[i].HistoryBndFluxes.Enqueue(edgeFlux);
                    }

                    dt = RungeKuttaScheme.Perform(dt);

                    DGCoordinates.CopyTo(upDGC, 0);

                    // Save history: DG coordiantes
                    for (int i = 0; i < subgridList.Count; i++) {
                        localABevolve[i].HistoryDGCoordinate.Enqueue(OrderValuesBySubgrid(subgridList[i], upDGC));
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
        /// Returns the update times of the boundary cells of a particular cluster.
        /// Needed for the flux interpolation when updating the cells of smaller clusters.
        /// </summary>
        /// <param name="clusterId">Id of the cluster</param>
        /// <param name="cell">Boundary cell</param>
        /// <returns>Array of update times of the boundary cells of a cluster</returns>
        protected override double[] GetBoundaryCellTimeHistory(int clusterId, int cell) {
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
        /// Checks the current and previous clustering for changes.
        /// </summary>
        /// <returns>True, if clustering has changed. False, if clustering has not changed.</returns>
        private bool CheckForNewClustering() {
            bool localResult = false;   // false = no reclustering needed

            if (subgridList.Count != localABevolve.Length)
                localResult = true;
            else {
                for (int i = 0; i < subgridList.Count; i++) {
                    if (!subgridList[i].VolumeMask.Equals(localABevolve[i].sgrd.VolumeMask)) {
                        localResult = true;
                    }
                }
            }

            bool globalResult;
            unsafe {
                int localResultAsInt = localResult ? 1 : 0;
                int globalResultAsInt;
                csMPI.Raw.Allreduce((IntPtr)(&localResultAsInt), (IntPtr)(&globalResultAsInt), numOfSubgrids, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.LOR, csMPI.Raw._COMM.WORLD);
                globalResult = globalResultAsInt == 1 ? true : false;
            }

            return globalResult;
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

        private MultidimensionalArray GetNUnitMetric() {
            if (m_Time < ChangeMetricTime)
                return MetricOne;
            else
                return MetricTwo;
        }

    }
}
