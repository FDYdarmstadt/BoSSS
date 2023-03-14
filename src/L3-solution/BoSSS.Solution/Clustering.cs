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


using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using static BoSSS.Foundation.XDG.CellAgglomerator;
using BoSSS.Foundation.Comm;

namespace BoSSS.Solution.Utils {

    /// <summary>
    /// Class for a cell clustering that devides the grid into sub-grids
    /// </summary>
    public class Clusterer {

        public class Clustering {

            public List<SubGrid> Clusters {
                get;
                private set;
            }

            public int NumberOfClusters {
                get {
                    return Clusters.Count();
                }
            }

            public SubGrid SubGrid {
                get;
                private set;
            }

            public List<int> SubStepsInitial {
                get;
            }

            public Clustering(List<SubGrid> clusters, SubGrid subGrid, List<int> subStepsInitial = null) {
                this.Clusters = clusters;
                this.SubGrid = subGrid;
                this.SubStepsInitial = subStepsInitial;
            }
        }

        /// <summary>
        /// Information about the grid
        /// </summary>
        private IGridData gridData;

        public int MaxSubSteps {
            get;
            private set;
        }

        public bool Restrict {
            get;
            private set;
        }

        private bool consoleOutput;

        private CellAgglomerator cellAgglomerator;

        /// <summary>
        /// Constructor for the grid clustering
        /// </summary>
        public Clusterer(IGridData gridData, int maxNumOfSubSteps, bool consoleOutput = false, CellAgglomerator cellAgglomerator = null) {
            this.gridData = gridData;
            this.MaxSubSteps = maxNumOfSubSteps;
            if (this.MaxSubSteps != 0) {
                this.Restrict = true;
            }
            this.consoleOutput = consoleOutput;
            this.cellAgglomerator = cellAgglomerator;
        }

        /// <summary>
        /// Creates the sub-grids of the clustering
        /// </summary>     
        /// <returns>A list of sub-grids</returns>
        public Clustering CreateClustering(int numOfClusters, IList<TimeStepConstraint> timeStepConstraints, SubGrid subGrid = null) {
            //using (var tr = new ilPSP.Tracing.FuncTrace()) {
                if (subGrid == null) {
                    subGrid = new SubGrid(CellMask.GetFullMask(gridData));
                }

                // Attention: numOfCells can equal all local cells or only the local cells of a subgrid,
                // e.g., the fluid cells in an IBM simulation
                int numOfCells = subGrid.LocalNoOfCells;

                MultidimensionalArray cellMetric = GetStableTimestepSize(subGrid, timeStepConstraints);
                MultidimensionalArray means = CreateInitialMeans(cellMetric, numOfClusters);
                Kmeans Kmean = new Kmeans(cellMetric.To1DArray(), numOfClusters, means.To1DArray());

                // The corresponding sub-grid IDs
                int[] subGridCellToClusterMap = Kmean.Cluster();
                int[] noOfCellsPerCluster = Kmean.ClusterCount;

                unsafe {
                    int[] globalCC = new int[numOfClusters];
                    // send = means[]
                    // receive = globalMeans[]
                    fixed (int* pSend = &noOfCellsPerCluster[0], pRcv = &globalCC[0]) {
                        csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), numOfClusters, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
                    }
                    noOfCellsPerCluster = globalCC;
                }

                int counter = numOfClusters;
                for (int i = 0; i < numOfClusters; i++) {
                    if (noOfCellsPerCluster[i] == 0) {
                        if (consoleOutput) {
                            Console.WriteLine("Sub-grid/Cluster " + (i + 1) + ", with mean value " + Kmean.Means[i] + ", is empty and not used anymore!");
                        }
                        counter--;
                    }
                }

                // Generating BitArray for all Subgrids, even for those which are empty, i.e ClusterCount == 0
                BitArray[] baMatrix = new BitArray[numOfClusters];
                for (int i = 0; i < numOfClusters; i++) {
                    //baMatrix[i] = new BitArray(gridData.iLogicalCells.NoOfCells);
                    baMatrix[i] = new BitArray(gridData.iLogicalCells.NoOfLocalUpdatedCells);
                }

                // Filling the BitArrays
                for (int i = 0; i < numOfCells; i++) {
                    baMatrix[subGridCellToClusterMap[i]][subGrid.SubgridIndex2LocalCellIndex[i]] = true;
                }

                // IBM source cells are assigned to the cluster of the corresponding target cells
                // This code is only excuted in IBM simulation runs
                if (this.cellAgglomerator != null) {
                    // MPI exchange in order to get cellToClusterMap (local + external cells)
                    int JE = gridData.iLogicalCells.Count;
                    int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
                    int[] cellToClusterMap = new int[JE];

                    int JSub = subGrid.LocalNoOfCells;
                    int[] jSub2j = subGrid.SubgridIndex2LocalCellIndex;
                    for (int jsub = 0; jsub < JSub; jsub++) {
                        cellToClusterMap[jSub2j[jsub]] = subGridCellToClusterMap[jsub];
                    }
                    cellToClusterMap.MPIExchange(gridData);

                    foreach (AgglomerationPair aggPair in this.cellAgglomerator.AggInfo.AgglomerationPairs) {
                        // AgglomerationPairs can contain combinations where jCellSource is on one MPI rank
                        // and the corresponding target cell is on another MPI rank. These duplications have to be eleminated.
                        if (aggPair.jCellSource < J) {
                            // Assign source cell to the cluster of the corresponding target cell
                            int clusterOfTargetCell = cellToClusterMap[aggPair.jCellTarget];
                            baMatrix[clusterOfTargetCell][aggPair.jCellSource] = true;

                            // Delete source cell from other clusters
                            for (int j = 0; j < numOfClusters; j++) {
                                if (clusterOfTargetCell != j) {
                                    baMatrix[j][aggPair.jCellSource] = false;
                                }
                            }
                        }
                    }
                }

                // Generating the sub-grids
                List<SubGrid> clusters = new List<SubGrid>(counter);
                for (int i = 0; i < numOfClusters; i++) {
                    // Generating only the sub-grids which are not empty
                    if (noOfCellsPerCluster[i] != 0) {
                        BitArray ba = baMatrix[i];
                        clusters.Add(new SubGrid(new CellMask(gridData, ba)));
                    }
                }

                if (consoleOutput) {
                    for (int i = 0; i < clusters.Count; i++) {
                        Console.WriteLine("CreateClustering:\t id=" + i + " ->\t\telements=" + clusters[i].GlobalNoOfCells);
                    }
                }

                return new Clustering(clusters, subGrid);
            }
        //}

        /// <summary>
        /// Creates an array with an tanh spaced distribution of the mean
        /// values between maximum and minimum value of a given cell metric, 
        /// e.g., minimal distance between two nodes in a cell <see cref="GridData.CellData.h_min"/>
        /// </summary>
        /// <returns>Double[] with the length of the number of given sub-grids></returns>
        private MultidimensionalArray CreateInitialMeans(MultidimensionalArray cellMetric, int numOfClusters) {
            System.Diagnostics.Debug.Assert(
                cellMetric.Storage.All(d => double.IsNaN(d) == false),
                "Cell metrics contains fucked up entries");

            double h_min = cellMetric.Min();
            double h_max = cellMetric.Max();
            if (h_max == double.MaxValue) {// Occurs for cells that are no REAL fluid cells (void cells, that are considered as fluid cells by the IBM tracker)
                h_max = FindSecondLargestElement(cellMetric.To1DArray());
            }

            //Console.WriteLine("Clustering: Create tanh spaced means");

            // Getting global h_min and h_max
            ilPSP.MPICollectiveWatchDog.Watch();
            h_min = h_min.MPIMin();
            h_max = h_max.MPIMax();

            if (h_min == h_max)
                h_max += 0.1 * h_max; // Dirty hack for IBM cases with equidistant grids

            // Tanh Spacing, which yields to more cell cluster for smaller cells
            var means = Grid1D.TanhSpacing(h_min, h_max, numOfClusters, 4.0, true).Reverse().ToArray();

            // Equidistant spacing, in general not the best choice
            //means = GenericBlas.Linspace(h_min, h_max, NumOfSgrd).Reverse().ToArray();

            return MultidimensionalArray.CreateWrapper(means, numOfClusters);
        }

        /// <summary>
        /// Checks for changes between two clusterings
        /// </summary>
        /// <returns>True, if clustering has changed. False, if clustering has not changed.</returns>
        public bool CheckForNewClustering(Clustering oldClustering, Clustering newClustering) {
            bool localResult = false;   // false = no reclustering needed

            if (newClustering.NumberOfClusters != oldClustering.NumberOfClusters) {
                localResult = true;
            } else {
                for (int i = 0; i < newClustering.NumberOfClusters; i++) {
                    if (!newClustering.Clusters[i].VolumeMask.Equals(oldClustering.Clusters[i].VolumeMask)) {
                        localResult = true;
                    }
                }
            }

            bool globalResult;
            unsafe {
                int localResultAsInt = localResult ? 1 : 0;
                int globalResultAsInt;
                csMPI.Raw.Allreduce((IntPtr)(&localResultAsInt), (IntPtr)(&globalResultAsInt), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.LOR, csMPI.Raw._COMM.WORLD);
                globalResult = globalResultAsInt == 1 ? true : false;
            }

            return globalResult;
        }

        /// <summary>
        /// Returns a cell metric value in every cell
        /// </summary>
        /// <returns>Cell metric as <see cref="MultidimensionalArray"/></returns>
        private MultidimensionalArray GetStableTimestepSize(SubGrid subGrid, IList<TimeStepConstraint> timeStepConstraints) {
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(subGrid.LocalNoOfCells);

            for (int subGridCell = 0; subGridCell < subGrid.LocalNoOfCells; subGridCell++) {
                int localCellIndex = subGrid.SubgridIndex2LocalCellIndex[subGridCell];
                //cellMetric[subGridCell] = timeStepConstraints.Min(c => c.GetLocalStepSize(localCellIndex, 1));    // cell metric based on smallest time step constraint
                //cellMetric[subGridCell] = 1.0 / timeStepConstraints.Sum(c => 1.0 / c.GetLocalStepSize(localCellIndex, 1));  // cell metric based on harmonic sum of time step constraints

                List<double> result = new List<double>();
                foreach (TimeStepConstraint constraint in timeStepConstraints) {
                    result.Add(constraint.GetLocalStepSize(localCellIndex, 1));
                }
                if (result.All(c => c >= double.MaxValue)) {                    // For IBM source cells: All timeStepConstraints return double.MaxValue --> No influence on clustering FUNKTIONIERT NICHT!!!!!!!!!!
                    cellMetric[subGridCell] = double.MaxValue;
                } else {
                    cellMetric[subGridCell] = 1.0 / result.Sum(c => 1.0 / c);  // cell metric based on harmonic sum of time step constraints
                }
            }

            return cellMetric;
        }

        public Clustering TuneClustering(Clustering clustering, double time, IList<TimeStepConstraint> timeStepConstraints) {
            //using (var tr = new ilPSP.Tracing.FuncTrace()) {

                // Calculate cluster time step sizes and sub-steps
                var result = GetPerCluster_dtMin_SubSteps(clustering, timeStepConstraints, 1.0e-2);
                List<int> subSteps = result.Item2;

                // Combine clusters with same number of sub-steps
                List<SubGrid> newClusters = new List<SubGrid>();
                List<int> newSubSteps = new List<int>();
                newClusters.Add(clustering.Clusters.First());
                newSubSteps.Add(subSteps.First());

                for (int i = 1; i < clustering.NumberOfClusters; i++) {
                    if (subSteps[i] == newSubSteps.Last()) {
                        // Combine both clusters and remove the previous one
                        SubGrid combinedSubGrid = new SubGrid(newClusters.Last().VolumeMask.Union(clustering.Clusters[i].VolumeMask));
                        newClusters.RemoveAt(newClusters.Count - 1);    // LATER: Implement this without removing old clusters, just adding, when combination is finished
                        newClusters.Add(combinedSubGrid);
                        if (consoleOutput) {
                            Console.WriteLine("TuneClustering: Clustering leads to clusters which are too similar. They are combined.");
                        }
                        //newSubSteps.Add(subSteps[i]);
                        //i++;
                    } else {
                        newClusters.Add(clustering.Clusters[i]);
                        newSubSteps.Add(subSteps[i]);
                    }
                }
                if (consoleOutput) {
                    for (int i = 0; i < newClusters.Count; i++) {
                        Console.WriteLine("TuneClustering:\t\t id=" + i + " -> sub-steps=" + newSubSteps[i] + "\telements=" + newClusters[i].GlobalNoOfCells);
                    }
                }

                return new Clustering(newClusters, clustering.SubGrid, newSubSteps);
            }
        //}

        public (double[], List<int>) GetPerCluster_dtHarmonicSum_SubSteps(Clustering clustering, double time, IList<TimeStepConstraint> timeStepConstraints, double eps) {
            double[] clusterDts = new double[clustering.NumberOfClusters];

            for (int i = 0; i < clustering.NumberOfClusters; i++) {
                // Use "harmonic sum" of step - sizes, see
                // WatkinsAsthanaJameson2016 for the reasoning
                double dt = 1.0 / timeStepConstraints.Sum(
                        c => 1.0 / c.GetGloballyAdmissibleStepSize(clustering.Clusters[i]));
                if (dt == 0.0) {
                    throw new ArgumentException(
                        "Time-step size is exactly zero.");
                } else if (double.IsNaN(dt)) {
                    throw new ArgumentException(
                        "Could not determine stable time-step size in sub-grid " + i + ". This indicates illegal values in some cells.");
                }

                // Restrict timesteps
                dt = Math.Min(dt, timeStepConstraints.First().Endtime - time);
                dt = Math.Min(Math.Max(dt, timeStepConstraints.First().dtMin), timeStepConstraints.First().dtMax);

                clusterDts[i] = dt;
            }

            List<int> subSteps = CalculateSubSteps(clusterDts, eps);

            for (int i = 0; i < clusterDts.Length; i++) {
                clusterDts[i] = clusterDts[0] / subSteps[i];
            }

            if (subSteps.Last() > this.MaxSubSteps && this.Restrict) {
                //#if DEBUG
                List<int> oldSubSteps = subSteps;
                double[] oldClusterDts = clusterDts;
                //#endif
                (clusterDts, subSteps) = RestrictDtsAndSubSteps(clusterDts, subSteps);
                if (consoleOutput) {
                    Console.WriteLine("### RESTRICTION OF SUB-STEPS ### (dt min)");
                    for (int i = 0; i < subSteps.Count; i++) {
                        Console.WriteLine("RestrictDtsAndSubSteps:\t id={0} -> sub-steps={1}\tdt={2:0.#######E-00} -> substeps={3}\tdt={4:0.#######E-00}", i, oldSubSteps[i], oldClusterDts[i], subSteps[i], clusterDts[i]);
                    }
                }
            }

            return (clusterDts, subSteps);
        }

        private (double[], List<int>) GetPerCluster_dtMin_SubSteps(Clustering clustering, IList<TimeStepConstraint> timeStepConstraints, double eps) {
            // Get smallest time step size of every cluster --> loop over all clusters
            // Currently: CFLFraction is not taken into account
            double[] sendDtMin = new double[clustering.NumberOfClusters];
            double[] rcvDtMin = new double[clustering.NumberOfClusters];

            MultidimensionalArray cellMetric = GetStableTimestepSize(clustering.SubGrid, timeStepConstraints);
            for (int i = 0; i < clustering.NumberOfClusters; i++) {
                double dtMin = double.MaxValue;
                CellMask volumeMask = clustering.Clusters[i].VolumeMask;
                foreach (Chunk c in volumeMask) {
                    int JE = c.JE;
                    for (int j = c.i0; j < JE; j++) {
                        dtMin = Math.Min(cellMetric[clustering.SubGrid.LocalCellIndex2SubgridIndex[j]], dtMin);
                    }
                }
                sendDtMin[i] = dtMin;
            }

            // MPI Allreduce necessary to exchange the smallest time step size of each cluster on each processor
            unsafe {
                fixed (double* pSend = sendDtMin, pRcv = rcvDtMin) {
                    csMPI.Raw.Allreduce((IntPtr)(pSend), (IntPtr)(pRcv), clustering.NumberOfClusters, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }
            }

            // Take CFLFraction into account for testing
            //for (int i = 0; i < rcvHmin.Length; i++) {
            //    rcvHmin[i] *= 0.3;
            //}

            List<int> subSteps = CalculateSubSteps(rcvDtMin, eps);

            if (subSteps.Last() > this.MaxSubSteps && this.Restrict) {
                //#if DEBUG
                List<int> oldSubSteps = subSteps;
                double[] oldRcvDtMin = rcvDtMin;
                //#endif
                (rcvDtMin, subSteps) = RestrictDtsAndSubSteps(rcvDtMin, subSteps);
                if (consoleOutput) {
                    Console.WriteLine("### RESTRICTION OF SUB-STEPS ### (dt min)");
                    for (int i = 0; i < subSteps.Count; i++) {
                        Console.WriteLine("RestrictDtsAndSubSteps:\t id={0} -> sub-steps={1}\tdt={2:0.#######E-00} -> substeps={3}\tdt={4:0.#######E-00}", i, oldSubSteps[i], oldRcvDtMin[i], subSteps[i], rcvDtMin[i]);
                    }
                }
            }

            return (rcvDtMin, subSteps);
        }

        private List<int> CalculateSubSteps(double[] timeStepSizes, double eps = 1.0e-1) {
            List<int> subSteps = new List<int>();

            for (int i = 0; i < timeStepSizes.Length; i++) {
                int result = RoundToInt(timeStepSizes[0] / timeStepSizes[i], eps);  // eps was 1.0e-1 
                Debug.Assert(result > 0 || result < 1e4, String.Format("Number of sub-steps in cluster {0} is {1}. Does this make sense?", i, result));
                subSteps.Add(result);
            }

            return subSteps;
        }

        /// <summary>
        /// CAUTION: Check if it is properly working in IBM simulations!!! Maybe, there is also a parallel bug
        /// Not needed right now, as IBM + LTS works also for a great difference in sub-steps between clusters
        /// </summary>
        /// <param name="clusterDts"></param>
        /// <param name="subSteps"></param>
        /// <returns></returns>
        private (double[], List<int>) RestrictDtsAndSubSteps(double[] clusterDts, List<int> subSteps) {
            // Restrict sub-steps
            List<int> restrictedSubSteps = new List<int>(subSteps);
            restrictedSubSteps[0] = (int)Math.Ceiling((subSteps.Last() / (double)this.MaxSubSteps));
            restrictedSubSteps[restrictedSubSteps.Count - 1] = subSteps.Last();

            for (int i = 1; i < (restrictedSubSteps.Count - 1); i++) {  // Leave first and last entry untouched
                if (subSteps[i] < restrictedSubSteps[i - 1]) {
                    restrictedSubSteps[i] = restrictedSubSteps[i - 1];
                }
            }

            // Restrict cluster time step sizes
            double[] restrictedClusterDts = new double[restrictedSubSteps.Count];
            for (int i = 0; i < restrictedClusterDts.Length; i++) {
                restrictedClusterDts[i] = clusterDts[0] / restrictedSubSteps[i];
            }

            List<int> newSubSteps = new List<int>();
            //for (int i = 0; i < restrictedClusterDts.Length; i++) {
            //    newSubSteps.Add((int)Math.Ceiling(restrictedClusterDts[0] / restrictedClusterDts[i]));
            //}
            newSubSteps = CalculateSubSteps(restrictedClusterDts, eps: 1.0e-1);

            return (restrictedClusterDts, newSubSteps);
        }

        public int RoundToInt(double number, double eps) {
            // Accounting for roundoff errors
            int result;

            if (number > Math.Floor(number) + eps) {
                result = (int)Math.Ceiling(number);
            } else {
                result = (int)Math.Floor(number);
            }

            return result;
        }

        private double FindSecondLargestElement(double[] inputElements) {
            Array.Sort(inputElements);
            Array.Reverse(inputElements);

            double largest = double.MinValue;
            double second = double.MinValue;
            for (int i = 0; i < inputElements.Length; i++)
                if (inputElements[i] > largest) {
                    second = largest;
                    largest = inputElements[i];
                } else if (inputElements[i] > second) {
                    second = inputElements[i];
                }

            return second;
        }
    }
}