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
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

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

            public bool RestrictDtsAndSubSteps {
                get;
                private set;
            }

            public int MaxSubSteps {
                get;
                private set;
            }

            public Clustering(List<SubGrid> clusters, SubGrid subGrid, bool restrictDtsAndSubSteps, List<int> subStepsInitial = null, int maxSubSteps = 50) {
                this.Clusters = clusters;
                this.SubGrid = subGrid;
                this.RestrictDtsAndSubSteps = restrictDtsAndSubSteps;
                this.SubStepsInitial = subStepsInitial;
                if (this.RestrictDtsAndSubSteps) {
                    this.MaxSubSteps = 50;
                }
            }
        }

        /// <summary>
        /// Information about the grid
        /// </summary>
        private IGridData gridData;

        /// <summary>
        /// Constructor for the grid clustering
        /// </summary>
        /// <param name="gridData">Information about the grid</param>
        public Clusterer(IGridData gridData) {
            this.gridData = gridData;
        }

        /// <summary>
        /// Creates the sub-grids of the clustering
        /// </summary>     
        /// <param name="numOfClusters">Number of clusters</param>
        /// <returns>A list of sub-grids</returns>
        public Clustering CreateClustering(int numOfClusters, IList<TimeStepConstraint> timeStepConstraints, bool restrictDtsAndSubSteps, SubGrid subGrid = null) {
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
#if DEBUG
                    Console.WriteLine("Sub-grid/Cluster " + (i + 1) + ", with mean value " + Kmean.Means[i] + ", is empty and not used anymore!");
#endif
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

            // Generating the sub-grids
            List<SubGrid> clusters = new List<SubGrid>(counter);
            for (int i = 0; i < numOfClusters; i++) {
                // Generating only the sub-grids which are not empty
                if (noOfCellsPerCluster[i] != 0) {
                    BitArray ba = baMatrix[i];
                    clusters.Add(new SubGrid(new CellMask(gridData, ba)));
                }
            }
#if DEBUG
            for (int i = 0; i < clusters.Count; i++) {
                Console.WriteLine("CreateClustering: id=" + i + " -> elements=" + clusters[i].GlobalNoOfCells);
            }
#endif
            return new Clustering(clusters, subGrid, restrictDtsAndSubSteps);
        }

        /// <summary>
        /// Creates an array with an tanh spaced distribution of the mean
        /// values between maximum and minimum value of a given cell metric, 
        /// e.g., minimal distance between two nodes in a cell <see cref="GridData.CellData.h_min"/>
        /// </summary>
        /// <param name="cellMetric">Given cell metric</param>
        /// <returns>Double[] with the length of the number of given sub-grids></returns>
        private MultidimensionalArray CreateInitialMeans(MultidimensionalArray cellMetric, int numOfClusters) {
            System.Diagnostics.Debug.Assert(
                cellMetric.Storage.All(d => double.IsNaN(d) == false),
                "Cell metrics contains fucked up entries");

            double h_min = cellMetric.Min();
            double h_max = cellMetric.Max();
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
                cellMetric[subGridCell] = 1.0 / timeStepConstraints.Sum(c => 1.0 / c.GetLocalStepSize(localCellIndex, 1));  // cell metric based on harmonic sum of time step constraints
            }

            return cellMetric;
        }

        public Clustering TuneClustering(Clustering clustering, double time, IList<TimeStepConstraint> timeStepConstraints) {
            // Calculate cluster time step sizes and sub-steps
            var result = GetPerCluster_dtMin_SubSteps(clustering, timeStepConstraints, 1.0e-2);
            List<int> numOfSubSteps = result.Item2;

            // Combine clusters with same number of sub-steps
            List<SubGrid> newClusters = new List<SubGrid>();
            List<int> newSubSteps = new List<int>();

            for (int i = 0; i < clustering.NumberOfClusters; i++) {
                if (i < clustering.NumberOfClusters - 1 && numOfSubSteps[i] == numOfSubSteps[i + 1]) {
                    // Combine both clusters and remove the previous one
                    SubGrid combinedSubGrid = new SubGrid(clustering.Clusters[i].VolumeMask.Union(clustering.Clusters[i + 1].VolumeMask));
                    newClusters.Add(combinedSubGrid);
#if DEBUG
                    Console.WriteLine("TuneClustering: Clustering leads to sub-grids which are too similar, i.e. they have the same number of local time steps. They are combined.");
#endif
                    newSubSteps.Add(numOfSubSteps[i]);
                    i++;
                } else {
                    newClusters.Add(clustering.Clusters[i]);
                    newSubSteps.Add(numOfSubSteps[i]);
                }
            }
#if DEBUG
            for (int i = 0; i < newClusters.Count; i++) {
                Console.WriteLine("TuneClustering:\t id=" + i + " -> sub-steps=" + newSubSteps[i] + " and elements=" + newClusters[i].GlobalNoOfCells);
            }
#endif
            return new Clustering(newClusters, clustering.SubGrid, clustering.RestrictDtsAndSubSteps, newSubSteps);
        }

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

            if (subSteps.Last() > clustering.MaxSubSteps && clustering.RestrictDtsAndSubSteps) {
                (clusterDts, subSteps) = RestrictDtsAndSubSteps(clusterDts, subSteps, clustering.MaxSubSteps);
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

            if (subSteps.Last() > clustering.MaxSubSteps && clustering.RestrictDtsAndSubSteps) {
                (rcvDtMin, subSteps) = RestrictDtsAndSubSteps(rcvDtMin, subSteps, clustering.MaxSubSteps);
            }

            return (rcvDtMin, subSteps);
        }

        private List<int> CalculateSubSteps(double[] timeStepSizes, double eps = 1.0e-1) {
            List<int> subSteps = new List<int>();

            for (int i = 0; i < timeStepSizes.Length; i++) {
                subSteps.Add(RoundToInt(timeStepSizes[0] / timeStepSizes[i], eps));    // eps was 1.0e-1 
            }

            return subSteps;
        }

        private (double[], List<int>) RestrictDtsAndSubSteps(double[] clusterDts, List<int> subSteps, int maxSubSteps) {
            // Restrict sub-steps
            List<int> restrictedSubSteps = new List<int>(subSteps);
            restrictedSubSteps[0] = (int)Math.Ceiling((subSteps.Last() / (double)maxSubSteps));
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
            for (int i = 0; i < restrictedClusterDts.Length; i++) {
                newSubSteps.Add((int)Math.Ceiling(restrictedClusterDts[0] / restrictedClusterDts[i]));
            }
#if DEBUG
            Console.WriteLine("### RESTRICTION OF SUB-STEPS ###");
            for (int i = 0; i < restrictedClusterDts.Length; i++) {
                Console.WriteLine("RestrictDtsAndSubSteps:\t id=" + i + " -> sub-steps=" + newSubSteps[i]);
            }
#endif
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
    }
}