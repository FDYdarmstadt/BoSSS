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

            public Clustering(List<SubGrid> clusters, SubGrid subGrid) {
                this.Clusters = clusters;
                this.SubGrid = subGrid;
            }
        }

        /// <summary>
        /// Information about the grid
        /// </summary>
        private IGridData gridData;

        /// <summary>
        /// The time step constraints that the Clustering is based on
        /// </summary>
        private IList<TimeStepConstraint> timeStepConstraints;

        /// <summary>
        /// Constructor for the grid clustering
        /// </summary>
        /// <param name="gridData">Information about the grid</param>
        /// <param name="timeStepConstraints">Time step constraings used as cell metric for the clustering</param>
        public Clusterer(IGridData gridData, IList<TimeStepConstraint> timeStepConstraints) {
            this.gridData = gridData;
            this.timeStepConstraints = timeStepConstraints;
        }

        /// <summary>
        /// Creates the sub-grids of the clustering
        /// </summary>     
        /// <param name="numOfClusters">Number of clusters</param>
        /// <returns>A list of sub-grids</returns>
        public Clustering CreateClustering(int numOfClusters, SubGrid subGrid = null) {
            if (subGrid == null) {
                subGrid = new SubGrid(CellMask.GetFullMask(gridData));
            }

            // Attention: numOfCells can equal all local cells or only the local cells of a subgrid,
            // e.g., the fluid cells in an IBM simulation
            int numOfCells = subGrid.LocalNoOfCells;

            MultidimensionalArray cellMetric = GetStableTimestepSize(subGrid);
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
                    System.Console.WriteLine("Sub-grid/Cluster " + (i + 1) + ", with mean value " + Kmean.Means[i] + ", is empty and not used anymore!");
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

            return new Clustering(clusters, subGrid);
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
        /// <param name="oldClustering">A clustering which should be compared to</param>
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
        public MultidimensionalArray GetStableTimestepSize(SubGrid subGrid){
            MultidimensionalArray cellMetric = MultidimensionalArray.Create(subGrid.LocalNoOfCells);

            for (int subGridCell = 0; subGridCell < subGrid.LocalNoOfCells; subGridCell++) {
                int localCellIndex = subGrid.SubgridIndex2LocalCellIndex[subGridCell];
                cellMetric[subGridCell] = this.timeStepConstraints.Min(c => c.GetLocalStepSize(localCellIndex, 1));
            }

            return cellMetric;
        }
    }
}