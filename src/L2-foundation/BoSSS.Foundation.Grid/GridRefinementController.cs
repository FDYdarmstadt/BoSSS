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

using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Basic algorithms for refinement and coarsening.
    /// </summary>
    public class GridRefinementController {

        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on a refinement indicator. If the calculation should work parallel it is recommend to use <see cref="ComputeGridChange(GridData, List{Tuple{int, CellMask}}, out List{int}, out List{int[]})"/>.
        /// </summary>
        /// <param name="currentGrid">
        /// Current grid.
        /// </param>
        /// <param name="levelIndicator">
        /// Mapping from (local cell index, current refinement level) to desired refinement level for the respective cell,
        /// see <see cref="Cell.RefinementLevel"/>.
        /// </param>
        /// <param name="cellsToRefine">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="cellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CutCells">
        /// If not null, a mask of cells in which coarsening is forbidden (usually cut-cells);
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public static bool ComputeGridChange(GridData currentGrid, CellMask CutCells, Func<int, int, int> levelIndicator, out List<int> cellsToRefine, out List<int[]> cellsToCoarsen) {
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);
            int[] globalDesiredLevel = GetGlobalDesiredLevel(currentGrid, levelIndicator, globalCellNeigbourship);
            cellsToRefine = GetCellsToRefine(currentGrid, globalDesiredLevel);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, CutCells, globalDesiredLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            cellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = (cellsToRefine.Count() == 0 && cellsToCoarsen.Count() == 0) ? false : true;
            return (anyChangeInGrid);
        }

        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on the max refinement level provided by the calling solver. This method is fully parallized.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="cellsToRefine">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="cellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CellsMaxRefineLevel">
        /// All cells, sorted in CellMask with their desired maximum level of refinement. Due to its nature as a list it is possible to define multiple different mask.
        /// ATTENTION: If to many CellMasks are defined the refinement algorithm might be slow.
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public static bool ComputeGridChange(GridData currentGrid, List<Tuple<int, CellMask>> CellsMaxRefineLevel, out List<int> cellsToRefine, out List<int[]> cellsToCoarsen) {
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);

            int[] CellsWithMaxRefineLevel = GetAllCellsWithMaxRefineLevel(currentGrid, CellsMaxRefineLevel);
            int[] levelIndicator = GetGlobalLevelIndicator(currentGrid, CellsWithMaxRefineLevel);
            int[] globalDesiredLevel = GetGlobalDesiredLevel(currentGrid, levelIndicator, globalCellNeigbourship);
            cellsToRefine = GetCellsToRefine(currentGrid, globalDesiredLevel);

            // to be refactored
            CellMask AllRelevantCells = CellMask.GetEmptyMask(currentGrid);
            if (CellsMaxRefineLevel.Count() > 0)
                AllRelevantCells = CellsMaxRefineLevel[0].Item2;
            for (int i = 1; i < CellsMaxRefineLevel.Count(); i++) {
                AllRelevantCells.Union(CellsMaxRefineLevel[i].Item2);
            }

            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, AllRelevantCells, globalDesiredLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            cellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = (cellsToRefine.Count() == 0 && cellsToCoarsen.Count() == 0) ? false : true;
            return (anyChangeInGrid);
        }

        /// <summary>
        /// Writes the max desired level of the specified cells into an int-array (mpi global).
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="CellsMaxRefineLevel">
        /// All cells, sorted in CellMask with their desired maximum level of refinement. Due to its nature as a list it is possible to define multiple different mask.
        /// ATTENTION: If to many CellMasks are defined the refinement algorithm might be slow.
        /// </param>
        private static int[] GetAllCellsWithMaxRefineLevel(GridData currentGrid, List<Tuple<int, CellMask>> CellsMaxRefineLevel) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int[] i0 = cellPartitioning.GetI0s();
            int globalJ = cellPartitioning.TotalLength;
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            int[] localCellsMaxRefineLvl = new int[localJ];

            for (int i = 0; i < CellsMaxRefineLevel.Count(); i++) {
                for (int j = 0; j < localJ; j++) {
                    if (CellsMaxRefineLevel[i].Item2.Contains(j) && localCellsMaxRefineLvl[j] < CellsMaxRefineLevel[i].Item1) {
                        localCellsMaxRefineLvl[j] = CellsMaxRefineLevel[i].Item1;
                    }
                }
            }

            int[][] exchangeCellsMaxRefineLvl = localCellsMaxRefineLvl.MPIGatherO(0);
            exchangeCellsMaxRefineLvl = exchangeCellsMaxRefineLvl.MPIBroadcast(0);

            int[] globalCellsMaxRefineLvl = new int[globalJ];
            for (int m = 0; m < exchangeCellsMaxRefineLvl.Length; m++) {
                for (int j = 0; j < exchangeCellsMaxRefineLvl[m].Length; j++) {
                    globalCellsMaxRefineLvl[j + i0[m]] = exchangeCellsMaxRefineLvl[m][j];
                }
            }

            return globalCellsMaxRefineLvl;
        }

        /// <summary>
        /// Computes the global cell neighbourship of all cell. Returns an jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        private static int[][] GetGlobalCellNeigbourship(GridData currentGrid) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            int[] i0 = cellPartitioning.GetI0s();
            int local_i0 = cellPartitioning.i0;

            long[] externalCellsGlobalIndices = currentGrid.iParallel.GlobalIndicesExternalCells;
            int[][] localCellNeighbourship = currentGrid.Cells.CellNeighbours;
            for (int j = 0; j < localCellNeighbourship.Length; j++) {
                for (int i = 0; i < localCellNeighbourship[j].Length; i++) {
                    if (localCellNeighbourship[j][i] < localJ)
                        localCellNeighbourship[j][i] = localCellNeighbourship[j][i] + local_i0;
                    else
                        localCellNeighbourship[j][i] = (int)externalCellsGlobalIndices[localCellNeighbourship[j][i] - localJ];
                }
            }

            int[][][] exchangeCellNeighbourship = localCellNeighbourship.MPIGatherO(0);
            exchangeCellNeighbourship = exchangeCellNeighbourship.MPIBroadcast(0);

            int[][] globalCellNeigbourship = new int[globalJ][];
            for (int m = 0; m < currentGrid.MpiSize; m++) {
                for (int j = 0; j < exchangeCellNeighbourship[m].Length; j++) {
                    globalCellNeigbourship[j + i0[m]] = exchangeCellNeighbourship[m][j];
                }
            }
            return globalCellNeigbourship;
        }

        /// <summary>
        /// Computes the level indicator for each cell (mpi global). The level indicator is then used in <see cref="GetGlobalDesiredLevel(GridData, int[], int[][])"/> to calculate the actual desired level,
        /// which is not allowed to be larger than currentLevel + 1.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="CellsWithMaxRefineLevel">
        /// Int-array with the length of the global no of cells. Contains the desired max level of each cell.
        /// </param>
        private static int[] GetGlobalLevelIndicator(GridData currentGrid, int[] CellsWithMaxRefineLevel) {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int globalJ = cellPartitioning.TotalLength;
            int[] levelIndicator = new int[globalJ];
            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);
            for (int j = 0; j < oldJ; j++) {
                int globalIndex = j + i0;
                if (CellsWithMaxRefineLevel[globalIndex] != 0 && globalCurrentLevel[globalIndex] < CellsWithMaxRefineLevel[globalIndex]) {
                    levelIndicator[globalIndex] = globalCurrentLevel[globalIndex] + 1;
                }
            }

            for (int j = 0; j < globalJ; j++) {
                levelIndicator[j] = levelIndicator[j].MPIMax();
            }

            return levelIndicator;
        }

        /// <summary>
        /// Calculates the desired level for each global cell. Note that the desired level can only increase by 1 compared to the current level of the cell.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="levelIndicator">
        /// The level indicator func provided by the calling solver.
        /// </param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        private static int[] GetGlobalDesiredLevel(GridData currentGrid, Func<int, int, int> levelIndicator, int[][] globalNeighbourship) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int i0 = cellPartitioning.i0;
            int[] globalDesiredLevel = new int[globalJ];

            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);

            for (int globalCellIndex = 0; globalCellIndex < globalJ; globalCellIndex++) {
                int localCellIndex = globalCellIndex - i0;

                int currentLevel_j = globalCurrentLevel[globalCellIndex];

                int desiredLevel_j = (localCellIndex < J && localCellIndex >= 0) ? levelIndicator(localCellIndex, currentLevel_j) : 0;
                desiredLevel_j = desiredLevel_j.MPIMax();

                if (globalDesiredLevel[globalCellIndex] < desiredLevel_j) {
                    globalDesiredLevel[globalCellIndex] = desiredLevel_j;
                    RefineNeighboursRecursive(currentLevel_j, globalDesiredLevel, globalCellIndex, desiredLevel_j - 1, globalNeighbourship, globalCurrentLevel);
                }
            }
            return globalDesiredLevel;
        }

        /// <summary>
        /// Calculates the desired level for each global cell. Note that the desired level can only increase by 1 compared to the current level of the cell.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="levelIndicator">
        /// The level indicator array provided by <see cref="GetGlobalLevelIndicator(GridData, int[])"/>
        /// </param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        private static int[] GetGlobalDesiredLevel(GridData currentGrid, int[] levelIndicator, int[][] globalNeighbourship) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;

            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);
            int[] globalDesiredLevel = new int[globalJ];

            for (int j = 0; j < globalJ; j++) {
                int currentLevel_j = globalCurrentLevel[j];
                int desiredLevel_j;
                desiredLevel_j = levelIndicator[j];

                if (globalDesiredLevel[j] < desiredLevel_j) {
                    globalDesiredLevel[j] = desiredLevel_j;
                    RefineNeighboursRecursive(currentLevel_j, globalDesiredLevel, j, desiredLevel_j - 1, globalNeighbourship, globalCurrentLevel);
                }
            }
            return globalDesiredLevel;
        }

        /// <summary>
        /// Recursive computation for the desired level of each global cell.
        /// </summary>
        /// <param name="currentLevel">
        /// The current level of the cell with the index globalCellIndex.
        /// </param>
        /// <param name="DesiredLevel">
        /// Int-array of the desired level of each global cell.
        /// </param>
        /// <param name="globalCellIndex">
        /// The global index of the current cell.
        /// </param>
        /// <param name="DesiredLevelNeigh">
        /// The desired level of the neighbours of the current cell.
        /// </param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        /// <param name="globalCurrentLevel"></param>
        static void RefineNeighboursRecursive(int currentLevel, int[] DesiredLevel, int globalCellIndex, int DesiredLevelNeigh, int[][] globalNeighbourship, int[] globalCurrentLevel) {
            if (DesiredLevelNeigh <= 0)
                return;

            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                int jNeigh = globalNeighbourship[globalCellIndex][j];
                currentLevel = globalCurrentLevel[jNeigh];
                if (currentLevel < DesiredLevelNeigh && DesiredLevel[jNeigh] < DesiredLevelNeigh) {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(currentLevel, DesiredLevel, jNeigh, DesiredLevelNeigh - 1, globalNeighbourship, globalCurrentLevel);
                }
            }
        }

        /// <summary>
        /// Gets the current refinement level for each global cell.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        private static int[] GetGlobalCurrentLevel(GridData currentGrid) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int[] localCurrentLevel = new int[J];
            int[] globalCurrentLevel = new int[globalJ];
            int[] i0 = cellPartitioning.GetI0s();

            for (int j = 0; j < J; j++) {
                localCurrentLevel[j] = currentGrid.Cells.GetCell(j).RefinementLevel;
            }

            int[][] exchangeGlobalCurrentLevel = localCurrentLevel.MPIGatherO(0);
            exchangeGlobalCurrentLevel = exchangeGlobalCurrentLevel.MPIBroadcast(0);

            for (int m = 0; m < currentGrid.MpiSize; m++) {
                for (int j = 0; j < exchangeGlobalCurrentLevel[m].Length; j++) {
                    globalCurrentLevel[j + i0[m]] = exchangeGlobalCurrentLevel[m][j];
                }
            }

            return globalCurrentLevel;
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a int-list.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        private static List<int> GetCellsToRefine(GridData currentGrid, int[] globalDesiredLevel) {
            int i0 = currentGrid.CellPartitioning.i0;
            List<int> cellToRefine = new List<int>();
            for (int j = 0; j < currentGrid.Cells.NoOfLocalUpdatedCells; j++) {
                int ActualLevel_j = currentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = globalDesiredLevel[i0 + j];
                if (ActualLevel_j < DesiredLevel_j)
                    cellToRefine.Add(j);
            }
            return cellToRefine;
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a int-list.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="cellsNotOK2Coarsen">
        /// A CellMask of all cells which should never be coarsend, e.g. all cut cells.
        /// </param>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        /// <param name="globalCellNeigbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        private static BitArray GetCellsOk2Coarsen(GridData currentGrid, CellMask cellsNotOK2Coarsen, int[] globalDesiredLevel, int[][] globalCellNeigbourship) {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            BitArray oK2Coarsen = new BitArray(oldJ);
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int myRank = currentGrid.MpiRank;
            int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;

            for (int globalCellIndex = i0; globalCellIndex < i0 + oldJ; globalCellIndex++) {
                int localCellIndex = globalCellIndex - i0;
                int ActualLevel_j = currentGrid.Cells.GetCell(localCellIndex).RefinementLevel;

                if (ActualLevel_j > globalDesiredLevel[globalCellIndex]
                    && globalDesiredLevel[globalCellIndex] >= globalCellNeigbourship[globalCellIndex].Select(neighbourIndex => globalDesiredLevel[neighbourIndex]).Max() - 1) {
                    oK2Coarsen[localCellIndex] = true;
                }
            }

            if (cellsNotOK2Coarsen != null) {
                foreach (int j in cellsNotOK2Coarsen.ItemEnum) {
                    oK2Coarsen[j] = false;
                }
            }

            return oK2Coarsen;
        }

        /// <summary>
        /// Gets all cells to be coarsend. This is not mpi-parallel, because coarsening is not allowed over process boundaries.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="oK2Coarsen">
        /// A BitArray of all cells which should be coarsend.
        /// </param>
        static int[][] FindCoarseningClusters(BitArray oK2Coarsen, GridData currentGrid) {
            int JE = currentGrid.Cells.Count;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;

            int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;
            List<Cell> temp = new List<Cell>();
            List<int> tempCoarseningCluster = new List<int>();

            int RecursionDepht = -1;
            if (currentGrid.SpatialDimension == 1)
                RecursionDepht = 1;
            else if (currentGrid.SpatialDimension == 2)
                RecursionDepht = 2;
            else if (currentGrid.SpatialDimension == 3)
                RecursionDepht = 4;
            else
                throw new NotSupportedException();

            int[][] CoarseningCluster = new int[JE][];

            BitArray marker = new BitArray(JE);
            for (int j = 0; j < J; j++) {
                if (marker[j])
                    continue;
                if (!oK2Coarsen[j])
                    continue;

                Cell currentCell = currentGrid.Cells.GetCell(j);
                int currentRefinmentLevel = currentCell.RefinementLevel;

                if (currentRefinmentLevel == 0) {
                    marker[j] = true;
                    continue;
                }

                temp.Clear();
                temp.Add(currentCell);

                tempCoarseningCluster.Clear();
                tempCoarseningCluster.Add(j);

                int searchClusterID = currentCell.CoarseningClusterID;
                int currentClusterSize = currentCell.CoarseningClusterSize;
                Debug.Assert(searchClusterID > 0);

                bool complete = false;
                FindCoarseiningClusterRecursive(currentGrid, j, RecursionDepht, cellNeighbours, marker, currentGrid.Cells.GetCell, currentRefinmentLevel, searchClusterID, currentClusterSize, tempCoarseningCluster, ref complete);
                foreach (int jC in tempCoarseningCluster) {
                    marker[jC] = true;
                }

                if (complete) {
                    if (!tempCoarseningCluster.Any(jC => oK2Coarsen[jC] == false)) {
                        int[] CC = tempCoarseningCluster.ToArray();
                        foreach (int jC in CC) {
                            Debug.Assert(CoarseningCluster[jC] == null);
                            CoarseningCluster[jC] = CC;
                        }
                    }
                }
            }
            return CoarseningCluster;
        }

        /// <summary>
        /// Recursive computing of all coarsening clusters.
        /// </summary>
        /// <param name="j">
        /// The local index of the current cell.
        /// </param>
        /// <param name="MaxRecursionDeph">
        /// The max recursion depht, depending on the spatial dimension of the grid.
        /// </param>
        /// <param name="CellNeighbours">
        /// </param>
        /// <param name="marker">
        /// Returns true for all cells not ok to coarsen.
        /// </param>
        /// <param name="GetCell">
        /// A func to find the current cell in the current grid.
        /// </param>
        /// <param name="currentRefinementLevel">
        /// The refinement level of the current cell.
        /// </param>
        /// <param name="searchClusterID">
        /// The ID of the current cluster.
        /// </param>
        /// <param name="clusterSize">
        /// The size of the current cluster.
        /// </param>
        /// <param name="coarseningCluster">
        /// The current coarsening cluster
        /// </param>
        /// <param name="complete">
        /// Bool is true if recursion has found all cells of the current cluster. 
        /// </param>
        private static void FindCoarseiningClusterRecursive(GridData currentGrid, int j, int MaxRecursionDeph,
            int[][] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int currentRefinementLevel, int searchClusterID, int clusterSize,
            List<int> coarseningCluster, ref bool complete) {
            if (!coarseningCluster.Contains(j))
                throw new Exception("Error in coarsening algortihm: Coarsening cluster does not contain a cell with the local ID: " + j);

            int J = currentGrid.Cells.NoOfLocalUpdatedCells;

            int[] Neighs = CellNeighbours[j];
            foreach (int neighbourCellIndex in Neighs) {
                if (neighbourCellIndex > J)
                    continue;

                if (marker[neighbourCellIndex] == true)
                    continue;

                if (coarseningCluster.Contains(neighbourCellIndex))
                    continue;

                Cell neighbourCell = GetCell(neighbourCellIndex);
                if (neighbourCell.RefinementLevel != currentRefinementLevel)
                    continue;

                if (neighbourCell.CoarseningClusterID != searchClusterID)
                    continue;

                coarseningCluster.Add(neighbourCellIndex);

                if (coarseningCluster.Count == clusterSize) {
                    complete = true;
#if DEBUG
                    foreach (int j1 in coarseningCluster)
                    {
                        Cell Cell_j1 = GetCell(j1);
                        Debug.Assert(Cell_j1.RefinementLevel > 0);
                        Debug.Assert(Cell_j1.RefinementLevel == currentRefinementLevel);
                        Debug.Assert(Cell_j1.CoarseningClusterID == searchClusterID);
                        Debug.Assert(Cell_j1.CoarseningClusterSize == clusterSize);
                    }
#endif
                    return;
                }

                if (MaxRecursionDeph > 0)
                    FindCoarseiningClusterRecursive(currentGrid, neighbourCellIndex, MaxRecursionDeph - 1,
                        CellNeighbours, marker, GetCell, currentRefinementLevel, searchClusterID, clusterSize, coarseningCluster, ref complete);

                if (complete) {
                    return;
                }
            }
        }

        /// <summary>
        /// Gets all coarsening clusters and writes them to a list.
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="coarseningClusters">
        /// All coarsening clusters (1st index) with their respective cells (2nd index)
        /// </param>
        private static List<int[]> GetCoarseningCells(GridData currentGrid, int[][] coarseningClusters) {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            List<int[]> coarseningCells = new List<int[]>();
            for (int j = 0; j < oldJ; j++) {
                if (coarseningClusters[j] != null) {
                    Debug.Assert(coarseningClusters[j].Contains(j));
                    if (j == coarseningClusters[j].Min()) {
                        coarseningCells.Add(coarseningClusters[j]);
                    }
                }
            }
            return coarseningCells;
        }
    }
}

