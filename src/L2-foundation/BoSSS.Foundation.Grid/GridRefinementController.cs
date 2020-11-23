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

        private readonly GridData CurrentGrid;
        private readonly Partitioning CellPartitioning;
        private readonly int LocalNumberOfCells;
        private readonly int myI0;
        private readonly int GlobalNumberOfCells;
        private readonly CellMask CutCells;
        private readonly BitArray CellsNotOK2Coarsen;
        private readonly bool EnsureHighestLevelAtPeriodicBoundary;

        /// <summary>
        /// Constructor for the grid refinement controller, defines input cells
        /// </summary>
        /// <param name="currentGrid">
        /// Current grid.
        /// </param>
        /// <param name="cutCells">
        /// Cut cells will have always the max refinement level. Null is a valid input if no level-set is used.
        /// </param>
        /// <param name="cellsNotOK2Coarsen">
        /// Cells which are not allowed to be coarsend. It is not necessary to include cut cells here, as they are handled by the cutCells CellMask.
        /// </param>
        public GridRefinementController(GridData currentGrid, CellMask cutCells, CellMask cellsNotOK2Coarsen = null, bool EnsureHighestLevelAtPeriodicBoundary = false) {
            CurrentGrid = currentGrid;
            CellPartitioning = CurrentGrid.CellPartitioning;
            LocalNumberOfCells = CurrentGrid.Cells.NoOfLocalUpdatedCells;
            myI0 = CellPartitioning.i0;
            GlobalNumberOfCells = CellPartitioning.TotalLength;
            CutCells = cutCells;
            if (CutCells == null)
                CutCells = CellMask.GetEmptyMask(CurrentGrid);
            if (cellsNotOK2Coarsen == null)
                cellsNotOK2Coarsen = CellMask.GetEmptyMask(CurrentGrid);
            CellsNotOK2Coarsen = cellsNotOK2Coarsen.Union(CutCells).GetBitMask();
            this.EnsureHighestLevelAtPeriodicBoundary = EnsureHighestLevelAtPeriodicBoundary;
        }

        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on a refinement indicator. If the calculation should work parallel it is recommend to use <see cref="ComputeGridChange(GridData, List{Tuple{int, CellMask}}, out List{int}, out List{int[]})"/>.
        /// </summary>
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
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public bool ComputeGridChange(Func<int, int, int> levelIndicator, out List<int> cellsToRefine, out List<int[]> cellsToCoarsen) {
            //throw new Exception("Legacy method. Will be deleted at 01.11.2020");
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship();
            int[] globalDesiredLevel = GetGlobalDesiredLevel(levelIndicator, globalCellNeigbourship);
            cellsToRefine = GetCellsToRefine(globalDesiredLevel);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(globalDesiredLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen);
            cellsToCoarsen = GetCoarseningCells(coarseningClusters);

            bool anyChangeInGrid = (cellsToRefine.Count() == 0 && cellsToCoarsen.Count() == 0) ? false : true;
            return (anyChangeInGrid);
        }

        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on the max refinement level provided by the calling solver. This method is fully parallized.
        /// </summary>
        /// <param name="cellsToRefine">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="cellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CellDesiredLevel">
        /// The refinement level of each local cell.
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public bool ComputeGridChange(int[] CellDesiredLevel, out List<int> cellsToRefine, out List<int[]> cellsToCoarsen) {
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship();
            BitArray cutCellsWithNeighbours = GetGlobalNearBand(globalCellNeigbourship);
            int[] globalDesiredLevel = GetGlobalDesiredLevel(cutCellsWithNeighbours, CellDesiredLevel, globalCellNeigbourship);
            int[] globalRefinementLevel = GetGlobalRefinementLevel(globalDesiredLevel, globalCellNeigbourship);
            cellsToRefine = GetCellsToRefine(globalRefinementLevel);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(globalRefinementLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen);
            cellsToCoarsen = GetCoarseningCells(coarseningClusters);

            bool anyChangeInGrid = (cellsToRefine.Count() == 0 && cellsToCoarsen.Count() == 0) ? false : true;
            return anyChangeInGrid.MPIOr();
        }



        /// <summary>
        /// Computes the global cell neighbourship of all cell. Returns an jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </summary>
        private int[][] GetGlobalCellNeigbourship() {
            int[] i0 = CellPartitioning.GetI0s();

            long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;
            int[][] localCellNeighbourship = new int[LocalNumberOfCells][];
            for (int j = 0; j < localCellNeighbourship.Length; j++) {
                Tuple<int, int, int>[] cellNeighbours = CurrentGrid.GetCellNeighboursViaEdges(j);
                localCellNeighbourship[j] = new int[cellNeighbours.Length];
                for (int i = 0; i < cellNeighbours.Length; i++) {
                    localCellNeighbourship[j][i] = cellNeighbours[i].Item1;
                }
                for (int i = 0; i < localCellNeighbourship[j].Length; i++) {
                    if (localCellNeighbourship[j][i] < LocalNumberOfCells)
                        localCellNeighbourship[j][i] = localCellNeighbourship[j][i] + myI0;
                    else
                        localCellNeighbourship[j][i] = (int)externalCellsGlobalIndices[localCellNeighbourship[j][i] - LocalNumberOfCells];
                }
            }

            int[][][] exchangeCellNeighbourship = localCellNeighbourship.MPIGatherO(0);
            exchangeCellNeighbourship = exchangeCellNeighbourship.MPIBroadcast(0);

            int[][] globalCellNeigbourship = new int[GlobalNumberOfCells][];
            for (int m = 0; m < CurrentGrid.MpiSize; m++) {
                for (int j = 0; j < exchangeCellNeighbourship[m].Length; j++) {
                    globalCellNeigbourship[j + i0[m]] = exchangeCellNeighbourship[m][j];
                }
            }
            return globalCellNeigbourship;
        }

        /// <summary>
        /// Returns the cutcells + neighbours on a global level.
        /// </summary>
        /// <param name="globalCellNeighbourship">
        /// </param>
        private BitArray GetGlobalNearBand(int[][] globalCellNeighbourship) {
            if (CutCells == null)
                return null;
            BitArray localCutCells = CutCells.GetBitMask();
            BitArray globalCutCells = new BitArray(GlobalNumberOfCells);
            for (int j = 0; j < LocalNumberOfCells; j++) {
                if (localCutCells[j]) {
                    globalCutCells[j + myI0] = true;
                    for (int i = 0; i < globalCellNeighbourship[j + myI0].Length; i++) {
                        globalCutCells[globalCellNeighbourship[j + myI0][i]] = true;
                    }
                }
            }
            for (int j = 0; j < globalCutCells.Length; j++) {
                globalCutCells[j] = globalCutCells[j].MPIOr();
            }
            return globalCutCells;
        }

        /// <summary>
        /// Writes the max desired level of the specified cells into an int-array (mpi global).
        /// </summary>
        /// <param name="cutCellsWithNeighbours"></param>
        /// <param name="CellRefinementLevel"></param>
        private int[] GetGlobalDesiredLevel(BitArray cutCellsWithNeighbours, int[] CellRefinementLevel, int[][] globalCellNeigbourship) {
            int[] i0 = CellPartitioning.GetI0s();
            int myRank = CurrentGrid.MpiRank;
            int[] globalDesiredLevel = new int[GlobalNumberOfCells];
            int levelSetMaxLevel = 1;

            for (int j = 0; j < LocalNumberOfCells; j++) {
                if (CellRefinementLevel[j] > 0 && globalDesiredLevel[j + i0[myRank]] <= CellRefinementLevel[j]) {
                    globalDesiredLevel[j + i0[myRank]] = CellRefinementLevel[j];
                    if (globalDesiredLevel[j + i0[myRank]] > levelSetMaxLevel)
                        levelSetMaxLevel = globalDesiredLevel[j + i0[myRank]];
                }
            }
            levelSetMaxLevel = levelSetMaxLevel.MPIMax();

            int noOfPeriodicCells = 0;
            if (EnsureHighestLevelAtPeriodicBoundary) {
                for (int j = 0; j < LocalNumberOfCells; j++) {
                    CellFaceTag[] cellFaceTags = CurrentGrid.Cells.GetCell(j).CellFaceTags;
                    if (cellFaceTags != null) {
                        for (int c = 0; c < cellFaceTags.Length; c++) {
                            if (cellFaceTags[c].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                                globalDesiredLevel[j + i0[myRank]] = levelSetMaxLevel;
                                noOfPeriodicCells += 1;
                                break;
                            }
                        }
                    }
                }
            }
            int noOfCutCells = 0;
            if (cutCellsWithNeighbours != null) {
                for (int j = 0; j < GlobalNumberOfCells; j++) {
                    if (cutCellsWithNeighbours[j]) {
                        globalDesiredLevel[j] = levelSetMaxLevel;
                        noOfCutCells += 1;
                    }
                }
            }
            return globalDesiredLevel.MPIMax();
        }

        /// <summary>
        /// Computes the level indicator for each cell (mpi global). 
        /// </summary>
        /// <param name="GlobalDesiredLevel">
        /// Int-array with the length of the global no of cells. Contains the desired max level of each cell.
        /// </param>
        /// <param name="globalNeighbourship"></param>
        private int[] GetGlobalRefinementLevel(int[] GlobalDesiredLevel, int[][] globalNeighbourship) {
            int[] globalRefinementLevel = new int[GlobalNumberOfCells];
            int[] globalCurrentLevel = GetGlobalCurrentLevel();
            Debug.Assert(globalRefinementLevel.Length == globalCurrentLevel.Length);

            for (int j = 0; j < LocalNumberOfCells; j++) {
                int globalCellIndex = j + myI0;
                if (globalRefinementLevel[globalCellIndex] < GlobalDesiredLevel[globalCellIndex] && globalCurrentLevel[globalCellIndex] <= GlobalDesiredLevel[globalCellIndex]) {
                    if (globalCurrentLevel[globalCellIndex] == GlobalDesiredLevel[globalCellIndex])
                        globalRefinementLevel[globalCellIndex] = globalCurrentLevel[globalCellIndex];
                    else
                        globalRefinementLevel[globalCellIndex] = globalCurrentLevel[globalCellIndex] + 1;
                    GetRefinementLevelRecursive(globalCellIndex, globalRefinementLevel[globalCellIndex] - 1, globalNeighbourship, globalRefinementLevel, globalCurrentLevel);
                }
                else if (globalRefinementLevel[globalCellIndex] <= globalCurrentLevel[globalCellIndex] - 1 && globalCurrentLevel[globalCellIndex] > 0)
                    globalRefinementLevel[globalCellIndex] = globalCurrentLevel[globalCellIndex] - 1;
            }
            return globalRefinementLevel.MPIMax();
        }

        /// <summary>
        /// Recursive computation for the desired level of each global cell.
        /// </summary>
        /// <param name="globalCellIndex">
        /// The global index of the current cell.
        /// </param>
        /// <param name="desiredLevel"></param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        /// <param name="globalRefinementLevel"></param>
        private void GetRefinementLevelRecursive(int globalCellIndex, int desiredLevel, int[][] globalNeighbourship, int[] globalRefinementLevel, int[] globalCurrentLevel) {
            if (desiredLevel <= 0)
                return;
            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                int jNeigh = globalNeighbourship[globalCellIndex][j];
                if (globalRefinementLevel[globalCellIndex] < desiredLevel && globalCurrentLevel[globalCellIndex] <= desiredLevel) {
                    if (globalCurrentLevel[globalCellIndex] == desiredLevel)
                        globalRefinementLevel[globalCellIndex] = globalCurrentLevel[globalCellIndex];
                    else
                        globalRefinementLevel[globalCellIndex] = globalCurrentLevel[globalCellIndex] + 1;
                    GetRefinementLevelRecursive(jNeigh, globalRefinementLevel[globalCellIndex] - 1, globalNeighbourship, globalRefinementLevel, globalCurrentLevel);
                }
            }
        }

        /// <summary>
        /// Calculates the desired level for each global cell. Note that the desired level can only increase by 1 compared to the current level of the cell.
        /// </summary>
        /// <param name="levelIndicator">
        /// The level indicator func provided by the calling solver.
        /// </param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        private int[] GetGlobalDesiredLevel(Func<int, int, int> levelIndicator, int[][] globalNeighbourship) {
            int[] globalDesiredLevel = new int[GlobalNumberOfCells];

            int[] globalCurrentLevel = GetGlobalCurrentLevel();

            for (int globalCellIndex = 0; globalCellIndex < GlobalNumberOfCells; globalCellIndex++) {
                int localCellIndex = globalCellIndex - myI0;

                int currentLevel_j = globalCurrentLevel[globalCellIndex];

                int desiredLevel_j = (localCellIndex < LocalNumberOfCells && localCellIndex >= 0) ? levelIndicator(localCellIndex, currentLevel_j) : 0;
                desiredLevel_j = desiredLevel_j.MPIMax();

                if (globalDesiredLevel[globalCellIndex] < desiredLevel_j) {
                    globalDesiredLevel[globalCellIndex] = desiredLevel_j;
                    RefineNeighboursRecursive(currentLevel_j, globalDesiredLevel, globalCellIndex, desiredLevel_j - 1, globalNeighbourship, globalCurrentLevel);
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
        private void RefineNeighboursRecursive(int currentLevel, int[] DesiredLevel, int globalCellIndex, int DesiredLevelNeigh, int[][] globalNeighbourship, int[] globalCurrentLevel) {
            if (DesiredLevelNeigh <= 0)
                return;
            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                int jNeigh = globalNeighbourship[globalCellIndex][j];
                if (currentLevel < DesiredLevelNeigh && DesiredLevel[jNeigh] < DesiredLevelNeigh) {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(currentLevel, DesiredLevel, jNeigh, DesiredLevelNeigh - 1, globalNeighbourship, globalCurrentLevel);
                }
            }
        }

        /// <summary>
        /// Gets the current refinement level for each global cell.
        /// </summary>
        private int[] GetGlobalCurrentLevel() {
            int[] localCurrentLevel = new int[LocalNumberOfCells];
            int[] globalCurrentLevel = new int[GlobalNumberOfCells];
            int[] i0 = CellPartitioning.GetI0s();

            for (int j = 0; j < LocalNumberOfCells; j++) {
                localCurrentLevel[j] = CurrentGrid.Cells.GetCell(j).RefinementLevel;
            }
            int[][] exchangeGlobalCurrentLevel = localCurrentLevel.MPIGatherO(0);
            exchangeGlobalCurrentLevel = exchangeGlobalCurrentLevel.MPIBroadcast(0);

            for (int m = 0; m < CurrentGrid.MpiSize; m++) {
                for (int j = 0; j < exchangeGlobalCurrentLevel[m].Length; j++) {
                    globalCurrentLevel[j + i0[m]] = exchangeGlobalCurrentLevel[m][j];
                }
            }

            return globalCurrentLevel;
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a int-list.
        /// </summary>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        private List<int> GetCellsToRefine(int[] globalDesiredLevel) {
            int i0 = CurrentGrid.CellPartitioning.i0;
            List<int> cellToRefine = new List<int>();
            for (int j = 0; j < CurrentGrid.Cells.NoOfLocalUpdatedCells; j++) {
                int ActualLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = globalDesiredLevel[i0 + j];
                if (ActualLevel_j < DesiredLevel_j)
                    cellToRefine.Add(j);
            }
            return cellToRefine;
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a int-list.
        /// </summary>
        /// </param>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        /// <param name="globalCellNeigbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        private BitArray GetCellsOk2Coarsen(int[] globalDesiredLevel, int[][] globalCellNeigbourship) {
            BitArray oK2Coarsen = new BitArray(LocalNumberOfCells);
            int myRank = CurrentGrid.MpiRank;
            int[][] cellNeighbours = CurrentGrid.Cells.CellNeighbours;

            for (int globalCellIndex = myI0; globalCellIndex < myI0 + LocalNumberOfCells; globalCellIndex++) {
                int localCellIndex = globalCellIndex - myI0;
                int ActualLevel_j = CurrentGrid.Cells.GetCell(localCellIndex).RefinementLevel;

                if (ActualLevel_j > globalDesiredLevel[globalCellIndex] && ActualLevel_j > globalCellNeigbourship[globalCellIndex].Select(neighbourIndex => globalDesiredLevel[neighbourIndex]).Max()) {
                    oK2Coarsen[localCellIndex] = true;
                }
            }

            if (CellsNotOK2Coarsen != null) {
                for (int j = 0; j < LocalNumberOfCells; j++) {
                    if (CellsNotOK2Coarsen[j])
                        oK2Coarsen[j] = false;
                }
            }

            return oK2Coarsen;
        }

        /// <summary>
        /// Gets all cells to be coarsend. This is not mpi-parallel, because coarsening is not allowed over process boundaries.
        /// </summary>
        /// <param name="oK2Coarsen">
        /// A BitArray of all cells which should be coarsend.
        /// </param>
        private int[][] FindCoarseningClusters(BitArray oK2Coarsen) {
            int JE = CurrentGrid.Cells.Count;

            //int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;
            List<Cell> temp = new List<Cell>();
            List<int> tempCoarseningCluster = new List<int>();

            int RecursionDepth = -1;
            if (CurrentGrid.SpatialDimension == 1)
                RecursionDepth = 1;
            else if (CurrentGrid.SpatialDimension == 2)
                RecursionDepth = 2;
            else if (CurrentGrid.SpatialDimension == 3)
                RecursionDepth = 4;
            else
                throw new NotSupportedException();

            int[][] CoarseningCluster = new int[JE][];

            BitArray marker = new BitArray(JE);
            for (int j = 0; j < LocalNumberOfCells; j++) {
                if (marker[j])
                    continue;
                if (!oK2Coarsen[j])
                    continue;

                Cell currentCell = CurrentGrid.Cells.GetCell(j);
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
                CurrentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighborsEdges, out _);
                CurrentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] CellNeighborsVertices, out _);
                List<int> cellNeighbours = new List<int>();
                cellNeighbours.AddRange(CellNeighborsEdges);
                cellNeighbours.AddRange(CellNeighborsVertices);
                bool complete = false;
                FindCoarseiningClusterRecursive(j, RecursionDepth, cellNeighbours.ToArray(), marker, CurrentGrid.Cells.GetCell, currentRefinmentLevel, searchClusterID, currentClusterSize, tempCoarseningCluster, ref complete);
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
        private void FindCoarseiningClusterRecursive(int j, int MaxRecursionDeph,
            int[] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int currentRefinementLevel, int searchClusterID, int clusterSize,
            List<int> coarseningCluster, ref bool complete) {
            if (!coarseningCluster.Contains(j))
                throw new Exception("Error in coarsening algortihm: Coarsening cluster does not contain a cell with the local ID: " + j);

            foreach (int neighbourCellIndex in CellNeighbours) {

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
                    FindCoarseiningClusterRecursive(neighbourCellIndex, MaxRecursionDeph - 1,
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
        private List<int[]> GetCoarseningCells(int[][] coarseningClusters) {
            List<int[]> coarseningCells = new List<int[]>();
            for (int j = 0; j < LocalNumberOfCells; j++) {
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

