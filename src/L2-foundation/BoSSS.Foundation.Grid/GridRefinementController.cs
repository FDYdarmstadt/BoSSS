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
        /// <param name="CutCells">A bit array of all cut cells on the local process</param>
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
        public static bool ComputeGridChange(GridData currentGrid, BitArray CutCells, List<Tuple<int, BitArray>> CellsMaxRefineLevel, out List<int> cellsToRefine, out List<int[]> cellsToCoarsen) {
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);
            BitArray cutCellsWithNeighbours = GetCutCellNeighbours(currentGrid, CutCells, globalCellNeigbourship);
            int[] CellsWithMaxRefineLevel = GetAllCellsWithMaxRefineLevel(currentGrid, cutCellsWithNeighbours, CellsMaxRefineLevel);
            int[] levelIndicator = GetGlobalLevelIndicator(currentGrid, CellsWithMaxRefineLevel, globalCellNeigbourship);
            cellsToRefine = GetCellsToRefine(currentGrid, levelIndicator);

            BitArray notOk2Coarsen = GetCellsNotOk2Coarsen(currentGrid, CellsMaxRefineLevel);
            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, notOk2Coarsen, levelIndicator, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            cellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = (cellsToRefine.Count() == 0 && cellsToCoarsen.Count() == 0) ? false : true;
            bool[] exchangeGridChange = anyChangeInGrid.MPIGatherO(0);
            exchangeGridChange = exchangeGridChange.MPIBroadcast(0);
            for (int m = 0; m < exchangeGridChange.Length; m++) {
                if (exchangeGridChange[m])
                    anyChangeInGrid = true;
            }
            return (anyChangeInGrid);
        }

        private static BitArray GetCutCellNeighbours(GridData currentGrid, BitArray localCutCells, int[][] globalCellNeighbourship) {
            if (localCutCells == null)
                return null;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;
            int[] i0 = cellPartitioning.GetI0s();
            BitArray[] exchangeCutCells = localCutCells.MPIGatherO(0);
            exchangeCutCells = exchangeCutCells.MPIBroadcast(0);
            BitArray globalCutCells = new BitArray(globalJ);   
            for (int m = 0; m < exchangeCutCells.Length; m++) {
                for (int j = 0; j < exchangeCutCells[m].Length; j++) {
                    globalCutCells[j + i0[m]] = exchangeCutCells[m][j];
                    if (globalCutCells[j + i0[m]])
                        for (int i = 0; i < globalCellNeighbourship[j + i0[m]].Length; i++) {
                            globalCutCells[globalCellNeighbourship[j + i0[m]][i]] = true;
                        }
                }
            }
            return globalCutCells;
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
        private static int[] GetAllCellsWithMaxRefineLevel(GridData currentGrid, BitArray cutCells, List<Tuple<int, BitArray>> CellsMaxRefineLevel) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int[] i0 = cellPartitioning.GetI0s();
            int globalJ = cellPartitioning.TotalLength;
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            int[] localCellsMaxRefineLvl = new int[localJ];
            int levelSetMaxLevel = 0;

            for (int i = 0; i < CellsMaxRefineLevel.Count(); i++) {
                BitArray currentArray = CellsMaxRefineLevel[i].Item2;
                for (int j = 0; j < localJ; j++) {
                    if (currentArray[j] && localCellsMaxRefineLvl[j] <= CellsMaxRefineLevel[i].Item1) {
                        localCellsMaxRefineLvl[j] = CellsMaxRefineLevel[i].Item1;
                        if (localCellsMaxRefineLvl[j] > levelSetMaxLevel)
                            levelSetMaxLevel = localCellsMaxRefineLvl[j];
                    }
                }
            }
            levelSetMaxLevel = levelSetMaxLevel.MPIMax();
            int[][] exchangeCellsMaxRefineLvl = localCellsMaxRefineLvl.MPIGatherO(0);
            exchangeCellsMaxRefineLvl = exchangeCellsMaxRefineLvl.MPIBroadcast(0);
            int[] globalCellsMaxRefineLvl = new int[globalJ];          
                for (int m = 0; m < exchangeCellsMaxRefineLvl.Length; m++) {
                    for (int j = 0; j < exchangeCellsMaxRefineLvl[m].Length; j++) {
                        globalCellsMaxRefineLvl[j + i0[m]] = exchangeCellsMaxRefineLvl[m][j];
                    }
                }
            if (cutCells != null) {
                for (int j = 0; j < globalJ; j++) {
                    if (cutCells[j]) {
                        globalCellsMaxRefineLvl[j] = levelSetMaxLevel;
                    }
                }
            }
            return globalCellsMaxRefineLvl;
        }

        private static BitArray GetCellsNotOk2Coarsen(GridData currentGrid, List<Tuple<int, BitArray>> CellsMaxRefineLevel) {
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            BitArray localCellsNotOk2Coarsen = new BitArray(localJ);
            bool anyCellsNotOk = false;
            for (int i = 0; i < CellsMaxRefineLevel.Count(); i++) {
                if (CellsMaxRefineLevel[i].Item1 == -1) {
                    anyCellsNotOk = true;
                    BitArray currentArray = CellsMaxRefineLevel[i].Item2;
                    for (int j = 0; j < localJ; j++) {
                        if (!localCellsNotOk2Coarsen[j])
                            localCellsNotOk2Coarsen[j] = currentArray[j];
                    }
                }
            }
            if (anyCellsNotOk)
                return localCellsNotOk2Coarsen;
            else
                return null;
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
            int[][] localCellNeighbourship = new int[localJ][];
            for (int j = 0; j < localCellNeighbourship.Length; j++) {
                currentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] CellNeighbours, out _);
                localCellNeighbourship[j] = CellNeighbours;
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
        /// Computes the level indicator for each cell (mpi global). 
        /// </summary>
        /// <param name="currentGrid">
        /// </param>
        /// <param name="CellsWithMaxRefineLevel">
        /// Int-array with the length of the global no of cells. Contains the desired max level of each cell.
        /// </param>
        /// <param name="globalNeighbourship"></param>
        private static int[] GetGlobalLevelIndicator(GridData currentGrid, int[] CellsWithMaxRefineLevel, int[][] globalNeighbourship) {
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int globalJ = cellPartitioning.TotalLength;
            int[] levelIndicator = new int[globalJ];
            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);

            for (int j = 0; j < localJ; j++) {
                int globalIndex = j + i0;
                if (levelIndicator[globalIndex] < CellsWithMaxRefineLevel[globalIndex] && globalCurrentLevel[globalIndex] <= CellsWithMaxRefineLevel[globalIndex]) {
                    if (globalCurrentLevel[globalIndex] == CellsWithMaxRefineLevel[globalIndex])
                        levelIndicator[globalIndex] = globalCurrentLevel[globalIndex];
                    else 
                        levelIndicator[globalIndex] = globalCurrentLevel[globalIndex] + 1;
                    GetLevelIndicatiorRecursive(globalIndex, levelIndicator[globalIndex] - 1, globalNeighbourship, levelIndicator);
                }
                else if (levelIndicator[globalIndex] <= globalCurrentLevel[globalIndex] - 1 && globalCurrentLevel[globalIndex] > 0)
                    levelIndicator[globalIndex] = globalCurrentLevel[globalIndex] - 1;
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
        /// Recursive computation for the desired level of each global cell.
        /// </summary>
        /// <param name="globalCellIndex">
        /// The global index of the current cell.
        /// </param>
        /// <param name="LevelIndNeighbour"></param>
        /// <param name="globalNeighbourship">
        /// Jaggerd int-array where the first index refers to the current cell and the second one to the neighbour cells.
        /// </param>
        /// <param name="levelIndicator"></param>
        static void GetLevelIndicatiorRecursive(int globalCellIndex, int LevelIndNeighbour, int[][] globalNeighbourship, int[] levelIndicator) {
            if (LevelIndNeighbour <= 0)
                return;
            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                int jNeigh = globalNeighbourship[globalCellIndex][j];
                if (levelIndicator[jNeigh] < LevelIndNeighbour) {
                    levelIndicator[jNeigh] = LevelIndNeighbour;
                    GetLevelIndicatiorRecursive(jNeigh, LevelIndNeighbour - 1, globalNeighbourship, levelIndicator);
                }
            }
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
        private static BitArray GetCellsOk2Coarsen(GridData currentGrid, BitArray cellsNotOK2Coarsen, int[] globalDesiredLevel, int[][] globalCellNeigbourship) {
            int localJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            BitArray oK2Coarsen = new BitArray(localJ);
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int myRank = currentGrid.MpiRank;
            int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;

            for (int globalCellIndex = i0; globalCellIndex < i0 + localJ; globalCellIndex++) {
                int localCellIndex = globalCellIndex - i0;
                int ActualLevel_j = currentGrid.Cells.GetCell(localCellIndex).RefinementLevel;

                if (ActualLevel_j > globalDesiredLevel[globalCellIndex] && ActualLevel_j > globalCellNeigbourship[globalCellIndex].Select(neighbourIndex => globalDesiredLevel[neighbourIndex]).Max()) {
                    oK2Coarsen[localCellIndex] = true;
                }
            }

            if (cellsNotOK2Coarsen != null) {
                for (int j = 0; j < localJ; j++) {
                    if (cellsNotOK2Coarsen[j])
                        oK2Coarsen[j] = false;
                }
            }

            return oK2Coarsen;
        }

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

                if (ActualLevel_j > globalDesiredLevel[globalCellIndex] && ActualLevel_j > globalCellNeigbourship[globalCellIndex].Select(neighbourIndex => globalDesiredLevel[neighbourIndex]).Max()) {
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

            //int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;
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
                currentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighborsEdges, out _);
                currentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] CellNeighborsVertices, out _);
                List<int> cellNeighbours = new List<int>();
                cellNeighbours.AddRange(CellNeighborsEdges);
                cellNeighbours.AddRange(CellNeighborsVertices);
                bool complete = false;
                FindCoarseiningClusterRecursive(currentGrid, j, RecursionDepht, cellNeighbours.ToArray(), marker, currentGrid.Cells.GetCell, currentRefinmentLevel, searchClusterID, currentClusterSize, tempCoarseningCluster, ref complete);
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
            int[] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int currentRefinementLevel, int searchClusterID, int clusterSize,
            List<int> coarseningCluster, ref bool complete) {
            if (!coarseningCluster.Contains(j))
                throw new Exception("Error in coarsening algortihm: Coarsening cluster does not contain a cell with the local ID: " + j);

            int J = currentGrid.Cells.NoOfLocalUpdatedCells;

            //int[] Neighs = CellNeighbours[j];
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

