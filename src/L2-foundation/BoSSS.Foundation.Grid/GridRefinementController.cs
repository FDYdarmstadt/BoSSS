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

using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Basic algorithms for refinement and coarsening.
    /// </summary>
    public class GridRefinementController {
        
        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on a refinement indicator.
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
        /// <param name="clustersOfCellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CutCells">
        /// If not null, a mask of cells in which coarsening is forbidden (usually cut-cells);
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public static bool ComputeGridChange(GridData currentGrid, CellMask CutCells, Func<int, int, int> levelIndicator, out List<int> cellsToRefine, out List<int[]> clustersOfCellsToCoarsen)
        {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;

            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);
            int[] globalDesiredLevel = new int[globalJ];

            FindRefiningCells(currentGrid, levelIndicator, globalDesiredLevel, globalCellNeigbourship);
            int[] localDesiredLevel = GetLocalDesiredLevel(currentGrid, globalDesiredLevel);
            cellsToRefine = GetRefiningCells(currentGrid, localDesiredLevel);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, CutCells, globalDesiredLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            clustersOfCellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = CheckForRefineOrCoarsening(cellsToRefine, clustersOfCellsToCoarsen);

            return (anyChangeInGrid);
        }

        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on the cut cells and the maximum refinementLevel.
        /// </summary>
        /// <param name="currentGrid">
        /// Current grid.
        /// </param>
        /// <param name="maxRefinementLevel">
        /// predefined maximum refinement level
        /// </param>
        /// <param name="cellsToRefine">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="clustersOfCellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CutCells">
        /// If not null, a mask of cells in which coarsening is forbidden (usually cut-cells);
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public static bool ComputeGridChange(GridData currentGrid, CellMask CutCells, int maxRefinementLevel, out List<int> cellsToRefine, out List<int[]> clustersOfCellsToCoarsen)
        { 
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;

            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);
            int[] globalDesiredLevel = new int[globalJ];
            int[] levelIndicator = GetLevelIndicator(currentGrid, CutCells, maxRefinementLevel, globalCellNeigbourship);

            FindRefiningCells(currentGrid, levelIndicator, globalDesiredLevel, globalCellNeigbourship);
            int[] localDesiredLevel = GetLocalDesiredLevel(currentGrid, globalDesiredLevel);
            cellsToRefine = GetRefiningCells(currentGrid, localDesiredLevel);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, CutCells, globalDesiredLevel, globalCellNeigbourship);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            clustersOfCellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = CheckForRefineOrCoarsening(cellsToRefine, clustersOfCellsToCoarsen);

            return (anyChangeInGrid);
        }

        private static int[] GetLevelIndicator(GridData currentGrid, CellMask CutCells, int maxRefinementLevel, int[][] globalCellNeighbourship)
        {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int globalJ = cellPartitioning.TotalLength;
            int[] levelIndicator = new int[globalJ];
            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);
            for (int j = 0; j < oldJ; j++)
            {
                int globalIndex = j + i0;
                if (CutCells.Contains(j) && globalCurrentLevel[globalIndex] < maxRefinementLevel)
                {
                    levelIndicator[globalIndex] = globalCurrentLevel[globalIndex] + 1;
                    GetNeighbourLevelIndicatorRecursive(globalCellNeighbourship, levelIndicator, globalIndex, globalCurrentLevel[globalIndex]);
                }  
            }

            for (int j = 0; j < globalJ; j++)
            {
                levelIndicator[j] = levelIndicator[j].MPIMax();
            }

            return levelIndicator;
        }

        private static void GetNeighbourLevelIndicatorRecursive(int[][] globalCellNeighbourship, int[] levelIndicator, int globalIndex, int refinementLevel)
        {
            foreach (int jNeigh in globalCellNeighbourship[globalIndex])
            {
                if (levelIndicator[jNeigh] < refinementLevel)
                    levelIndicator[jNeigh] = refinementLevel;
                if (refinementLevel > 1)
                    GetNeighbourLevelIndicatorRecursive(globalCellNeighbourship, levelIndicator, jNeigh, refinementLevel - 1);
            }
        }

        private static int[] GetLocalDesiredLevel(GridData currentGrid, int[] globalDesiredLevel)
        {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            int myRank = currentGrid.MpiRank;
            int[] localDesiredLevel = new int[oldJ];
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int[] i0 = cellPartitioning.GetI0s();
            for (int j = 0; j < oldJ; j++)
            {
                localDesiredLevel[j] = globalDesiredLevel[i0[myRank] + j];
            }

            return localDesiredLevel;
        }

        private static int[][] GetGlobalCellNeigbourship(GridData currentGrid)
        {
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;
            int[] i0 = cellPartitioning.GetI0s();
            int[][] localCellNeighbourship = new int[J][];
            for (int j = 0; j < J; j++)
            {
                currentGrid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] neighbourCells, out _);
                int[] globalIndexNeighbourCells = neighbourCells.CloneAs();
                for (int i = 0; i < neighbourCells.Length; i++)
                {
                    globalIndexNeighbourCells[i] = (int)currentGrid.Parallel.GetGlobalCellIndex(neighbourCells[i]);
                }
                localCellNeighbourship[j] = globalIndexNeighbourCells;
            }
            int[][][] exchangeCellNeighbourship = localCellNeighbourship.MPIGatherO(0);
            exchangeCellNeighbourship = exchangeCellNeighbourship.MPIBroadcast(0);
            int[][] globalCellNeigbourship = new int[globalJ][];
            for(int m = 0; m < currentGrid.MpiSize; m++)
            {
                for (int j = 0; j < exchangeCellNeighbourship[m].Length; j++)
                {
                    globalCellNeigbourship[j + i0[m]] = exchangeCellNeighbourship[m][j];
                }
            }
            return globalCellNeigbourship;
        }

        private static void FindRefiningCells(GridData currentGrid, Func<int, int, int> LevelIndicator, int[] desiredLevel, int[][] globalNeighbourship)
        {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int[] globalLevelIndicator = new int[globalJ];
            int[] localLevelIndicator = new int[J];
            int[] i0 = cellPartitioning.GetI0s();
            
            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);

            for (int j = 0; j < globalJ; j++)
            {
                if (j == 15)
                    Console.WriteLine("jsngj");
                int CurrentLevel_j = globalCurrentLevel[j];
                int desiredLevel_j;
                desiredLevel_j = GetGlobalLevelIndicator(currentGrid, j, LevelIndicator);
                desiredLevel_j = desiredLevel_j.MPIMax();

                if (desiredLevel[j] < desiredLevel_j)
                {
                    desiredLevel[j] = desiredLevel_j;
                    RefineNeighboursRecursive(CurrentLevel_j, desiredLevel, j, desiredLevel_j - 1, globalNeighbourship);
                }
            }
        }

        private static void FindRefiningCells(GridData currentGrid, int[] LevelIndicator, int[] desiredLevel, int[][] globalNeighbourship)
        {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int[] globalLevelIndicator = new int[globalJ];
            int[] localLevelIndicator = new int[J];
            int[] i0 = cellPartitioning.GetI0s();

            int[] globalCurrentLevel = GetGlobalCurrentLevel(currentGrid);

            for (int j = 0; j < globalJ; j++)
            {
                if (j == 15)
                    Console.WriteLine("jsngj");
                int CurrentLevel_j = globalCurrentLevel[j];
                int desiredLevel_j;
                desiredLevel_j = LevelIndicator[j];

                if (desiredLevel[j] < desiredLevel_j)
                {
                    desiredLevel[j] = desiredLevel_j;
                    RefineNeighboursRecursive(CurrentLevel_j, desiredLevel, j, desiredLevel_j - 1, globalNeighbourship);
                }
            }
        }

        private static int[] GetGlobalCurrentLevel(GridData currentGrid)
        {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int[] localCurrentLevel = new int[J];
            int[] globalCurrentLevel = new int[globalJ];
            int[] i0 = cellPartitioning.GetI0s();

            for (int j = 0; j < J; j++)
            {
                localCurrentLevel[j] = currentGrid.Cells.GetCell(j).RefinementLevel;
            }

            int[][] exchangeGlobalCurrentLevel = localCurrentLevel.MPIGatherO(0);
            exchangeGlobalCurrentLevel = exchangeGlobalCurrentLevel.MPIBroadcast(0);

            for (int m = 0; m < currentGrid.MpiSize; m++)
            {
                for(int j = 0; j < exchangeGlobalCurrentLevel[m].Length; j++)
                {
                    globalCurrentLevel[j + i0[m]] = exchangeGlobalCurrentLevel[m][j];
                }
            }

            return globalCurrentLevel;
        }

        private static int GetGlobalLevelIndicator(GridData currentGrid, int globalCellIndex, Func<int, int, int> LevelIndicator)
        {
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int i0 = currentGrid.CellPartitioning.i0;
            int localCellIndex = globalCellIndex - i0;

            int currentLevel_j = (localCellIndex < J && localCellIndex >= 0) ? currentGrid.Cells.GetCell(localCellIndex).RefinementLevel : 0;
            int localLevelIndicator = (localCellIndex < J && localCellIndex >= 0) ? LevelIndicator(localCellIndex, currentLevel_j) : 0;
            return localLevelIndicator.MPIMax();
        }

        static void RefineNeighboursRecursive(int currentLevel, int[] DesiredLevel, int j, int DesiredLevelNeigh, int[][] globalNeighbourship)
        {
            if (DesiredLevelNeigh <= 0)
                return;

            foreach (int jNeigh in globalNeighbourship[j])
            {
                if (jNeigh == 15)
                    Console.WriteLine("jsngj");
                if (currentLevel < DesiredLevelNeigh && DesiredLevel[jNeigh] < DesiredLevelNeigh)
                {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(currentLevel, DesiredLevel, jNeigh, DesiredLevelNeigh - 1, globalNeighbourship);
                }
            }
        }

        private static List<int> GetRefiningCells(GridData CurrentGrid, int[] DesiredLevel)
        {
            int oldJ = CurrentGrid.Cells.NoOfLocalUpdatedCells;
            List<int> CellsToRefineList = new List<int>();
            for (int j = 0; j < oldJ; j++)
            {
                int ActualLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = DesiredLevel[j];

                if (ActualLevel_j < DesiredLevel_j)
                    CellsToRefineList.Add(j);
            }
            return CellsToRefineList;
        }

        private static BitArray GetCellsOk2Coarsen(GridData currentGrid, CellMask cutCells, int[] globalDesiredLevel, int[][] globalCellNeigbourship)
        {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            BitArray oK2Coarsen = new BitArray(oldJ);
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int i0 = cellPartitioning.i0;
            int myRank = currentGrid.MpiRank;
            int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;

            for (int globalCellIndex = i0; globalCellIndex < i0 + oldJ; globalCellIndex++)
            {
                int localCellIndex = globalCellIndex - i0;
                int ActualLevel_j = currentGrid.Cells.GetCell(localCellIndex).RefinementLevel;

                if (ActualLevel_j > globalDesiredLevel[globalCellIndex]
                    && globalDesiredLevel[globalCellIndex] >= globalCellNeigbourship[globalCellIndex].Select(neighbourIndex => globalDesiredLevel[neighbourIndex]).Max() - 1)
                {
                    oK2Coarsen[localCellIndex] = true;
                }
            }

            if (cutCells != null)
            {
                foreach (int j in cutCells.ItemEnum)
                {
                    oK2Coarsen[j] = false;
                }
            }

            return oK2Coarsen;
        }

        static int[][] FindCoarseningClusters(BitArray oK2Coarsen, GridData currentGrid)
        {
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
            for (int j = 0; j < J; j++)
            {
                if (marker[j])
                    continue;
                if (!oK2Coarsen[j])
                    continue;

                Cell currentCell = currentGrid.Cells.GetCell(j);
                int currentRefinmentLevel = currentCell.RefinementLevel;

                if (currentRefinmentLevel == 0)
                {
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
                FindCoarseiningClusterRecursive(j, RecursionDepht, cellNeighbours, marker, currentGrid.Cells.GetCell, currentRefinmentLevel, searchClusterID, currentClusterSize, tempCoarseningCluster, ref complete);
                foreach (int jC in tempCoarseningCluster)
                {
                    marker[jC] = true;
                }

                if (complete)
                {
                    if (!tempCoarseningCluster.Any(jC => oK2Coarsen[jC] == false))
                    {
                        int[] CC = tempCoarseningCluster.ToArray();
                        foreach (int jC in CC)
                        {
                            Debug.Assert(CoarseningCluster[jC] == null);
                            CoarseningCluster[jC] = CC;
                        }
                    }
                    else
                    {
                        //Console.WriteLine("Not ok to coarsen.");
                    }
                }
            }
            return CoarseningCluster;
        }

        private static void FindCoarseiningClusterRecursive(int j, int MaxRecursionDeph,
            int[][] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int currentRefinementLevel, int searchClusterID, int ClusterSize,
            List<int> CoarseiningCluster, ref bool complete)
        {
            if (!CoarseiningCluster.Contains(j))
                throw new Exception("Error in coarsening algortihm: Coarsening cluster does not contain a cell with the local ID: " + j);

            int[] Neighs = CellNeighbours[j];
            foreach (int neighbourCellIndex in Neighs)
            {
                if (marker[neighbourCellIndex] == true)
                    continue;

                if (CoarseiningCluster.Contains(neighbourCellIndex))
                    continue;

                Cell neighbourCell = GetCell(neighbourCellIndex);
                if (neighbourCell.RefinementLevel != currentRefinementLevel)
                    continue;

                if (neighbourCell.CoarseningClusterID != searchClusterID)
                    continue;

                CoarseiningCluster.Add(neighbourCellIndex);

                if (CoarseiningCluster.Count == ClusterSize)
                {
                    complete = true;
#if DEBUG
                    foreach (int j1 in CoarseiningCluster)
                    {
                        Cell Cell_j1 = GetCell(j1);
                        Debug.Assert(Cell_j1.RefinementLevel > 0);
                        Debug.Assert(Cell_j1.RefinementLevel == currentRefinementLevel);
                        Debug.Assert(Cell_j1.CoarseningClusterID == searchClusterID);
                        Debug.Assert(Cell_j1.CoarseningClusterSize == ClusterSize);
                    }
#endif
                    return;
                }

                if (MaxRecursionDeph > 0)
                    FindCoarseiningClusterRecursive(neighbourCellIndex, MaxRecursionDeph - 1,
                        CellNeighbours, marker, GetCell, currentRefinementLevel, searchClusterID, ClusterSize, CoarseiningCluster, ref complete);

                if (complete)
                {
                    return;
                }
            }
        }

        private static List<int[]> GetCoarseningCells(GridData currentGrid, int[][] coarseningClusters)
        {
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            List<int[]> coarseningCells = new List<int[]>();
            for (int j = 0; j < oldJ; j++)
            {
                if (coarseningClusters[j] != null)
                {
                    Debug.Assert(coarseningClusters[j].Contains(j));
                    if (j == coarseningClusters[j].Min())
                    {
                        coarseningCells.Add(coarseningClusters[j]);
                    }
                }
            }
            return coarseningCells;
        }

        private static bool CheckForRefineOrCoarsening(List<int> cellsToRefine, List<int[]> cellsToCoarsen)
        {
            bool isRefining = true;
            if (cellsToRefine.Count == 0 && cellsToCoarsen.Count == 0)
                isRefining = false;
            return isRefining;
        }
    }
}

