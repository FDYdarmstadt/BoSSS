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
            int oldJE = currentGrid.Cells.NoOfExternalCells + currentGrid.Cells.NoOfLocalUpdatedCells;
            int oldJ = currentGrid.Cells.NoOfLocalUpdatedCells;
            int myRank = currentGrid.MpiRank;
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int[] i0 = cellPartitioning.GetI0s();
            int globalJ = cellPartitioning.TotalLength;
            Debugger.Launch();
            int[][] globalCellNeigbourship = GetGlobalCellNeigbourship(currentGrid);
            int[] globalDesiredLevel = new int[globalJ];
            FindRefiningCells(currentGrid, levelIndicator, globalDesiredLevel, globalCellNeigbourship);
            int[] localDesiredLevel = GetLocalDesiredLevel(currentGrid, globalDesiredLevel);
            cellsToRefine = GetRefiningCells(currentGrid, localDesiredLevel);

            globalDesiredLevel.MPIExchange(currentGrid);

            BitArray oK2Coarsen = GetCellsOk2Coarsen(currentGrid, CutCells, globalDesiredLevel);
            int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen, currentGrid);
            clustersOfCellsToCoarsen = GetCoarseningCells(currentGrid, coarseningClusters);

            bool anyChangeInGrid = CheckForRefineOrCoarsening(cellsToRefine, clustersOfCellsToCoarsen);

            return (anyChangeInGrid);
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
            long[] externalCellsGlobalIndices = currentGrid.iParallel.GlobalIndicesExternalCells;
            int[][] localCellNeighbourship = new int[J][];
            Cell cl = currentGrid.Cells.GetCell(0);
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

        private static void FindRefiningCells(GridData currentGrid, Func<int, int, int> LevelIndicator, int[] DesiredLevel, int[][] globalNeighbourship)
        {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int globalJ = cellPartitioning.TotalLength;
            int[] globalLevelIndicator = new int[globalJ];
            int[] localLevelIndicator = new int[J];
            int[] i0 = cellPartitioning.GetI0s();
            for (int j = 0; j < globalJ; j++)
            {
                int CurrentLevel_j = currentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = GetGlobalLevelIndicator(currentGrid, j, LevelIndicator);

                if (DesiredLevel[j] < DesiredLevel_j)
                {
                    DesiredLevel[j] = DesiredLevel_j;
                    RefineNeighboursRecursive(CurrentLevel_j, DesiredLevel, j, DesiredLevel_j - 1, globalNeighbourship);
                }
            }
        }

        private static int GetGlobalLevelIndicator(GridData currentGrid, int currentCellIndex, Func<int, int, int> LevelIndicator)
        {
            int J = currentGrid.Cells.NoOfLocalUpdatedCells;
            int currentLevel_j = currentCellIndex < J ? currentGrid.Cells.GetCell(currentCellIndex).RefinementLevel : 0;
            int localLevelIndicator = currentCellIndex < J ? LevelIndicator(currentCellIndex, currentLevel_j) : 0;
            return localLevelIndicator.MPIMax();
        }

        static void RefineNeighboursRecursive(int currentLevel, int[] DesiredLevel, int j, int DesiredLevelNeigh, int[][] globalNeighbourship)
        {
            if (DesiredLevelNeigh <= 0)
                return;

            foreach (int jNeigh in globalNeighbourship[j])
            {
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

        private static BitArray GetCellsOk2Coarsen(GridData CurrentGrid, CellMask CutCells, int[] DesiredLevel)
        {
            int oldJ = CurrentGrid.Cells.NoOfLocalUpdatedCells;
            BitArray Ok2Coarsen = new BitArray(oldJ);
            for (int j = 0; j < oldJ; j++)
            {
                int ActualLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = DesiredLevel[j];

                int[][] CellNeighbours = CurrentGrid.Cells.CellNeighbours;

                if (ActualLevel_j > DesiredLevel_j && DesiredLevel_j >= CellNeighbours[j].Select(cn => DesiredLevel[cn]).Max() - 1)
                {
                    Ok2Coarsen[j] = true;
                }
            }
            if (CutCells != null)
            {
                foreach (int j in CutCells.ItemEnum)
                {
                    Ok2Coarsen[j] = false;
                }
            }
            return Ok2Coarsen;
        }

        static int[][] FindCoarseningClusters(BitArray Ok2Coarsen, GridData CurrentGrid)
        {
            int JE = CurrentGrid.Cells.Count;
            int J = CurrentGrid.Cells.NoOfLocalUpdatedCells;

            //if (Ok2Coarsen.Length != JE)
            //    throw new ArgumentException();

            int[][] CellNeighbours = CurrentGrid.Cells.CellNeighbours;

            List<Cell> temp = new List<Cell>();
            List<int> tempCC = new List<int>();

            int RecursionDepht = -1;
            if (CurrentGrid.SpatialDimension == 1)
                RecursionDepht = 1;
            else if (CurrentGrid.SpatialDimension == 2)
                RecursionDepht = 2;
            else if (CurrentGrid.SpatialDimension == 3)
                RecursionDepht = 4;
            else
                throw new NotSupportedException();

            int[][] CoarseningCluster = new int[JE][];

            BitArray marker = new BitArray(JE);
            for (int j = 0; j < J; j++)
            {
                if (marker[j])
                    continue;
                if (!Ok2Coarsen[j])
                    continue;

                Cell Cell_j = CurrentGrid.Cells.GetCell(j);
                int Level = Cell_j.RefinementLevel;

                if (Level == 0)
                {
                    marker[j] = true;
                    continue;
                }

                temp.Clear();
                temp.Add(Cell_j);

                tempCC.Clear();
                tempCC.Add(j);

                int SearchID = Cell_j.CoarseningClusterID;
                int ClusterSize = Cell_j.CoarseningClusterSize;
                Debug.Assert(SearchID > 0);

                bool complete = false;
                FindCoarseiningClusterRecursive(j, RecursionDepht, CellNeighbours, marker, CurrentGrid.Cells.GetCell, Level, SearchID, ClusterSize, tempCC, ref complete);
                foreach (int jC in tempCC)
                {
                    marker[jC] = true;
                }

                if (complete)
                {
                    if (!tempCC.Any(jC => Ok2Coarsen[jC] == false))
                    {
                        int[] CC = tempCC.ToArray();
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
            int[][] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int Level, int SearchID, int ClusterSize,
            List<int> CoarseiningCluster, ref bool complete)
        {
            Debug.Assert(CoarseiningCluster.Contains(j) == true);

            int[] Neighs = CellNeighbours[j];
            foreach (var jNeigh in Neighs)
            {
                if (marker[jNeigh] == true)
                    continue;
                if (CoarseiningCluster.Contains(jNeigh))
                    continue;

                Cell Cell_neigh = GetCell(jNeigh);
                if (Cell_neigh.RefinementLevel != Level)
                    continue;

                if (Cell_neigh.CoarseningClusterID != SearchID)
                    continue;

                CoarseiningCluster.Add(jNeigh);

                if (CoarseiningCluster.Count == ClusterSize)
                {
                    complete = true;
#if DEBUG
                    foreach (int j1 in CoarseiningCluster)
                    {
                        Cell Cell_j1 = GetCell(j1);
                        Debug.Assert(Cell_j1.RefinementLevel > 0);
                        Debug.Assert(Cell_j1.RefinementLevel == Level);
                        Debug.Assert(Cell_j1.CoarseningClusterID == SearchID);
                        Debug.Assert(Cell_j1.CoarseningClusterSize == ClusterSize);
                    }
#endif
                    return;
                }

                if (MaxRecursionDeph > 0)
                    FindCoarseiningClusterRecursive(jNeigh, MaxRecursionDeph - 1,
                        CellNeighbours, marker, GetCell, Level, SearchID, ClusterSize, CoarseiningCluster, ref complete);

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

