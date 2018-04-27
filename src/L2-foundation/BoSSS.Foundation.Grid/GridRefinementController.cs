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
        /// <param name="CurrentGrid">
        /// Current grid.
        /// </param>
        /// <param name="LevelIndicator">
        /// Mapping from (local cell index, current refinement level) to desired refinement level for the respective cell,
        /// see <see cref="Cell.RefinementLevel"/>.
        /// </param>
        /// <param name="CellsToRefineList">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="Coarsening">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells.
        /// </param>
        /// <param name="CutCells">
        /// If not null, a mask of cells in which coarsening is forbidden (usually cut-cells);
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public static bool ComputeGridChange(GridData CurrentGrid, CellMask CutCells, Func<int, int, int> LevelIndicator, out List<int> CellsToRefineList, out List<int[]> Coarsening) {
            int oldJ = CurrentGrid.Cells.NoOfLocalUpdatedCells;

            bool NoRefinement = true;
            int[] DesiredLevel = new int[oldJ];
            for (int j = 0; j < oldJ; j++) {

                int CurrentLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = LevelIndicator(j, CurrentLevel_j);

                if (DesiredLevel[j] < DesiredLevel_j) {
                    DesiredLevel[j] = DesiredLevel_j;
                    NoRefinement = false;
                    RefineNeighboursRecursive(CurrentGrid, DesiredLevel, j, DesiredLevel_j - 1);
                }
            }

            BitArray Ok2Coarsen = new BitArray(oldJ);
            for (int j = 0; j < oldJ; j++) {
                int ActualLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                int DesiredLevel_j = DesiredLevel[j];

                int[][] CellNeighbours = CurrentGrid.Cells.CellNeighbours;

                if (ActualLevel_j > DesiredLevel_j && DesiredLevel_j >= CellNeighbours[j].Select(cn => DesiredLevel[cn]).Max() - 1 ) {
                    Ok2Coarsen[j] = true;
                }
            }
            if (CutCells != null) {
                foreach (int j in CutCells.ItemEnum) {
                    Ok2Coarsen[j] = false;
                }
            }


            int[][] CClusters = FindCoarseningClusters(Ok2Coarsen, CurrentGrid);

            Coarsening = new List<int[]>();
            int NoOfCellsToCoarsen = 0;
            for (int j = 0; j < oldJ; j++) {

                if (CClusters[j] != null) {
                    NoOfCellsToCoarsen++;

                    Debug.Assert(CClusters[j].Contains(j));
                    if (j == CClusters[j].Min()) {
                        Coarsening.Add(CClusters[j]);
                    }
                }
            }

            CellsToRefineList = new List<int>();
            if ((!NoRefinement) || (Coarsening.Count > 0)) {

                for (int j = 0; j < oldJ; j++) {
                    int ActualLevel_j = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                    int DesiredLevel_j = DesiredLevel[j];

                    if (ActualLevel_j < DesiredLevel_j) 
                        CellsToRefineList.Add(j);
                }
            }

            // If any cells which should refined are members of CutCells
            if (CellsToRefineList.Count == 0 && Coarsening.Count == 0)
                NoRefinement = true;


            return (!NoRefinement);
        }



        static int[][] FindCoarseningClusters(BitArray Ok2Coarsen, GridData CurrentGrid) {
            int JE = CurrentGrid.Cells.NoOfCells;
            int J = CurrentGrid.Cells.NoOfLocalUpdatedCells;

            if (Ok2Coarsen.Length != JE)
                throw new ArgumentException();

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
            for (int j = 0; j < J; j++) {
                if (marker[j])
                    continue;
                if (!Ok2Coarsen[j])
                    continue;

                Cell Cell_j = CurrentGrid.Cells.GetCell(j);
                int Level = Cell_j.RefinementLevel;

                if (Level == 0) {
                    marker[j] = true;
                    continue;
                }

                temp.Clear();
                temp.Add(Cell_j);

                tempCC.Clear();
                tempCC.Add(j);

                //SearchGids = new long[Cell_j.CoarseningPeers.Length + 1];
                //Array.Copy(Cell_j.CoarseningPeers, SearchGids, SearchGids.Length - 1);
                //SearchGids[SearchGids.Length - 1] = Cell_j.GlobalID;
                int SearchID = Cell_j.CoarseningClusterID;
                int ClusterSize = Cell_j.CoarseningClusterSize;
                Debug.Assert(SearchID > 0);


                bool complete = false;
                FindCoarseiningClusterRecursive(j, RecursionDepht, CellNeighbours, marker, CurrentGrid.Cells.GetCell, Level, SearchID, ClusterSize, tempCC, ref complete);
                foreach (int jC in tempCC) {
                    marker[jC] = true;
                }

                if (complete) {
                    if (!tempCC.Any(jC => Ok2Coarsen[jC] == false)) {
                        int[] CC = tempCC.ToArray();
                        foreach (int jC in CC) {
                            Debug.Assert(CoarseningCluster[jC] == null);
                            CoarseningCluster[jC] = CC;
                        }
                    } else {
                        //Console.WriteLine("Not ok to coarsen.");
                    }
                }
            }


            return CoarseningCluster;
        }

        static void FindCoarseiningClusterRecursive(int j, int MaxRecursionDeph,
            int[][] CellNeighbours, BitArray marker, Func<int, Cell> GetCell, int Level, int SearchID, int ClusterSize,
            List<int> CoarseiningCluster, ref bool complete) {

            Debug.Assert(CoarseiningCluster.Contains(j) == true);

            int[] Neighs = CellNeighbours[j];
            foreach (var jNeigh in Neighs) {
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

                if (CoarseiningCluster.Count == ClusterSize) {
                    complete = true;
#if DEBUG
                    foreach (int j1 in CoarseiningCluster) {
                        Cell Cell_j1 = GetCell(j1);
                        Debug.Assert(Cell_j1.RefinementLevel > 0);
                        //Debug.Assert(Cell_j1.CoarseningPeers != null);
                        //Debug.Assert(Cell_j1.CoarseningPeers.Length == CoarseiningCluster.Count - 1);
                        Debug.Assert(Cell_j1.RefinementLevel == Level);
                        Debug.Assert(Cell_j1.CoarseningClusterID == SearchID);
                        Debug.Assert(Cell_j1.CoarseningClusterSize == ClusterSize);
                        //Debug.Assert(Array.IndexOf(SearchGids, Cell_j1.GlobalID) >= 0);
                        //foreach(long gid_cp in Cell_j1.CoarseningPeers) {
                        //    Debug.Assert(Array.IndexOf(SearchGids, gid_cp) >= 0);
                        //}
                    }
#endif
                    return;
                }

                if (MaxRecursionDeph > 0)
                    FindCoarseiningClusterRecursive(jNeigh, MaxRecursionDeph - 1,
                        CellNeighbours, marker, GetCell, Level, SearchID, ClusterSize, CoarseiningCluster, ref complete);

                if (complete) {
                    return;
                }
            }

        }

        static void RefineNeighboursRecursive(GridData gdat, int[] DesiredLevel, int j, int DesiredLevelNeigh) {
            if (DesiredLevelNeigh <= 0)
                return;

            foreach (var jNeigh in gdat.Cells.CellNeighbours[j]) {
                var cl = gdat.Cells.GetCell(jNeigh);
                if (cl.RefinementLevel < DesiredLevelNeigh && DesiredLevel[jNeigh] < DesiredLevelNeigh) {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(gdat, DesiredLevel, jNeigh, DesiredLevelNeigh - 1);
                }
            }
        }

    }
}

