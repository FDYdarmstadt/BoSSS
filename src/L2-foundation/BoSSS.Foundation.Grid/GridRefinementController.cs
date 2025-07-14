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
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
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
        
        
        private readonly CellMask CutCells;
        private readonly BitArray CellsNotOK2Coarsen;
        private readonly bool EnsureHighestLevelAtPeriodicBoundary;

        /// <summary>
        /// Constructor for the grid refinement controller, defines input cells
        /// </summary>
        /// <param name="CurrentGrid">
        /// Current grid.
        /// </param>
        /// <param name="CutCells">
        /// Cut cells will have always the max refinement level. Null is a valid input if no level-set is used.
        /// </param>
        /// <param name="cellsNotOK2Coarsen">
        /// Cells which are not allowed to be coarsened. It is not necessary to include cut cells here, as they are handled by the cutCells CellMask.
        /// </param>
        /// <param name="EnsureHighestLevelAtPeriodicBoundary"></param>
        public GridRefinementController(GridData CurrentGrid, CellMask CutCells, CellMask cellsNotOK2Coarsen = null, bool EnsureHighestLevelAtPeriodicBoundary = false) {
            this.CurrentGrid = CurrentGrid;
            CellPartitioning = this.CurrentGrid.CellPartitioning;
            
            this.CutCells = CutCells;
            if (this.CutCells == null)
                this.CutCells = CellMask.GetEmptyMask(this.CurrentGrid);

            if (cellsNotOK2Coarsen == null)
                cellsNotOK2Coarsen = CellMask.GetEmptyMask(this.CurrentGrid);
            CellsNotOK2Coarsen = cellsNotOK2Coarsen.Union(this.CutCells).GetBitMask();

            this.EnsureHighestLevelAtPeriodicBoundary = EnsureHighestLevelAtPeriodicBoundary;
        }
        
        /// <summary>
        /// Computes refinement and coarsening lists 
        /// (inputs for <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>),
        /// based on the max refinement level provided by the calling solver. This method is fully parallelized.
        /// </summary>
        /// <param name="localCellsToRefine">
        /// Output, local indices of cells which should be refined.
        /// </param>
        /// <param name="localCellsToCoarsen">
        /// Output, clusters of cells (identified by local cell indices) which can be combined into coarser cells,
        /// input for the second argument of 
        /// <see cref="GridData.Adapt(IEnumerable{int}, IEnumerable{int[]}, out GridCorrelation)"/>.
        /// 
        /// </param>
        /// <param name="LocalDesiredLevel">
        /// The refinement level of each local cell.
        /// </param>
        /// <returns>
        /// True if any refinement or coarsening of the current grid should be performed; otherwise false.
        /// </returns>
        public bool ComputeGridChange(int[] LocalDesiredLevel, out List<int> localCellsToRefine, out List<int[]> localCellsToCoarsen) {
            using(new FuncTrace()) {
                // init
                long[][] globalCellNeigbourship = GetGlobalCellNeigbourship();
                BitArray cutCellsWithNeighbours = GetGlobalNearBand(globalCellNeigbourship);

                // refinement
                int[] globalRefinementLevel = GetGlobalRefinementLevel(globalCellNeigbourship, cutCellsWithNeighbours, LocalDesiredLevel, out int changes);
                VerifyPeriodicRefinement(globalRefinementLevel);
                localCellsToRefine = GetCellsToRefine(globalRefinementLevel);

                // coarsening of those cells where refinement is not needed (anymore)
                localCellsToCoarsen = GetCellsToCoarsen(globalRefinementLevel, globalCellNeigbourship, cutCellsWithNeighbours);
                localCellsToCoarsen = FixPeriodicCoarsening(localCellsToCoarsen);
                VerifyPeriodicCoarsening(localCellsToCoarsen);

                // return 
                bool anyChangeInGrid = (localCellsToRefine.Count() == 0 && localCellsToCoarsen.Count() == 0) ? false : true;
                return anyChangeInGrid.MPIOr();
            }

        }

        /// <summary>
        /// Removes a coarsening cluster if it occurs only on one side of a periodic edge
        /// </summary>
        List<int[]> FixPeriodicCoarsening(List<int[]> localCellsToCoarsen) {
            using(var tr = new FuncTrace()) {
                if(this.CurrentGrid.Grid.PeriodicTrafo.Count <= 0)
                    return localCellsToCoarsen; // no need to check periodicity


                int Nupdate = CellPartitioning.LocalLength;
                long myI0 = CellPartitioning.i0;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;


                long[] allCellsToCoarsen;
                Partitioning allCellsToCoarsenPartition;
                int[] localCell2CoarseningCluster = new int[Nupdate];
                {
                    localCell2CoarseningCluster.SetAll(-10000001);

                    List<long> _allCellsToCoarsen = new List<long>();
                    int icc = 0;
                    foreach(var CoarsenCluster in localCellsToCoarsen) {
                        foreach(int j in CoarsenCluster) {
                            if(j < 0 || j >= Nupdate)
                                throw new ApplicationException("unable to coarsen external cell");

                            _allCellsToCoarsen.Add(j + myI0);
                            localCell2CoarseningCluster[j] = icc;
                        }
                        icc++;
                    }

                    allCellsToCoarsen = _allCellsToCoarsen.ToArray().MPIAllGatherv();
                    allCellsToCoarsenPartition = new Partitioning(_allCellsToCoarsen.Count);
                }

                BitArray allCellsToCoarsen_mask = new BitArray(checked((int)GlobalNumberOfCells));
                foreach(long jGlob in allCellsToCoarsen) {
                    allCellsToCoarsen_mask[checked((int)jGlob)] = true;
                }


                int[,] E2C = this.CurrentGrid.Edges.CellIndices;
                byte[] tags = this.CurrentGrid.Edges.EdgeTags;
                long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;


                int NoOfEdges = tags.Length;
                bool anyRemoved = false; ;
                for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) { // loop over local edges...
                    byte tag = tags[iEdge];
                    if(tag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                        // found some periodic edge
                        int jCell0_loc = E2C[iEdge, 0];
                        int jCell1_loc = E2C[iEdge, 1];

                        if(jCell0_loc >= Nupdate)
                            throw new ApplicationException("at an periodic edge, the 0-th cell is supposed to be local");
                        long jCell0_glb = jCell0_loc + myI0;
                        long jCell1_glb = jCell1_loc < Nupdate ? jCell1_loc + myI0 : externalCellsGlobalIndices[jCell1_loc - Nupdate];


                        if(allCellsToCoarsen_mask[checked((int)jCell0_glb)] != allCellsToCoarsen_mask[checked((int)jCell1_glb)]) {
                            // the respective coarsening cluster has to be removed...
                            if(allCellsToCoarsen_mask[checked((int)jCell0_glb)] && jCell0_loc < Nupdate) {
                                localCellsToCoarsen[localCell2CoarseningCluster[jCell0_loc]] = null;
                            }
                            if(allCellsToCoarsen_mask[checked((int)jCell1_glb)] && jCell1_loc < Nupdate) {
                                localCellsToCoarsen[localCell2CoarseningCluster[jCell1_loc]] = null;
                            }
                            anyRemoved = true;
                        }
                    }
                }
                tr.InfoToConsole = true;
                
                // return
                // ======

                tr.Info("Found hanging nodes on periodic edges? " + anyRemoved);
                if(anyRemoved) {
                    var ret = new List<int[]>();
                    foreach(var rii in localCellsToCoarsen) {
                        if(rii != null)
                            ret.Add(rii);
                    }
                    return ret;
                } else {
                    return localCellsToCoarsen;
                }



            }
        }


        /// <summary>
        /// verifies that the <paramref name="localCellsToCoarsen"/> is equal on both sides of periodic edges.
        /// </summary>
        void VerifyPeriodicCoarsening(List<int[]> localCellsToCoarsen) {
            using(var tr = new FuncTrace()) {
                if(this.CurrentGrid.Grid.PeriodicTrafo.Count <= 0)
                    return; // no need to check periodicity


                int Nupdate = CellPartitioning.LocalLength;
                long myI0 = CellPartitioning.i0;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;


                long[] allCellsToCoarsen;
                {
                    List<long> _allCellsToCoarsen = new List<long>();
                    foreach(var CoarsenCluster in localCellsToCoarsen) {
                        foreach(int j in CoarsenCluster) {
                            if(j < 0 || j >= Nupdate)
                                throw new ApplicationException("unable to coarsen external cell");

                            _allCellsToCoarsen.Add(j + myI0);

                        }
                    }

                    allCellsToCoarsen = _allCellsToCoarsen.ToArray().MPIAllGatherv();
                }

                BitArray allCellsToCoarsen_mask = new BitArray(checked((int)GlobalNumberOfCells));
                foreach(long jGlob in allCellsToCoarsen) {
                    allCellsToCoarsen_mask[checked((int)jGlob)] = true;
                }


                int[,] E2C = this.CurrentGrid.Edges.CellIndices;
                byte[] tags = this.CurrentGrid.Edges.EdgeTags;
                long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;


                int NoOfEdges = tags.Length;
                for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                    byte tag = tags[iEdge];
                    if(tag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                        int jCell0_loc = E2C[iEdge, 0];
                        int jCell1_loc = E2C[iEdge, 1];

                        if(jCell0_loc >= Nupdate)
                            throw new ApplicationException("at an periodic edge, the 0-th cell is supposed to be local");
                        long jCell0_glb = jCell0_loc + myI0;
                        long jCell1_glb = jCell1_loc < Nupdate ? jCell1_loc + myI0 : externalCellsGlobalIndices[jCell1_loc - Nupdate];


                        if(allCellsToCoarsen_mask[checked((int)jCell0_glb)] != allCellsToCoarsen_mask[checked((int)jCell1_glb)]) {
                            throw new ApplicationException("internal error -- illegal coarsening detected across periodic edge.");
                        }



                    }
                }
                tr.InfoToConsole = true;
                tr.Info("checked coarsening across periodic edges");
            }
        }


        /// <summary>
        /// verifies that the <paramref name="globalRefinementLevel"/> is equal on both sides of periodic edges.
        /// </summary>
        void VerifyPeriodicRefinement(int[] globalRefinementLevel) {
            using(var tr = new FuncTrace()) {
                if(this.CurrentGrid.Grid.PeriodicTrafo.Count <= 0)
                    return; // no need to check periodicity


                int[,] E2C = this.CurrentGrid.Edges.CellIndices;
                byte[] tags = this.CurrentGrid.Edges.EdgeTags;
                long myI0 = CellPartitioning.i0;
                int Nupdate = CellPartitioning.LocalLength;
                long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;


                int NoOfEdges = tags.Length;
                for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                    byte tag = tags[iEdge];
                    if(tag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                        int jCell0_loc = E2C[iEdge, 0];
                        int jCell1_loc = E2C[iEdge, 1];

                        if(jCell0_loc >= Nupdate)
                            throw new ApplicationException("at an periodic edge, the 0-th cell is supposed to be local");
                        long jCell0_glb = jCell0_loc + myI0;
                        long jCell1_glb = jCell1_loc < Nupdate ? jCell1_loc + myI0 : externalCellsGlobalIndices[jCell1_loc - Nupdate];

                        if(globalRefinementLevel[jCell0_glb] != globalRefinementLevel[jCell1_glb]) {
                            throw new ApplicationException("internal error -- different global levels requested for periodic neighbors.");
                        }
                    }
                }

                tr.Info("checked periodic neighborship");
            }
        }





        /// <summary>
        /// Computes the global cell neighborship of all cells. 
        /// Returns an jagged array where the 
        /// - first index: global index of the current cell 
        /// - second index: enumeration over neighbors.
        /// - content: global index of neighbor cells.
        /// </summary>
        private long[][] GetGlobalCellNeigbourship() {
            using(var tr = new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
                long myI0 = CellPartitioning.i0;


                long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;
                long[][] globalCellNeigbourship_local = new long[LocalNumberOfCells][];


                for(int j = 0; j < LocalNumberOfCells; j++) {
                    //long globalIndex = j + myI0;
                    // we use GetCellNeighboursViaEdges(j) to also find neigbours at periodic boundaries
                    var cellNeighbours = CurrentGrid.GetCellNeighboursViaEdges(j, OmmitPeriodic:false);
                    int NN = cellNeighbours.Length;
                    globalCellNeigbourship_local[j] = new long[NN];

                    for(int i = 0; i < NN; i++) {
                        globalCellNeigbourship_local[j][i] = cellNeighbours[i].Item1;
                    }

                    // translate local neighbour index into global index
                    for(int i = 0; i < NN; i++) {
                        long n_ji = globalCellNeigbourship_local[j][i];

                        if(n_ji < LocalNumberOfCells)
                            globalCellNeigbourship_local[j][i] = n_ji + myI0;
                        else
                            globalCellNeigbourship_local[j][i] = externalCellsGlobalIndices[n_ji - LocalNumberOfCells];
                    }

                }

                return globalCellNeigbourship_local.MPI_AllGaterv();
            }
        }


        /// <summary>
        /// Returns the cutcells + neighbours on a global level.
        /// </summary>
        /// <param name="globalCellNeighbourship">
        /// </param>
        private BitArray GetGlobalNearBand(long[][] globalCellNeighbourship) {
            using(new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
                long myI0 = CellPartitioning.i0;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;

                if(CutCells == null)
                    return new BitArray(checked((int)GlobalNumberOfCells));

                BitArray localCutCells = CutCells.GetBitMask();
                BitArray globalCutCells = new BitArray(checked((int)GlobalNumberOfCells));
                for(int j = 0; j < LocalNumberOfCells; j++) {
                    if(localCutCells[j]) {
                        int globalIndex = checked((int)(j + myI0));
                        globalCutCells[globalIndex] = true;
                        for(int i = 0; i < globalCellNeighbourship[globalIndex].Length; i++) {
                            globalCutCells[(int)globalCellNeighbourship[globalIndex][i]] = true;
                        }
                    }
                }

                globalCutCells.MPIOr();

                return globalCutCells;
            }
        }
        

        /// <summary>
        /// Computes the level indicator for each cell (MPI global). 
        /// </summary>
        private int[] GetGlobalRefinementLevel(long[][] globalNeighbourship, BitArray cutCellsWithNeighbours, int[] LocalDesiredLevel, out int changes) {
            using(new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
                long myI0 = CellPartitioning.i0;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;


                changes = 0;
                int[] globalDesiredLevel = GetGlobalDesiredLevel(cutCellsWithNeighbours, LocalDesiredLevel, ref changes);
                int[] globalCurrentLevel = GetGlobalCurrentLevel();
                int[] globalRefinementLevel = new int[GlobalNumberOfCells];
                Debug.Assert(globalRefinementLevel.Length == globalCurrentLevel.Length);

                for(long j = 0; j < LocalNumberOfCells; j++) {
                    long globalCellIndex = j + myI0;
                    if(globalRefinementLevel[globalCellIndex] < globalDesiredLevel[globalCellIndex]) {
                        globalRefinementLevel[globalCellIndex] = globalDesiredLevel[globalCellIndex];
                        GetRefinementLevelRecursive(globalCellIndex, globalRefinementLevel[globalCellIndex] - 1, globalNeighbourship, globalRefinementLevel);
                    }
                }
                return globalRefinementLevel.MPIMax();
            }
        }


        /// <summary>
        /// Writes the max desired level of the specified cells into an int-array (mpi global).
        /// </summary>
        private int[] GetGlobalDesiredLevel(BitArray cutCellsWithNeighbours, int[] CellRefinementLevel, ref int changes) {
            int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
            long myI0 = CellPartitioning.i0;
            long GlobalNumberOfCells = CellPartitioning.TotalLength;


            using(new FuncTrace()) {
                int[] globalDesiredLevel = new int[GlobalNumberOfCells];
                int levelSetMaxLevel = CellRefinementLevel.Max().MPIMax();

                for(int j = 0; j < LocalNumberOfCells; j++) {
                    int globalIndex = checked((int)(j + myI0));
                    if(CellRefinementLevel[j] > 0 && globalDesiredLevel[globalIndex] <= CellRefinementLevel[j]) {
                        globalDesiredLevel[globalIndex] = CellRefinementLevel[j];
                    }

                    if(cutCellsWithNeighbours[globalIndex]) {
                        globalDesiredLevel[globalIndex] = levelSetMaxLevel;
                    }

                    if(EnsureHighestLevelAtPeriodicBoundary) {
                        CellFaceTag[] cellFaceTags = CurrentGrid.Cells.GetCell(j).CellFaceTags;
                        if(cellFaceTags != null) {
                            for(int c = 0; c < cellFaceTags.Length; c++) {
                                if(cellFaceTags[c].EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                                    globalDesiredLevel[globalIndex] = levelSetMaxLevel;
                                    break;
                                }
                            }
                        }
                    }
                }

                return globalDesiredLevel.MPIMax();
            }
        }


        /// <summary>
        /// Gets the current refinement level for each global cell.
        /// </summary>
        private int[] GetGlobalCurrentLevel() {
            using(new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;


                int[] localCurrentLevel = new int[LocalNumberOfCells];
                long[] i0 = CellPartitioning.GetI0s();

                for(int j = 0; j < LocalNumberOfCells; j++) {
                    localCurrentLevel[j] = CurrentGrid.Cells.GetCell(j).RefinementLevel;
                }
                
                int[] globalCurrentLevel = localCurrentLevel.MPIAllGatherv();
                return globalCurrentLevel;
            }
        }

        /// <summary>
        /// Recursive computation for the desired level of each global cell,
        /// i.e., ensures that if some cell has refinement level n, that all neighbors have at least n-1.
        /// </summary>
        /// <param name="globalCellIndex">
        /// The global index of the current cell.
        /// </param>
        /// <param name="desiredLevel"></param>
        /// <param name="globalNeighbourship">
        /// Jagged int-array where the first index refers to the current cell and the second one to the neighbor cells.
        /// </param>
        /// <param name="globalRefinementLevel"></param>
        private void GetRefinementLevelRecursive(long globalCellIndex, int desiredLevel, long[][] globalNeighbourship, int[] globalRefinementLevel) {
            if (desiredLevel <= 0)
                return;
            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                long globalCellIndexNeighbour = globalNeighbourship[globalCellIndex][j];
                if (globalRefinementLevel[globalCellIndexNeighbour] < desiredLevel){
                    globalRefinementLevel[globalCellIndexNeighbour] = desiredLevel;
                    GetRefinementLevelRecursive(globalCellIndexNeighbour, desiredLevel - 1, globalNeighbourship, globalRefinementLevel);
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
        private void RefineNeighboursRecursive(int currentLevel, int[] DesiredLevel, long globalCellIndex, int DesiredLevelNeigh, long[][] globalNeighbourship, int[] globalCurrentLevel) {
            if (DesiredLevelNeigh <= 0)
                return;
            for (int j = 0; j < globalNeighbourship[globalCellIndex].Length; j++) {
                long jNeigh = globalNeighbourship[globalCellIndex][j];
                if (currentLevel <= DesiredLevelNeigh && DesiredLevel[jNeigh] <= DesiredLevelNeigh) {
                    DesiredLevel[jNeigh] = DesiredLevelNeigh;
                    RefineNeighboursRecursive(currentLevel, DesiredLevel, jNeigh, DesiredLevelNeigh - 1, globalNeighbourship, globalCurrentLevel);
                }
            }
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a local int-list.
        /// </summary>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        private List<int> GetCellsToRefine(int[] globalDesiredLevel) {
            long i0 = CurrentGrid.CellPartitioning.i0;
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
        /// Returns local cells to coarsen
        /// </summary>
        private List<int[]> GetCellsToCoarsen(int[] globalRefinementLevel, long[][] globalCellNeigbourship, BitArray cutCellsWithNeighbours) {
            using(new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;

                BitArray oK2Coarsen = GetCellsOk2Coarsen(globalRefinementLevel, globalCellNeigbourship, cutCellsWithNeighbours);
                int[][] coarseningClusters = FindCoarseningClusters(oK2Coarsen);
                List<int[]> coarseningCells = new List<int[]>();
                for(int j = 0; j < LocalNumberOfCells; j++) {
                    if(coarseningClusters[j] != null) {
                        Debug.Assert(coarseningClusters[j].Contains(j));
                        if(j == coarseningClusters[j].Min()) {
                            coarseningCells.Add(coarseningClusters[j]);
                        }
                    }
                }
                return coarseningCells;
            }
        }

        /// <summary>
        /// Gets all cells to refine and writes them to a int-list.
        /// </summary>
        /// <param name="globalDesiredLevel">
        /// The desired level of all global cells.
        /// </param>
        /// <param name="globalCellNeigbourship">
        /// Jagged int-array where the first index refers to the current cell and the second one to the neighbor cells.
        /// </param>
        /// <param name="cutCellsWithNeighbours"></param>
        private BitArray GetCellsOk2Coarsen(int[] globalDesiredLevel, long[][] globalCellNeigbourship, BitArray cutCellsWithNeighbours) {
            int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
            long myI0 = CellPartitioning.i0;
            long GlobalNumberOfCells = CellPartitioning.TotalLength;


            BitArray oK2Coarsen = new BitArray(checked((int)GlobalNumberOfCells));
            int[] globalCurrentLevel = GetGlobalCurrentLevel();

            for (long globalCellIndex = myI0; globalCellIndex < myI0 + LocalNumberOfCells; globalCellIndex++) {
                int currentLevel = globalCurrentLevel[globalCellIndex];

                int maxNeighbourDesiredLevel = 0;
                int maxNeighbourCurrentLevel = 0;
                for(int i = 0; i < globalCellNeigbourship[globalCellIndex].Length; i++) {
                    if (globalDesiredLevel[globalCellNeigbourship[globalCellIndex][i]] > maxNeighbourDesiredLevel)
                        maxNeighbourDesiredLevel = globalDesiredLevel[globalCellNeigbourship[globalCellIndex][i]];
                    int neighbourCurrentLevel = globalCurrentLevel[globalCellNeigbourship[globalCellIndex][i]];
                    if (neighbourCurrentLevel > maxNeighbourCurrentLevel)
                        maxNeighbourCurrentLevel = neighbourCurrentLevel;
                }

                if (currentLevel > globalDesiredLevel[globalCellIndex] 
                    && currentLevel > maxNeighbourDesiredLevel 
                    && currentLevel > maxNeighbourCurrentLevel - 1 
                    && !cutCellsWithNeighbours[checked((int)globalCellIndex)]) {
                    oK2Coarsen[checked((int)globalCellIndex)] = true;
                }
            }

            BitArray globalCellsNotOK2Coarsen = GetGlobalCellsNotOK2Coarsen();

            for (int j = 0; j < globalCellsNotOK2Coarsen.Length; j++) {
                if (globalCellsNotOK2Coarsen[j]) {
                    oK2Coarsen[j] = false;
                }
            }

            var oK2Coarsen_b4 = oK2Coarsen.CloneAs();
            oK2Coarsen.MPIAnd();

            bool changed = false;
            for(int j = 0; j < oK2Coarsen.Length; j++) {
                if(oK2Coarsen[j] != oK2Coarsen_b4[j])
                    changed = true;
            }
            changed = changed.MPIOr();
            Console.WriteLine("-----------------------------------------  coarse changed " + changed);


            return oK2Coarsen;
        }


        /// <summary>
        /// Globalize cellsNotOk2Coarsen
        /// </summary>
        private BitArray GetGlobalCellsNotOK2Coarsen() {
            int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
            long myI0 = CellPartitioning.i0;
            long GlobalNumberOfCells = CellPartitioning.TotalLength;


            BitArray globalCellsNotOK2Coarsen = new BitArray(checked((int)GlobalNumberOfCells));
            for (int j = 0; j < LocalNumberOfCells; j++) {
                if (CellsNotOK2Coarsen[j]) {
                    int globalIndex = checked((int)(j + myI0));
                    globalCellsNotOK2Coarsen[globalIndex] = true;
                }
            }

            //for (int j = 0; j < globalCellsNotOK2Coarsen.Length; j++) {
            //    globalCellsNotOK2Coarsen[j] = globalCellsNotOK2Coarsen[j].MPIOr();
            //}
            globalCellsNotOK2Coarsen.MPIOr();

            return globalCellsNotOK2Coarsen;
        }

        /// <summary>
        /// Gets all cells to be coarsened. This operates process local, no coarsening allowed over process boundaries.
        /// </summary>
        /// <param name="oK2CoarsenGlobal">
        /// A BitArray of all cells which should be coarsened.
        /// </param>
        private int[][] FindCoarseningClusters(BitArray oK2CoarsenGlobal) {
            using(new FuncTrace()) {
                int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;
                long myI0 = CellPartitioning.i0;
                long GlobalNumberOfCells = CellPartitioning.TotalLength;

                int JE = CurrentGrid.Cells.Count;
                BitArray oK2Coarsen = new BitArray(LocalNumberOfCells);
                for(long j = myI0; j < myI0 + LocalNumberOfCells; j++) {
                    oK2Coarsen[checked((int)(j - myI0))] = oK2CoarsenGlobal[checked((int)j)];
                }
                //int[][] cellNeighbours = currentGrid.Cells.CellNeighbours;
                List<Cell> temp = new List<Cell>();
                List<int> tempCoarseningCluster = new List<int>();

                int RecursionDepth = -1;
                if(CurrentGrid.SpatialDimension == 1)
                    RecursionDepth = 1;
                else if(CurrentGrid.SpatialDimension == 2)
                    RecursionDepth = 2;
                else if(CurrentGrid.SpatialDimension == 3)
                    RecursionDepth = 4;
                else
                    throw new NotSupportedException();

                int[][] CoarseningCluster = new int[JE][];

                BitArray marker = new BitArray(JE); // true for alls cells which should **not** be coarsened
                for(int j = 0; j < LocalNumberOfCells; j++) {
                    if(marker[j]) 
                        continue; // coarsening not allowed
                    if(!oK2Coarsen[j]) 
                        continue; // coarsening not allowed

                    Cell currentCell = CurrentGrid.Cells.GetCell(j);
                    int currentRefinmentLevel = currentCell.RefinementLevel;

                    if(currentRefinmentLevel == 0) {
                        marker[j] = true; 
                        continue; // coarsening not possible - already on coarsest level
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
                    foreach(int jC in tempCoarseningCluster) {
                        marker[jC] = true;
                    }

                    if(complete) {
                        if(!tempCoarseningCluster.Any(jC => oK2Coarsen[jC] == false)) {
                            int[] CC = tempCoarseningCluster.ToArray();
                            foreach(int jC in CC) {
                                Debug.Assert(CoarseningCluster[jC] == null);
                                CoarseningCluster[jC] = CC;
                            }
                        }
                    }
                }
                return CoarseningCluster;
            }
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
        /// Returns true for all cells **not** ok to coarsen.
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
            int[] CellNeighbours, 
            BitArray marker, Func<int, Cell> GetCell, int currentRefinementLevel, int searchClusterID, int clusterSize,
            List<int> coarseningCluster, ref bool complete) {
            if (!coarseningCluster.Contains(j))
                throw new Exception("Error in coarsening algorithm: Coarsening cluster does not contain a cell with the local ID: " + j);
            int LocalNumberOfCells = this.CurrentGrid.Cells.NoOfLocalUpdatedCells;


            foreach(int neighbourCellIndex in CellNeighbours) {

                Cell neighbourCell = GetCell(neighbourCellIndex);

                bool ContinueCrit = (marker[neighbourCellIndex] == true)
                    || coarseningCluster.Contains(neighbourCellIndex)
                    || neighbourCell.RefinementLevel != currentRefinementLevel
                    || neighbourCell.CoarseningClusterID != searchClusterID
                    || neighbourCellIndex >= LocalNumberOfCells;

                if (ContinueCrit) continue;

                //if (marker[neighbourCellIndex] == true)
                //    continue;

                    //if (coarseningCluster.Contains(neighbourCellIndex))
                    //    continue;

                    //if (neighbourCell.RefinementLevel != currentRefinementLevel)
                    //    continue;

                    //if (neighbourCell.CoarseningClusterID != searchClusterID)
                    //    continue;

                    //if (neighbourCellIndex > LocalNumberOfCells)
                    //    continue;

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

       


    }
}

