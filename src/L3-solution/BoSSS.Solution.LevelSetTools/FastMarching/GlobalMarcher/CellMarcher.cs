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

using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.FastMarcher;
using BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher;
using ilPSP;
using MPI.Wrappers;

namespace BoSSS.Solution.LevelSetTools.FastMarching.GlobalMarcher {

    class CellMarcher {

        GridData gridDat;
        ILocalSolver localSolver;
        
        /// <summary>
        /// Fast marching solver. Initializes a Domain by fast marching. 
        /// Each cell must be initialized locally with a <paramref name="LocalSolver"/>.
        /// </summary>
        /// <param name="LevelSetBasis"></param>
        /// <param name="LocalSolver"> A solver that initializes only one cell</param>
        public CellMarcher(Basis LevelSetBasis, ILocalSolver LocalSolver) {
            gridDat = (GridData)(LevelSetBasis.GridDat);
            localSolver = LocalSolver;
        }

        /// <summary>
        /// Firstorder Reinit of <paramref name="Phi"/> on <paramref name="ReinitField"/> field.
        /// The order in which each cell is initialized is determined by fast marching. 
        /// Locally, each cell is initialized by the localSolver specified in the constructor.
        /// </summary>
        /// <param name="Phi">Field to reinitialize</param>
        /// <param name="Accepted">Start values</param>
        /// <param name="ReinitField">Specific Domain, e.g. whole domain or nearField</param>
        public void Reinit(SinglePhaseField Phi, CellMask Accepted, CellMask ReinitField) {

            //Build Marcher 
            IFastMarchingQueue<IMarchingNode> Heap = new MarchingHeap(this.gridDat.Cells.Count);
            Fastmarcher Solver = new Fastmarcher(Heap);

            //GetNeighbours for external accepted cells
            int[][] externalCellNeighbours = GetNeighboursForExternalCells(gridDat, Accepted);

            //Initialize Graph for Marching and build initial accepted nodes
            MarchingCell.Initialize(localSolver, Phi, gridDat, ReinitField, externalCellNeighbours);
            MarchingCell[] AcceptedCells = MarchingCell.BuildInitialAcceptedCells(Accepted);

            //Solve
            Solver.march(AcceptedCells);

        }


        internal int[][] GetNeighboursForExternalCells(GridData gridDat, CellMask Accepted) {

            int LocalLength = gridDat.Cells.NoOfLocalUpdatedCells;
            int LocalLengthExt = gridDat.Cells.Count; 
            int NoExternal = gridDat.Cells.NoOfExternalCells;

            int[][] neighbours = new int[NoExternal][];

            long[] externalCellsGlobalIndices = gridDat.iParallel.GlobalIndicesExternalCells;
            int[][] globNeigh = GetGlobalCellNeigbourship(gridDat);
            for (int j = LocalLength; j < LocalLengthExt; j++) {
                int jGlob = (int)gridDat.iParallel.GlobalIndicesExternalCells[j - LocalLength];
                int[] globNeighExt = globNeigh[jGlob];
                List<int> neighboursList = new List<int>();
                for (int i = 0; i < globNeighExt.Length; i++) {
                    int neighIndGlob = globNeighExt[i];
                    if (gridDat.CellPartitioning.IsInLocalRange(neighIndGlob))
                        neighboursList.Add(gridDat.CellPartitioning.TransformIndexToLocal(neighIndGlob));
                    else {
                        for (int k = 0; k < globNeighExt.Length; k++) {
                            if (externalCellsGlobalIndices[k] == neighIndGlob)
                                neighboursList.Add(k + LocalLength);
                        }
                    }
                }
                neighbours[j - LocalLength] = neighboursList.ToArray();
            }

            return neighbours;
        }


        private int[][] GetGlobalCellNeigbourship(GridData CurrentGrid) {

            // intermediate solution  

            long[] externalCellsGlobalIndices = CurrentGrid.iParallel.GlobalIndicesExternalCells;
            int GlobalNumberOfCells = CurrentGrid.CellPartitioning.TotalLength;
            int[][] globalCellNeigbourship = new int[GlobalNumberOfCells][];
            int LocalNumberOfCells = CurrentGrid.CellPartitioning.LocalLength;

            for (int j = 0; j < LocalNumberOfCells; j++) {
                int globalIndex = j + CurrentGrid.CellPartitioning.i0;
                // we use GetCellNeighboursViaEdges(j) to also find neigbours at periodic boundaries
                Tuple<int, int, int>[] cellNeighbours = CurrentGrid.GetCellNeighboursViaEdges(j);
                globalCellNeigbourship[globalIndex] = new int[cellNeighbours.Length];
                for (int i = 0; i < cellNeighbours.Length; i++) {
                    int neighIndLoc = cellNeighbours[i].Item1;
                    if (neighIndLoc < LocalNumberOfCells)
                        globalCellNeigbourship[globalIndex][i] = neighIndLoc + CurrentGrid.CellPartitioning.i0;
                    else
                        globalCellNeigbourship[globalIndex][i] = (int)externalCellsGlobalIndices[neighIndLoc - LocalNumberOfCells];
                }

                // translate local neighbour index into global index
                //for (int i = 0; i < globalCellNeigbourship[globalIndex].Length; i++) {
                //    if (globalCellNeigbourship[globalIndex][i] < LocalNumberOfCells)
                //        globalCellNeigbourship[globalIndex][i] = globalCellNeigbourship[globalIndex][i] + CurrentGrid.CellPartitioning.i0;
                //    else
                //        globalCellNeigbourship[globalIndex][i] = (int)externalCellsGlobalIndices[globalCellNeigbourship[globalIndex][i] - LocalNumberOfCells];
                //}
            }

            for (int globalIndex = 0; globalIndex < GlobalNumberOfCells; globalIndex++) {
                int processID = !globalCellNeigbourship[globalIndex].IsNullOrEmpty() ? CurrentGrid.MpiRank : 0;
                processID = processID.MPIMax();
                globalCellNeigbourship[globalIndex] = globalCellNeigbourship[globalIndex].MPIBroadcast(processID);
            }

            return globalCellNeigbourship;
        }
    }
}
