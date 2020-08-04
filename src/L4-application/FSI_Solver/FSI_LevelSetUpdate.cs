/* =======================================================================
Copyright 2019 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using BoSSS.Application.FSI_Solver;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace FSI_Solver {
    internal class FSI_LevelSetUpdate {

        internal FSI_LevelSetUpdate(IGridData GridData, double MinGridLength) {
            this.GridData = GridData;
            this.MinGridLength = MinGridLength;
        }

        private readonly IGridData GridData;
        private readonly double MinGridLength;

        /// <summary>
        /// Initialization of coloured cells based on particle geometry.
        /// </summary>
        internal int[] InitializeColoring(LevelSetTracker LsTrk, Particle[] Particles, double MaxGridLength) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            int[] coloredCells = new int[J];
            for (int p = 0; p < Particles.Length; p++) {
                for (int j = 0; j < J; j++) {
                    if (Particles[p].Contains(new Vector(CellCenters[j, 0], CellCenters[j, 1]), MaxGridLength))
                        coloredCells[j] = p + 1;
                }
            }
            RecolorCellsOfNeighborParticles(coloredCells, (GridData)GridData);
            return coloredCells;
        }

        /// <summary>
        /// Update of all coloured cells.
        /// </summary>
        internal int[] UpdateColoring(LevelSetTracker LsTrk) {
            int[] coloredCells = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int[] coloredCellsExchange = coloredCells.CloneAs();
            coloredCellsExchange.MPIExchange(GridData);
            ColorNeighborCells(coloredCells, coloredCellsExchange);
            RecolorCellsOfNeighborParticles(coloredCells, (GridData)GridData);
            return coloredCells;
        }

        /// <summary>
        /// searchs for all particles with the same color
        /// </summary>
        /// <param name="ParticleColor">
        /// a list of all particles with their specific color
        /// </param>
        /// <param name="CurrentColor">
        /// the color which is associated with the Cell Mask
        /// </param
        internal int[] FindParticlesWithSameColor(int[] ParticleColor, int CurrentColor) {
            List<int> particleIDsWithSameColor = new List<int>();
            for (int i = 0; i < ParticleColor.Length; i++) {
                if (ParticleColor[i] == CurrentColor)
                    particleIDsWithSameColor.Add(i);
            }
            return particleIDsWithSameColor.ToArray();
        }
    
        /// <summary>
        /// Method find the color of all particles combined with MPI communication.
        /// </summary>
        /// <param name="CellColor">
        /// all cells with their color
        /// </param>
        /// <param name="GridData">
        /// IGridData
        /// </param>
        /// <param name="Particles">
        /// A list of all particles.
        /// </param>
        internal int[] DetermineGlobalParticleColor(IGridData GridData, int[] CellColor, List<Particle> Particles) {
            List<int[]> ColoredCellsSorted = ColoredCellsFindAndSort(CellColor);
            int[] ParticleColorArray = FindParticleColor(GridData, Particles, ColoredCellsSorted);
            int[] GlobalParticleColor = ParticleColorArray.MPIMax();
            return GlobalParticleColor;
        }

        /// <summary>
        /// Method to color new ghost particles with the same color as their masters.
        /// </summary>
        /// <param name="globalParticleColor">
        /// Colors of all particles
        /// </param>
        /// <param name="ParticleList">
        /// A list of all particles
        /// </param>
        internal void DetermineParticleColorOfGhostParticles(int[] globalParticleColor, List<Particle> ParticleList) {
            for (int i = 0; i < globalParticleColor.Length; i++) {
                if (globalParticleColor[i] == 0) {
                    int masterID = ParticleList[i].MasterGhostIDs[0] - 1;
                    globalParticleColor[i] = globalParticleColor[masterID];
                }
            }
        }

        /// <summary>
        /// Method to check whether the current color is used on this process.
        /// </summary>
        /// <param name="CellColor">
        /// all cells with their color
        /// </param>
        /// <param name="CurrentColor">
        /// The current color to be checked.
        /// </param>
        /// <param name="NoOfLocalCells">
        /// Number of locally stored cells.
        /// </param>
        internal bool MPIProcessContainsCurrentColor(int[] CellColor, int CurrentColor, int NoOfLocalCells) {
            for (int j = 0; j < NoOfLocalCells; j++) {
                if (CellColor[j] == CurrentColor && CurrentColor != 0)
                    return true;
            }
            return false;
        }

        /// <summary>
        /// Transforms the int[] CellColor into an BitArray.
        /// </summary>
        /// <param name="CellColor">
        /// all cells with their color
        /// </param>
        /// <param name="CurrentColor">
        /// The current color to be checked.
        /// </param>
        /// <param name="NoOfLocalCells">
        /// Number of locally stored cells.
        /// </param>
        internal BitArray CreateBitArrayFromColoredCells(int[] CellColor, int CurrentColor, int NoOfLocalCells) {
            BitArray coloredCells = new BitArray(NoOfLocalCells);
            for (int j = 0; j < NoOfLocalCells; j++) {
                if (CellColor[j] == CurrentColor && CurrentColor != 0) {
                    coloredCells[j] = true;
                }
            }
            return coloredCells;
        }

        /// <summary>
        /// method to find the color of all particle
        /// </summary>
        /// <param name="gridData">
        /// the grid that this mask will be associated with
        /// </param>
        /// <param name="ColoredCellsSorted">
        /// a list of all cells sorted by color
        /// </param>
        /// <param name="Particles">
        /// a list of all particles
        /// </param>
        private int[] FindParticleColor(IGridData gridData, List<Particle> Particles, List<int[]> ColoredCellsSorted) {
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            List<int> colorOfParticle = new List<int>();
            for (int p = 0; p < Particles.Count; p++) {
                int colorOfCurrentParticle = 0;
                for (int i = 0; i < ColoredCellsSorted.Count; i++) {
                    if (ColoredCellsSorted[i][0] < J) {
                        Vector center = new Vector(gridData.iLogicalCells.GetCenter(ColoredCellsSorted[i][0]));
                        if (ColoredCellsSorted[i][1] != 0 && Particles[p].Contains(center, 1.5 * MinGridLength)) {
                            colorOfCurrentParticle = ColoredCellsSorted[i][1];
                            break;
                        }
                    }
                }
                colorOfParticle.Add(colorOfCurrentParticle);
            }
            return colorOfParticle.ToArray();
        }

        /// <summary>
        /// method to sort all cells by their color
        /// </summary>
        /// <param name="CellColor">
        /// all cells with their color
        /// </param>
        private List<int[]> ColoredCellsFindAndSort(int[] CellColor) {
            List<int[]> ColoredCellsSorted = new List<int[]>();
            int ListIndex;
            for (int CellID = 0; CellID < CellColor.Length; CellID++) {
                if (CellColor[CellID] == 0)
                    continue;
                ListIndex = 0;
                if (ColoredCellsSorted.Count != 0) {
                    while (ListIndex < ColoredCellsSorted.Count && CellColor[CellID] >= ColoredCellsSorted[ListIndex][1]) {
                        ListIndex += 1;
                    }
                }
                ColoredCellsSorted.Insert(ListIndex, new int[] { CellID, CellColor[CellID] });
            }
            return ColoredCellsSorted;
        }

        private void ColorNeighborCells(int[] coloredCells, int[] coloredCellsExchange) {
            int neighbourSearchDepth = 2;
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int k = 0; k < neighbourSearchDepth; k++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    GridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (coloredCellsExchange[CellNeighbors[i]] != 0 && coloredCellsExchange[j] == 0) {
                            coloredCells[j] = coloredCellsExchange[CellNeighbors[i]];
                        }
                    }
                }
                coloredCellsExchange = coloredCells.CloneAs();
                coloredCellsExchange.MPIExchange(GridData);
            }
        }

        private int[] FindCellsToRecolor(int[] coloredCells, GridData currentGrid) {
            int[] globalCellColor = GetGlobalCellColor(coloredCells, currentGrid);
            int[][] globalCellNeighbourship = GetGlobalCellNeigbourship(currentGrid);
            int maxColor = globalCellColor.Max().MPIMax();
            int[] newColor = new int[maxColor + 1];
            //for (int i = 0; i < globalCellColor.Length; i++) {
            //    for (int j = 0; j < globalCellNeighbourship[i].Length; j++) {
            //        if (globalCellColor[i] != globalCellColor[globalCellNeighbourship[i][j]] && globalCellColor[globalCellNeighbourship[i][j]] != 0) {
            //            if (newColor[globalCellColor[globalCellNeighbourship[i][j]]] != 0) {
            //                if (newColor[globalCellColor[i]] != 0) {// && newColor[globalCellColor[i]] > globalCellColor[globalCellNeighbourship[i][j]]) {
            //                    for (int k = newColor.Length - 1; k > 0; k--) {
            //                        if (k == newColor[globalCellColor[i]]) {
            //                            RecolorAlreadyRecoloredCellsRecursive(newColor, newColor[globalCellColor[i]], newColor[globalCellColor[globalCellNeighbourship[i][j]]]);
            //                            newColor[k] = newColor[globalCellColor[globalCellNeighbourship[i][j]]];
            //                        }
            //                    }
            //                }
            //                newColor[globalCellColor[i]] = newColor[globalCellColor[globalCellNeighbourship[i][j]]];
            //            }
            //            else {
            //                if (newColor[globalCellColor[i]] != 0) {// && newColor[globalCellColor[i]] > globalCellColor[globalCellNeighbourship[i][j]]) {
            //                    for (int k = newColor.Length - 1; k > 0; k--) {
            //                        if (k == newColor[globalCellColor[i]]) {
            //                            RecolorAlreadyRecoloredCellsRecursive(newColor, newColor[globalCellColor[i]], globalCellColor[globalCellNeighbourship[i][j]]);
            //                            newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
            //                        }
            //                    }
            //                }
            //                newColor[globalCellColor[i]] = globalCellColor[globalCellNeighbourship[i][j]];
            //            }
            //        }
            //    }
            //}
            for (int i = 0; i < globalCellColor.Length; i++) {
                for (int j = 0; j < globalCellNeighbourship[i].Length; j++) {
                    if (globalCellColor[i] != globalCellColor[globalCellNeighbourship[i][j]] && globalCellColor[globalCellNeighbourship[i][j]] != 0 && globalCellColor[i] != 0) {
                        if(newColor[globalCellColor[i]] != 0) {
                            for(int k = 1; k < newColor.Length; k++) {
                                if(newColor[k] == newColor[globalCellColor[i]]) {
                                    newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
                                }
                                if(k == newColor[globalCellColor[i]]) {
                                    newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
                                }
                            }
                        }
                        newColor[globalCellColor[i]] = globalCellColor[globalCellNeighbourship[i][j]];
                        if(newColor[globalCellColor[globalCellNeighbourship[i][j]]] != 0) {
                            for (int k = 1; k < newColor.Length; k++) {
                                if (newColor[k] == newColor[globalCellColor[globalCellNeighbourship[i][j]]]) {
                                    newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
                                }
                                if (k == newColor[globalCellColor[globalCellNeighbourship[i][j]]]) {
                                    newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
                                }
                            }
                        }
                        newColor[globalCellColor[globalCellNeighbourship[i][j]]] = globalCellColor[globalCellNeighbourship[i][j]];
                    }
                }
            }
            return newColor;
        }

        private void RecolorAlreadyRecoloredCellsRecursive(int[] newColorArray, int currentColor, int newColor) {
            if (newColor == newColorArray[currentColor])
                return;
            if (newColorArray[currentColor] != 0){// && currentColor > newColor) {
                for (int k = newColorArray.Length - 1; k > 0; k--) {
                    if (k == newColorArray[currentColor]) {
                        newColorArray[k] = newColor;
                        RecolorAlreadyRecoloredCellsRecursive(newColorArray, k, newColor);
                    }
                }
            }
        }

        private int[] GetGlobalCellColor(int[] localCellColor, GridData currentGrid) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;
            int[] i0 = cellPartitioning.GetI0s();
            int[] globalCellColor = new int[globalJ];

            int[][] exchangeCellColor = localCellColor.MPIGatherO(0);
            exchangeCellColor = exchangeCellColor.MPIBroadcast(0);

            for (int m = 0; m < currentGrid.MpiSize; m++) {
                for (int j = i0[m]; j < i0[m+1]; j++) {
                    globalCellColor[j] = exchangeCellColor[m][j - i0[m]];
                }
            }
            return globalCellColor;
        }

        private int[][] GetGlobalCellNeigbourship(GridData currentGrid) {
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

        private void RecolorCellsOfNeighborParticles(int[] coloredCells, GridData currentGrid) {
            int[] newColor = FindCellsToRecolor(coloredCells, currentGrid);
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int i = 1; i < newColor.Length; i++) {
                if (newColor[i] != 0) {
                    for (int j = 0; j < J; j++) {
                        if (coloredCells[j] == i) {
                            coloredCells[j] = newColor[i];
                        }
                    }
                }
            }
        }
    }
}
