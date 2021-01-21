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
using System.Diagnostics;
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
            int JE = GridData.iLogicalCells.NoOfExternalCells + GridData.iLogicalCells.NoOfLocalUpdatedCells;
            MultidimensionalArray CellCenters = LsTrk.GridDat.Cells.CellCenter;
            int[] coloredCells = new int[JE];
            for (int p = 0; p < Particles.Length; p++) {
                for (int j = 0; j < GridData.iLogicalCells.NoOfLocalUpdatedCells; j++) {
                    if (Particles[p].Contains(new Vector(CellCenters[j, 0], CellCenters[j, 1]), MaxGridLength))
                        coloredCells[j] = p + 1;
                }
            }
            coloredCells.MPIExchange(GridData);
            RecolorNeighbouringCells(coloredCells);
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
            RecolorNeighbouringCells(coloredCells);
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
                    int masterID = ParticleList[i].MasterDuplicateIDs[0] - 1;
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
        internal bool CurrentProcessContainsCurrentColor(int[] CellColor, int CurrentColor, int NoOfLocalCells) {
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
            int noOfLocalCells = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < noOfLocalCells; j++) {
                Tuple<int, int, int>[] cellNeigbours = GridData.GetCellNeighboursViaEdges(j);
                for (int i = 0; i < cellNeigbours.Length; i++) {
                    if (coloredCellsExchange[cellNeigbours[i].Item1] != 0 && coloredCellsExchange[j] == 0) {
                        coloredCells[j] = coloredCellsExchange[cellNeigbours[i].Item1];
                    }
                }
            }
            coloredCellsExchange = coloredCells.CloneAs();
            coloredCellsExchange.MPIExchange(GridData);
        }

        private void RecolorNeighbouringCells(int[] ColoredCells) {
            int[][] globalColorNeighbours = GlobalColorNeighbours(ColoredCells);
            int[] recolorWith = new int[globalColorNeighbours.Length];
            for (int i = globalColorNeighbours.Length - 1; i > 0; i--) {
                if (globalColorNeighbours[i].Length > 0) {
                    if (globalColorNeighbours[i].Min() < i) {
                        recolorWith[i] = globalColorNeighbours[i].Min();
                        for (int j = 0; j < globalColorNeighbours[i].Length; j++) {
                            recolorWith[globalColorNeighbours[i][j]] = globalColorNeighbours[i].Min();
                            RecolorRecursive(globalColorNeighbours[i][0], globalColorNeighbours[i][j], recolorWith, globalColorNeighbours);
                        }
                    } else {
                        recolorWith[i] = i;
                        for (int j = 0; j < globalColorNeighbours[i].Length; j++) {
                            recolorWith[globalColorNeighbours[i][j]] = i;
                            RecolorRecursive(i, globalColorNeighbours[i][j], recolorWith, globalColorNeighbours);
                        }
                    }
                } else {
                    recolorWith[i] = i;
                }
            }

            for (int i = 0; i < ColoredCells.Length; i++) {
                if (ColoredCells[i] != 0) {
                    ColoredCells[i] = recolorWith[ColoredCells[i]];
                }
            }
        }

        private void RecolorRecursive(int currentColor, int neighbourColor, int[] recolorWith, int[][] GlobalColorNeighbours) {
            if (currentColor == 0)
                return;
            for (int j = 0; j < GlobalColorNeighbours[neighbourColor].Length; j++) {
                if (recolorWith[GlobalColorNeighbours[neighbourColor][j]] != currentColor) {
                    recolorWith[GlobalColorNeighbours[neighbourColor][j]] = currentColor;
                    RecolorRecursive(currentColor, GlobalColorNeighbours[neighbourColor][j], recolorWith, GlobalColorNeighbours);
                }
            }
        }

        private int[][] GlobalColorNeighbours(int[] coloredCells) {
            List<int>[] localColorNeighbours = LocalColorNeighbours(coloredCells, coloredCells.Max().MPIMax());
            int[][] globalColorNeighbours = new int[localColorNeighbours.Length][];
            for (int i = 1; i < localColorNeighbours.Length; i++) {
                List<int>[] exchangeVariable = localColorNeighbours[i].MPIAllGatherO();
                for (int j = 0; j < exchangeVariable.Length; j++)
                    localColorNeighbours[i].AddRange(exchangeVariable[j]);
                globalColorNeighbours[i] = localColorNeighbours[i].ToArray();
            }
            return globalColorNeighbours;
        }

        private List<int>[] LocalColorNeighbours(int[] coloredCells, int GlobalMaxColor) {
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            List<int>[] colorNeighbours = new List<int>[GlobalMaxColor + 1];
            for (int i = 0; i < colorNeighbours.Length; i++)
                colorNeighbours[i] = new List<int>();
            for (int j = 0; j < J; j++) {
                int currentColor = coloredCells[j];
                if (currentColor != 0) {
                    Tuple<int, int, int>[] cellNeigbours = GridData.GetCellNeighboursViaEdges(j);
                    for (int i = 0; i < cellNeigbours.Length; i++) {
                        int neighbourColor = coloredCells[cellNeigbours[i].Item1];
                        if (neighbourColor != 0 && neighbourColor != currentColor) {
                            if (!colorNeighbours[currentColor].Contains(neighbourColor))
                                colorNeighbours[currentColor].Add(neighbourColor);
                            if (!colorNeighbours[neighbourColor].Contains(currentColor))
                                colorNeighbours[neighbourColor].Add(currentColor);
                        }
                    }
                }
            }
            return colorNeighbours;
        }
    }
}
