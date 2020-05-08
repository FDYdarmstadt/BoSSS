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

        internal FSI_LevelSetUpdate(IGridData gridData, double minGridLength) {
            this.gridData = gridData;
            this.minGridLength = minGridLength;
        }

        private readonly IGridData gridData;
        private readonly double minGridLength;


        /// <summary>
        /// compiles a cell mask from all cells with a specific color
        /// </summary>
        /// <param name="gridData">
        /// the grid that this mask will be associated with
        /// </param>
        /// <param name="ColoredCellsSorted">
        /// a list of all cells sorted by color
        /// </param>
        /// <param name="CurrentColor">
        /// the color which is associated with the Cell Mask
        /// </param>
        /// /// <param name="J">
        /// the number of locally updated cells
        /// </param>
        internal CellMask CellsOneColor(IGridData gridData, int[][] ColoredCellsSorted, int CurrentColor, int J) {
            int[] ListID = BinarySearchWithNeighbors(ColoredCellsSorted, CurrentColor);
            BitArray ColoredCells = new BitArray(J);
            for (int i = 0; i < ListID.Length; i++) {
                if (ColoredCellsSorted[ListID[i]][0] < J)
                    ColoredCells[ColoredCellsSorted[ListID[i]][0]] = true;
            }
            CellMask ColoredCellMask = new CellMask(gridData, ColoredCells);
            return ColoredCellMask;
        }

        /// <summary>
        /// Binary search algorithm, output is the target and its neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="Target">
        /// The element to be found
        /// </param>
        internal int[] BinarySearchWithNeighbors(int[][] SortedList, int Target) {
            int L = 0;
            int R = SortedList.Length - 1;
            List<int> targetID = new List<int>();
            while (L <= R) { // find a random cell with the current color
                int currentID = (L + R) / 2;
                Console.WriteLine("L " + L + " R " + R + " color " + SortedList[currentID][1]);
                if (SortedList[currentID][1] == Target) {
                    targetID.Add(currentID);
                    break;
                }
                else if (SortedList[currentID][1] > Target)
                    R = currentID - 1;
                else
                    L = currentID + 1;
            }
            Console.WriteLine("TargetID " + targetID.Count() +  "Target " +Target + "dafsd " + SortedList.Length);
            while (targetID[0] > 0) {
                if (SortedList[targetID[0] - 1][1] == Target) { // find all cells with the current color on the left side of the previous found cells
                    targetID.Insert(0, targetID[0] - 1);
                }
                else
                    break;
            }
            while (targetID.Last() < SortedList.Length - 1){
                if (SortedList[targetID.Last() + 1][1] == Target) { // find all cells with the current color on the right side of the previous found cells
                    targetID.Add(targetID.Last() + 1);
                }
                else
                    break;
            }
            return targetID.ToArray();
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
        internal int[] FindParticlesOneColor(int[] ParticleColor, int CurrentColor) {
            List<int> temp = new List<int>();
            for (int i = 0; i < ParticleColor.Length; i++) {
                if (ParticleColor[i] == CurrentColor) {
                    temp.Add(i);
                }
            }
            return temp.ToArray();
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
            int NoOfParticles = ParticleColorArray.Length;
            int[] GlobalParticleColor = new int[NoOfParticles];
            double[] StateBuffer = new double[NoOfParticles];
            for (int i = 0; i < NoOfParticles; i++) {
                StateBuffer[i] = Convert.ToDouble(ParticleColorArray[i]);
            }
            double[] GlobalStateBuffer = StateBuffer.MPIMax();
            for (int i = 0; i < NoOfParticles; i++) {
                GlobalParticleColor[i] = Convert.ToInt32(GlobalStateBuffer[i]);
            }
            return GlobalParticleColor;
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
        public int[] FindParticleColor(IGridData gridData, List<Particle> Particles, List<int[]> ColoredCellsSorted) {
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            List<int> CurrentColor = new List<int>();
            for (int p = 0; p < Particles.Count; p++) {
                double h_min = minGridLength > Particles[p].GetLengthScales().Min()
                    ? 1.5 * minGridLength
                    : 0;
                int temp = 0;
                for (int i = 0; i < ColoredCellsSorted.Count; i++) {
                    if (ColoredCellsSorted[i][0] < J) {
                        Vector center = new Vector(gridData.iLogicalCells.GetCenter(ColoredCellsSorted[i][0]));
                        if (ColoredCellsSorted[i][1] != 0 && Particles[p].Contains(center, h_min)) {
                            temp = ColoredCellsSorted[i][1];
                            break;
                        }
                    }

                }
                CurrentColor.Add(temp);
            }
            return CurrentColor.ToArray();
        }

        /// <summary>
        /// method to sort all cells by their color
        /// </summary>
        /// <param name="CellColor">
        /// all cells with their color
        /// </param>
        internal List<int[]> ColoredCellsFindAndSort(int[] CellColor) {
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

        internal void ColorNeighborCells(int[] coloredCells, int[] coloredCellsExchange) {
            int neighbourSearchDepth = 1;
            int noOfLocalCells = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int k = 0; k < neighbourSearchDepth; k++) {
                for (int j = 0; j < noOfLocalCells; j++) {
                    gridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (coloredCellsExchange[CellNeighbors[i]] != 0 && coloredCellsExchange[j] == 0) {
                            coloredCells[j] = coloredCellsExchange[CellNeighbors[i]];
                        }
                    }
                }
                coloredCellsExchange = coloredCells.CloneAs();
                coloredCellsExchange.MPIExchange(gridData);
            }
        }

        internal int[] FindCellsToRecolor(int[] coloredCells, GridData currentGrid) {
            //Debugger.Launch();
            int[] globalCellColor = GetGlobalCellColor(coloredCells, currentGrid);
            int[][] globalCellNeighbourship = GetGlobalCellNeigbourship(currentGrid);
            int maxColor = globalCellColor.Max().MPIMax();
            int[] newColor = new int[maxColor + 1];
            for (int i = 0; i < globalCellColor.Length; i++) {
                for (int j = 0; j < globalCellNeighbourship[i].Length; j++) {
                    if (globalCellColor[i] > globalCellColor[globalCellNeighbourship[i][j]] && globalCellColor[globalCellNeighbourship[i][j]] != 0) {
                        if (newColor[globalCellColor[i]] != 0 && newColor[globalCellColor[i]] > globalCellColor[globalCellNeighbourship[i][j]]) {
                            for (int k = newColor.Length - 1; k > 0; k--) {
                                if (k == newColor[globalCellColor[i]]) {
                                    RecolorAlreadyRecoloredCellsRecursive(newColor, newColor[globalCellColor[i]], globalCellColor[globalCellNeighbourship[i][j]]);
                                    newColor[k] = globalCellColor[globalCellNeighbourship[i][j]];
                                }
                            }
                        }
                        newColor[globalCellColor[i]] = globalCellColor[globalCellNeighbourship[i][j]];
                    }
                }
            }
            return newColor;
        }

        internal void RecolorAlreadyRecoloredCellsRecursive(int[] newColorArray, int currentColor, int newColor) {
            if (newColor == newColorArray[currentColor])
                return;
            if (newColorArray[currentColor] != 0 && currentColor > newColor) {
                for (int k = newColorArray.Length - 1; k > 0; k--) {
                    if (k == newColorArray[currentColor]) {
                        RecolorAlreadyRecoloredCellsRecursive(newColorArray, k, newColor);
                        newColorArray[k] = newColor;
                    }
                }
            }
        }

        private int[] GetGlobalCellColor(int[] localCellColor, GridData currentGrid) {
            Partitioning cellPartitioning = currentGrid.CellPartitioning;
            int globalJ = cellPartitioning.TotalLength;
            int[] i0 = cellPartitioning.GetI0s();
            int local_i0 = cellPartitioning.i0;
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

        internal void RecolorCellsOfNeighborParticles(int[] coloredCells, GridData currentGrid) {
            int[] newColor = FindCellsToRecolor(coloredCells, currentGrid);
            int J = gridData.iLogicalCells.NoOfLocalUpdatedCells;
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

        private int[][] FindCellsToRecolor(int[] coloredCells, int[] coloredCellsExchange, int MPISize) {
            int maxColor = coloredCells.Max().MPIMax();
            int noOfLocalCells = gridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[][] colorToRecolorWith = new int[maxColor + 1][];
            for (int k = 0; k < colorToRecolorWith.Length; k++) {
                colorToRecolorWith[k] = new int[2];
            }
            for (int j = 0; j < noOfLocalCells; j++) {
                if (coloredCells[j] != 0) {
                    gridData.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaEdges, out int[] CellNeighbors, out _);
                    for (int i = 0; i < CellNeighbors.Length; i++) {
                        if (coloredCellsExchange[CellNeighbors[i]] != coloredCells[j] && coloredCellsExchange[CellNeighbors[i]] > 0) {
                            if (coloredCellsExchange[CellNeighbors[i]] < coloredCells[j] || colorToRecolorWith[coloredCells[j]][1] > coloredCellsExchange[CellNeighbors[i]]) {
                                colorToRecolorWith[coloredCells[j]][0] = coloredCells[j];
                                colorToRecolorWith[coloredCells[j]][1] = coloredCellsExchange[CellNeighbors[i]];
                                if(colorToRecolorWith[coloredCells[j]][1] != 0) {
                                    for (int k = 0; k < colorToRecolorWith[k].Length; k++) {
                                        if(colorToRecolorWith[k][0] == colorToRecolorWith[coloredCells[j]][1])
                                            colorToRecolorWith[k][1] = coloredCellsExchange[CellNeighbors[i]];
                                    }
                                }
                            }
                            if (coloredCellsExchange[CellNeighbors[i]] > coloredCells[j]) {
                                if (colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][0] == 0 || colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][1] > coloredCells[j]) {
                                    colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][0] = coloredCellsExchange[CellNeighbors[i]];
                                    colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][1] = coloredCells[j];
                                    if (colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][1] != 0) {
                                        for (int k = 0; k < colorToRecolorWith[k].Length; k++) {
                                            if (colorToRecolorWith[k][0] == colorToRecolorWith[coloredCellsExchange[CellNeighbors[i]]][1])
                                                colorToRecolorWith[k][1] = coloredCellsExchange[CellNeighbors[i]];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Communicate
            // -----------------------------
            int[][][] GlobalColorToRecolorWith = colorToRecolorWith.MPIGatherO(0);
            GlobalColorToRecolorWith = GlobalColorToRecolorWith.MPIBroadcast(0);
            for (int m = 0; m < MPISize; m++) {
                for (int i = 0; i < maxColor + 1; i++) {
                    if (GlobalColorToRecolorWith[0][i][1] == 0 || GlobalColorToRecolorWith[0][i][1] > GlobalColorToRecolorWith[m][i][1] && GlobalColorToRecolorWith[m][i][1] != 0) {
                        GlobalColorToRecolorWith[0][i][0] = GlobalColorToRecolorWith[m][i][0];
                        GlobalColorToRecolorWith[0][i][1] = GlobalColorToRecolorWith[m][i][1];
                    }
                }
            }
            colorToRecolorWith = GlobalColorToRecolorWith[0];
            return colorToRecolorWith;
        }
    }
}
