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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;

namespace FSI_Solver {
    class FSI_LevelSetUpdate {

        internal FSI_LevelSetUpdate(LevelSetTracker levelSetTracker) {
            m_LevelSetTracker = levelSetTracker;
        }

        private readonly LevelSetTracker m_LevelSetTracker;

        /// <summary>
        /// compiles a cell mask from all cells with a specific color and their direct neighbors
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
        internal CellMask CellsOneColor(IGridData gridData, List<int[]> ColoredCellsSorted, int CurrentColor, int J, bool FindNeighbours = true) {
            int[] CellIDCurrentColor = FindCellIDs(ColoredCellsSorted, CurrentColor);
            BitArray ColoredCells = new BitArray(J);
            for (int i = 0; i < CellIDCurrentColor.Length; i++) {
                if (CellIDCurrentColor[i] < J)
                    ColoredCells[CellIDCurrentColor[i]] = true;
            }
            CellMask ColoredCellMask = new CellMask(gridData, ColoredCells);

            if (FindNeighbours) {
                CellMask ColoredCellMaskNeighbour = ColoredCellMask.AllNeighbourCells();
                ColoredCellMask = ColoredCellMask.Union(ColoredCellMaskNeighbour);
            }
            return ColoredCellMask;
        }

        /// <summary>
        /// searchs for all cell ID with the current color
        /// </summary>
        /// <param name="ColoredCellsSorted">
        /// a list of all cells sorted by color
        /// </param>
        /// <param name="CurrentColor">
        /// the color which is associated with the Cell Mask
        /// </param>
        private int[] FindCellIDs(List<int[]> ColoredCellsSorted, int CurrentColor) {
            List<int> ColorList = new List<int>();
            for (int i = 0; i < ColoredCellsSorted.Count(); i++) {
                ColorList.Add(ColoredCellsSorted[i][1]);
            }
            int[] ListID = BinarySearchWithNeighbors(ColorList, CurrentColor);
            int[] CellID = new int[ListID.Length];
            for (int i = 0; i < ListID.Length; i++) {
                CellID[i] = ColoredCellsSorted[ListID[i]][0];
            }
            return CellID;
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
        internal int[] BinarySearchWithNeighbors(List<int> SortedList, int Target) {
            int StartIndex = BinarySearch(SortedList, Target);
            if (StartIndex == 0)
                return LeftSearch(SortedList, StartIndex);
            else if (StartIndex == SortedList.Count() - 1)
                return RightSearch(SortedList, StartIndex);
            else {
                int[] Right = RightSearch(SortedList, StartIndex);
                int[] Left = LeftSearch(SortedList, StartIndex);
                return Left.Union(Right).ToArray();
            }

        }

        /// <summary>
        /// Binary search algorithm
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="Target">
        /// The element to be found
        /// </param>
        int BinarySearch(List<int> SortedList, int Target) {
            int L = 0;
            int R = SortedList.Count() - 1;
            int TargetIndex = 0;

            while (L <= R) {
                TargetIndex = (L + R) / 2;
                if (SortedList[TargetIndex] < Target)
                    L = TargetIndex + 1;
                else if (SortedList[TargetIndex] > Target)
                    R = TargetIndex - 1;
                else
                    break;
            }

            return TargetIndex;
        }

        /// <summary>
        /// Search for left neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="StartIndex"></param>
        private int[] LeftSearch(List<int> SortedList, int StartIndex) {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() + 1 < SortedList.Count() && SortedList[temp.Last() + 1] == SortedList[temp.Last()]) {
                temp.Add(temp.Last() + 1);
            }
            return temp.ToArray();
        }

        /// <summary>
        /// Search for right neighbours
        /// </summary>
        /// <param name="SortedList">
        /// A sorted list of all elements
        /// </param>
        /// <param name="StartIndex"></param>
        private int[] RightSearch(List<int> SortedList, int StartIndex) {
            List<int> temp = new List<int> { StartIndex };
            while (temp.Last() - 1 > 0 && SortedList[temp.Last() - 1] == SortedList[temp.Last()]) {
                temp.Add(temp.Last() - 1);
            }
            return temp.ToArray();
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
        /// searchs for all particles with the same color
        /// </summary>
        /// <param name="ParticleColor">
        /// a list of all particles with their specific color
        /// </param>
        /// <param name="CurrentColor">
        /// the color which is associated with the Cell Mask
        /// </param
        internal List<Particle> GetParticleListOneColor(List<Particle> AllParticles, int[] ParticleColor, int CurrentColor) {
            List<Particle> temp = new List<Particle>();
            for (int i = 0; i < ParticleColor.Length; i++) {
                if (ParticleColor[i] == CurrentColor) {
                    temp.Add(AllParticles[i]);
                }
            }
            return temp;
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
                double h_min = m_LevelSetTracker.GridDat.Cells.h_minGlobal > Particles[p].GetLengthScales().Min()
                    ? 1.5 * m_LevelSetTracker.GridDat.Cells.h_minGlobal
                    : 2 * m_LevelSetTracker.GridDat.Cells.h_minGlobal;
                double[] ParticlePos = Particles[p].Motion.position[0];
                double Upperedge = ParticlePos[1] + h_min;
                double Loweredge = ParticlePos[1] - h_min;
                double Leftedge = ParticlePos[0] - h_min;
                double Rightedge = ParticlePos[0] + h_min;
                int temp = 0;
                for (int i = 0; i < ColoredCellsSorted.Count; i++) {
                    if (ColoredCellsSorted[i][0] < J) {
                        double[] center = gridData.iLogicalCells.GetCenter(ColoredCellsSorted[i][0]);
                        if (center[0] > Leftedge && center[0] < Rightedge && center[1] > Loweredge && center[1] < Upperedge && ColoredCellsSorted[i][1] != 0) {
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
                ListIndex = 0;
                if (ColoredCellsSorted.Count != 0) {
                    while (ListIndex < ColoredCellsSorted.Count && CellColor[CellID] >= ColoredCellsSorted[ListIndex][1]) {
                        ListIndex += 1;
                    }
                }
                int[] temp = new int[2];
                temp[0] = CellID;
                temp[1] = CellColor[CellID];
                ColoredCellsSorted.Insert(ListIndex, temp);
            }
            return ColoredCellsSorted;
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
    }
}
