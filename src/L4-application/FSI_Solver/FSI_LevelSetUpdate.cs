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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class FSI_LevelSetUpdate
    {
        private FSI_Auxillary Auxillary = new FSI_Auxillary();

        /// <summary>
        /// Cell Mask of all cells with a specific color "CurrentColor"
        /// </summary>
        internal CellMask CellsOneColor(IGridData gridData, List<int[]> ColoredCellsSorted, int CurrentColor, int J)
        {
            int[] CellIDCurrentColor = FindCellIDs(ColoredCellsSorted, CurrentColor);
            BitArray ColoredCells = new BitArray(J);
            for (int i = 0; i < CellIDCurrentColor.Length; i++)
            {
                if (CellIDCurrentColor[i] < J)
                    ColoredCells[CellIDCurrentColor[i]] = true;
            }
            CellMask ColoredCellMask = new CellMask(gridData, ColoredCells);
            CellMask ColoredCellMaskNeighbour = ColoredCellMask.AllNeighbourCells();
            return ColoredCellMask.Union(ColoredCellMaskNeighbour);
        }

        private int[] FindCellIDs(List<int[]> ColoredCellsSorted, int CurrentColor)
        {
            List<int> ColorList = new List<int>();
            for (int i = 0; i < ColoredCellsSorted.Count(); i++)
            {
                ColorList.Add(ColoredCellsSorted[i][1]);
            }
            int[] ListID = Auxillary.BinarySearchWithNeighbors(ColorList, CurrentColor);
            int[] CellID = new int[ListID.Length];
            for (int i = 0; i < ListID.Length; i++)
            {
                CellID[i] = ColoredCellsSorted[ListID[i]][0];
            }
            return CellID;
        }

        internal int[] FindParticlesOneColor(int[] ParticleColor, int CurrentColor)
        {
            List<int> temp = new List<int>();
            for (int i = 0; i < ParticleColor.Length; i++)
            {
                if (ParticleColor[i] == CurrentColor)
                {
                    temp.Add(i);
                }
            }
            return temp.ToArray();
        }

        // Method to find the Color of a specific particle
        // ===================================================
        public int[] FindParticleColor(IGridData GrdDat, List<Particle> Particles, List<int[]> ColoredCellsSorted)
        {
            List<int> CurrentColor = new List<int>();
            for (int p = 0; p < Particles.Count; p++)
            {
                double[] ParticleScales = Particles[p].GetLengthScales();
                double Lengthscale = ParticleScales.Min();
                double[] ParticlePos = Particles[p].Position[0];
                double Upperedge = ParticlePos[1] + Lengthscale;
                double Loweredge = ParticlePos[1] - Lengthscale;
                double Leftedge = ParticlePos[0] - Lengthscale;
                double Rightedge = ParticlePos[0] + Lengthscale;
                int temp = 0;
                for (int i = 0; i < ColoredCellsSorted.Count; i++)
                {
                    if (Math.Sqrt(GrdDat.iGeomCells.GetCellVolume(ColoredCellsSorted[i][0])) > Lengthscale)
                        throw new ArithmeticException("Hmin of the cells is larger than the particles. Please use a finer grid (or grid refinement).");

                    double[] center = GrdDat.iLogicalCells.GetCenter(ColoredCellsSorted[i][0]);
                    if (center[0] > Leftedge && center[0] < Rightedge && center[1] > Loweredge && center[1] < Upperedge)
                    {
                        temp = ColoredCellsSorted[i][1];
                        break;
                    }
                }
                CurrentColor.Add(temp);
            }
            return CurrentColor.ToArray();
        }

        internal List<int[]> ColoredCellsFindAndSort(int[] CellColor, bool RemoveUncoloredCells)
        {
            List<int[]> ColoredCellsSorted = new List<int[]>();
            int ListIndex = 0;
            for (int CellID = 0; CellID < CellColor.Length; CellID++)
            {
                if (RemoveUncoloredCells == false || (CellColor[CellID] != 0 && RemoveUncoloredCells == true))
                {
                    ListIndex = 0;
                    if (ColoredCellsSorted.Count != 0)
                    {
                        while (ListIndex < ColoredCellsSorted.Count && CellColor[CellID] >= ColoredCellsSorted[ListIndex][1])
                        {
                            ListIndex += 1;
                        }
                    }
                    int[] temp = new int[2];
                    temp[0] = CellID;
                    temp[1] = CellColor[CellID];
                    ColoredCellsSorted.Insert(ListIndex, temp);
                }
            }
            return ColoredCellsSorted;
        }
    }
}
