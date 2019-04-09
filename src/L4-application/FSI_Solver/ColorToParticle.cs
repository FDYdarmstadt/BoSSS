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

using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Application.FSI_Solver;
using BoSSS.Foundation.Grid;
using System;

namespace FSI_Solver
{
    class ColorToParticle
    {
        private Auxillary Quicksort = new Auxillary();

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

        // Method to find all Cells of a specific color
        // ===================================================
        public int[] FindCellsWithColor(int Color, List<int[]> ColoredCellsSorted)
        {
            List<int> Cells = new List<int>();
            for (int i = 0; i < ColoredCellsSorted.Count; i++)
            {
                if (ColoredCellsSorted[i][1] == Color)
                    Cells.Add(ColoredCellsSorted[i][0]);
            }
            return Cells.ToArray();
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

        internal int[] FindMulitpleParticlesInColor(int[] ParticleColors, int CurrentColor)
        {
            int[] NumberOfParticlesInColor = new int[2];
            int[] _particleColors = ParticleColors;
            Auxillary.Main(0, ParticleColors.Length - 1, ref _particleColors);
            for (int i = 0; i < ParticleColors.Length; i++)
            {
                while (ParticleColors[i] == CurrentColor)
                    NumberOfParticlesInColor[0] += 1;

                if (NumberOfParticlesInColor[0] != 0 && ParticleColors[i] > CurrentColor)
                    break;
            }
            NumberOfParticlesInColor[1] = CurrentColor;
            return NumberOfParticlesInColor;
        }
    }
}
