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
        // The main method to assign the color to the particle
        // ===================================================
        public void AssignToParticle(LevelSetTracker LsTrk, List<Particle> Particles, int[] CellColor)
        {
            // Step 1
            // Define variables
            // ================
            int MaxColor = CellColor.Max();
            int NoOfParticles = Particles.Count();
            List<int[]> ColoredCellsSorted = new List<int[]>();
            List<int[]> ColorsWithMultipleParticles = new List<int[]>();

            // Step 2 
            // Find all colored cells and sort them
            // ====================================
            ColoredCellsSorted = ColoredCellsFindAndSort(CellColor, RemoveUncoloredCells: true);

            // Step 3
            // Save all colored cells for later usage
            // ======================================
            // int[,] ColoredCells = new int[2, ColoredCellsSorted.Count()];
            // for (int i = 0; i < ColoredCellsSorted.Count(); i++)
            // {
            //    ColoredCells[0, i] = ColoredCellsSorted[i][0];
            //    ColoredCells[1, i] = ColoredCellsSorted[i][1];
            // }

            // Step 4
            // Assign all colors to a particle
            // ===============================
            for (int i = ColoredCellsSorted.Count - 1; i >= 0; i--)
            {
                for (int p = 0; p < Particles.Count; p++)
                {
                    if (LsTrk.GridDat.Cells.IsInCell(Particles[p].Position[0], i))
                    {
                        bool SameNoOfColorParticles = CheckForSameNumber(MaxColor, NoOfParticles);
                        int CurrentColor = ColoredCellsSorted[i][1];

                        // Assign color and cells to particle
                        // ========================
                        int ListPos = ColoredCellsSorted.Count - 1;
                        while (ListPos >= 0 && CurrentColor == ColoredCellsSorted[ListPos][1])
                        {
                            Particles[p].ParticleColoredCells.Add(ColoredCellsSorted[ListPos]);

                            // Remove all cells with current colors if there are no other particles with the same color left (saves some operations)
                            // =====================================================================================================================
                            if (SameNoOfColorParticles)
                            {
                                ColoredCellsSorted.RemoveAt(ListPos);
                                MaxColor -= 1;
                            }
                            ListPos -= 1;
                        }

                        // Save colors with more than one particle in a seperate list
                        // ==========================================================
                        ColorsWithMultipleParticles = CheckForMultipleParticlesPerColor(ColorsWithMultipleParticles, CurrentColor, p, SameNoOfColorParticles);
                        NoOfParticles -= 1;
                    }
                }
            }
        }

        // Method to find the Color of a specific particle
        // ===================================================
        public int[] FindParticleColor(IGridData GrdDat, List<Particle> Particles, List<int[]> ColoredCellsSorted)
        {
            int[] CurrentColor = new int[Particles.Count];
            for (int p = 0; p < Particles.Count; p++)
            {
                double[] ParticleScales = Particles[p].GetLengthScales();
                double Lengthscale = ParticleScales.Min();
                double[] ParticlePos = Particles[p].Position[0];
                double Upperedge = ParticlePos[1] + Lengthscale;
                double Loweredge = ParticlePos[1] - Lengthscale;
                double Leftedge = ParticlePos[0] - Lengthscale;
                double Rightedge = ParticlePos[0] + Lengthscale;
                for (int i = 0; i < ColoredCellsSorted.Count; i++)
                {
                    if (Math.Sqrt(GrdDat.iGeomCells.GetCellVolume(ColoredCellsSorted[i][0])) > Lengthscale)
                        throw new ArithmeticException("Hmin of the cells is larger than the particles. Please use a finer grid (or grid refinement).");

                    double[] center = GrdDat.iLogicalCells.GetCenter(ColoredCellsSorted[i][0]);
                    if (center[0] > Leftedge && center[0] < Rightedge && center[1] > Loweredge && center[1] < Upperedge)
                        CurrentColor[p] = ColoredCellsSorted[i][1];
                }
            }
            return CurrentColor;
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

        internal int[] FindNeighborCells(List<int[]> ParticleColoredCells, LevelSetTracker LsTrk)
        {
            List<int> ColoredAndNeighborCells = new List<int>();
            for (int i = 0; i < ParticleColoredCells.Count; i++)
            {
                ColoredAndNeighborCells.Add(ParticleColoredCells[i][0]);
            }
            
            for (int i = 0; i < ParticleColoredCells.Count(); i++)
            {
                LsTrk.GridDat.GetCellNeighbours(ParticleColoredCells[i][0], GetCellNeighbours_Mode.ViaEdges, out int[] NeighborsOfCurrentCell, out int[] ConectingEntities);
                for (int j = NeighborsOfCurrentCell.Length - 1; j >= 0; j--)
                {
                    int ListIndex = ParticleColoredCells.Count();
                    while (ListIndex > 0 && NeighborsOfCurrentCell[j] < ParticleColoredCells[ListIndex - 1][0])
                    {
                        ListIndex -= 1;
                    }
                    if (ListIndex == 0 || ParticleColoredCells[ListIndex - 1][0] != NeighborsOfCurrentCell[j])
                    {
                        int[] temp = new int[2];
                        temp[0] = NeighborsOfCurrentCell[0];
                        temp[1] = ParticleColoredCells[i][1];
                        ColoredAndNeighborCells.Insert(ListIndex, NeighborsOfCurrentCell[0]);
                    }
                }
            }
            return ColoredAndNeighborCells.ToArray();
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

        // Compare the number of particles with the number of colors
        // =========================================================
        private static bool CheckForSameNumber(int MaxColor, int NoOfParticles)
        {
            return NoOfParticles - MaxColor == 0 ? true : false;
        }

        // There may be cases where multiple particles have the same color and need a special treatment.
        // With this method we save those particles in a seperate list
        // =============================================================================================
        private static List<int[]> CheckForMultipleParticlesPerColor(List<int[]> ColorsWithMultipleParticles, int CurrentColor, int CurrentParticle, bool SameNoOfColorParticles = true)
        {
            int Length = ColorsWithMultipleParticles.Count();
            if(SameNoOfColorParticles == false || Length == 0 || CurrentColor == ColorsWithMultipleParticles[Length - 1][0])
            {
                int[] temp = new int[2];
                temp[0] = CurrentColor;
                temp[1] = CurrentParticle;
                ColorsWithMultipleParticles.Add(temp);
            }
            else if(ColorsWithMultipleParticles[Length - 1][0] != ColorsWithMultipleParticles[Length - 2][0])
            {
                ColorsWithMultipleParticles.RemoveAt(Length - 1);
            }
            return ColorsWithMultipleParticles;
        }

        void ParticleToColor()
        {

        }
    }
}
