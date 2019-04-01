using BoSSS.Application.FSI_Solver;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FSI_Solver
{
    class ColorToParticle
    {        
        // Compare the number of particles with the number of colors
        // =========================================================
        static bool CheckForSameNumber(int MaxColor, int NoOfParticles)
        {
            if (NoOfParticles - MaxColor == 0)
            {
                return true;
            }
            else
                return false;
        }

        static List<int[]> RemoveCellsWithCurrentColor(List<int[]> ColoredCells, int CurrentColor)
        {
            int ListID = ColoredCells.Count() - 1;
            while (ColoredCells.Count() > 0 && ColoredCells[ListID][1] == CurrentColor)
            {
                ColoredCells.RemoveAt(ListID);
                ListID -= 1;
            }
            return ColoredCells;
        }

        List<int[]> AddToList(List<int[]> ColorsWithMultipleParticles, int CurrentColor, int CurrentParticle)
        {
            int[] temp = new int[2];
            temp[0] = CurrentColor;
            temp[1] = CurrentParticle;
            ColorsWithMultipleParticles.Add(temp);
            return ColorsWithMultipleParticles;
        }

        List<int[]> CheckForMultipleParticlesPerColor(List<int[]> ColorsWithMultipleParticles, int CurrentColor, int CurrentParticle, bool SameNoOfColorParticles = true)
        {
            int Length = ColorsWithMultipleParticles.Count();
            if(SameNoOfColorParticles == false || Length == 0 || CurrentColor == ColorsWithMultipleParticles[Length - 1][0])
            {
                ColorsWithMultipleParticles = AddToList(ColorsWithMultipleParticles, CurrentColor, CurrentParticle);
            }
            else if(ColorsWithMultipleParticles[Length - 1][0] != ColorsWithMultipleParticles[Length - 2][0])
            {
                ColorsWithMultipleParticles.RemoveAt(Length - 1);
            }
            return ColorsWithMultipleParticles;
        }

        void ColoredCellsToParticle(LevelSetTracker LsTrk, List<Particle> Particles)
        {
            // Step 1
            // Define variables
            // ================
            int[] ParticleColor = LsTrk.Regions.ColorMap4Spc[LsTrk.GetSpeciesId("B")];
            int MaxColor = ParticleColor.Max();
            int NoOfParticles = Particles.Count();
            List<int[]> ColoredCellsSorted = new List<int[]>();
            int ListIndex = 0;
            List<int[]> ColorsWithMultipleParticles = new List<int[]>();

            // Step 2 
            // Find all colored cells and sort them
            // ====================================
            for (int CellID = 0; CellID < ParticleColor.Length; CellID++)
            {
                if (ParticleColor[CellID] != 0)
                {
                    ListIndex = 0;
                    if (ColoredCellsSorted.Count != 0)
                    {
                        while (ListIndex < ColoredCellsSorted.Count && ParticleColor[CellID] >= ColoredCellsSorted[ListIndex][1])
                        {
                            ListIndex += 1;
                        }
                    }
                    int[] temp = new int[2];
                    temp[0] = CellID;
                    temp[1] = ParticleColor[CellID];
                    ColoredCellsSorted.Insert(ListIndex, temp);
                }
            }

            // Step 3
            // Assign all colors to a particle
            // ===============================
            for (int i = ColoredCellsSorted.Count - 1; i >= 0; i--)
            {
                for (int p = 0; p < Particles.Count; p++)
                {
                    if (LsTrk.GridDat.Cells.IsInCell(Particles[p].Position[0], i))
                    {
                        bool SameNoOfColorParticles = Sa
                        int CurrentColor = ColoredCellsSorted[i][1];
                        Particles[p].ParticleColor = CurrentColor;

                        if (CheckForSameNumber(MaxColor, NoOfParticles))
                        {
                            ColoredCellsSorted = RemoveCellsWithCurrentColor(ColoredCellsSorted, CurrentColor);
                            MaxColor -= 1;
                            
                        }
                        ColorsWithMultipleParticles = CheckForMultipleParticlesPerColor(ColorsWithMultipleParticles, CurrentColor, p, );
                        NoOfParticles -= 1;
                    }
                }
            }
        }

        void ParticleToColor()
        {

        }
    }
}
