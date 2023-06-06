using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.LoadBalancing {



    /// <summary>
    /// Employs the state of the level-set-tracker (cut-cell, double-cut-cell, near-cell, etc. as a classification criterion<see cref="CutStateClassifier.CellTypeFlags"/>)
    /// </summary>
    [Serializable]
    public class CutStateClassifier : TrackerStateClassifier {
        
        override public int[] ClassifyCells(IApplication app) {

            var lsTrk = app.LsTrk;
            int J = app.GridData.iLogicalCells.NoOfLocalUpdatedCells;

            bool CountAlsoVoidSpecies = false;

            var ret = new int[J];
            if (lsTrk != null) {
                var rg = lsTrk.Regions;
                int NoOfLevSets = lsTrk.NoOfLevelSets;

                int MinDistLevSet(int j, out int NoOfCuts) {
                    int minDist = int.MaxValue;
                    NoOfCuts = 0;
                    for (int iLs = 0; iLs < NoOfLevSets; iLs++) {
                        int dist_iLs = Math.Abs(rg.GetLevelSetDistance(iLs, j));
                        if (dist_iLs == 0)
                            NoOfCuts++;
                        minDist = Math.Min(minDist, dist_iLs);
                    }
                    return minDist;
                }


                SpeciesId voidId = default(SpeciesId);
                SpeciesId[] AllNonVoidSpecies;

                if (!VoidSpecies.IsEmptyOrWhite() && lsTrk.SpeciesNames.Contains(VoidSpecies)) {
                    voidId = lsTrk.GetSpeciesId(this.VoidSpecies);
                    var _AllNonVoidSpecies = lsTrk.SpeciesIdS.ToList();
                    _AllNonVoidSpecies.Remove(voidId);
                    AllNonVoidSpecies = _AllNonVoidSpecies.ToArray();
                } else {
                    CountAlsoVoidSpecies = true;
                    AllNonVoidSpecies = lsTrk.SpeciesIdS.ToArray();
                }

                bool ContainsAnyNonVoid(int j) {
                    foreach (var spc in AllNonVoidSpecies)
                        if (rg.IsSpeciesPresentInCell(spc, j))
                            return true;

                    return false;
                }

                // Assign default CellTypeFlag for each cell
                ret.SetAll((int)CellTypeFlags.Ordinary);


                // Classify cells w.r.t. their distance (this is not geometric  distance) check GetLevelSetDistance
                int NearWidth = lsTrk.NearRegionWidth;

                for (int j = 0; j < J; j++) {
                    int dist_j = MinDistLevSet(j, out int NoOfCuts_j);
                    if (dist_j == 0)
                        ret[j] |= (int)CellTypeFlags.Cut;
                    if (NoOfCuts_j > 1)
                        ret[j] |= (int)CellTypeFlags.DoubleCut;

                    if (dist_j > 0 && dist_j <= NearWidth)
                        ret[j] |= (int)CellTypeFlags.Near;

                    if (CountAlsoVoidSpecies == false) {
                        bool isPresent = rg.IsSpeciesPresentInCell(voidId, j);
                        bool nixAnderes = ContainsAnyNonVoid(j);

                        if (nixAnderes && isPresent)
                            ret[j] |= (int)CellTypeFlags.Void;
                    }

                }
            } else {
                Console.WriteLine("LevelSet Tracker is not set. Load balancing may not work as intended!!");
            }

            return ret;
        }



        /// <summary>
        /// Constants/Flags enumerating different cell types
        /// </summary>
        public enum CellTypeFlags {

            /// <summary>
            /// Default
            /// </summary>
            Ordinary = 0x1,

            /// <summary>
            /// At least one Level-Set is zero
            /// </summary>
            Cut = 0x2,

            /// <summary>
            /// More than one Level-Set is zero
            /// </summary>
            DoubleCut = 0x4,

            /// <summary>
            /// Cell is not cut, but next to some level-set
            /// </summary>
            Near = 0x01,

            /// <summary>
            /// cell contains only the void species <see cref="TrackerStateClassifier.VoidSpecies"/>, if specified
            /// </summary>
            Void = 0x0,

        }

        public override object Clone() {
            return new CutStateClassifier() {
                VoidSpecies = this.VoidSpecies
            };
        }

    }

}