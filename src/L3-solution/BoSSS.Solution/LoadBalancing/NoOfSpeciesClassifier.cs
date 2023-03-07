using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.LoadBalancing {



    /// <summary>
    /// Employs the number of species present in a cell as a classification criterion
    /// </summary>
    [Serializable]
    public class NoOfSpeciesClassifier : TrackerStateClassifier {

        /// <summary>
        /// - false: multiple species are only reported in cut-cells
        /// - true: multiple species are reported in near- and cut-cells
        /// </summary>
        public bool ConsiderAlsoNearCells = false;

        public override int[] ClassifyCells(IApplication app) {
        
            var lsTrk = app.LsTrk;
            int J = app.GridData.iLogicalCells.NoOfLocalUpdatedCells;

            bool CountAlsoVoidSpecies = false;

            var ret = new int[J];
            if (lsTrk != null) {
                var rg = lsTrk.Regions;
                int NoOfLevSets = lsTrk.NoOfLevelSets;

                SpeciesId voidId = default(SpeciesId);
                SpeciesId[] AllNonVoidSpecies;
                if (!VoidSpecies.IsEmptyOrWhite()) {
                    voidId = lsTrk.GetSpeciesId(this.VoidSpecies);
                    var _AllNonVoidSpecies = lsTrk.SpeciesIdS.ToList();
                    _AllNonVoidSpecies.Remove(voidId);
                    AllNonVoidSpecies = _AllNonVoidSpecies.ToArray();
                } else {
                    CountAlsoVoidSpecies = true;
                    AllNonVoidSpecies = lsTrk.SpeciesIdS.ToArray();
                }

                int CountNonVoid(int j) {
                    int r = 0;
                    foreach (var spc in AllNonVoidSpecies)
                        if (rg.IsSpeciesPresentInCell(spc, j))
                            r++;

                    return r;
                }

                int MinDistLevSet(int j) {
                    int dist = Math.Abs(rg.GetLevelSetDistance(0, j));
                    for (int iLs = 1; iLs < NoOfLevSets; iLs++)
                        dist = Math.Min(dist, Math.Abs(rg.GetLevelSetDistance(0, j)));
                    return dist;
                }


                if (ConsiderAlsoNearCells == false && CountAlsoVoidSpecies == false) {
                    for (int j = 0; j < J; j++)
                        ret[j] = rg.GetNoOfSpecies(j);
                } else if (ConsiderAlsoNearCells == false && CountAlsoVoidSpecies == true) {
                    for (int j = 0; j < J; j++)
                        ret[j] = CountNonVoid(j) + (rg.IsSpeciesPresentInCell(voidId, j) ? 1 : 0);

                } else if (ConsiderAlsoNearCells == true && CountAlsoVoidSpecies == false) {
                    for (int j = 0; j < J; j++)
                        ret[j] = CountNonVoid(j);

                } else if (ConsiderAlsoNearCells == true && CountAlsoVoidSpecies == true) {
                    for (int j = 0; j < J; j++)
                        ret[j] = rg.GetNoOfSpecies(j);

                } else {
                    throw new NotImplementedException();
                }


                for (int j = 0; j < J; j++) {



                    if (ConsiderAlsoNearCells) {
                        ret[j] = rg.GetNoOfSpecies(j);
                    } else {
                        if (MinDistLevSet(j) == 0)
                            ret[j] = rg.GetNoOfSpecies(j);
                        else
                            ret[j] = 1;
                    }

                    bool j_hasVoid;
                    if (CountAlsoVoidSpecies == false) {
                        j_hasVoid = rg.IsSpeciesPresentInCell(voidId, j);
                    } else {
                        j_hasVoid = false;
                    }

                }
            } else {
                ret.SetAll(1);
            }

            return ret;
        }

      


    }


   

}