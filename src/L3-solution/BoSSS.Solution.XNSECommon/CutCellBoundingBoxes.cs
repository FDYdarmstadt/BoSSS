/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using ilPSP;
using BoSSS.Platform;
using ilPSP.Utils;

namespace BoSSS.Solution.XNSECommon {

    [Obsolete]
    public class CutCellBoundingBoxes {

        BoundingBox[][] BoundingBoxes;


        public CutCellBoundingBoxes(LevelSetTracker __lsTrk) {
            this.lsTrk = __lsTrk;
            int J = lsTrk.GridDat.Cells.Count;
            int D = lsTrk.GridDat.SpatialDimension;
            this.BoundingBoxes = new BoundingBox[J][];
            this.UpdateBoundingBoxes();
        }

        LevelSetTracker lsTrk;

        public BoundingBox GetSpeciesBB(SpeciesId spid, int jCell) {
            int spIdx = lsTrk.Regions.GetSpeciesIndex(spid, jCell);
            if(spIdx >= 0) {
                bool isPresent = lsTrk.Regions.IsSpeciesPresentInCell(spid, jCell);
                if(isPresent) {
                    int iii = this.BoundingBoxes[jCell].Length > 1 ? spIdx : 0;
                    BoundingBox BB = this.BoundingBoxes[jCell][iii];
                    return BB;
                } else {
                    return null;
                }
            } else {
                return null;
            }
        }

        public double Get_hminBB(SpeciesId spid, int jCell) {
            var BB = GetSpeciesBB(spid, jCell);
            if(BB == null) {
                return 0.0;
            } else {
                return BB.h_min;
            }
        }


        void UpdateBoundingBoxes() {

            var grd = lsTrk.GridDat;
            int D = grd.SpatialDimension;
            int J = grd.Cells.Count;
            var KrefS = grd.Grid.RefElements;
            var VolQRs = KrefS.Select(Kref => Kref.GetBruteForceQuadRule(6, 2)).ToArray();
            var BruteForceVolFam = VolQRs;


            // cache some vals
            SubGrid allCut = lsTrk.Regions.GetCutCellSubGrid();
            int JsubCC = allCut.LocalNoOfCells;
            int NoOfLevSets = lsTrk.LevelSets.Count();
            MultidimensionalArray[] LevSetVals = new MultidimensionalArray[NoOfLevSets];
            MultidimensionalArray[] GlobalNodes = BruteForceVolFam.Select(nscc => MultidimensionalArray.Create(nscc.Nodes.NoOfNodes, D)).ToArray();

            int[] m_CutCellSubGrid_SubgridIndex2LocalCellIndex = lsTrk.Regions.GetCutCellSubGrid().SubgridIndex2LocalCellIndex;

            // set BoundingBox'es for cut cells
            // ================================
            
            double[] pt = new double[D];
            for(int jsub = 0; jsub < JsubCC; jsub++) { // loop over all local cells in subgrid ...
                int j_cell = m_CutCellSubGrid_SubgridIndex2LocalCellIndex[jsub];
                int iKref = grd.Cells.GetRefElementIndex(j_cell);
                int NoOfQuadNodes = BruteForceVolFam[iKref].NoOfNodes;

                // loop over level sets (evaluation) ...
                for(int levSetIdx = 0; levSetIdx < NoOfLevSets; levSetIdx++) {
                    LevSetVals[levSetIdx] = lsTrk.DataHistories[levSetIdx].Current.GetLevSetValues(BruteForceVolFam[iKref].Nodes, j_cell, 1);
                }

                grd.TransformLocal2Global(BruteForceVolFam[iKref].Nodes, GlobalNodes[iKref], j_cell);
                                
                ReducedRegionCode rrc;
                int NoOfSpec = lsTrk.Regions.GetNoOfSpecies(j_cell, out rrc);

                this.BoundingBoxes[j_cell] = NoOfSpec.ForLoop(iSpc => new BoundingBox(D));
                                 
                for(int m = 0; m < NoOfQuadNodes; m++) {
                    GlobalNodes[iKref].GetRow(m, pt);

                    double val0 = LevSetVals[0][0, m];
                    double val1 = 0;
                    if(NoOfLevSets > 1)
                        val1 = LevSetVals[1][0, m];
                    double val2 = 0;
                    if(NoOfLevSets > 2)
                        val2 = LevSetVals[2][0, m];
                    double val3 = 0;
                    if(NoOfLevSets > 3)
                        val3 = LevSetVals[3][0, m];

                    LevelSetSignCode LevSetCode = LevelSetSignCode.ComputeLevelSetBytecode(val0, val1, val2, val3);
                    int iSpec = lsTrk.GetSpeciesIndex(rrc, LevSetCode);

                    this.BoundingBoxes[j_cell][iSpec].AddPoint(pt);
                }
            }

            // set Boundingboxes for un-cut cells
            // ==================================
            for(int j = 0; j < J; j++) {
                if(this.BoundingBoxes[j] == null) {
                    ReducedRegionCode rrc;
                    int NoOfSpc = lsTrk.Regions.GetNoOfSpecies(j, out rrc);

                    

                    BoundingBox BB = new BoundingBox(D);
                    grd.Cells.GetCellBoundingBox(j, BB);
                    this.BoundingBoxes[j] = new BoundingBox[] { BB };
                }
            }
            
        }

    }
}
