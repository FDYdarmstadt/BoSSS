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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// refinement on cells which contain the Contact-Line 
    /// only taking boundary edges into account, for contactlines between multiple levelsets <see cref="AMRatInnerContactLine"/>
    /// </summary>
    [Serializable]
    public class AMRatContactLine : AMRLevelIndicatorWithLevelset {

        public int[] levelSets; // null all level set intersections, otherwise only between those specified
        public byte[] edgeTags; // null all boundaries, otherwise those specified
        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            CellMask ccMask = CellMask.GetEmptyMask(GridData);
            // cut cell masks for level sets
            if(levelSets != null) {
                foreach(int ls in levelSets) {
                    ccMask = ccMask.Union(this.LsTrk.Regions.GetCutCellMask4LevSet(ls));
                }
            } else {
                ccMask = ccMask.Union(this.LsTrk.Regions.GetCutCellMask());
            }
            
            CellMask bMask = CellMask.GetEmptyMask(GridData);
            // boundary masks
            if (edgeTags != null) {
                foreach (byte tag in edgeTags) {

                    BitArray boundaryCells = new BitArray(J);

                    var Edges = GridData.BoundaryEdges;
                    int E = Edges.NoOfItemsLocally;
                    // loop over all Edges
                    foreach(Chunk ch in Edges) {
                        foreach(int e in ch.Elements) {
                            int Cel = GridData.Edges.CellIndices[e, 0];
                            if (edgeTags.Contains(GridData.Edges.EdgeTags[e])) {
                                boundaryCells[Cel] = true;
                            }
                        }
                    }
                    bMask = bMask.Union(new CellMask(GridData, boundaryCells));
                }
            } else {
                bMask = bMask.Union(GridData.BoundaryCells.VolumeMask);
            }

            // get the contactline cell
            CellMask clMask = bMask.Intersect(ccMask);

            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                if (clMask.Contains(j) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!clMask.Contains(j) && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }
    }

    /// <summary>
    /// AMR at levelset/levelset intersections (inner contactlines)
    /// </summary>
    public class AMRatInnerContactLine : AMRLevelIndicatorWithLevelset {

        public int[] levelSets; // null all level set intersections, otherwise only between those specified.
                                // If only one is selected, intersection between this one and other levelsets are taken into account,
                                // but intersections not involving this are not
        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            CellMask clMask = CellMask.GetEmptyMask(GridData);

            if (levelSets == null) {
                levelSets = Enumerable.Range(0, this.LsTrk.LevelSets.Count).ToArray();
            }

            if (levelSets.Length > 1) {
                for (int i = 0; i < levelSets.Length; i++) {
                    CellMask ccMask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSets.First());
                    foreach (int ls in levelSets.Skip(i + 1)) {
                        clMask = clMask.Union(ccMask.Intersect(this.LsTrk.Regions.GetCutCellMask4LevSet(ls)));
                    }
                }
            } else {
                // only one level set specified, refine all intersections, this one has with any other
                CellMask ccMask = this.LsTrk.Regions.GetCutCellMask4LevSet(levelSets[0]);
                levelSets = Enumerable.Range(0, this.LsTrk.LevelSets.Count).Except(levelSets).ToArray();
                foreach (int ls in levelSets) {
                    clMask = clMask.Union(ccMask.Intersect(this.LsTrk.Regions.GetCutCellMask4LevSet(ls)));
                }
            }


            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                if (clMask.Contains(j) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!clMask.Contains(j) && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }
    }
}
