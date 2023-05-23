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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution {
    public class AMRLevelIndicatorLibrary {

        /// <summary>
        /// refinement on cells which are on the specified boundary
        /// </summary>
        [Serializable]
        public class AMRonBoundary : AMRLevelIndicator {

            
            private byte[] m_EdgeTags;
            
            public AMRonBoundary(params byte[] EdgeTags) {
                m_EdgeTags = EdgeTags;
            }

            public override int[] DesiredCellChanges() {

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];
                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for (int j = 0; j < J; j++) {
                    int[] edges = GridData.Cells.Cells2Edges[j];
                    int currentLevel = cells[j].RefinementLevel;

                    bool refine = false;
                    foreach (int edge in edges) {
                        if (m_EdgeTags.Contains(GridData.Edges.EdgeTags[Math.Abs(edge) - 1]))
                            refine = true;
                    }

                    if (refine && currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if (!refine && currentLevel > 0) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }

                return levels;
            }

            public override bool Equals(object obj) {
                if (!base.Equals(obj))
                    return false;
                if (!m_EdgeTags.SetEquals((obj as AMRonBoundary)?.m_EdgeTags))
                    return false;
                return true;
            }

            public override int GetHashCode() {
                return base.GetHashCode();
            }
        }

        /// <summary>
        /// refinement on cells which are on the specified boundary
        /// </summary>
        [Serializable]
        public class AMRInBoundingBox : AMRLevelIndicator {

            private BoundingBox bb;
            public AMRInBoundingBox(BoundingBox _bb) {
                bb = _bb;
            }

            public override int[] DesiredCellChanges() {

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];
                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for (int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;


                    bool refine = false;
                    foreach (var cell in cells) {
                        if (bb.Contains(GridData.Cells.CellCenter.ExtractSubArrayShallow(j, -1).To1DArray()))
                            refine = true;
                    }

                    if (refine && currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if (!refine && currentLevel > 0) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }

                return levels;
            }

            public override bool Equals(object obj) {
                if (!base.Equals(obj))
                    return false;
                if (!bb.Equals((obj as AMRInBoundingBox)?.bb))
                    return false;
                return true;
            }

            public override int GetHashCode() {
                return base.GetHashCode();
            }
        }

        /// <summary>
        /// refinement everywhere
        /// </summary>
        [Serializable]
        public class AMReveryWhere : AMRLevelIndicator {
            public AMReveryWhere() {
            }

            public override int[] DesiredCellChanges() {

                int J = GridData.CellPartitioning.LocalLength;
                int[] levels = new int[J];
                int cellsToRefine = 0;
                int cellsToCoarse = 0;
                Cell[] cells = GridData.Grid.Cells;
                for (int j = 0; j < J; j++) {
                    int currentLevel = cells[j].RefinementLevel;

                    if (currentLevel < maxRefinementLevel) {
                        levels[j] = 1;
                        cellsToRefine++;
                    } else if (currentLevel > maxRefinementLevel) {
                        levels[j] = -1;
                        cellsToCoarse++;
                    }
                }

                return levels;
            }
        }
    }
}
