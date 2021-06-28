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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// refinement on cells which are inside the narrow band 
    /// (cut-cells and neighboring cells sharing at least one point)
    /// </summary>
    [Serializable]
    public class AMRonNarrowband : AMRLevelIndicatorWithLevelset {


        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            CellMask band = this.LsTrk.Regions.GetNearFieldMask(1);

            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                if (band.Contains(j) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!band.Contains(j) && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }
    }
}
