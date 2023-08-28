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
using ilPSP;
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

        public int levelSet = -1; // level set this Indicator should be active on
        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            CellMask band;
            if (levelSet == -1) {
                band = this.LsTrk.Regions.GetNearFieldMask(1);
            } else {
                band = this.LsTrk.Regions.GetNearMask4LevSet(levelSet, 1);
            }

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


        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            var other = obj as AMRonNarrowband;
            if (other == null)
                return false;
            if (other.levelSet != this.levelSet)
                return false;
            return true;
        }

        public override int GetHashCode() {
            return base.GetHashCode();
        }
    }

    public class AMRForRigidObject : AMRLevelIndicatorWithLevelset {

        public AMRForRigidObject(List<Func<Vector, double, bool>> ContainsFunction, double LengthScale) {
            AllContainsFunctions = ContainsFunction;
            this.LengthScale = LengthScale;
        }

        private readonly List<Func<Vector, double, bool>> AllContainsFunctions;
        private readonly double LengthScale;

        bool ContainsFunction(double[] X) {
            bool containsFunction = false;
            int i = 0;
            while(!containsFunction && i < AllContainsFunctions.Count()) {
                containsFunction = AllContainsFunctions[i](X, LengthScale);
                i += 1;
            }
            return containsFunction;
        }

        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;
                Vector cellCenter = new Vector(GridData.iGeomCells.GetCenter(j));
                if (ContainsFunction(cellCenter) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!ContainsFunction(cellCenter) && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }
    }
}
