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

using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;


namespace BoSSS.Application.XNSERO_Solver {

    /// <summary>
    /// refinement on cells which are inside a test range around the rigid object
    /// Checks for the veticies of the cells in the narrwo band
    /// </summary>
    [DataContract]
    [Serializable]
    public class AMRaroundRigidObject : AMRLevelIndicatorWithLevelset {

        [DataMember]
        List<Particle> ParticlesToCheck;

        /// <summary>
        /// if zero, contains returns only true if cell is inside the rigid object
        /// </summary>
        [DataMember]
        double TestRange;   

        public AMRaroundRigidObject(List<Particle> ParticlesToAdapt, double Range = 0.0) {
            ParticlesToCheck = ParticlesToAdapt;
            TestRange = Range;
        }


        bool IsInRange(int jCell) {
            int[] jCellVerticies = this.GridData.Cells.CellVertices[jCell];
            foreach (int vertId in jCellVerticies) {
                double[] vCoord = this.GridData.Vertices.Coordinates.ExtractSubArrayShallow(new int[] { vertId, -1}).To1DArray();
                foreach (Particle p in ParticlesToCheck) {
                    if (p.Contains(new Vector(vCoord), TestRange) && !p.Contains(new Vector(vCoord), 0.0)) {
                        return true;
                    }
                }
            }
            return false;
        }


        public override int[] DesiredCellChanges() {

            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            //CellMask CutCells = this.LsTrk.Regions.GetNearMask4LevSet(1, 1); // narrow band for rigid object level set

            int cellsToRefine = 0;
            int cellsToCoarse = 0;
            Cell[] cells = GridData.Grid.Cells;
            for (int j = 0; j < J; j++) {
                int currentLevel = cells[j].RefinementLevel;

                if (IsInRange(j) && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                    cellsToRefine++;
                } else if (!IsInRange(j) && currentLevel > 0) {
                    levels[j] = -1;
                    cellsToCoarse++;
                }
            }

            return levels;
        }

        public override bool Equals(object obj) {
            if (!base.Equals(obj))
                return false;
            var other = obj as AMRaroundRigidObject;
            if (other == null)
                return false;
            return true;
        }

        public override int GetHashCode() {
            return base.GetHashCode();
        }
    }
}
