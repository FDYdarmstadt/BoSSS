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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using System;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// refinement of cells using predefined limits of a variable
    /// </summary>
    [Serializable]
    public class AMR_BasedOnVariableLimits : AMRLevelIndicatorWithLevelset {

        [DataMember]
        private string FieldName;
        [DataMember]
        private double[] VariableLimits;

        /// <summary>
        /// Empty constructor for serialization
        /// </summary> 
        private AMR_BasedOnVariableLimits() { }


        public AMR_BasedOnVariableLimits(string _FieldName, double[] _VariableLimits, int _maxRefinementLevel) {
            FieldName = _FieldName;
            if(_VariableLimits.Length != 2)
                throw new Exception("A minimum and maximum value is needed");
            VariableLimits = _VariableLimits;
            maxRefinementLevel = _maxRefinementLevel;
        }

        protected DGField FieldForRefinement {
            get {              
                return ((XNSEC)this.SolverMain).CurrentStateVector.Fields.Where(f => f.Identification == FieldName).SingleOrDefault();
            }
        }

        public override int[] DesiredCellChanges() {
            int J = GridData.CellPartitioning.LocalLength;
            int[] levels = new int[J];

            var field = this.FieldForRefinement;
            if(field == null)
                return levels;

            double eps = 1e-1 * 0;
            double maxAcceptedValue = VariableLimits[1] + eps;
            double minAcceptedValue = VariableLimits[0] - eps;
            Cell[] cells = GridData.Grid.Cells;
            for(int j = 0; j < J; j++) {
                double localMinVal, localMaxVal;
                field.GetExtremalValuesInCell(out localMinVal, out localMaxVal, j);

                int currentLevel = cells[j].RefinementLevel;
                if(localMaxVal > maxAcceptedValue && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                }

                if(localMinVal < minAcceptedValue && currentLevel < maxRefinementLevel) {
                    levels[j] = 1;
                }
            }
            Console.WriteLine("Number of refined cells  where temperature is higher than 5: " + levels.Sum() + ".\n");
            return levels;
        }
    }
}