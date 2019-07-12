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

using CNS.IBM;
using System;

namespace CNS.LoadBalancing {

    /// <summary>
    /// Standard cells are "0", cut cells are "1", void cells are "2"
    /// </summary>
    public class IBMCellClassifier : ICellClassifier {

        /// <summary>
        /// Standard cells are "0", cut cells are "1", void cells are "2"
        /// </summary>
        /// <param name="program"></param>
        /// <returns></returns>
        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
            IBMControl ibmControl = program.Control as IBMControl;
            if (speciesMap == null || ibmControl == null) {
                throw new Exception("IBM classifier only valid for IBM runs");
            }

            // Pure fluid cells are "0"
            int J = program.gridData.iLogicalCells.NoOfLocalUpdatedCells;
            int[] cellToPerformanceClassMap = new int[J];

            // Cut cells are "1"
            foreach (int j in speciesMap.Tracker.Regions.GetCutCellMask().ItemEnum) {
                cellToPerformanceClassMap[j] = 1;
            }

            // Void cells are "2"
            foreach (int j in speciesMap.Tracker.Regions.GetSpeciesMask(ibmControl.VoidSpeciesName).ItemEnum) {
                cellToPerformanceClassMap[j] = 2;
            }

            int noOfClasses = 3;
            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}
