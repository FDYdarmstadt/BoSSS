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
    /// Fluid cells are "0", void cells are "1"
    /// </summary>
    public class IBMCellClassifierTwo : ICellClassifier {

        /// <summary>
        /// Fluid cells are "0", void cells are "1"
        /// </summary>
        /// <param name="program"></param>
        /// <returns></returns>
        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            ImmersedSpeciesMap speciesMap = program.SpeciesMap as ImmersedSpeciesMap;
            IBMControl ibmControl = program.Control as IBMControl;
            if (speciesMap == null || ibmControl == null) {
                throw new Exception("IBM classifier only valid for IBM runs");
            }

            // Fluid and cut cells are "0"
            int[] cellToPerformanceClassMap = new int[program.Grid.NoOfUpdateCells];

            // Void cells are "1"
            foreach (int j in speciesMap.Tracker.Regions.GetSpeciesMask(ibmControl.VoidSpeciesName).ItemEnum) {
                cellToPerformanceClassMap[j] = 1;
            }

            int noOfClasses = 2;
            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}
