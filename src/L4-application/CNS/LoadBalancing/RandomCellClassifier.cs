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

namespace CNS.LoadBalancing {

    /// <summary>
    /// Random classification of cells for testing purposes
    /// </summary>
    public class RandomCellClassifier : ICellClassifier {

        private int noOfClasses;

        public RandomCellClassifier(int noOfClasses) {
            this.noOfClasses = noOfClasses;
        }

        public (int noOfClasses, int[] cellToPerformanceClassMap) ClassifyCells(IProgram<CNSControl> program) {
            // Use seed to make runs deterministic, without having the same
            // random numbers for all cores
            Random rand = new Random(program.ResLogger.TimeStep + program.MPIRank);

            int[] cellToPerformanceClassMap = new int[program.GridData.iLogicalCells.NoOfLocalUpdatedCells];
            for (int i = 0; i < cellToPerformanceClassMap.Length; i++) {
                cellToPerformanceClassMap[i] = rand.Next(0, noOfClasses);
            }

            return (noOfClasses, cellToPerformanceClassMap);
        }
    }
}
