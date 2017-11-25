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

using BoSSS.Solution;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS.LoadBalancing {

    public class LTSCellCostEstimator : ICellCostEstimator {

        private int[] cellToCostMap;

        private int selectedPerformanceClass;

        public LTSCellCostEstimator(int selectedPerformanceClass) {
            this.selectedPerformanceClass = selectedPerformanceClass;
        }

        public double EstimatedLocalCost {
            get;
            private set;
        }

        public int CurrentPerformanceClassCount {
            get;
            private set;
        }

        public int[] GetEstimatedCellCosts() {
            return cellToCostMap;
        }

        public void UpdateEstimates(int performanceClassCount, int[] cellToPerformanceClassMap) {
            CurrentPerformanceClassCount = performanceClassCount;

            // One balance constraint per cluster
            cellToCostMap = new int[cellToPerformanceClassMap.Length];
            cellToCostMap.SetAll(1);
            for (int j = 0; j < cellToPerformanceClassMap.Length; j++) {
                if (cellToPerformanceClassMap[j] == selectedPerformanceClass) {
                    cellToCostMap[j] = 10;
                }
            }

            EstimatedLocalCost = cellToCostMap.Sum();
        }

        public static IEnumerable<Func<IApplication<AppControl>, int, ICellCostEstimator>> Factory(int numberOfClusters) {
            for (int i = 0; i < numberOfClusters; i++) {
                int temp = i; // Avoid delegate creation from capturing variable $i
                yield return (app, classCount) => new LTSCellCostEstimator(temp);
            }
        }
    }
}
