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

using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Solution {

    /// <summary>
    /// Cell cost estimate with a fixed mapping between a performance class and
    /// and cost value
    /// </summary>
    public class StaticCellCostEstimator : ICellCostEstimator {

        /// <summary>
        /// <see cref="ICellCostEstimator"/>
        /// </summary>
        public int PerformanceClassCount {
            get;
            private set;
        }

        /// <summary>
        /// <see cref="ICellCostEstimator"/>
        /// </summary>
        public double EstimatedLocalCost {
            get;
            private set;
        }

        private int[] performanceClassToCostMap;

        private int[] cellToCostMap;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="performanceClassToCostMap">
        /// The fixed mapping between performance classes and their costs
        /// </param>
        public StaticCellCostEstimator(int[] performanceClassToCostMap) {
            this.performanceClassToCostMap = performanceClassToCostMap;
            this.PerformanceClassCount = performanceClassToCostMap.Length;
        }

        /// <summary>
        /// <see cref="ICellCostEstimator"/>
        /// </summary>
        /// <p<param name="performanceClassCount"></param>
        /// <param name="cellToPerformanceClassMap"></param>
        public void UpdateEstimates(int performanceClassCount, int[] cellToPerformanceClassMap) {
            cellToCostMap = new int[cellToPerformanceClassMap.Length];
            for (int j = 0; j < cellToPerformanceClassMap.Length; j++) {
                int performanceClass = cellToPerformanceClassMap[j];
                cellToCostMap[j] = performanceClassToCostMap[performanceClass];
            }

            EstimatedLocalCost = cellToCostMap.Sum();
        }

        /// <summary>
        /// <see cref="ICellCostEstimator"/>
        /// </summary>
        /// <returns></returns>
        public int[] GetEstimatedCellCosts() {
            return cellToCostMap;
        }
    }
}
