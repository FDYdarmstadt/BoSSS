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
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Configuration;
using System.Diagnostics;
using System.Linq;

namespace CNS.LoadBalancing {

    /// <summary>
    /// A cell cost estimator based on aritificial viscosity
    /// </summary>
    public class ArtificialViscosityCellCostEstimator : ICellCostEstimator {

        private int[] cellToCostMap;

        private int performanceClass;

        /// <summary>
        /// <see cref="ICellCostEstimator"/>
        /// </summary>
        public int CurrentPerformanceClassCount {
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

        /// <summary>
        /// A cell cost estimator based on aritificial viscosity
        /// Performance class 0: Cells without AV
        /// Performance class 1: Cells with AV
        /// </summary>
        public ArtificialViscosityCellCostEstimator(int performanceClass) {
            this.performanceClass = performanceClass;
        }

        public void UpdateEstimates(int performanceClassCount, int[] cellToPerformanceClassMap) {
            CurrentPerformanceClassCount = performanceClassCount;

            // One balance constraint per cluster
            cellToCostMap = new int[cellToPerformanceClassMap.Length];
            cellToCostMap.SetAll(1);
            for (int j = 0; j < cellToPerformanceClassMap.Length; j++) {
                if (cellToPerformanceClassMap[j] == this.performanceClass) {
                    cellToCostMap[j] = 10;
                }
            }

            EstimatedLocalCost = cellToCostMap.Sum();
        }

        /// <summary>
        /// AV cells are ten times more expensive than non-AV cells
        /// </summary>
        /// <returns></returns>
        public static Func<IApplication<AppControl>, int, ICellCostEstimator> GetStaticCostBasedEstimator() {
            return (p, i) => new StaticCellCostEstimator(
                new int[] { 1, 10 });
        }

        /// <summary>
        /// Create two <see cref="ArtificialViscosityCellCostEstimator"/> for non-AV and AV cells
        /// with converse costs: (1, 10) and (10, 1), cell costs do not matter in this case
        /// (can be arbitrary)
        /// </summary>
        /// <returns></returns>
        public static IEnumerable<Func<IApplication<AppControl>, int, ICellCostEstimator>> GetMultiBalanceConstraintsBasedEstimators() {
            int noOfPerformancesClasses = 2; // Cells with AV + cells without AV
            for (int i = 0; i < noOfPerformancesClasses; i++) {
                int temp = i; // Avoid delegate creation from capturing variable $i
                yield return (app, classCount) => new ArtificialViscosityCellCostEstimator(temp);
            }
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