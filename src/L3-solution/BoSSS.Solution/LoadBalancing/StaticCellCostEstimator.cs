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

using ilPSP;
using System;
using System.Linq;

namespace BoSSS.Solution.LoadBalancing {

    /// <summary>
    /// Cell cost estimate with a fixed mapping between a performance class and
    /// and cost value
    /// </summary>
    [Serializable]
    public class StaticCellCostEstimator : CellTypeBasedEstimator {

        /// <summary>
        /// Serialization constructor
        /// </summary>
        private StaticCellCostEstimator() {
            
        }


        private int[] performanceClassToCostMap;

        [NonSerialized]
        private int[] cellToCostMap;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="performanceClassToCostMap">
        /// The fixed mapping between performance classes and their costs
        /// </param>
        public StaticCellCostEstimator(int[] performanceClassToCostMap) {
            this.performanceClassToCostMap = performanceClassToCostMap;
        }

        /// <summary>
        /// <see cref="ICellCostEstimator.UpdateEstimates"/>
        /// </summary>
        override public void UpdateEstimates(IApplication app) {
            int J = m_app.GridData.CellPartitioning.LocalLength;

            var cellToPerformanceClassMap = base.CellClassifier.ClassifyCells(app);

            cellToCostMap = new int[J];
            for (int j = 0; j < J; j++) {
                int performanceClass = cellToPerformanceClassMap[j];
                cellToCostMap[j] = performanceClassToCostMap[performanceClass];
            }

            
        }

        /// <summary>
        /// <see cref="ICellCostEstimator.GetEstimatedCellCosts"/>
        /// </summary>
        override public int[][] GetEstimatedCellCosts() {
            return new int[][] { cellToCostMap };
        }
               

        public override object Clone() {
            return new StaticCellCostEstimator() {
                CellClassifier = this.CellClassifier.CloneAs(),
                performanceClassToCostMap = this.performanceClassToCostMap.CloneAs()
            };
        }

    }
}
