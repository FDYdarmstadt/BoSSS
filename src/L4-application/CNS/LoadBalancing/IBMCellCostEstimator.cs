using BoSSS.Solution;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNS.LoadBalancing {

    public class IBMCellCostEstimator : ICellCostEstimator {

        private int[] cellToCostMap;

        private int selectedCellType;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="selectedCellType">
        /// See <see cref="IBMCellClassifier"/>:
        /// 0: Fluid
        /// 1: Cut
        /// 2: Void
        /// </param>
        public IBMCellCostEstimator(int selectedCellType) {
            this.selectedCellType = selectedCellType;
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
                if (cellToPerformanceClassMap[j] == selectedCellType) {
                    cellToCostMap[j] = 10;
                }
            }

            EstimatedLocalCost = cellToCostMap.Sum();
        }

        public static Func<IApplication<AppControl>, int, ICellCostEstimator> GetStaticCostBasedEstimator() {
            return (p, i) => new StaticCellCostEstimator(
                new int[] { 10, 100, 1 });
        }

        public static IEnumerable<Func<IApplication<AppControl>, int, ICellCostEstimator>> GetMultiBalanceConstraintedBasedEstimators() {
            int noOfCellTypes = 3; // Fluid + Cut + Void
            for (int i = 0; i < noOfCellTypes; i++) {
                int temp = i; // Avoid delegate creation from capturing variable $i
                yield return (app, classCount) => new IBMCellCostEstimator(temp);
            }
        }
    }
}
