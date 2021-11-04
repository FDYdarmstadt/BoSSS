using BoSSS.Solution;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Loadbalancing {
    public class XNSECellCostEstimator : ICellCostEstimator {
        private int[] cellToCostMap;

        private int selectedCellType;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="selectedCellType">
        /// See <see cref="IBMCellClassifier"/>:
        /// - 0: Void (Solid)
        /// - 1: Fluid
        /// - 2: cut cells
        /// </param>
        public XNSECellCostEstimator(int selectedCellType) {
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
                    if (selectedCellType == 2)
                        cellToCostMap[j] = (int)Math.Pow((double)selectedCellType, 10) - selectedCellType;
                    else
                        cellToCostMap[j] = 10;
                }
            }

            EstimatedLocalCost = cellToCostMap.Sum();
        }

        public static IEnumerable<Func<IApplication, int, ICellCostEstimator>> Factory() {
            int noOfCellTypes = 3; // Fluid + Cut + Void
            for (int i = 0; i < noOfCellTypes; i++) {
                int temp = i; // Avoid delegate creation from capturing variable $i
                yield return (app, classCount) => new XNSECellCostEstimator(temp);
            }
        }

    }
}
