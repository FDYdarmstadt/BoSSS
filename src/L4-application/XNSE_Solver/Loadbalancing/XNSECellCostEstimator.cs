using BoSSS.Solution;
using BoSSS.Solution.LoadBalancing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Loadbalancing {


    [Serializable]
    public class XNSECellCostEstimator : CellTypeBasedEstimator {

        [NonSerialized]
        private int[][] cellToCostMaps;


        /// <summary>
        /// 
        /// </summary>
        public XNSECellCostEstimator() {
            base.CellClassifier = new CutStateClassifier();



            
        }

        

        public override int[][] GetEstimatedCellCosts() {
            if (cellToCostMaps == null)
                UpdateEstimates(m_app);

            return cellToCostMaps;
        }

        /// <summary>
        /// returns multiple sets to cell-weights for a multi-constrained optimization
        /// </summary>
        public override void UpdateEstimates(IApplication app) {
            int J = app.GridData.CellPartitioning.LocalLength;
            var _cellToCostMapS = new List<int[]>();


            var cellToPerformanceClassMap = base.CellClassifier.ClassifyCells(app);

            var typeToWgt = new[] {
                //(CutStateClassifier.CellTypeFlags.Void, 0), As void does not impact weights, it is unnecessary to include it in calculations.
                (CutStateClassifier.CellTypeFlags.Ordinary, 10),
                (CutStateClassifier.CellTypeFlags.Cut, (int)Math.Pow(2, 10))
            };

                foreach(var tt in typeToWgt) {
                    // One balance constraint per cluster
                    var cellToCostMap = new int[J];
                    _cellToCostMapS.Add(cellToCostMap);
                    int cellFlag = (int) tt.Item1;
                    int wgt = tt.Item2;

                    for (int j = 0; j < J; j++) { // loop over cells...
                        if ((cellToPerformanceClassMap[j] & cellFlag) != 0) {
                            cellToCostMap[j] = wgt;
                    }
                }
            }

            cellToCostMaps = _cellToCostMapS.ToArray();
        }

        //public static IEnumerable<ICellCostEstimator> Factory() {
        //    int noOfCellTypes = 3; // Fluid + Cut + Void
        //    for (int i = 0; i < noOfCellTypes; i++) {
        //        int temp = i; // Avoid delegate creation from capturing variable $i
        //        yield return new XNSECellCostEstimator(temp);
        //    }
        //}

        public override object Clone() {
            throw new NotImplementedException();
        }

    }
}
