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
            return cellToCostMaps;
        }

        /// <summary>
        /// returns multiple sets to cell-weights for a multi-constrained optimization
        /// </summary>
        public override void UpdateEstimates(IApplication app) {
            int J = app.GridData.CellPartitioning.LocalLength;
            var _cellToCostMapS = new List<int[]>();

            var cellToPerformanceClassMap = base.CellClassifier.ClassifyCells(app);

            //int NoOfAddtitionalConstraints = 

            var blabla = new[] {
                (CutStateClassifier.CellTypeFlags.Ordinary, 10),
                (CutStateClassifier.CellTypeFlags.Cut, (int)Math.Pow(2, 10))
            };
            int NoOfAddtitionalConstraints = blabla.Length;

            foreach(var tt in blabla) {
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
