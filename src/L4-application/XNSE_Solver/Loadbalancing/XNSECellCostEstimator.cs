using BoSSS.Solution;
using BoSSS.Solution.LoadBalancing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Loadbalancing {

    /// <summary>
    /// The default cell cost estimator for XNSE;
    /// The primary purpose of this class is to assign multiple weights/costs to each cell, so that 
    /// <see cref="LoadBalancer.GetNewPartitioning"/>
    /// can compute a cell partitioning (i.e., which cell is assigned to which MPI process).
    /// 
    /// Why multiple weights?
    /// The general idea is to have multiple classes (aka. clusters) of cells.
    /// E.g., cluster 0 are ordinary, un-cut cells and cluster 1 are cut-cells.
    /// Then, <see cref="LoadBalancer.GetNewPartitioning"/> tries not to balance the total weight (sum over all weights) 
    /// across the MPI processors, but it also tries to balance the weight in each cluster
    /// (i.e., the sum of all weights over all cells, **for each cluster**, is roughly the same for each MPI process.)
    /// </summary>
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
        /// returns the weight for the individual cell type
        /// </summary>
        public int FindWeightFor(CutStateClassifier.CellTypeFlags CellType) {
            return TypeToWgt.Where(p => p.CellType == CellType).First().Weight;
        }

        (CutStateClassifier.CellTypeFlags CellType, int Weight)[] typeToWgt;

        /// <summary>
        /// returns tuple array for cell-weights for a multi-constrained optimization
        /// 
        /// Currently only treats 
        /// </summary>
        public (CutStateClassifier.CellTypeFlags CellType, int Weight)[] TypeToWgt {
            get {
                typeToWgt = typeToWgt ?? new[] {
                //(CutStateClassifier.CellTypeFlags.Void, 0), As void does not impact weights, it is unnecessary to include it in calculations.
                (CutStateClassifier.CellTypeFlags.Ordinary, 10),
                (CutStateClassifier.CellTypeFlags.Cut, (int)Math.Pow(2, 10))
            };
                return typeToWgt;
            }
            set { typeToWgt = value; }
        }

        /// <summary>
        /// returns multiple sets to cell-weights for a multi-constrained optimization
        /// </summary>
        public override void UpdateEstimates(IApplication app) {
            int J = app.GridData.CellPartitioning.LocalLength;
            var _cellToCostMapS = new List<int[]>();


            var cellToPerformanceClassMap = base.CellClassifier.ClassifyCells(app);


                foreach(var tt in TypeToWgt) {
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
