using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature {

    

    /// <summary>
    /// integration metric for $D-1$-dimensional integral over the level-set surface in each cell, i.e., 
    /// \[
    ///    \oint_{K_j \cap \mathfrak{I}}  f \mathrm{dS} .
    /// \]
    /// </summary>
    public class LevelSetIntegrationMetric : IIntegrationMetric {


        public LevelSetIntegrationMetric(LevelSetTracker.LevelSetData levelSetData) { 
            m_levelSetData = levelSetData;
        }

        readonly LevelSetTracker.LevelSetData m_levelSetData;

        /// <summary>
        /// For the integral over the level-set-surface, the metric depends on the node (if the element transformation applies some shearing transformation),
        /// therefore always true
        /// </summary>
        public bool AlwaysUsePerNodeScaling => true;


        /// <summary>
        /// not used, because of reasons given for <see cref="AlwaysUsePerNodeScaling"/>
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int jCell, int L) {
            throw new NotSupportedException("shall not be called");
            
        }

        /// <summary>
        ///
        /// </summary>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int jCell0, int L) {
            if(!object.ReferenceEquals(gridData, m_levelSetData.GridDat))
                throw new ArgumentException();
            var metrics = m_levelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(qr.Nodes, jCell0, L);

            int D = gridData.SpatialDimension;

            int NoOfNodes = qr.NoOfNodes;
            Debug.Assert(metrics.Dimension == 2);
            Debug.Assert(metrics.GetLength(0) == L);
            Debug.Assert(metrics.GetLength(1) == qr.NoOfNodes);

            if(gridData.iGeomCells.IsCellAffineLinear(jCell0)) {
                //
                // note: implemented in the refactoring to reproduce the original implementation,
                // without any mathematical consideration put in place
                //
                
                var cellJacDet = gridData.iGeomCells.JacobiDet;

                for(int i = 0; i < L; i++) {
                    for(int k = 0; k < NoOfNodes; k++) {
                        metrics[i, k] = cellJacDet[i + jCell0] / metrics[i, k];

                        //if(D == 3)
                        //    Console.WriteLine(metrics[i, k] + " " + (gridData as GridData).JacobianDeterminat.GetValue_Cell(qr.Nodes,jCell0 + i)[0]);
                    }
                }


            } else {
                throw new NotSupportedException("todo");
            }


            return metrics;
            //return gridData.JacobianDeterminat.GetValue_Cell(qr.Nodes, jCell0, L);
        }
    }
}
