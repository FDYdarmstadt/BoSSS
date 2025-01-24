using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace BoSSS.Foundation.XDG.Quadrature {

    /// <summary>
    /// integration metric for $D-2$-dimensional integral the boundary of surface elements for each edge, i.e., for an edge $e_i$
    /// \[
    ///    \int_{e_j \cap \mathfrak{I}}  f \mathrm{dl} .
    /// \]
    /// </summary>
    internal class SurfaceElementEdgeIntegrationMetric : IIntegrationMetric {


        public SurfaceElementEdgeIntegrationMetric(LevelSetTracker.LevelSetData levelSetData) {
            m_levelSetData = levelSetData;
        }

        readonly LevelSetTracker.LevelSetData m_levelSetData;

        /// <summary>
        /// For the integral over the surface element boundary, 
        /// the metric depends on the node (if the element transformation applies some shearing transformation),
        /// but only in 3D
        /// </summary>
        public bool AlwaysUsePerNodeScaling {
            get {
                if(m_levelSetData.GridDat.SpatialDimension == 2)
                    return false;
                if(m_levelSetData.GridDat.SpatialDimension == 3)
                    return true;
                throw new NotSupportedException();
            }
        }


        /// <summary>
        /// in 2D, the boundary of surface elements are just points -- one needs a count measure, where the integration metric is always 1
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int jCell, int L) {
            if(m_levelSetData.GridDat.SpatialDimension == 2) {
                var ret = MultidimensionalArray.Create(L);
                ret.SetAll(1.0);
                return ret;
            } else { 
                throw new NotSupportedException();
            }

        }

        /// <summary>
        ///
        /// </summary>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int jCell0, int L) {
            if(!object.ReferenceEquals(gridData, m_levelSetData.GridDat))
                throw new ArgumentException();

            throw new NotImplementedException("todo");


            //var metrics = m_levelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(qr.Nodes, jCell0, L);

            //int NoOfNodes = qr.NoOfNodes;
            //Debug.Assert(metrics.Dimension == 2);
            //Debug.Assert(metrics.GetLength(0) == L);
            //Debug.Assert(metrics.GetLength(1) == qr.NoOfNodes);

            //if(gridData.iGeomCells.IsCellAffineLinear(jCell0)) {
            //    //
            //    // note: implemented in the refactoring to reproduce the original implementation,
            //    // without any mathematical consideration put in place
            //    //

            //    var cellJacDet = gridData.iGeomCells.JacobiDet;

            //    for(int i = 0; i < L; i++) {
            //        for(int k = 0; k < NoOfNodes; k++) {
            //            metrics[i, k] = cellJacDet[i + jCell0] / metrics[i, k];
            //        }
            //    }


            //} else {
            //    throw new NotSupportedException("todo");
            //}


            //return metrics;
            //return gridData.JacobianDeterminat.GetValue_Cell(qr.Nodes, jCell0, L);
        }
    }
}


