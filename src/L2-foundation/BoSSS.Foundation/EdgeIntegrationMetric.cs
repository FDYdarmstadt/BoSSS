using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Foundation.Quadrature {


    /// <summary>
    /// integration metric for area integrals on edges
    /// </summary>
    public class EdgeIntegrationMetric : IIntegrationMetric {

        /// <summary>
        /// For edges, the scaling is constant within the edge for linear elements, therefore always false;
        /// c.f. <see cref="IIntegrationMetric.AlwaysUsePerNodeScaling"/>;
        /// </summary>
        public bool AlwaysUsePerNodeScaling => false;


        /// <summary>
        /// the square root of the Gramian determinant
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int i0, int L) {
            var R = gridData.iGeomEdges.SqrtGramian.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + L - 1 });
            return R;
        }
        /// <summary>
        /// the transformation metric for nonlinear edges.
        /// </summary>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int i0, int L) {
            return gridData.iGeomEdges.NormalsCache.GetIntegrationMetric(qr.Nodes, i0, L);
        }
    }
}
