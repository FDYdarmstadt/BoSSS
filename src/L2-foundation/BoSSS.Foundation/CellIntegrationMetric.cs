using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Foundation.Quadrature {


    /// <summary>
    /// integration metric for volume integrals
    /// </summary>
    public class CellIntegrationMetric : IIntegrationMetric {



        /// <summary>
        /// the absolute value of the Jacobian determinate, for the 
        /// transformation that translates cell reference coordinates into physical coordinates.
        /// </summary>
        public MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int i0, int L) {
            return gridData.iGeomCells.JacobiDet.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + L - 1 });
        }

        /// <summary>
        /// the absolute value of the Jacobian determinate, for the 
        /// transformation that translates cell reference coordinates into physical coordinates.
        /// </summary>
        public MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int jCell0, int L) {
            return gridData.JacobianDeterminat.GetValue_Cell(qr.Nodes, jCell0, L);
        }
    }
}
