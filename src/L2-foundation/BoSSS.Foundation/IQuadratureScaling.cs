using BoSSS.Foundation.Grid;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Quadrature rules (<see cref="QuadRule"/>) are defined in the reference coordinate system of cells or edges.
    /// To obtain a quadrature rule in the physical domain, the integrand must be applied with an integration metric/scaling factor,
    /// see also https://en.wikipedia.org/wiki/Integral_transform 
    /// 
    /// This interface is required e.g., because not all quadrature rules which are defined in cell coordinates
    /// are necessarily also volume integrals, therefore the quadrature rule must provide its scaling factors
    /// for the quadrature kernel (<see cref="CellQuadrature"/>, <see cref="EdgeQuadrature"/>, <see cref="CellBoundaryQuadrature{TQuadRule}"/>).
    /// One example are integrals over embedded manifolds, e.g. level-set-surface integrals; Those are defined 
    /// in cell coordinates, but require different scaling factors than volume integrals.
    /// </summary>
    public interface IQuadratureScaling {

        /// <summary>
        /// Provides integration metric factors for linear edges or cells, i.e., 
        /// where the transformation from the reference to the physical domain is linear, resp. affine-linear.
        /// For such edges or cells, the metric factor is constant within a cell.
        /// </summary>
        /// <param name="i0">Index of first cell or edge for which to obtain the scaling factors.</param>
        /// <param name="L">Number of cells or edges for which to obtain the scaling factors.</param>
        /// <param name="qr">Quadrature rule in reference element for edges/cells <paramref name="i0"/> to <paramref name="i0"/>+<paramref name="L"/>-1</param>
        /// <param name="gridData">current grid</param>
        /// <returns>
        /// Scaling factors for each cell or edge
        /// - index: cell or edge index, from 0 (including) to <paramref name="L"/> (excluding)
        /// </returns>        
        MultidimensionalArray GetScalingsForLinearElements(IGridData gridData, QuadRule qr, int i0, int L);

        /// <summary>
        /// Provides integration metric factors for edges or cells, 
        /// where the transformation from the reference to the physical domain is nonlinear (and also not affine-linear).
        /// For such edges or cells, the metric factor might change for each node in rule <paramref name="qr"/>.
        /// </summary>
        /// <param name="i0">Index of first cell or edge for which to obtain the scaling factors.</param>
        /// <param name="L">Number of cells or edges for which to obtain the scaling factors.</param>
        /// <param name="qr">Quadrature rule in reference element for edges/cells <paramref name="i0"/> to <paramref name="i0"/>+<paramref name="L"/>-1</param>
        /// <param name="gridData">current grid</param>
        /// <returns>
        /// Scaling factors for each node in each cell or edge
        /// - 1st index: cell or edge index, from 0 (including) to <paramref name="L"/> (excluding)
        /// - 2nd index: node index, from 0 (including) to number of nodes in quadrature rule (excluding)
        /// </returns>
        MultidimensionalArray GetScalingsForNonlinElements(IGridData gridData, QuadRule qr, int i0, int L);
    }
}
