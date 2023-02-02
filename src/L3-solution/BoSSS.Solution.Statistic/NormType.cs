using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.Statistic {
    
    
    /// <summary>
    /// Switch between various norms in which spatial convergence can be measured;
    /// There are two general concepts for computing 
    /// the norm between two DG fields on different meshes:
    /// - Embedded: in a comparison of two DG fields, the one on the coarser mesh is projected onto the finer mesh.
    ///   It is required that the coarser mesh is embedded in the finer one, which makes the projection exact (i.e. an injection).
    ///   Then the norm (of the difference) between the two DG fields is computed on the finer mesh.
    /// - Approximate: if meshes are not embedded, one can still compute the norm w.r.t. quadrature nodes on the fine mesh.
    ///   The DG field on the coarser mesh is evaluated at these points, but one might run into aliasing in-accuracies 
    ///   if the meshes are not perfectly embedded.
    ///   **This approach is recommended for curved cells meshes.**
    /// </summary>
    public enum NormType {

        /// <summary>
        /// Norm by <see cref="DGFieldComparisonEmbedded.ComputeErrors_L2(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>; 
        /// very accurate, but requires geometrically embedded meshes
        /// </summary>
        L2_embedded,

        /// <summary>
        /// Norm by <see cref="DGFieldComparisonEmbedded.ComputeErrors_L2noMean(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>; 
        /// very accurate, but requires geometrically embedded meshes
        /// </summary>
        L2noMean_embedded,

        /// <summary>
        /// Norm by <see cref="DGFieldComparisonEmbedded.ComputeErrors_H1(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>; 
        /// very accurate, but requires geometrically embedded meshes
        /// </summary>
        H1_embedded,

        /// <summary>
        /// Norm computed by <see cref="DGFieldComparisonNonEmb.ComputeErrors_L2(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>
        /// less accurate than <see cref="L2_embedded"/>, but provides results on arbitrary meshes
        /// </summary>
        L2_approximate,

        /// <summary>
        /// Approximate comparison, the mean value (average) of the DG fields is ignored.
        /// This is typically used for comparing pressure in incompressible simulations.
        /// 
        /// Norm computed by <see cref="DGFieldComparisonNonEmb.ComputeErrors_L2noMean(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>.
        /// </summary>
        L2noMean_approximate,


        /// <summary>
        /// Approximate comparison in the H1-Norm (Sobolev-Norm, i.e. L2-norm plus L2-norm of the gradient).
        /// 
        /// Norm computed by <see cref="DGFieldComparisonNonEmb.ComputeErrors_H1(IList{IEnumerable{Foundation.DGField}}, out double[], out Dictionary{string, long[]}, out Dictionary{string, double[]})"/>
        /// </summary>
        H1_approximate
    }
}
