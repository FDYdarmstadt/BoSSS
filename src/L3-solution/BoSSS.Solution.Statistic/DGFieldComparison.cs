using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.Statistic {
    
    /// <summary>
    /// Utility functions to compare DG fields from different, but geometrically embedded, grids.
    /// </summary>
    public class DGFieldComparison {

        /// <summary>
        /// Driver routine to 
        /// computes norms between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the L2 error 
        /// (for each field specified in <paramref name="fields"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="fields"/>).
        /// </param>
        /// <param name="nt">
        /// Selection of the appropriate norm (the same norm for all fields)
        /// </param>
        public static void ComputeErrors(IList<IEnumerable<DGField>> fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors,
            NormType nt) {
            using(var tr = new FuncTrace()) {
                tr.Info($"ComputeErrors {fields[0].Select(f => f.Identification).ToConcatString("", ", ", "")} in norm {nt}...");


                switch (nt) {
                    case NormType.H1_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_H1(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    case NormType.L2_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_L2(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    case NormType.L2noMean_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_L2noMean(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    case NormType.H1_embedded:
                    DGFieldComparisonEmbedded.ComputeErrors_H1(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    case NormType.L2_embedded:
                    DGFieldComparisonEmbedded.ComputeErrors_L2(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    case NormType.L2noMean_embedded:
                    DGFieldComparisonEmbedded.ComputeErrors_L2noMean(fields, out GridRes, out __DOFs, out Errors);
                    return;

                    default:
                    throw new NotImplementedException();
                }
            }
        }


        /// <summary>
        /// Driver routine to 
        /// computes compute norms between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="FieldsToCompare">
        /// Identification (<see cref="DGField.Identification"/>) of the fields which should be compared.
        /// </param>
        /// <param name="timestepS">
        /// A collection of solutions on different grid resolutions.
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the L2 error 
        /// (for each field specified in <paramref name="FieldsToCompare"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="FieldsToCompare"/>).
        /// </param>
        /// <param name="timestepIds">
        /// on exit, the timestep id which correlate with the resolutions <paramref name="GridRes"/>
        /// (remarks: <paramref name="timestepIds"/> may be re-sorted internally according to grid resolution).
        /// </param>
        /// <param name="nt">
        /// Selection of the appropriate norm
        /// </param>        
        public static void ComputeErrors(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors, out Guid[] timestepIds,
            NormType nt) {
            using(var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                tr.Info($"ComputeErrors {FieldsToCompare.ToConcatString("",", ","")} in norm {nt}...");
                switch (nt) {
                    case NormType.H1_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_H1(FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
                    return;

                    case NormType.L2_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_L2(FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
                    return;

                    case NormType.L2noMean_approximate:
                    DGFieldComparisonNonEmb.ComputeErrors_L2noMean(FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
                    return;

                    case NormType.L2_embedded:
                    DGFieldComparisonEmbedded.ComputeErrors_L2(FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
                    return;

                    case NormType.L2noMean_embedded:
                    DGFieldComparisonEmbedded.ComputeErrors_L2noMean(FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
                    return;

                    default:
                    throw new NotImplementedException();
                }
            }
        }
    }
}
