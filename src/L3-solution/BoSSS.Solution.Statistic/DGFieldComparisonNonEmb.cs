using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;

namespace BoSSS.Solution.Statistic {

    /// <summary>
    /// Utility functions to compare DG fields from different, non-embedded, grids.
    /// </summary>
    public static class DGFieldComparisonNonEmb {

        /// <summary>
        /// Computes L2 norms between DG fields on different grid resolutions, i.e. for a 
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
        public static void ComputeErrors_L2(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
          out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors, out Guid[] timestepIds) {
            
            ComputeErrors(Enumerable.Repeat(NormFactory(NormType.L2_approximate), FieldsToCompare.Count()).ToArray(), 
                FieldsToCompare, timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
        }

        /// <summary>
        /// Computes L2 norms, ignoring the mean value (useful e.g. for pressure) between DG fields on different grid resolutions, i.e. for a 
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
        public static void ComputeErrors_L2noMean(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
          out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors, out Guid[] timestepIds) {

            

            ComputeErrors(Enumerable.Repeat(NormFactory(NormType.L2noMean_approximate), FieldsToCompare.Count()).ToArray(), FieldsToCompare,
                timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
        }

        /// <summary>
        /// Computes H1 norms between DG fields on different grid resolutions, i.e. for a 
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
        /// On exit, the H1 error 
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
        public static void ComputeErrors_H1(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
          out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors, out Guid[] timestepIds) {

            


            ComputeErrors(Enumerable.Repeat(NormFactory(NormType.H1_approximate), FieldsToCompare.Count()).ToArray(), FieldsToCompare,
                timestepS, out GridRes, out __DOFs, out Errors, out timestepIds);
        }


        /// <summary>
        /// Computes L2 norms between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="__fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the L2 error 
        /// (for each field specified in <paramref name="__fields"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="__fields"/>).
        /// </param>
        public static void ComputeErrors_L2(IList<IEnumerable<DGField>> __fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(__fields, Enumerable.Repeat(NormType.L2_approximate, __fields.First().Count()).ToArray(),
                out GridRes, out __DOFs, out Errors);
        }

        /// <summary>
        /// Computes L2 norms, ignoring the mean value (useful e.g. for pressure) between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="__fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the L2 error 
        /// (for each field specified in <paramref name="__fields"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="__fields"/>).
        /// </param>
        public static void ComputeErrors_L2noMean(IList<IEnumerable<DGField>> __fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(__fields, Enumerable.Repeat(NormType.L2noMean_approximate, __fields.First().Count()).ToArray(), 
                out GridRes, out __DOFs, out Errors);
        }

        /// <summary>
        /// Computes H1 norms between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="__fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the H1 error 
        /// (for each field specified in <paramref name="__fields"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="__fields"/>).
        /// </param>
        public static void ComputeErrors_H1(IList<IEnumerable<DGField>> __fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(__fields, Enumerable.Repeat(NormType.H1_approximate, __fields.First().Count()).ToArray(),
                out GridRes, out __DOFs, out Errors);
        }


        /// <summary>
        /// Computes H1 norms between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="__fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="normTypes">
        /// norm selection for eacjh field; index correlates with second index into <paramref name="__fields"/>
        /// </param>
        /// <param name="GridRes">
        /// On exit, the resolution of the different grids.
        /// </param>
        /// <param name="Errors">
        /// On exit, the H1 error 
        /// (for each field specified in <paramref name="__fields"/>)
        /// in comparison to the solution on the finest grid.
        /// </param>
        /// <param name="__DOFs">
        /// On exit, the number of degrees-of-freedom 
        /// (for each field specified in <paramref name="__fields"/>).
        /// </param>
        public static void ComputeErrors(IList<IEnumerable<DGField>> __fields, NormType[] normTypes,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            


            ComputeErrors(normTypes.Select(nt => NormFactory(nt)).ToArray(), __fields, out GridRes, out __DOFs, out Errors);
        }

        static double DistFunc_L2nomean(ConventionalDGField coarse, ConventionalDGField fine) {
            return coarse.L2Distance(fine, IgnoreMeanValue: true);
        }

        static double DistFunc_L2(ConventionalDGField coarse, ConventionalDGField fine) {
            return coarse.L2Distance(fine);
        }

        static double DistFunc_H1(ConventionalDGField coarse, ConventionalDGField fine) {
            return coarse.H1Distance(fine);
        }

        static Func<ConventionalDGField, ConventionalDGField, double> NormFactory(NormType nt) {
            switch (nt) {
                case NormType.H1_approximate:
                    return DistFunc_H1;

                case NormType.L2_approximate:
                    return DistFunc_L2;

                case NormType.L2noMean_approximate:
                    return DistFunc_L2nomean;

                default:
                    throw new ArgumentOutOfRangeException();
            }
        }


        /*
        /// <summary>
        /// computation based on time-steps
        /// </summary>
        static void ComputeErrors(Func<ConventionalDGField, ConventionalDGField, double> distFunc,
            IEnumerable<string> FieldsToCompare,
            IEnumerable<ITimestepInfo> timestepS,
            out double[] GridRes,
            out Dictionary<string, long[]> __DOFs,
            out Dictionary<string, double[]> Errors,
            out Guid[] timestepIds) {
            int L = FieldsToCompare.Count();
            
            ComputeErrors(Enumerable.Repeat(distFunc, L).ToArray(),
                FieldsToCompare,
                timestepS,
                out GridRes,
                out __DOFs,
                out Errors,
                out timestepIds);

        }
        */


        /// <summary>
        /// computation based on time-steps
        /// </summary>
        static void ComputeErrors(Func<ConventionalDGField, ConventionalDGField, double>[] distFuncS,
          IEnumerable<string> FieldsToCompare,
          IEnumerable<ITimestepInfo> timestepS,
          out double[] GridRes,
          out Dictionary<string, long[]> __DOFs,
          out Dictionary<string, double[]> Errors,
          out Guid[] timestepIds) {
            using (var tr = new FuncTrace()) {

                if (FieldsToCompare == null || FieldsToCompare.Count() <= 0)
                    throw new ArgumentException("empty list of field names.");
                if (timestepS == null || timestepS.Count() < 1)
                    throw new ArgumentException("requiring at least two different solutions.");

                // load the DG-Fields
                List<IEnumerable<DGField>> fields = new List<IEnumerable<DGField>>(); // 1st index: grid / 2nd index: enumeration
                int i = 1;
                foreach (var timestep in timestepS) {
                    tr.Info(string.Format("Loading timestep {0} of {1}, ({2})...", i, timestepS.Count(), timestep.ID));
                    fields.Add(timestep.Fields.Where(f => FieldsToCompare.Contains(f.Identification)).ToArray());
                    i++;
                    tr.Info(string.Format("done (Grid has {0} cells).", fields.Last().First().GridDat.CellPartitioning.TotalLength));
                }


                // sort according to grid resolution
                {
                    var s = fields.OrderBy(f => f.First().GridDat.CellPartitioning.TotalLength).ToArray();
                    var orgfields = fields.ToArray();
                    fields.Clear();
                    fields.AddRange(s);
                    s = null;

                    // filter equal grids:
                    while (fields.Count >= 2
                        && (fields[fields.Count - 1].First().GridDat.CellPartitioning.TotalLength
                        == fields[fields.Count - 2].First().GridDat.CellPartitioning.TotalLength)) {
                        fields.RemoveAt(fields.Count - 2);
                    }

                    // extract timestep Id's
                    timestepIds = new Guid[fields.Count];
                    for (int z = 0; z < timestepIds.Length; z++) {
                        int idx = orgfields.IndexOf(fields[z], (f1, f2) => object.ReferenceEquals(f1, f2));
                        timestepIds[z] = timestepS.ElementAt(idx).ID;
                    }
                }

                ComputeErrors(distFuncS, fields, out GridRes, out __DOFs, out Errors);

            }
        }

        /*
        /// <summary>
        /// computation based on DG-fields
        /// </summary>
        static void ComputeErrors(Func<ConventionalDGField, ConventionalDGField, double> distFunc,
            IList<IEnumerable<DGField>> __fields,
            out double[] GridRes,
            out Dictionary<string, long[]> __DOFs,
            out Dictionary<string, double[]> Errors) {
            int L = __fields.First().Count();

            ComputeErrors(Enumerable.Repeat(distFunc, L).ToArray(),
                __fields,
                out GridRes,
                out __DOFs,
                out Errors);

        }
        */

        /// <summary>
        /// computation based on DG-fields
        /// </summary>
        static void ComputeErrors(Func<ConventionalDGField, ConventionalDGField, double>[] distFuncS,
            IList<IEnumerable<DGField>> __fields,
            out double[] GridRes,
            out Dictionary<string, long[]> __DOFs,
            out Dictionary<string, double[]> Errors) {
            using (var tr = new FuncTrace()) {
                tr.InfoToConsole = false;
                if (__fields == null || __fields.Count() <= 2)
                    throw new ArgumentException("expecting at least solutions on two different meshes.");




                // load the DG-Fields
                List<IEnumerable<DGField>> fields = new List<IEnumerable<DGField>>(__fields); // 1st index: grid / 2nd index: enumeration
                tr.Info("Computing errors for fields: " + fields[0].Select(f => f.Identification).ToConcatString("", ", ", ";"));

                // sort according to grid resolution
                {
                    var s = fields.OrderBy(f => f.First().GridDat.CellPartitioning.TotalLength).ToArray();
                    var orgfields = fields.ToArray();
                    fields.Clear();
                    fields.AddRange(s);
                    s = null;

                    // filter equal grids:
                    while (fields.Count >= 2
                        && (fields[fields.Count - 1].First().GridDat.CellPartitioning.TotalLength
                        == fields[fields.Count - 2].First().GridDat.CellPartitioning.TotalLength)) {
                        fields.RemoveAt(fields.Count - 2);
                    }
                }

                string[] FieldsToCompare = fields.First().Select(dgf => dgf.Identification).ToArray();
                if (FieldsToCompare.Length != distFuncS.Length)
                    throw new ArgumentException("mismatch between number of fields and distance functions.");
                foreach (var flds in fields.Skip(1)) {
                    if (!FieldsToCompare.SetEquals(flds.Select(dgf => dgf.Identification))) {
                        throw new ArgumentException("DG Field identifications must match on all mesh levels in order to be correlated.");
                    }
                }
                tr.Info("Comparison Names: " + FieldsToCompare.ToConcatString("", ", ", ";"));


                // grids and resolution
                GridData[] gDataS = fields.Select(fc => GridHelper.ExtractGridData(fc.First().GridDat)).ToArray();
                GridRes = gDataS.Take(gDataS.Length - 1).Select(gd => gd.Cells.h_minGlobal).ToArray();

                // compute the errors
                Errors = new Dictionary<string, double[]>();
                __DOFs = new Dictionary<string, long[]>();
                for (int iField = 0; iField < FieldsToCompare.Length; iField++) {

                    string Identification = FieldsToCompare[iField];

                    double[] L2Error = new double[gDataS.Length - 1];
                    long[] dof = new long[gDataS.Length - 1];

                    for (int iLevel = 0; iLevel < gDataS.Length - 1; iLevel++) {
                        //Console.WriteLine("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel);
                        tr.Info(string.Format("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel));

                        ConventionalDGField fine = (ConventionalDGField)(fields.Last().Single(fi => fi.Identification == Identification));
                        ConventionalDGField coarse = (ConventionalDGField)(fields.ElementAt(iLevel).Single(fi => fi.Identification == Identification));

                        L2Error[iLevel] = distFuncS[iField](coarse, fine);
                        dof[iLevel] = coarse.Mapping.TotalLength;

                        //Console.WriteLine("done (Error is {0:0.####E-00}).", L2Error[iLevel]);
                        tr.Info(string.Format("done '{0}' on level {1} (Error is {2:0.####E-00}, dof {3}, h is {4}).", Identification, iLevel, L2Error[iLevel], dof[iLevel], GridRes[iLevel]));

                    }

                    Errors.Add(Identification, L2Error);
                    __DOFs.Add(Identification, dof);
                }

            }
        }

        //*/
    }
}
