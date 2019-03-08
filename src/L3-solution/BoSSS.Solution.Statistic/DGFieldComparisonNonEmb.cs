using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;


namespace BoSSS.Solution.Statistic {
    /*
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
        /// <param name="L2Errors">
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
        public static void ComputeErrors(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
          out double[] GridRes, out Dictionary<string, int[]> __DOFs, out Dictionary<string, double[]> L2Errors, out Guid[] timestepIds) {
            using (var tr = new FuncTrace()) {
                
                if (FieldsToCompare == null || FieldsToCompare.Count() <= 0)
                    throw new ArgumentException("empty list of field names.");
                if (timestepS == null || timestepS.Count() < 1)
                    throw new ArgumentException("requiring at least two different solutions.");

                // load the DG-Fields
                List<IEnumerable<DGField>> fields = new List<IEnumerable<DGField>>();
                int i = 1;
                foreach (var timestep in timestepS) {
                    //Console.WriteLine("Loading timestep {0} of {1}, ({2})...", i, timestepS.Count(), timestep.ID);
                    fields.Add(timestep.Fields);
                    i++;
                    //Console.WriteLine("done (Grid has {0} cells).", fields.Last().First().GridDat.CellPartitioning.TotalLength);
                }


                // sort according to grid resolution
                {
                    var s = fields.OrderBy(f => f.First().GridDat.CellPartitioning.TotalLength).ToArray();
                    var orgfields = fields.ToArray();
                    fields.Clear();
                    fields.AddRange(s);
                    s = null;

                    // filter equal grids:
                    while(fields.Count >= 2 
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


                // Grids 
                GridData[] gDataS = fields.Select(fc => (GridData)(fc.First().GridDat)).ToArray();


                
            }
        }
        
    }*/
}
