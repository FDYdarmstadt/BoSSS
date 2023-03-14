﻿/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Tracing;

namespace BoSSS.Solution.Statistic {
    
    /// <summary>
    /// Utility functions to compare DG fields from different, but geometrically embedded, grids.
    /// </summary>
    public static class DGFieldComparisonEmbedded {

        static Func<DGField, double> NormTypeFactory(NormType nt) {
            switch(nt) {
                case NormType.H1_embedded:
                    return ((DGField f) => f.H1Norm());

                case NormType.L2_embedded:
                    return ((DGField f) => f.L2Norm());

                case NormType.L2noMean_embedded:
                    return ((DGField f) => f.L2Norm_IgnoreMean());

                default:
                    throw new ArgumentOutOfRangeException();
            }
        }

        static Func<DGField, double>[] NormTypeFactory(IList<IEnumerable<DGField>> fields, NormType nt) {
            int L = fields.First().Count();
            var r = NormTypeFactory(nt);
            var rr = new Func<DGField, double>[L];
            rr.SetAll(r);
            return rr;
        }



        /// <summary>
        /// Computes L2 norms between DG fields on different grid resolutions, i.e. for a 
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
        public static void ComputeErrors_L2(IList<IEnumerable<DGField>> fields,
        out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(NormTypeFactory(fields, NormType.L2_embedded),
                fields, out GridRes, out __DOFs, out Errors);
        }

        /// <summary>
        /// Computes H1 norms (aka. Sobolev norm) between DG fields on different grid resolutions, i.e. for a 
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
        public static void ComputeErrors_H1(IList<IEnumerable<DGField>> fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(NormTypeFactory(fields, NormType.H1_embedded),
                fields, out GridRes, out __DOFs, out Errors);


        }

        /// <summary>
        /// Computes H1 norms (aka. Sobolev norm) between DG fields on different grid resolutions, i.e. for a 
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
        public static void ComputeErrors_L2noMean(IList<IEnumerable<DGField>> fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            ComputeErrors(NormTypeFactory(fields, NormType.L2noMean_embedded),
                fields, out GridRes, out __DOFs, out Errors);
        }


        /// <summary>
        /// Computes norms (specified by <paramref name="normTypes"/>) between DG fields on different grid resolutions, i.e. for a 
        /// convergence study, where the solution on the finest grid is assumed to be exact.
        /// </summary>
        /// <param name="fields">
        /// - outer enumeration: sequence of meshes;
        /// - inner enumeration: a set of fields on the same mesh level
        /// </param>
        /// <param name="normTypes">
        /// norm selection for eacjh field; index correlates with second index into <paramref name="fields"/>
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
        public static void ComputeErrors(IList<IEnumerable<DGField>> fields, NormType[] normTypes,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors) {

            if (fields.First().Count() != normTypes.Length)
                throw new ArgumentException("mismatch between number of fields and number of specified norms");

            ComputeErrors(normTypes.Select(nt => NormTypeFactory(nt)).ToArray(),
                fields, out GridRes, out __DOFs, out Errors);
        }



        static void ComputeErrors(
            Func<DGField, double>[] NormFuncS,
            IList<IEnumerable<DGField>> fields,
            out double[] GridRes, out Dictionary<string, long[]> __DOFs, out Dictionary<string, double[]> Errors,
            Func<ilPSP.Vector, bool> SelectionFunc = null) {
            using(var tr = new FuncTrace()) {

                // Grids and coarse-to-fine -- mappings.
                GridData[] gDataS = fields.Select(fc => (GridData)(fc.First().GridDat)).ToArray();

                int[][] Fine2CoarseMapS = new int[gDataS.Length - 1][]; // 1st index: level; 2n index: cell index on finest level
                for(int iLevel = 0; iLevel < Fine2CoarseMapS.Length; iLevel++) {
                    ComputeFine2CoarseMap(gDataS.Last(), gDataS[iLevel], out Fine2CoarseMapS[iLevel]);
                }

                int[] FieldsToCompare = fields[0].Count().ForLoop(i => i); // haCK
                string[] IdentificationS = fields[0].Select(f => f.Identification).ToArray(); // more hack

                // extrapolate to fine grid
                Dictionary<int, List<DGField>> injectedFields = new Dictionary<int, List<DGField>>();
                Dictionary<string, List<int>> DOFs = new Dictionary<string, List<int>>();

                foreach(int index in FieldsToCompare) {
                 
                    List<DGField> fields_Identification = new List<DGField>(); // fields for different resolutions
                    List<int> dofs_Idenitification = new List<int>();

                    DGField finestSolution = fields.Last().ElementAt(index);

                    for(int iLevel = 0; iLevel < gDataS.Length - 1; iLevel++) {
                        //Console.WriteLine("Injecting '{0}' from level {1} to finest grid...", Identification, iLevel);
                        tr.Info(string.Format("Injecting '{0}' from level {1} to finest grid...", index, iLevel));


                        DGField coarseSolution = fields[iLevel].ElementAt(index);

                        if(finestSolution.GetType() != coarseSolution.GetType())
                            throw new NotSupportedException();
                        if(coarseSolution.Basis.Degree != finestSolution.Basis.Degree)
                            throw new NotSupportedException();

                        if(finestSolution is XDGField) {
                            XDGField _coarseSolution = (XDGField)coarseSolution;
                            XDGField _finestSolution = (XDGField)finestSolution;
                            XDGField injectedSolution = new XDGField(_finestSolution.Basis, index + "-inj-" + iLevel);

                            InjectXDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                            fields_Identification.Add(injectedSolution);
                            dofs_Idenitification.Add(coarseSolution.Mapping.GetTotalNoOfDOFs());
                        } else if(finestSolution is SinglePhaseField) {
                            SinglePhaseField _coarseSolution = (SinglePhaseField)coarseSolution;
                            SinglePhaseField _finestSolution = (SinglePhaseField)finestSolution;
                            SinglePhaseField injectedSolution = new SinglePhaseField(_finestSolution.Basis, index + "-inj-" + iLevel);

                            InjectDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                            fields_Identification.Add(injectedSolution);
                            dofs_Idenitification.Add(coarseSolution.Mapping.GetTotalNoOfDOFs());
                        } else {
                            throw new NotSupportedException("DG field type '" + finestSolution.GetType().FullName + "' not supported, Identification is '" + finestSolution.Identification + "'");
                        }

                        tr.Info(string.Format("done."));
                        //Console.WriteLine("done.");
                    }

                    fields_Identification.Add(finestSolution);
                    injectedFields.Add(index, fields_Identification);
                    DOFs.Add(IdentificationS[index], dofs_Idenitification);
                }
                __DOFs = new Dictionary<string, long[]>();
                foreach(var kv in DOFs) {
                    __DOFs.Add(kv.Key, kv.Value.Select((int iii) => (long)iii).ToArray());
                }


                // compute the errors
                Errors = new Dictionary<string, double[]>();
                foreach(int index in FieldsToCompare) {
                

                    double[] L2Error = new double[gDataS.Length - 1];

                    for(int iLevel = 0; iLevel < gDataS.Length - 1; iLevel++) {
                        //Console.WriteLine("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel);
                        tr.Info(string.Format("Computing L2 error of '{0}' on level {1} ...", index, iLevel));

                        DGField Error = injectedFields[index].Last().CloneAs();
                        DGField injSol = injectedFields[index].ElementAt(iLevel);
                        Error.Acc(-1.0, injSol);

                        L2Error[iLevel] = NormFuncS[index](Error);

                        //Console.WriteLine("done (Error is {0:0.####E-00}).", L2Error[iLevel]);
                        tr.Info(string.Format("done (Error is {0:0.####E-00}).", L2Error[iLevel]));
                    }

                    Errors.Add(IdentificationS[index], L2Error);
                }

                GridRes = gDataS.Take(gDataS.Length - 1).Select(gd => gd.Cells.h_minGlobal).ToArray();

            }
        }


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
        /// <param name="SelectionFunc">
        /// if specified, all cells where this evaluates as false are ignored. 
        /// </param>
        public static void ComputeErrors_L2noMean(IEnumerable<string> FieldsToCompare,
            IEnumerable<ITimestepInfo> timestepS,
            out double[] GridRes,
            out Dictionary<string, long[]> __DOFs,
            out Dictionary<string, double[]> L2Errors, out Guid[] timestepIds, Func<ilPSP.Vector, bool> SelectionFunc = null) {
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


                // clear cut-out cells 
                if (SelectionFunc != null) {
                    foreach (var fenum in fields) {
                        foreach (var f in fenum) {
                            GridData grd = (GridData)f.GridDat;
                            //SelectionFunc = X => (X[1] < (3.0/9.0) && ( ((X[0] > 3.0 * (1.5 / 9.0)) && (X[0] < 6.0 * (1.5 / 9.0))) || ((X[0] < -3.0 * (1.5 / 9.0)) && (X[0] > -6.0 * (1.5 / 9.0)))));
                            CellMask cutout = CellMask.GetCellMask(grd, SelectionFunc);
                            Console.WriteLine("number of cut-out cells: {0}", cutout.NoOfItemsLocally);
                            f.Clear(cutout);
                        }
                    }
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

                // continue in other routine
                ComputeErrors_L2noMean(fields, out GridRes, out __DOFs, out L2Errors);
            }
        }


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
        /// <param name="SelectionFunc">
        /// if specified, all cells where this evaluates as false are ignored. 
        /// </param>
        public static void ComputeErrors_L2(IEnumerable<string> FieldsToCompare,
            IEnumerable<ITimestepInfo> timestepS,
            out double[] GridRes, 
            out Dictionary<string,long[]> __DOFs, 
            out Dictionary<string, double[]> L2Errors, out Guid[] timestepIds, Func<ilPSP.Vector, bool> SelectionFunc = null) {  
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


                // clear cut-out cells 
                if (SelectionFunc != null) {
                    foreach (var fenum in fields) {
                        foreach (var f in fenum) {
                            GridData grd = (GridData)f.GridDat;
                            //SelectionFunc = X => (X[1] < (3.0/9.0) && ( ((X[0] > 3.0 * (1.5 / 9.0)) && (X[0] < 6.0 * (1.5 / 9.0))) || ((X[0] < -3.0 * (1.5 / 9.0)) && (X[0] > -6.0 * (1.5 / 9.0)))));
                            CellMask cutout = CellMask.GetCellMask(grd, SelectionFunc);
                            Console.WriteLine("number of cut-out cells: {0}", cutout.NoOfItemsLocally);
                            f.Clear(cutout);
                        }
                    }
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

                // continue in other routine
                ComputeErrors_L2(fields, out GridRes, out __DOFs, out L2Errors);
            }
        }

        /// <summary>
        /// Temporal convergence, 
        /// i.e. all compared timesteps are on the same grid (and physical time) and we compare to the simulation using the smallest temporal resolution
        /// </summary>
        /// <param name="FieldsToCompare">
        /// Identification (<see cref="DGField.Identification"/>) of the fields which should be compared.
        /// </param>
        /// <param name="timestepS">
        /// A collection of solutions on different grid resolutions.
        /// </param>
        /// <param name="TempRes">
        /// On exit, various timestep sizes of the solutions.
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
        /// on exit, the timestep id which correlate with the individual <paramref name="timestepS"/>.
        /// (remarks: <paramref name="timestepIds"/> may be re-sorted internally according to resolution).
        /// </param>
        /// <param name="SelectionFunc">
        /// if specified, all cells where this evaluates as false are ignored. 
        /// </param>
        public static void ComputeErrorsTemporal(IEnumerable<string> FieldsToCompare, IEnumerable<ITimestepInfo> timestepS,
          out double[] TempRes, 
          out Dictionary<string, long[]> __DOFs, 
          out Dictionary<string, double[]> L2Errors, 
          out Guid[] timestepIds, Func<ilPSP.Vector, bool> SelectionFunc = null) {
            using (var tr = new FuncTrace()) {
                if (FieldsToCompare == null || FieldsToCompare.Count() <= 0)
                    throw new ArgumentException("empty list of field names.");
                if (timestepS == null || timestepS.Count() < 1)
                    throw new ArgumentException("requiring at least two different solutions.");

                // sort according to temporal resolution
                {
                    timestepS = timestepS.OrderBy(f => f.TimeStepNumber.First()).ToArray();

                    // filter unequal timesteps:
                    double phystime = timestepS.Last().PhysicalTime;
                    foreach (var timestep in timestepS)
                        if (timestep.PhysicalTime != phystime) ((List<ITimestepInfo>)timestepS).Remove(timestep);

                    // filter unequal grids
                    int cells = timestepS.Last().Grid.NumberOfCells;
                    foreach (var timestep in timestepS)
                        if (timestep.Grid.NumberOfCells != cells) ((List<ITimestepInfo>)timestepS).Remove(timestep);

                    // extract timestep Id's
                    timestepIds = new Guid[timestepS.Count()];
                    for (int z = 0; z < timestepIds.Length; z++) {
                        timestepIds[z] = timestepS.ElementAt(z).ID;
                    }
                }

                // load the DG-Fields
                List<IEnumerable<DGField>> fields = new List<IEnumerable<DGField>>();
                int i = 1;
                foreach (var timestep in timestepS) {
                    //Console.WriteLine("Loading timestep {0} of {1}, ({2})...", i, timestepS.Count(), timestep.ID);
                    fields.Add(timestep.Fields);
                    i++;
                    //Console.WriteLine("done (Grid has {0} cells).", fields.Last().First().GridDat.CellPartitioning.TotalLength);
                }


                // clear cut-out cells 
                if (SelectionFunc != null) {
                    foreach (var fenum in fields) {
                        foreach (var f in fenum) {
                            GridData grd = (GridData)f.GridDat;
                            //SelectionFunc = X => (X[1] < (3.0/9.0) && ( ((X[0] > 3.0 * (1.5 / 9.0)) && (X[0] < 6.0 * (1.5 / 9.0))) || ((X[0] < -3.0 * (1.5 / 9.0)) && (X[0] > -6.0 * (1.5 / 9.0)))));
                            CellMask cutout = CellMask.GetCellMask(grd, SelectionFunc);
                            Console.WriteLine("number of cut-out cells: {0}", cutout.NoOfItemsLocally);
                            f.Clear(cutout);
                        }
                    }
                }

                // extrapolate all solutions to same grid, this is necessary so references for grid and tracker are the same
                Dictionary<string, List<DGField>> injectedFields = new Dictionary<string, List<DGField>>();
                Dictionary<string, List<int>> DOFs = new Dictionary<string, List<int>>();

                int[][] Fine2CoarseMapS = new int[timestepS.Count() - 1][]; // 1st index: level; 2n index: cell index on finest level
                for (int iLevel = 0; iLevel < Fine2CoarseMapS.Length; iLevel++) {
                    ComputeFine2CoarseMap((GridData)fields.Last().First().GridDat, (GridData)fields.ElementAt(iLevel).First().GridDat, out Fine2CoarseMapS[iLevel]);
                }

                foreach (string Identification in FieldsToCompare) {
                    List<DGField> fields_Identification = new List<DGField>(); // fields for different resolutions
                    List<int> dofs_Idenitification = new List<int>();

                    DGField finestSolution = fields.Last().Single(f => f.Identification == Identification);

                    for (int iLevel = 0; iLevel < timestepS.Count() - 1; iLevel++) {
                        //Console.WriteLine("Injecting '{0}' from level {1} to finest grid...", Identification, iLevel);
                        tr.Info(string.Format("Injecting '{0}' from level {1} to finest grid...", Identification, iLevel));


                        DGField coarseSolution = fields[iLevel].Single(f => f.Identification == Identification);

                        if (finestSolution.GetType() != coarseSolution.GetType())
                            throw new NotSupportedException();
                        if (coarseSolution.Basis.Degree != finestSolution.Basis.Degree)
                            throw new NotSupportedException();

                        if (finestSolution is XDGField) {
                            XDGField _coarseSolution = (XDGField)coarseSolution;
                            XDGField _finestSolution = (XDGField)finestSolution;
                            XDGField injectedSolution = new XDGField(_finestSolution.Basis, Identification + "-inj-" + iLevel);

                            InjectXDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                            fields_Identification.Add(injectedSolution);
                            dofs_Idenitification.Add(coarseSolution.Mapping.GetTotalNoOfDOFs());
                        } else if (finestSolution is SinglePhaseField) {
                            SinglePhaseField _coarseSolution = (SinglePhaseField)coarseSolution;
                            SinglePhaseField _finestSolution = (SinglePhaseField)finestSolution;
                            SinglePhaseField injectedSolution = new SinglePhaseField(_finestSolution.Basis, Identification + "-inj-" + iLevel);

                            InjectDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                            fields_Identification.Add(injectedSolution);
                            dofs_Idenitification.Add(coarseSolution.Mapping.GetTotalNoOfDOFs());
                        } else {
                            throw new NotSupportedException("DG field type '" + finestSolution.GetType().FullName + "' not supported, Identification is '" + finestSolution.Identification + "'");
                        }

                        tr.Info(string.Format("done."));
                        //Console.WriteLine("done.");
                    }

                    fields_Identification.Add(finestSolution);
                    injectedFields.Add(Identification, fields_Identification);
                    DOFs.Add(Identification, dofs_Idenitification);
                }
                __DOFs = new Dictionary<string, long[]>();
                foreach (var kv in DOFs) {
                    __DOFs.Add(kv.Key, kv.Value.Select((int iii) => (long)iii).ToArray());
                }


                // compute the errors
                L2Errors = new Dictionary<string, double[]>();
                foreach (string Identification in FieldsToCompare) {

                    double[] L2Error = new double[timestepS.Count() - 1];

                    for (int iLevel = 0; iLevel < timestepS.Count() - 1; iLevel++) {
                        //Console.WriteLine("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel);
                        tr.Info(string.Format("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel));

                        DGField Error = injectedFields[Identification].Last().CloneAs();
                        DGField injSol = injectedFields[Identification].ElementAt(iLevel);
                        Error.Acc(-1.0, injSol);

                        L2Error[iLevel] = Error.L2Norm();

                        //Console.WriteLine("done (Error is {0:0.####E-00}).", L2Error[iLevel]);
                        tr.Info(string.Format("done (Error is {0:0.####E-00}).", L2Error[iLevel]));
                    }

                    L2Errors.Add(Identification, L2Error);
                }                           

                TempRes = timestepS.Take(timestepS.Count() - 1).Select(gd => gd.PhysicalTime / (gd.TimeStepNumber.First() - 1)).ToArray();
            }
        }

        /// <summary>
        /// Injects an XDG field from a coarse grid to a fine grid.
        /// </summary>
        public static XDGField InjectXDGField(int[] Fine2Coarse, XDGField injected, XDGField cors_field, CellMask subGrd = null) {
            
            var trk = injected.Basis.Tracker;

            foreach (var spc in trk.SpeciesIdS) {
                var grd = trk.Regions.GetSpeciesMask(spc);

                InjectDGField(Fine2Coarse,
                    injected.GetSpeciesShadowField(spc),
                    cors_field.GetSpeciesShadowField(spc),
                    subGrd == null ? grd : grd.Intersect(subGrd));
            }
            return injected;
        }

        /// <summary>
        /// Creates a mapping between the cells of a coarse grid and a fine gird.
        /// </summary>
        /// <param name="ctxFine">Fine resolution grid.</param>
        /// <param name="ctxCoarse">Coarse resolution grid.</param>
        /// <param name="Fine2CoarseMap">
        /// For each cell in 
        /// </param>
        public static void ComputeFine2CoarseMap(GridData ctxFine, GridData ctxCoarse, out int[] Fine2CoarseMap) {

            int Jf = ctxFine.Grid.NoOfUpdateCells;
            Fine2CoarseMap = new int[Jf];

            int D = ctxFine.Grid.SpatialDimension;
            MultidimensionalArray Centas = MultidimensionalArray.Create(Jf, 1, D);

            for(int j = 0; j < Jf; ) {
                int Len = ctxFine.Cells.GetNoOfSimilarConsecutiveCells(CellInfo.AllOn, j, Jf - j);
                ctxFine.TransformLocal2Global(ctxFine.Cells.GetRefElement(j).Center, j, Len, Centas, j);
                j += Len;
            }

            CellLocalization cl = new CellLocalization(ctxCoarse);
            int NoOfUnassigned;
            cl.LocalizePointsWithinGrid(Centas.ResizeShallow(Jf,D), Fine2CoarseMap, out NoOfUnassigned);

            if (NoOfUnassigned > 0)
                throw new ApplicationException("fucking scheisse.");
        }


        /// <summary>
        /// Injects a DG field from a coarser grid to a fine grid.
        /// </summary>
        public static void InjectDGField(int[] CellIdxFine2Coarse, ConventionalDGField onFineGrid, ConventionalDGField onCoarseGrid, CellMask subGrd = null) {
            if(subGrd == null)
                subGrd = CellMask.GetFullMask(onFineGrid.GridDat);
            if (onFineGrid.Basis.GridDat.iGeomCells.RefElements.Length != 1)
                throw new NotImplementedException("todo");
            var QR = onFineGrid.Basis.GridDat.iGeomCells.RefElements[0].GetQuadratureRule(onFineGrid.Basis.Degree * 2);
            if(!object.ReferenceEquals(subGrd.GridData, onFineGrid.GridDat))
                throw new ArgumentException();

            int N = onFineGrid.Basis.Length;
            int K = QR.NoOfNodes;
            int D = onFineGrid.Basis.GridDat.SpatialDimension;
            
            NodeSet xi = QR.Nodes;                                           // quadrature nodes in local coordsys. of cell to extrapolate TO
            //MultidimensionalArray eta = MultidimensionalArray.Create(K, D);  // quadrature nodes in local coordsys. of cell to extrapolate FROM
            MultidimensionalArray tmp = MultidimensionalArray.Create(K, D);  // quadrature nodes in global coordinates
            MultidimensionalArray a = MultidimensionalArray.Create(K, N);

            for (int n = 0; n < N; n++) {
                onFineGrid.Basis.Polynomials[0][n].Evaluate(a.ExtractSubArrayShallow(-1, n), xi);
                for (int k = 0; k < K; k++)
                    a[k, n] *= QR.Weights[k];
            }
            MultidimensionalArray v = MultidimensionalArray.Create(K, N);
            MultidimensionalArray[] v_ = new MultidimensionalArray[N];
            for (int n = 0; n < N; n++)
                v_[n] = v.ExtractSubArrayShallow(-1, n);

            double[,] eps = new double[N, N];
            double[] u2 = new double[N];
            double[] u1 = new double[N];

            //double[] ooScales = m_context.GridDat.OneOverSqrt_AbsDetTransformation;

            IGridData ctxFine = onFineGrid.Basis.GridDat;
            IGridData ctxCors = onCoarseGrid.Basis.GridDat;

            var ooScalesFine = ctxFine.ChefBasis.Scaling;
            var ooScalesCors = ctxCors.ChefBasis.Scaling;


            foreach (var chunk in subGrd) {
                for (int j = chunk.i0; j < chunk.JE; j++) {

                    int _1 = CellIdxFine2Coarse[j]; // cell to extrapolate FROM (on coarse grid)
                    int _2 = j;                     // cell to extrapolate TO   (on fine grid)

                    int iKref_1 = ctxCors.iGeomCells.GetRefElementIndex(_1);
                    int iKref_2 = ctxFine.iGeomCells.GetRefElementIndex(_2);

                    if (!ctxCors.iGeomCells.IsCellAffineLinear(_1))
                        throw new NotImplementedException("todo");
                    if (!ctxFine.iGeomCells.IsCellAffineLinear(_2))
                        throw new NotImplementedException("todo");
                    
                    
                    // get DG coordinates
                    onCoarseGrid.Coordinates.GetRow(_1, u1);


                    // transform quad. nodes from cell 2 (extrapolate to) to cell 1 (extrapolate FROM)
                    ctxFine.TransformLocal2Global(xi, tmp, _2);
                    NodeSet eta = new NodeSet(ctxFine.iGeomCells.GetRefElement(_2), K, D, false);
                    ctxCors.TransformGlobal2Local(tmp, eta, _1, null);
                    eta.LockForever();

                    // evaluate Polynomials of cell 1 (fine)
                    v.Clear();
                    for(int n = 0; n < N; n++) {
                        if(n < onCoarseGrid.Basis.Length)
                            onCoarseGrid.Basis.Polynomials[iKref_1][n].Evaluate(v_[n], eta);
                        else
                            v_[n].Clear();
                    }

                    // perform quadrature
                    double scale = ooScalesCors[_1] / ooScalesFine[_2];
                    for (int m = 0; m < N; m++) {
                        for (int n = 0; n < N; n++) {
                            double eps_m_n = 0;
                            for (int k = 0; k < K; k++) {
                                eps_m_n += a[k, m] * v[k, n];
                            }
                            eps[m, n] = eps_m_n * scale;
                        }
                    }

                    // new DG coordinates
                    for (int m = 0; m < N; m++) {
                        double u2_m = 0;
                        for (int n = 0; n < N; n++) {
                            u2_m += eps[m, n] * u1[n];
                        }
                        u2[m] = u2_m;
                    }

                    // set DG coordinates
                    onFineGrid.Coordinates.SetRow(_2, u2);
                }
            }
        }

    }
}
