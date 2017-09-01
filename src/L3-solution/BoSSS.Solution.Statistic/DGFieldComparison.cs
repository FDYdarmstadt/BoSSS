/* =======================================================================
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

namespace BoSSS.Solution.Statistic {
    
    /// <summary>
    /// Utility functions to compare DG fields from different grids.
    /// </summary>
    public static class DGFieldComparison {

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
        public static void ComputeErrors(string[] FieldsToCompare, ITimestepInfo[] timestepS,
          out double[] GridRes, out Dictionary<string, double[]> L2Errors) {

            // load the DG-Fields
            List<IEnumerable<DGField>> fields = new List<IEnumerable<DGField>>();
            int i = 1;
            foreach(var timestep in timestepS) {
                Console.WriteLine("Loading timestep {0} of {1}, ({2})...", i, timestepS.Length, timestep.ID);
                fields.Add(timestep.Fields);
                i++;
                Console.WriteLine("done (Grid has {0} cells).", fields.Last().First().GridDat.CellPartitioning.TotalLength);
            }

            {
                var s = fields.OrderBy(f => f.First().GridDat.CellPartitioning.TotalLength).ToArray();
                fields.Clear();
                fields.AddRange(s);
            }

            // Grids and coarse-to-fine -- mappings.
            GridData[] gDataS = fields.Select(fc => (GridData)(fc.First().GridDat)).ToArray();

            int[][] Fine2CoarseMapS = new int[gDataS.Length - 1][];
            for(int iLevel = 0; iLevel < Fine2CoarseMapS.Length; iLevel++) {
                ComputeFine2CoarseMap(gDataS.Last(), gDataS[iLevel], out Fine2CoarseMapS[iLevel]);
            }

            // extrapolate to fine grid
            Dictionary<string, List<DGField>> injectedFields = new Dictionary<string, List<DGField>>();
            foreach(string Identification in FieldsToCompare) {
                List<DGField> blabla = new List<DGField>();

                DGField finestSolution = fields.Last().Single(f => f.Identification == Identification);

                for(int iLevel = 0; iLevel < gDataS.Length - 1; iLevel++) {
                    Console.WriteLine("Injecting '{0}' from level {1} to finest grid...", Identification, iLevel);

                    DGField coarseSolution = fields[iLevel].Single(f => f.Identification == Identification);

                    if(finestSolution.GetType() != coarseSolution.GetType())
                        throw new NotSupportedException();
                    if(coarseSolution.Basis.Degree != finestSolution.Basis.Degree)
                        throw new NotSupportedException();

                    if(finestSolution is XDGField) {
                        XDGField _coarseSolution = (XDGField)coarseSolution;
                        XDGField _finestSolution = (XDGField)finestSolution;
                        XDGField injectedSolution = new XDGField(_finestSolution.Basis, Identification + "-inj-" + iLevel);

                        InjectXDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                        blabla.Add(injectedSolution);
                    } else if(finestSolution is SinglePhaseField) {
                        SinglePhaseField _coarseSolution = (SinglePhaseField)coarseSolution;
                        SinglePhaseField _finestSolution = (SinglePhaseField)finestSolution;
                        SinglePhaseField injectedSolution = new SinglePhaseField(_finestSolution.Basis, Identification + "-inj-" + iLevel);

                        InjectDGField(Fine2CoarseMapS[iLevel], injectedSolution, _coarseSolution);

                        blabla.Add(injectedSolution);
                    } else {
                        throw new NotSupportedException("DG field type '" + finestSolution.GetType().FullName + "' not supported, Identification is '" + finestSolution.Identification + "'");
                    }

                    Console.WriteLine("done.");
                }

                blabla.Add(finestSolution);
                injectedFields.Add(Identification, blabla);
            }

            // compute the errors
            L2Errors = new Dictionary<string, double[]>();
            foreach(string Identification in FieldsToCompare) {

                double[] L2Error = new double[gDataS.Length - 1];

                for(int iLevel = 0; iLevel < gDataS.Length - 1; iLevel++) {
                    Console.WriteLine("Computing L2 error of '{0}' on level {1} ...", Identification, iLevel);

                    DGField Error = injectedFields[Identification].Last().CloneAs();
                    DGField injSol = injectedFields[Identification].ElementAt(iLevel);
                    Error.Acc(-1.0, injSol);

                    L2Error[iLevel] = Error.L2Norm();

                    Console.WriteLine("done (Error is {0:0.####E-00}).", L2Error[iLevel]);
                }

                L2Errors.Add(Identification, L2Error);
            }

            GridRes = gDataS.Take(gDataS.Length - 1).Select(gd => gd.Cells.h_minGlobal).ToArray();

            // convert to table
            /*
            MultidimensionalArray RES = MultidimensionalArray.Create(GridRes.Length, FieldsToCompare.Length + 1);
            RES.SetColumn(0, GridRes);
            for(int ii = 0; ii < FieldsToCompare.Length; ii++) {
                RES.SetColumn(ii + 1, L2Errors[FieldsToCompare[ii]]);
            }


            using(var fs = new FileStream("res.csv", FileMode.Append)) {
                var stw = new StreamWriter(fs);
                stw.Write("GridRes ");
                for(int ii = 0; ii < FieldsToCompare.Length; ii++) {
                    stw.Write(FieldsToCompare[ii]);
                    if(ii + 1 < FieldsToCompare.Length)
                        stw.Write(" ");
                }

                stw.Flush();

                RES.WriteToStream(fs);

            }
            */

        }

        /// <summary>
        /// Injects an XDG field from a coarsr grid to a fine grid.
        /// </summary>
        public static XDGField InjectXDGField(int[] Fine2Coarse, XDGField injected, XDGField cors_field, CellMask subGrd = null) {
            
            var trk = injected.Basis.Tracker;

            foreach (var spc in trk.SpeciesIdS) {
                var grd = trk._Regions.GetSpeciesMask(spc);

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
        /// Injects a DG field from a coarsr grid to a fine grid.
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
                    NodeSet eta = new NodeSet(ctxFine.iGeomCells.GetRefElement(_2), K, D);
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
