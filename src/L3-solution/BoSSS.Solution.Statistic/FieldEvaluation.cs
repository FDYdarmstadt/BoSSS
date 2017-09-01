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
using BoSSS.Platform;
using System.Collections;
using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic {
    /// <summary>
    /// evaluation of a DG field at arbitrary points
    /// </summary>
    public class FieldEvaluation {

        GridData m_Context;

        /// <summary>
        /// ctor
        /// </summary>
        public FieldEvaluation(GridData ctx) {
            m_Context = ctx;
            CellLoc = new CellLocalization(ctx);
        }

        CellLocalization CellLoc;


        /// <summary>
        /// bounding box of the grid
        /// </summary>
        BoundingBox GridBoundingBox {
            get { return CellLoc.GridBB; }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="Flds">
        /// a list of <em>M</em> DG fields.
        /// </param>
        /// <param name="Points">
        /// 2-dimensional: <br/>
        ///  - 1st index: point index, from 0 (including) to <em>N</em> (excluding)<br/>
        ///  - 2nd index: spatial dimension/direction
        /// </param>
        /// <param name="Result">
        /// result of the evaluation of the DG fields <paramref name="Flds"/> at points <paramref name="Points"/>.<br/>
        ///  - 1st index: point index, from 0 (including) to <em>N</em> (excluding)<br/>
        ///  - 2nd index: DG field, from 0 (including) to <em>M</em> (excluding).<br/>
        /// the result will be accumulated.
        /// </param>
        /// <param name="UnlocatedPoints">
        /// Optional, can be null; otherwise, an array of length <em>N</em>.
        /// On exit, in the latter case, true for every point that is not within the grid.
        /// </param>
        /// <param name="LocalCellIndices">
        /// Optional, can be null; otherwise, an array of length <em>N</em>.
        /// On exit, in the latter case, the <em>n</em>-th entry contains the local index of the cell  
        /// where the <em>n</em>-th point is found.
        /// </param>
        /// <param name="beta">
        /// pre-scaling of <paramref name="Result"/> on entry, i.e. before accumulation
        /// </param>
        /// <param name="alpha">
        /// scaling of <paramref name="Flds"/> for evaluation
        /// </param>
        /// <returns>
        /// the number of points that are not within the grid;
        /// </returns>
        public int Evaluate(double alpha, IEnumerable<DGField> Flds, MultidimensionalArray Points, double beta, MultidimensionalArray Result, BitArray UnlocatedPoints = null, int[] LocalCellIndices = null) {

            // init / check args
            // =================

            int D = m_Context.Grid.SpatialDimension;
            int J = m_Context.Grid.NoOfUpdateCells;

            if (Points.Dimension != 2)
                throw new ArgumentException("must be 2-dimensional", "Points");
            int N = Points.GetLength(0);
            if (Result.Dimension != 2)
                throw new ArgumentException("must be 2-dimensional", "Result");
            if (Result.GetLength(0) != N)
                throw new ArgumentException("mismatch in 1st dimension", "Points,Result");
            DGField[] Fields = Flds.ToArray();
            int M = Fields.Length;
            if (Result.GetLength(1) != Fields.Length)
                throw new ArgumentException("2nd dimension of 'Result' and length of 'Flds' must be equal");

            if (UnlocatedPoints == null)
                UnlocatedPoints = new BitArray(N, false);
            else
                UnlocatedPoints.SetAll(false);

            if (LocalCellIndices != null) {
                if (LocalCellIndices.Length != N)
                    throw new ArgumentException("2nd dimension of 'Result' and length of 'LocalCellIndices' must be equal");
                for (int n = 0; n < N; n++)
                    LocalCellIndices[n] = int.MinValue;
            }

            if (UnlocatedPoints.Length != N)
                throw new ArgumentException("2nd dimension of 'Result' and length of 'UnlocatedPoints' must be equal");

            double[] pt = new double[D];

            Result.Scale(beta);

            // build point localization
            // ========================

            // filter points that are outside the grid bounding box
            int[] Perm2 = new int[N];
            {
                int cnt = 0;
                var bb = CellLoc.GridBB;
                for (int n = 0; n < N; n++) {
                    for (int d = 0; d < D; d++)
                        pt[d] = Points[n, d];

                    if (!bb.Contains(pt)) {
                        UnlocatedPoints[n] = true;
                    } else {
                        Perm2[cnt] = n;
                        cnt++;
                    }
                }

                Array.Resize(ref Perm2, cnt);
            }

            if (Perm2.Length <= 0) {
                UnlocatedPoints.SetAll(true);
                return N; // we are done
            }

            int _cnt = 0;
            MultidimensionalArray pts = MultidimensionalArray.Create(Perm2.Length, D);
            for (int n = 0; n < N; n++) {
                if (!UnlocatedPoints[n]) {
                    for (int d = 0; d < D; d++)
                        pts[_cnt, d] = Points[n, d];
                    _cnt++;
                }
            }
            int NU = N - Perm2.Length;

            UnlocatedPoints.SetAll(true);
            N = Perm2.Length;
            int[] Perm = new int[N];
            int NoOfUnassignedNodes = N;
            PointLocalization pl = new PointLocalization(pts, CellLoc.GridBB, Perm);
            pts = null; // not required anymore

            // localize Points / evaluate at
            // =============================

            MultidimensionalArray vertGlobalSupect = new MultidimensionalArray(2);
            MultidimensionalArray vertLocalSuspect = new MultidimensionalArray(3);
            NodeSet vertLocal = null;

            MultidimensionalArray[] fieldVal = new MultidimensionalArray[M];
            for (int m = 0; m < M; m++) {
                fieldVal[m] = new MultidimensionalArray(2);
            }

            //Console.WriteLine("debgcode akt.");
            //int[] bbsfound = new int[Points.GetLength(0)];
            BoundingBox CellTreeBB = new BoundingBox(D);
            BoundingBox CellBB = new BoundingBox(D);


            // loop over cells ...
            for (int j = 0; j < J; j++) {


                // code of the cell: all bounding boxes in the tree that 
                // share a point with the cell
                GeomBinTreeBranchCode bbcode; int bbBits;
                {
                    BoundingBoxCode __b = CellLoc.GetCellBoundingBoxCode(j);
                    bbcode = __b.Branch;
                    bbBits = (int)__b.SignificantBits;
                    CellLoc.GridBB.SubBoxFromCode(CellTreeBB, __b);
                }

                // cell bounding box is smaller than the bounding box in the tree
                m_Context.Cells.GetCellBoundingBox(j, CellBB);
                if (!CellTreeBB.Contains(CellBB)) // test
                    throw new ApplicationException("internal error: should not happen");
                CellBB.ExtendByFactor(0.001); // safety factor

                // determine all points in cell
                int iP0, Len;
                pl.GetPointsInBranch(bbcode, bbBits, out iP0, out Len);
                if (Len <= 0)
                    // no points in cell j
                    continue;

                // transform points to cell-local coordinates
                vertGlobalSupect.Allocate(Len, D);
                vertLocalSuspect.Allocate(1, Len, D);
                for (int n = 0; n < Len; n++)
                    for (int d = 0; d < D; d++)
                        vertGlobalSupect[n, d] = pl.Points[n + iP0, d];
                CellLoc.GrdDat.TransformGlobal2Local(vertGlobalSupect, vertLocalSuspect, j, 1, 0);

                // test
                //for (int n = 0; n < Len; n++) {
                //    for (int d = 0; d < D; d++)
                //        pt[d] = vertGlobalSupect[n, d];

                //    if (!CellBB.Contains(pt))
                //        Console.WriteLine("Huramnet");
                //}

                var splx = m_Context.Cells.GetRefElement(j);

                // test whether the points in the bounding box of cell j
                // are also in cell j
                bool[] tatsaechlich = new bool[Len];
                int Z = 0;
                for (int n = 0; n < Len; n++) {
                    int nPt = Perm2[Perm[n + iP0]];

                    if (!UnlocatedPoints[nPt]) {
                        // point was already located in/assigned to another cell
                        continue;
                    }

                    for (int d = 0; d < D; d++)
                        pt[d] = vertGlobalSupect[n, d];
                    if (!CellBB.Contains(pt)) {
                        continue; // cell bounding box is usually smaller than the bounding box in the tree ('bbcode')
                    }

                    for (int d = 0; d < D; d++)
                        pt[d] = vertLocalSuspect[0, n, d];
                    if (splx.IsWithin(pt, 1.0e-8)) {
                        UnlocatedPoints[nPt] = false;
                        NoOfUnassignedNodes--;

                        if (LocalCellIndices != null)
                            LocalCellIndices[nPt] = j;

                        Z++;
                        tatsaechlich[n] = true;
                    }
                }
                if (Z <= 0)
                    // no points in bb
                    continue;

                // collect all vertices that are really in cell j
                vertLocal = new NodeSet(m_Context.Cells.GetRefElement(j), Z, D);
                int z = 0;
                for (int n = 0; n < Len; n++) {
                    if (!tatsaechlich[n])
                        continue;

                    for (int d = 0; d < D; d++)
                        vertLocal[z, d] = vertLocalSuspect[0, n, d];
                    z++;
                }
                vertLocal.LockForever();

                // evaluate Velocity Field there
                for (int m = 0; m < M; m++) {
                    fieldVal[m].Allocate(1, Z);
                }

                for (int m = 0; m < M; m++) {
                    Fields[m].Evaluate(j, 1, vertLocal, fieldVal[m]);
                }

                
                // store result of evaluation
                z = 0;
                for (int n = 0; n < Len; n++) {
                    int nPt = Perm2[Perm[n + iP0]];

                    if (!tatsaechlich[n])
                        continue;

                    for (int m = 0; m < M; m++) {
                        Result[nPt, m] += alpha * fieldVal[m][0, z];
                    }

                    z++;
                }
            }

            // return
            // ======
            return (NoOfUnassignedNodes + NU);
        }
    }
}
