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
using ilPSP.Tracing;
using MPI.Wrappers;
using System.Diagnostics;
using ilPSP.Utils;

namespace BoSSS.Solution.Statistic {
    /// <summary>
    /// evaluation of a DG field at arbitrary points
    /// </summary>
    public class FieldEvaluation {

        GridData m_Context;

        /// <summary>
        /// ctor
        /// </summary>
        public FieldEvaluation(IGridData ctx) {
            m_Context = (Foundation.Grid.Classic.GridData)ctx;
            CellLoc = new CellLocalization(m_Context);
        }

        CellLocalization CellLoc;


        /// <summary>
        /// bounding box of the grid
        /// </summary>
        BoundingBox GridBoundingBox {
            get { return CellLoc.GridBB; }
        }


        /// <summary>
        /// Evaluation at a singe point
        /// </summary>
        public double Evaluate(Vector Point, DGField F) {
            var gdat = CellLoc.GrdDat;
            if(!object.ReferenceEquals(gdat, F.GridDat)) {
                throw new ArgumentException("Grid mismatch.");
            }
            int D = gdat.SpatialDimension;
            if(Point.Dim != D) {
                throw new ArgumentException("Spatial dimension mismatch.");
            }

            MultidimensionalArray _Point = MultidimensionalArray.Create(1, Point.Dim);
            _Point.SetRowPt(0, Point);

            int[] cellIdx = new int[1];
            CellLoc.LocalizePointsWithinGrid(_Point, cellIdx, out int NoOfUnassi);
            if(NoOfUnassi > 0) {
                throw new ArithmeticException("unable to locate point " + Point);
            }
            int j = cellIdx[0];

            NodeSet ns = new NodeSet(gdat.iGeomCells.GetRefElement(j), 1, D, false);
            bool[] NewtonConvergence = new bool[1];
            gdat.TransformGlobal2Local(_Point, ns, j, NewtonConvergence);
            if(NewtonConvergence[0] == false) {
                throw new ArithmeticException("unable to transform point " + Point + " to cell " + j + "; Newton method did not converged.");
            }
            ns.LockForever();

            var R = MultidimensionalArray.Create(1, 1);
            F.Evaluate(j, 1, ns, R);
            return R[0, 0];
        }

        /// <summary>
        /// Like <see cref="Evaluate(double, IEnumerable{DGField}, MultidimensionalArray, double, MultidimensionalArray, BitArray, int[])"/>,
        /// but with MPI-Exchange
        /// </summary>
        public int EvaluateParallel(double alpha, IEnumerable<DGField> Flds, MultidimensionalArray Points, double beta, MultidimensionalArray Result, BitArray UnlocatedPoints = null) {
            using(new FuncTrace()) {
                MPICollectiveWatchDog.Watch();
                

                int L = Points != null ? Points.NoOfRows : 0;
                int MPIsize = m_Context.MpiSize;
                int D = m_Context.SpatialDimension;

                if(UnlocatedPoints != null) {
                    if(UnlocatedPoints.Length != L)
                        throw new ArgumentException("Length mismatch");
                }

                // evaluate locally
                // ================
                var unlocated = new System.Collections.BitArray(L);
                int NoOfUnlocated;
                if(L > 0)
                    NoOfUnlocated = this.Evaluate(alpha, Flds, Points, beta, Result, unlocated);
                else
                    NoOfUnlocated = 0;

                // return, if there are no unlocalized points
                // ==========================================

                int TotNoOfUnlocated = NoOfUnlocated.MPISum();
                if(TotNoOfUnlocated <= 0) {
                    if(UnlocatedPoints != null)
                        UnlocatedPoints.SetAll(false);

                    return 0;
                }

                // copy unlocalized to separate array
                // ==================================
                double[,] localUnlocated = new double[NoOfUnlocated, D]; // MultidimensionalArray does not allow zero length -- so use double[,] instead
                int[] IndexToOrgIndex = new int[NoOfUnlocated];
                int u = 0;
                for(int i = 0; i < L; i++) {
                    if(unlocated[i]) {
                        localUnlocated.SetRowPt(u, Points.GetRowPt(i));
                        IndexToOrgIndex[u] = i;
                        u++;
                    }
                }
                Debug.Assert(u == NoOfUnlocated);

                // collect on all ranks -- this won't scale well, but it may work
                // ==============================================================
                MultidimensionalArray globalUnlocated;
                int[] WhoIsInterestedIn; // index: point index, corresponds with 'globalUnlocated' rows; content: rank which needs the result
                int[] OriginalIndex; // index: detto; content: index which the point had on the processor that sent it.
                double[][,] __globalUnlocated;
                int LL;
                {
                    __globalUnlocated = localUnlocated.MPIAllGatherO();
                    Debug.Assert(__globalUnlocated.Length == MPIsize);
                    Debug.Assert(__globalUnlocated.Select(aa => aa.GetLength(0)).Sum() == TotNoOfUnlocated);
                    Debug.Assert(__globalUnlocated[m_Context.MpiRank].GetLength(0) == NoOfUnlocated);
                    LL = TotNoOfUnlocated - NoOfUnlocated;
                    if(LL > 0)
                        globalUnlocated = MultidimensionalArray.Create(LL, D);
                    else
                        globalUnlocated = null;
                    WhoIsInterestedIn = new int[LL];
                    OriginalIndex = new int[LL];
                    int g = 0;
                    for(int r = 0; r < MPIsize; r++) { // concat all point arrays from all processors
                        if(r == m_Context.MpiRank)
                            continue;

                        double[,] __globalPart = __globalUnlocated[r];
                        int Lr = __globalPart.GetLength(0);
                        if(Lr > 0)
                            globalUnlocated.ExtractSubArrayShallow(new[] { g, 0 }, new[] { g + Lr - 1, D - 1 }).Acc2DArray(1.0, __globalPart);
                        for(int i = 0; i < Lr; i++) {
                            WhoIsInterestedIn[i + g] = r;
                            OriginalIndex[i + g] = i;
                        }
                        g += Lr;
                    }
                }

                // try to evaluate the so-far-unlocalized points
                // ---------------------------------------------
           
                var unlocated2 = new System.Collections.BitArray(LL);
                var Result2 = LL > 0 ? MultidimensionalArray.Create(LL, Flds.Count()) : null;
                int NoOfUnlocated2 = LL > 0 ? this.Evaluate(1.0, Flds, globalUnlocated, 0.0, Result2, unlocated2) : 0;

                // backward MPI sending
                // --------------------
                IDictionary<int, EvaluateParallelHelper> resultFromOtherProcs;
                {
                    var backSend = new Dictionary<int, EvaluateParallelHelper>();
                    for(int ll = 0; ll < LL; ll++) {
                        if(!unlocated2[ll]) {
                            int iTarget = WhoIsInterestedIn[ll];
                            Debug.Assert(iTarget != m_Context.MpiRank);
                            if(!backSend.TryGetValue(iTarget, out EvaluateParallelHelper eph)) {
                                eph = new EvaluateParallelHelper();
                                backSend.Add(iTarget, eph);
                            }

                            eph.OriginalIndices.Add(OriginalIndex[ll]);
                            eph.Results.Add(Result2.GetRow(ll));
                        }
                    }

                    resultFromOtherProcs = SerialisationMessenger.ExchangeData(backSend);
                }

                // fill the results from other processors
                // ======================================
                foreach(var res in resultFromOtherProcs.Values) {
                    int K = res.OriginalIndices.Count();
                    Debug.Assert(res.OriginalIndices.Count == res.Results.Count);

                    for(int k = 0; k < K; k++) {
                        int iOrg = IndexToOrgIndex[res.OriginalIndices[k]];
                        if(unlocated[iOrg] == true)
                            NoOfUnlocated--;
                        unlocated[iOrg] = false;

                        Result.AccRow(iOrg, alpha, res.Results[k]);
                    }
                }

                // Return
                // ======
                if(UnlocatedPoints != null) {
                    for(int l = 0; l < L; l++)
                        UnlocatedPoints[l] = unlocated[l];
                }

                return NoOfUnlocated;
            }

        }

        [Serializable]
        class EvaluateParallelHelper {
            public List<int> OriginalIndices = new List<int>();
            public List<double[]> Results = new List<double[]>();
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
        /// result of the evaluation of the DG fields <paramref name="Flds"/> at points <paramref name="Points"/>.
        ///  - 1st index: point index, from 0 (including) to <em>N</em> (excluding)
        ///  - 2nd index: DG field, from 0 (including) to <em>M</em> (excluding).
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
            for(int m = 0; m < M; m++) {
                if(!object.ReferenceEquals(this.m_Context, Fields[m].GridDat))
                    throw new ArgumentException($"Mismatch in grid: field #{m} is defined on a different grid.");
            }

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

            MultidimensionalArray vertGlobalSupect1 = new MultidimensionalArray(2);
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
                vertGlobalSupect1.Allocate(Len, D);
                double[] tempPt = new double[D];
                int tempLen = 0;
                int[] nPts = new int[Len];

                for (int n = 0; n < Len; n++) {
                    int nPt = Perm2[Perm[n + iP0]];
                    if (!UnlocatedPoints[nPt])
                        continue;

                    for (int d = 0; d < D; d++)
                        tempPt[d] = pl.Points[n + iP0, d];

                    if (CellBB.Contains(tempPt)) {
                        for (int d = 0; d < D; d++)
                            vertGlobalSupect1[tempLen, d] = tempPt[d];
                        nPts[tempLen] = nPt;
                        tempLen++;
                    }
                }
                Len = tempLen;
                if (Len <= 0)
                    // no points in cell j
                    continue;

                vertGlobalSupect.Allocate(Len, D);
                vertGlobalSupect.Set(vertGlobalSupect1.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { Len - 1, D - 1 }));

                vertLocalSuspect.Allocate(1, Len, D);
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
                    int nPt = nPts[n];

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
                vertLocal = new NodeSet(m_Context.Cells.GetRefElement(j), Z, D, false);
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
                    //int nPt = Perm2[Perm[n + iP0]];
                    int nPt = nPts[n];

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
