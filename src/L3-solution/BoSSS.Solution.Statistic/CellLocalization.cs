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
using BoSSS.Foundation.Grid;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Platform;
using ilPSP.Tracing;
using System.Collections;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.Statistic {

    /// <summary>
    /// This class localizes every cell of the grid in a geometrical binary tree
    /// (<see cref="CellMinCode"/> and <see cref="CellMaxCode"/>).
    /// </summary>
    public class CellLocalization {

        GridData m_GrdDat;

        /// <summary>
        /// <see cref="GridData"/> - object which was provided with the constructor
        /// </summary>
        public GridData GrdDat { get { return m_GrdDat; } }
        
        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="grdDat"></param>
        public CellLocalization(GridData grdDat) {
            m_GrdDat = grdDat;

            int J = grdDat.Cells.NoOfLocalUpdatedCells;
            int D = grdDat.SpatialDimension;

            CellMinCode = new GeomBinTreeBranchCode[J];
            CellMaxCode = new GeomBinTreeBranchCode[J];

            // 1st pass: find bounding box of grid
            // ===================================
            GridBB = m_GrdDat.LocalBoundingBox;
            GridBB.ExtendByFactor(0.01);
            
            // 2nd pass: localize all cells
            // ============================
            CellMinCode = new GeomBinTreeBranchCode[J];
            CellMaxCode = new GeomBinTreeBranchCode[J];
            BoundingBox CellBB = new BoundingBox(D);
            BoundingBox TreeBB = new BoundingBox(D);
            var CellBB_points = new double[][] { CellBB.Min, CellBB.Max};
            for (int j = 0; j < J; j++) {
                m_GrdDat.Cells.GetCellBoundingBox(j, CellBB);
                //grdDat.TransformLocal2Global(verticesLoc, vericesGlob, j, 1, 0);

                uint mincode = uint.MaxValue;
                uint maxcode = 0;
                foreach (var _pt in CellBB_points) {
                    GeomBinTreeBranchCode cd = GeomBinTreeBranchCode.CreateFormPoint(GridBB, _pt);
                    mincode = Math.Min(mincode, cd.Code);
                    maxcode = Math.Max(maxcode, cd.Code);
                }

                CellMinCode[j].Code = mincode;
                CellMaxCode[j].Code = maxcode;

                // test
                BoundingBoxCode currentlyGen = GetCellBoundingBoxCode(j);
                GridBB.SubBoxFromCode(TreeBB, currentlyGen);
                
                if (!TreeBB.Contains(CellBB))
                    throw new ApplicationException("internal error - should not occur");
            }
        }

        /// <summary>
        /// bounding box of grid (on current MPI-process)
        /// </summary>
        public BoundingBox GridBB;
        
        /// <summary>
        /// 
        /// <br/>
        /// index: local cell index
        /// </summary>
        public GeomBinTreeBranchCode[] CellMinCode;
        
        /// <summary>
        /// index: local cell index
        /// </summary>
        public GeomBinTreeBranchCode[] CellMaxCode;

        /// <summary>
        /// finds the smallest bounding box (in the geometric binary tree) which contains cell <paramref name="jCell"/>;
        /// </summary>
        /// <param name="jCell"></param>
        /// <return>
        /// the bounding box information (branch and significant bits) for the bounding box of cell <paramref name="jCell"/>;
        /// </return>
        public BoundingBoxCode GetCellBoundingBoxCode(int jCell) {
            BoundingBoxCode ret;
            int sigbits;
            ret.Branch = GeomBinTreeBranchCode.Combine(CellMinCode[jCell], CellMaxCode[jCell], out sigbits);
            ret.SignificantBits = (uint)sigbits;
            return ret;
        }


        /// <summary>
        /// for a cloud of points, this method finds the cells which contain the points
        /// </summary>
        /// <param name="_pts">
        /// Input: a cloud of points; 1st index: point index, 2nd index: spatial dimension;
        /// </param>
        /// <param name="LocalCellIdx">
        /// Output: for each point in <paramref name="_pts"/>, the the (local) index of the cell which contains the specific point
        /// </param>
        /// <param name="NoOfUnassigned">
        /// on exit, the number of points which cannot be assigned to one cell
        /// </param>
        public void LocalizePointsWithinGrid<T>(MultidimensionalArray _pts, T LocalCellIdx, out int NoOfUnassigned) where T : IList<int> {
            using (new FuncTrace()) {
                // init and check
                // ==============

                NoOfUnassigned = 0;

                GridData grd = this.GrdDat;
                var CellLoc = this;
                var splxS = grd.Grid.RefElements;

                int D = grd.SpatialDimension;
                int N = _pts.GetLength(0);
                if (_pts.Dimension != 2)
                    throw new ArgumentException();
                if (_pts.GetLength(0) != LocalCellIdx.Count)
                    throw new ArgumentException();
                if (_pts.GetLength(1) != D)
                    throw new ArgumentException();

                var UnlocatedPoints = new BitArray(N, false);

                int J = grd.Cells.NoOfLocalUpdatedCells;

                double[] pt = new double[D];


                // build point localization
                // ========================

                // filter points that are outside the grid bounding box
                int[] Perm2 = new int[N];
                {
                    int cnt = 0;
                    var bb = CellLoc.GridBB;
                    for (int n = 0; n < N; n++) {
                        _pts.ExtractVector(pt, 1, 0, D, n, 0);
                        //for (int d = 0; d < D; d++)
                        //    pt[d] = Points[n, d];

                        if (!bb.Contains(pt)) {
                            UnlocatedPoints[n] = true;
                            NoOfUnassigned++;
                        } else {
                            Perm2[cnt] = n;
                            cnt++;
                        }
                    }

                    Array.Resize(ref Perm2, cnt);
                }

                if (Perm2.Length <= 0) {
                    // all points are outside the bounding box of the grid !
                    LocalCellIdx.SetAll(int.MinValue);
                    return; // we are done
                }

                int _cnt = 0;
                MultidimensionalArray pts = MultidimensionalArray.Create(Perm2.Length, D);
                for (int n = 0; n < N; n++) {
                    if (!UnlocatedPoints[n]) {
                        for (int d = 0; d < D; d++)
                            pts[_cnt, d] = _pts[n, d];
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
                    grd.Cells.GetCellBoundingBox(j, CellBB);
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

                    var Kref = grd.Cells.GetRefElement(j);

                    // test whether the points in the bounding box of cell j
                    // are also in cell j
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
                        if (Kref.IsWithin(pt)) {
                            UnlocatedPoints[nPt] = false;
                            NoOfUnassignedNodes--;

                            LocalCellIdx[nPt] = j;
                        }
                    }
                }

                // return
                NoOfUnassigned += NoOfUnassignedNodes;
            }
        }
    }
}
