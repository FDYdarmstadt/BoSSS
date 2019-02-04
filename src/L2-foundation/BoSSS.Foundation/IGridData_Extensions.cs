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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Utility functions regarding <see cref="IGridData"/>-objects.
    /// </summary>
    public static class IGridData_Extensions {

        class Logical2Geom_Enum : IEnumerator<int> {
            public int j0;
            public int Len;
            public int[][] Log2Geom;

            int jCellLogical = -1;
            int jSub = int.MinValue;
            int[] CurrLog2Geom = null;

            public int Current {
                get {
                    if (jCellLogical < j0 || jCellLogical >= (j0 + Len))
                        throw new InvalidOperationException();
                    if (CurrLog2Geom == null) {
                        // +++++++++++++++++++++++++++++++++++++++++++
                        // logical and geometrical cells are identical
                        // +++++++++++++++++++++++++++++++++++++++++++
                        return jCellLogical;
                    } else {
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++
                        // logical cell consists of multiple geometrical cells
                        // +++++++++++++++++++++++++++++++++++++++++++++++++++

                        return CurrLog2Geom[jSub];
                    }
                }
            }

            object IEnumerator.Current {
                get {
                    return Current;
                }
            }

            public void Dispose() {
            }

            public bool MoveNext() {
                if (Len < 1)
                    return false;
                if (jCellLogical < j0) {
                    jCellLogical = j0;
                    if (Log2Geom != null) {
                        CurrLog2Geom = Log2Geom[jCellLogical];
                        jSub = 0;
                        Debug.Assert(CurrLog2Geom.Length >= 1);
                    }
                    return true;
                }

                if (CurrLog2Geom == null) {
                    // logical cells are identical to geometrical 
                    jCellLogical += 1;
                    return (jCellLogical < (j0 + Len));
                } else {
                    Debug.Assert(CurrLog2Geom.Length >= 1);
                    if (jSub < CurrLog2Geom.Length - 1) {
                        // advance within logical cell
                        jSub++;
                        return true;
                    } else {
                        // advance to next logical cell
                        Debug.Assert(jSub == CurrLog2Geom.Length - 1);
                        jCellLogical += 1;
                        if (jCellLogical >= j0 + Len) {
                            return false;
                        } else {
                            CurrLog2Geom = Log2Geom[jCellLogical];
                            jSub = 0;
                            Debug.Assert(jSub == CurrLog2Geom.Length - 1);
                            return true;
                        }
                    }
                }

            }

            public void Reset() {
                jCellLogical = -1;
                jSub = int.MinValue;
                CurrLog2Geom = null;
            }
        }

        class Logical2Geom_Enumable : IEnumerable<int> {
            public int j0;
            public int Len;
            public int[][] Log2Geom;

            public IEnumerator<int> GetEnumerator() {
                return new Logical2Geom_Enum() { j0 = j0, Len = Len, Log2Geom = Log2Geom };
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return new Logical2Geom_Enum() { j0 = j0, Len = Len, Log2Geom = Log2Geom };
            }
        }

        /// <summary>
        /// Returns an enumeration of geometrical cell indices (<see cref="IGeometricalCellsData"/>)
        /// for a range <paramref name="C"/> of logical cell indices.
        /// </summary>
        /// <param name="g">
        /// The grid object.
        /// </param>
        /// <param name="C">
        /// A range of logical cell indices.
        /// </param>
        /// <returns>
        /// An enumeration of geometrical cell indices.
        /// </returns>
        public static IEnumerable<int> GetGeometricCellIndices(this IGridData g, Chunk C) {
            if (C.i0 < 0)
                throw new ArgumentException();
            if (C.JE > g.iLogicalCells.Count)
                throw new ArgumentException();
            return new Logical2Geom_Enumable() { j0 = C.i0, Len = C.Len, Log2Geom = g.iLogicalCells.AggregateCellToParts };
        }


        /// <summary>
        /// Returns an enumeration of geometrical cell indices (<see cref="IGeometricalCellsData"/>)
        /// for a logical cell index <paramref name="j"/>.
        /// </summary>
        /// <param name="g">
        /// The grid object.
        /// </param>
        /// <param name="j">
        /// A logical cell index.
        /// </param>
        /// <returns>
        /// An enumeration of geometrical cell indices.
        /// </returns>
        public static IEnumerable<int> GetGeometricCellIndices(this IGridData g, int j) {
            if (j < 0)
                throw new ArgumentException();
            if (j >= g.iLogicalCells.Count)
                throw new ArgumentException();


            var enu = new Logical2Geom_Enumable() { j0 = j, Len = 1, Log2Geom = g.iLogicalCells.AggregateCellToParts };
#if DEBUG
            int[] geom2log = g.iGeomCells.GeomCell2LogicalCell;
            if (geom2log == null) {
                int cnt = 0;
                foreach (int jG in enu) {
                    Debug.Assert(jG == j);
                    cnt++;
                }
                Debug.Assert(cnt == 1);
            } else {
                int cnt = 0;
                foreach (int jG in enu) {
                    Debug.Assert(geom2log[jG] == j);
                    cnt++;
                }
                Debug.Assert(g.iLogicalCells.AggregateCellToParts != null && g.iLogicalCells.AggregateCellToParts[j] != null);
                Debug.Assert(g.iLogicalCells.AggregateCellToParts[j].Length == cnt);
            }

#endif
            return enu;
        }


        class Mask2GeomChunks_Enum : IEnumerator<Tuple<int, int>> {

            public int MaxVecLen = 1;
            public CellMask CM;
            public CellInfo ConsecutiveMask = CellInfo.Undefined;

            IEnumerator<Chunk> CMenum = null;
            Tuple<int, int> m_Current = null;
            //List<Tuple<int, int>> bisher = new List<Tuple<int, int>>();

            public Tuple<int, int> Current {
                get {
                    if (CMenum == null)
                        throw new InvalidOperationException();
                    else
                        return m_Current;
                }
            }

            object IEnumerator.Current {
                get {
                    return Current;
                }
            }

            public void Dispose() {
                CMenum.Dispose();
                CMenum = null;
                CM = null; // object can never be used again.
            }

            Chunk Current_Chunk;
            int Current_jLog; // current logical cell index
            int[] Current_L2G; // current logical to geometrical mapping
            int iPart;
            bool reachedEnd = false;

            public bool MoveNext() {
                if (reachedEnd)
                    return false;

                IGridData grd = CM.GridData;
                int[][] L2G = grd.iLogicalCells.AggregateCellToParts;

                //if (Counto == 3 && bisher.Count == 2)
                //    Debugger.Break();

                if (this.CMenum == null) {
                    CMenum = CM.GetEnumerator();
                    bool CMenumRet = CMenum.MoveNext();
                    if (CMenumRet == false) {
                        return false;
                    } else {
                        Current_Chunk = CMenum.Current;
                        Current_jLog = Current_Chunk.i0;
                        if (L2G != null) {
                            Current_L2G = L2G[Current_jLog];
                        } else {
                            Current_L2G = null;
                        }

                        if (Current_L2G != null) {
                            iPart = 0;
                            Debug.Assert(Current_L2G.Length > 0);
                        }
                    }
                }

                Debug.Assert(Current_Chunk.Len > 0);
                Debug.Assert(Current_Chunk.i0 >= 0);
                Debug.Assert(Current_Chunk.JE <= grd.iLogicalCells.NoOfLocalUpdatedCells);
                Debug.Assert(Current_jLog >= Current_Chunk.i0);
                Debug.Assert(Current_jLog <= Current_Chunk.JE);
                Debug.Assert(Current_L2G == null || iPart < Current_L2G.Length);
                Debug.Assert(Current_L2G == null || iPart >= 0);

                int i0;
                if (Current_L2G == null) {
                    i0 = Current_jLog; // logical and geometrical cell indices are identical.
                } else {
                    i0 = Current_L2G[iPart];
                }

                int MaxLen;
                if (i0 < grd.iGeomCells.Count)
                    MaxLen = grd.iGeomCells.GetNoOfSimilarConsecutiveCells(this.ConsecutiveMask, i0, this.MaxVecLen);
                else
                    MaxLen = 1;
                Debug.Assert(MaxLen > 0);

                int iE = i0;
                while (true) {
                    Debug.Assert(Current_Chunk.Len > 0);
                    Debug.Assert(Current_Chunk.i0 >= 0);
                    Debug.Assert(Current_Chunk.JE <= grd.iLogicalCells.NoOfLocalUpdatedCells);
                    Debug.Assert(Current_jLog >= Current_Chunk.i0);
                    Debug.Assert(Current_jLog <= Current_Chunk.JE);
                    Debug.Assert(Current_L2G == null || iPart < Current_L2G.Length);
                    Debug.Assert(Current_L2G == null || iPart >= 0);

                    if (iE - i0 >= MaxLen)
                        break;

                    if (Current_L2G == null) {
                        if (iE < Current_Chunk.JE) {
                            Current_jLog++;
                            // nop
                        } else {
                            // move to next chunk
                            if (CMenum.MoveNext()) {
                                Current_Chunk = CMenum.Current;
                                Current_jLog = Current_Chunk.i0;
                                Debug.Assert(Current_Chunk.i0 > iE, "strange overlap of chunks");

                                if (Current_Chunk.i0 - iE > 0)
                                    // unable to concat chunks
                                    break;

                            } else {
                                // end of cell mask reached; not possible to advance any further
                                reachedEnd = true;
                                break;
                            }

                        }
                    } else {


                        throw new NotImplementedException("todo");
                    }

                    iE++;
                }

                if (iE > i0) {
                    m_Current = new Tuple<int, int>(i0, iE - i0);
                    //bisher.Add(m_Current);
                    Debug.Assert(iE - i0 <= MaxLen);
                    return true;
                } else {
                    m_Current = null;
                    return false;
                }
            }

            public void Reset() {
                CMenum.Dispose();
                CMenum = null;
                reachedEnd = false;
            }
        }

        class Mask2GeomChunks_Enumable : IEnumerable<Tuple<int, int>> {
            public int MaxVecLen = 1;
            public CellMask CM;
            public CellInfo ConsecutiveMask = CellInfo.Undefined;

            public IEnumerator<Tuple<int, int>> GetEnumerator() {
                return new Mask2GeomChunks_Enum() { CM = CM, MaxVecLen = MaxVecLen, ConsecutiveMask = ConsecutiveMask };
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return new Mask2GeomChunks_Enum() { CM = CM, MaxVecLen = MaxVecLen, ConsecutiveMask = ConsecutiveMask };
            }
        }



        /// <summary>
        /// Returns an enumeration of geometrical cell chunks (first index an length) for a given cell mask.
        /// </summary>
        /// <param name="CM"></param>
        /// <param name="MaxVecLen"></param>
        /// <param name="ConsecutiveMask"></param>
        /// <returns></returns>
        public static IEnumerable<Tuple<int, int>> GetGeometricCellChunks(this CellMask CM, int MaxVecLen, CellInfo ConsecutiveMask = CellInfo.Undefined) {
            var ret = new Mask2GeomChunks_Enumable() { CM = CM, MaxVecLen = MaxVecLen, ConsecutiveMask = ConsecutiveMask };
#if DEBUG
            int JG = CM.GridData.iGeomCells.Count;
            BitArray test = new BitArray(JG);
            foreach (var t_i0_len in ret) {
                int j0 = t_i0_len.Item1;
                int Len = t_i0_len.Item2;
                var Flag_j0 = CM.GridData.iGeomCells.InfoFlags[j0] & ConsecutiveMask;
                for (int j = j0; j < j0 + Len; j++) {
                    Debug.Assert(test[j] == false); // each geometric cell is touched only once.
                    test[j] = true;
                    var Flag_j = CM.GridData.iGeomCells.InfoFlags[j0] & ConsecutiveMask;
                    Debug.Assert(Flag_j == Flag_j0); // each geometric cell has the same 
                }
            }

            BitArray CMmask = CM.GetBitMask();
            int[][] L2G = CM.GridData.iLogicalCells.AggregateCellToParts;
            for (int jL = 0; jL < CMmask.Count; jL++) {
                if (L2G == null || L2G[jL] == null) {
                    int jG = jL;
                    if (CMmask[jL] != test[jG]) {
                        var r = new Mask2GeomChunks_Enumable() { CM = CM, MaxVecLen = MaxVecLen, ConsecutiveMask = ConsecutiveMask };

                        foreach (var _t_i0_len in r) {
                            int j0 = _t_i0_len.Item1;
                            int Len = _t_i0_len.Item2;
                            for (int j = j0; j < j0 + Len; j++) {
                                Console.WriteLine(j);
                            }
                        }

                        Debugger.Break();
                    }
                    Debug.Assert(CMmask[jL] == test[jG]);
                } else {
                    foreach (int jG in L2G[jL])
                        Debug.Assert(CMmask[jL] == test[jG]);
                }
            }

#endif
            return ret;
        }

        /// <summary>
        /// Returns an enumeration of geometrical edge indices (<see cref="IGeometricalEdgeData"/>)
        /// for a range <paramref name="C"/> of logical edge indices.
        /// </summary>
        /// <param name="g">
        /// The grid object.
        /// </param>
        /// <param name="C">
        /// A range of logical edge indices.
        /// </param>
        /// <returns>
        /// An enumeration of geometrical cell indices.
        /// </returns>
        public static IEnumerable<int> GetGeometricEdgeIndices(this IGridData g, Chunk C) {
            if (C.i0 < 0)
                throw new ArgumentException();
            if (C.JE > g.iLogicalEdges.Count)
                throw new ArgumentException();
            return new Logical2Geom_Enumable() { j0 = C.i0, Len = C.Len, Log2Geom = g.iLogicalEdges.EdgeToParts };
        }

        /// <summary>
        /// Returns an enumeration of geometrical edge indices (<see cref="IGeometricalCellsData"/>)
        /// for a logical edge index <paramref name="e"/>.
        /// </summary>
        /// <param name="g">
        /// The grid object.
        /// </param>
        /// <param name="e">
        /// A logical edge index.
        /// </param>
        /// <returns>
        /// An enumeration of geometrical cell indices.
        /// </returns>
        public static IEnumerable<int> GetGeometricEdgeIndices(this IGridData g, int e) {
            if (e < 0)
                throw new ArgumentException();
            if (e >= g.iLogicalEdges.Count)
                throw new ArgumentException();
            return new Logical2Geom_Enumable() { j0 = e, Len = 1, Log2Geom = g.iLogicalEdges.EdgeToParts };
        }



        /// <summary>
        /// Find the cell which contains some point <paramref name="pt"/>;
        /// If <paramref name="pt"/> is not within any cell, the cell with its
        /// center nearest to <paramref name="pt"/> is returned and in this
        /// case, <paramref name="IsInside"/> is false.
        /// </summary>
        /// <param name="pt"></param>
        /// <param name="GlobalId">
        /// the Global ID of the found cell;
        /// </param>
        /// <param name="GlobalIndex">
        /// the Global index of the found cell;
        /// </param>
        /// <param name="IsInside">
        /// true, if <paramref name="pt"/> is within the cell identified by
        /// <paramref name="GlobalId"/>;
        /// otherwise, false;
        /// </param>
        /// <param name="OnThisProcess">
        /// If true, the cell <paramref name="GlobalId"/> is located on the current MPI process.
        /// </param>
        /// <param name="CM">
        /// optional cell mask to restrict the search region
        /// </param>
        /// <remarks>
        /// This operation is relatively costly, as it needs to perform a sweep
        /// over all cells, it should not be used for performance-critical
        /// tasks.<br/>
        /// This operation is MPI-collective, the output-values are equal on
        /// all MPI-processors.
        /// </remarks>
        static public void LocatePoint(this IGridData gdat, double[] pt, out long GlobalId, out long GlobalIndex, out bool IsInside, out bool OnThisProcess, CellMask CM = null) {
            using (new FuncTrace()) {
                if (pt.Length != gdat.SpatialDimension)
                    throw new ArgumentException("length must be equal to spatial dimension", "pt");
                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                int MpiRank = gdat.CellPartitioning.MpiRank;
                int MpiSize = gdat.CellPartitioning.MpiSize;

                int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
                int D = gdat.SpatialDimension;

                MultidimensionalArray center = MultidimensionalArray.Create(1, 1, D);  // cell center in global coordinates

                MultidimensionalArray _pt = MultidimensionalArray.Create(1, D);          // point to search for
                for (int d = 0; d < D; d++)
                    _pt[0, d] = pt[d];
                MultidimensionalArray _pt_local = MultidimensionalArray.Create(1, 1, D); // .. in cell-local coordinate
                double[] pt_local = new double[D];


                // sweep over locally updated cells ...
                // ====================================

                int jL_MinDistCel = 0; // logical cell with minimum distance
                int jG_MinDistCel = 0; // geometrical cell with minimum distance
                double MinDist = double.MaxValue;

                int j_Within = -1;
                if (CM == null)
                    CM = CellMask.GetFullMask(gdat);

                foreach (int jL in CM.ItemEnum) {
                    foreach (int j in gdat.GetGeometricCellIndices(jL)) {
                        // compute distance
                        // ================
                        {
                            var smplx = gdat.iGeomCells.GetRefElement(j);
                            gdat.TransformLocal2Global(smplx.Center, j, 1, center, 0);

                            double dist = 0;
                            for (int d = 0; d < D; d++) {
                                double del = pt[d] - center[0, 0, d];
                                dist += del * del;
                            }

                            dist = Math.Sqrt(dist);

                            if (dist < MinDist) {
                                MinDist = dist;
                                jL_MinDistCel = jL;
                                jG_MinDistCel = j;
                            }
                        }

                        // transform back into cell
                        // ========================
                        {
                            try {
                                gdat.TransformGlobal2Local(_pt, _pt_local, j, 1, 0);
                                for (int d = 0; d < D; d++)
                                    pt_local[d] = _pt_local[0, 0, d];

                                var smplx = gdat.iGeomCells.GetRefElement(j);

                                if (smplx.IsWithin(pt_local))
                                    j_Within = j;
                            } catch (ArithmeticException) {
                                // probably outside...

                            }
                        }
                    }
                }

                // ========================================
                // First case: point is inside of some cell
                // ========================================

                int lowestWithinRank = j_Within >= 0 ? MpiRank : int.MaxValue;
                lowestWithinRank = lowestWithinRank.MPIMin(); // lowest rank which found a cell that contains the point:
                //                                               this rank will be the official finder!

                if (lowestWithinRank < MpiSize) {
                    LocatPointHelper(gdat, lowestWithinRank, j_Within, out GlobalId, out GlobalIndex, out OnThisProcess);
                    IsInside = true;
                    return;
                }


                // ========================================
                // Second case: point is outside of all cells
                // ========================================

                double MinDistGlobal = MinDist.MPIMin();
                int lowestMinimumRank = MinDistGlobal == MinDist ? MpiRank : int.MaxValue;
                lowestMinimumRank = lowestMinimumRank.MPIMin(); // find minimum rank on which the global minimum was reached.
                if (lowestMinimumRank < 0 || lowestMinimumRank >= MpiSize)
                    throw new ApplicationException();

                LocatPointHelper(gdat, lowestMinimumRank, jL_MinDistCel, out GlobalId, out GlobalIndex, out OnThisProcess);
                IsInside = false;
                return;

            }
        }

        static void LocatPointHelper(IGridData gdat, int RootRank, int jL_MinDistCel, out long GlobalId, out long GlobalIndex, out bool OnThisProcess) {
            int j0 = gdat.CellPartitioning.i0;
            int MpiSize = gdat.CellPartitioning.MpiSize;
            int MpiRank = gdat.CellPartitioning.MpiRank;

            if (RootRank < 0 || RootRank >= MpiSize)
                throw new ArgumentException();

            if (RootRank == MpiRank) {
                OnThisProcess = true;
                GlobalIndex = jL_MinDistCel + j0;
                GlobalId = gdat.iLogicalCells.GetGlobalID(jL_MinDistCel);
            } else {
                OnThisProcess = false;
                GlobalIndex = long.MinValue;
                GlobalId = long.MinValue;
            }

            unsafe {
                long* buf = stackalloc long[2];
                buf[0] = GlobalId;
                buf[1] = GlobalIndex;

                csMPI.Raw.Bcast((IntPtr)buf, 2, csMPI.Raw._DATATYPE.LONG, RootRank, gdat.CellPartitioning.MPI_Comm);
                GlobalId = buf[0];
                GlobalIndex = buf[1];

            }
            return;

        }



        /// <summary>
        /// Returns an edge mask which contains all boundary cells.
        /// </summary>
        public static EdgeMask GetBoundaryEdgeMask(this IGridData gdat) {

            int E = gdat.iLogicalEdges.Count;
            BitArray boundaryEdges = new BitArray(E);
            BitArray boundaryCellsEdges = new BitArray(E);

            int[,] C2E = gdat.iLogicalEdges.CellIndices;

            // loop over all Edges
            for (int e = 0; e < E; e++) {
                int Cel1 = C2E[e, 0];
                int Cel2 = C2E[e, 1];

                if (Cel2 < 0) {
                    // edge is located on the computational domain boundary
                    boundaryEdges[e] = true;
                }
            }

            return new EdgeMask(gdat, boundaryEdges);
        }



        /// <summary>
        /// Returns a mask which contains all boundary edges
        /// </summary>
        static public EdgeMask GetBoundaryEdges(this IGridData gdat) {

            int E = gdat.iLogicalEdges.Count;
            BitArray boundaryEdges = new BitArray(E);


            int[,] CellIndices = gdat.iLogicalEdges.CellIndices;

            // loop over all Edges
            for (int e = 0; e < E; e++) {
                int Cel1 = CellIndices[e, 0];
                int Cel2 = CellIndices[e, 1];

                if (Cel2 < 0) {
                    // edge is located on the computational domain boundary
                    boundaryEdges[e] = true;

                }
            }

            return new EdgeMask(gdat, boundaryEdges, MaskType.Logical);
        }



        /// <summary>
        /// Finds all neighbor cells for a given cell; 
        /// </summary>
        /// <param name="jCell"></param>
        /// <returns>
        /// - 1st entry: local cell index of neighbor cell 
        /// - 2nd index: edge index of connecting edge 
        /// - 3rd index: whether the other cell is the in-cell (0) or the out-cell (1) of the edge
        /// </returns>
        /// <param name="OmmitPeriodic">
        /// If true, neighborship relations originating from periodic boundary conditions are ignored.
        /// </param>
        /// <param name="g">
        /// Grid object.
        /// </param>
        static public Tuple<int, int, int>[] GetCellNeighboursViaEdges(this IGridData g, int jCell, bool OmmitPeriodic = false) {
            List<Tuple<int, int, int>> ret = new List<Tuple<int, int, int>>();

            int[,] __Edges = g.iLogicalEdges.CellIndices;
            var Cell2Edges_jCell = g.iLogicalCells.Cells2Edges[jCell];
            int K = Cell2Edges_jCell.GetLength(0);
            byte[] __EdgeTags = g.iGeomEdges.EdgeTags;

            if (OmmitPeriodic && g.iLogicalEdges.EdgeToParts != null)
                throw new NotImplementedException("todo"); // problem: edge tags have geometric index.

            for (int k = 0; k < K; k++) {
                int q = Cell2Edges_jCell[k];
                int iEdge = Math.Abs(q) - 1;
                int iOther = q >= 0 ? 1 : 0;
                Debug.Assert(jCell == __Edges[iEdge, iOther == 1 ? 0 : 1]);

                int jCellOut = __Edges[iEdge, iOther];
                if (jCellOut >= 0 && (!OmmitPeriodic || __EdgeTags[iEdge] < GridCommons.FIRST_PERIODIC_BC_TAG)) {
                    ret.Add(new Tuple<int, int, int>(jCellOut, iEdge, iOther));
                }
            }

            return ret.ToArray();
        }

        /// <summary>
        /// Finds all neighbor cells for a given cell.
        /// </summary>
        /// <param name="jCell">
        /// a local cell index in the range of locally updated cells
        /// </param>
        /// <param name="CellNeighBours">
        /// a collection of neighbor cells (local indices);
        /// </param>
        /// <param name="ConectingEntities">
        /// a collection of all edges/all vertices of cell
        /// <paramref name="jCell"/>
        /// </param>
        /// <param name="mode">
        /// Interpretation of neighborhood: sharing a whole edge vs.
        /// sharing at least one vertex.
        /// </param>
        /// <param name="g">
        /// Grid object.
        /// </param>
        static public void GetCellNeighbours(
            this IGridData g,
            int jCell,
            GetCellNeighbours_Mode mode,
            out int[] CellNeighBours,
            out int[] ConectingEntities) {

            if (jCell < 0 || jCell >= g.iLogicalCells.NoOfLocalUpdatedCells)
                throw new ArgumentOutOfRangeException("localIndex", "must be between 0 (including) and 'NoOfLocalUpdatedCells' (excluding)");

            List<int> ret = new List<int>();
            List<int> ret2 = new List<int>();

            switch (mode) {
                case GetCellNeighbours_Mode.ViaEdges: {
                    var R = g.GetCellNeighboursViaEdges(jCell);
                    foreach (var rr in R) {
                        ret.Add(rr.Item1);
                        ret2.Add(rr.Item2);
                    }
                    break;
                }

                case GetCellNeighbours_Mode.ViaVertices: {
                    foreach (int jCellGeom in g.GetGeometricCellIndices(jCell)) {
                        var CellVtx_jCell = g.iGeomCells.CellVertices[jCellGeom];
                        int K = CellVtx_jCell.GetLength(0);
                        var VerticeIndex2Cell = g.iVertices.VerticeToCell;

                        for (int k = 0; k < K; k++) {
                            int iVtx = CellVtx_jCell[k];
                            ret2.Add(iVtx);

                            foreach (int jN in VerticeIndex2Cell[iVtx]) {
                                if (jN != jCell) {
                                    if (!ret.Contains(jN))
                                        ret.Add(jN);
                                }

                            }
                        }
                    }
                    break;
                }

                default:
                throw new NotImplementedException();
            }

            CellNeighBours = ret.ToArray();
            ConectingEntities = ret2.ToArray();
        }

        /// <summary>
        /// For edge number <paramref name="e"/>, whether it is conformal
        /// with cell 1
        /// </summary>
        static public bool IsEdgeConformalWithCell1(this IGeometricalEdgeData ge, int e) {
            return ((ge.Info[e] & EdgeInfo.Cell1_Nonconformal) == 0);
        }

        /// <summary>
        /// For edge number <paramref name="e"/>, whether it is conformal
        /// with cell 2
        /// </summary>
        static public bool IsEdgeConformalWithCell2(this IGeometricalEdgeData ge, int e) {
            return ((ge.Info[e] & EdgeInfo.Cell2_Nonconformal) == 0);
        }


        /// <summary>
        /// For edge number <paramref name="iedge"/>, whether it is
        /// conformal with the adjacent cells
        /// </summary>
        /// <param name="iedge"></param>
        /// <param name="InOrOut">
        /// 0 for IN-cell/first cell, see <see cref="IsEdgeConformalWithCell1"/><br/>
        /// 1 for OUT-cell/second cell, see <see cref="IsEdgeConformalWithCell2"/>
        /// </param>
        /// <returns></returns>
        static public bool IsEdgeConformal(this IGeometricalEdgeData ge, int iedge, int InOrOut) {
            switch (InOrOut) {
                case 0:
                return ge.IsEdgeConformalWithCell1(iedge);
                case 1:
                return ge.IsEdgeConformalWithCell2(iedge);
                default:
                throw new ArgumentOutOfRangeException();
            }

        }

        /// <summary>
        /// true, if edge <paramref name="e"/> is a boundary-edge.
        /// </summary>
        static public bool IsEdgeBoundaryEdge(this IGeometricalEdgeData ge, int e) {
            bool R = (ge.Info[e] & EdgeInfo.Boundary) != 0;
            return R;
        }


        ///// <summary>
        ///// true, if edge <paramref name="e"/> is affine-linear, false if not
        ///// </summary>
        //static public bool IsEdgeAffineLinear(this IGeometricalEdgeData ge, int e) {
        //    return (ge.Info[e] & EdgeInfo.EdgeIsAffineLinear) != 0;
        //}


        /// <summary>
        /// Returns a mask containing all cells which lie at the domain boundary
        /// </summary>
        static public CellMask GetBoundaryCells(this IGridData gdat) {
            int J = gdat.iLogicalCells.NoOfLocalUpdatedCells;
            BitArray boundaryCells = new BitArray(J);

            int E = gdat.iLogicalEdges.Count;


            int[,] CellIndices = gdat.iLogicalEdges.CellIndices;

            // loop over all Edges
            for (int e = 0; e < E; e++) {
                int Cel1 = CellIndices[e, 0];
                int Cel2 = CellIndices[e, 1];

                if (Cel2 < 0) {
                    // edge is located on the computational domain boundary


                    boundaryCells[Cel1] = true;
                }



            }

            return new CellMask(gdat, boundaryCells, MaskType.Logical);


        }


        /// <summary>
        /// computes a global time-step length ("delta t") according to the 
        /// Courant-Friedrichs-Lax - criterion, based on a velocity
        /// vector (<paramref name="velvect"/>) and the cell size
        /// this.<see cref="Cells"/>.<see cref="CellData.h_min"/>;
        /// </summary>
        /// <param name="velvect">
        /// components of a velocity vector
        /// </param>
        /// <param name="max">
        /// an upper maximum for the return value; This is useful if the velocity
        /// defined by <paramref name="velvect"/> is 0 or very small everywhere;
        /// </param>
        /// <param name="cm">
        /// optional restriction of domain.
        /// </param>
        /// <returns>
        /// the minimum (over all cells j in all processes) of <see cref="CellData.h_min"/>[j]
        /// over v, where v is the Euclidean norm of a vector build from 
        /// <paramref name="velvect"/>;
        /// This vector is evaluated at cell center and all cell vertices.
        /// The return value is the same on all processes;
        /// </returns>
        static public double ComputeCFLTime<T>(this IGridData __gdat, IEnumerable<T> velvect, double max, CellMask cm = null)
            where T : DGField //
        {
            using (var tr = new FuncTrace()) {
                GridData gdat = (GridData)__gdat;

                ilPSP.MPICollectiveWatchDog.Watch(MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                T[] _velvect = velvect.ToArray();

                if (cm == null)
                    cm = CellMask.GetFullMask(gdat);

                int D = gdat.SpatialDimension;
                var KrefS = gdat.Grid.RefElements;

                // find cfl number on this processor
                // ---------------------------------

                var m_CFL_EvalPoints = new NodeSet[KrefS.Length];
                for (int i = 0; i < KrefS.Length; i++) {
                    var Kref = KrefS[i];
                    int N = Kref.NoOfVertices + 1;

                    MultidimensionalArray vert = MultidimensionalArray.Create(N, D);
                    vert.SetSubArray(Kref.Vertices, new int[] { 0, 0 }, new int[] { N - 2, D - 1 });

                    m_CFL_EvalPoints[i] = new NodeSet(Kref, vert);
                }
                


                // evaluators an memory for result
                int VecMax = 1000;
                DGField[] evalers = new DGField[_velvect.Length];
                MultidimensionalArray[] fieldValues = new MultidimensionalArray[_velvect.Length];
                for (int i = 0; i < _velvect.Length; i++) {
                    evalers[i] = _velvect[i];
                    fieldValues[i] = MultidimensionalArray.Create(VecMax, m_CFL_EvalPoints[0].NoOfNodes);
                }

                var h_min = gdat.Cells.h_min;
                int K = _velvect.Length;
                double cflhere = max;

                //for (int j = 0; j < J; j += VectorSize) {
                foreach (Chunk chk in cm) {
                    int VectorSize = VecMax;
                    for (int j = chk.i0; j < chk.JE; j += VectorSize) {
                        if (j + VectorSize > chk.JE + 1)
                            VectorSize = chk.JE - j;
                        VectorSize = gdat.Cells.GetNoOfSimilarConsecutiveCells(CellInfo.RefElementIndex_Mask, j, VectorSize);


                        int iKref = gdat.Cells.GetRefElementIndex(j);
                        int N = m_CFL_EvalPoints[iKref].GetLength(0);

                        if (fieldValues[0].GetLength(0) != VectorSize) {
                            for (int i = 0; i < _velvect.Length; i++) {
                                fieldValues[i].Allocate(VectorSize, N);
                            }
                        }

                        for (int k = 0; k < K; k++)
                            evalers[k].Evaluate(j, VectorSize, m_CFL_EvalPoints[iKref], fieldValues[k], 0, 0.0);

                        // loop over cells ...
                        for (int jj = j; jj < j + VectorSize; jj++) {

                            // loop over nodes ...
                            for (int n = 0; n < N; n++) {
                                double velabs = 0;

                                // loop over velocity components ...
                                for (int k = 0; k < K; k++) {
                                    double v = fieldValues[k][jj - j, n];
                                    velabs += v * v;
                                }

                                velabs = Math.Sqrt(velabs);

                                double cfl = h_min[jj] / velabs;
                                cflhere = Math.Min(cfl, cflhere);
                            }
                        }
                    }

                }

                // find the minimum over all processes via MPI and return
                // ------------------------------------------------------
                double cfltotal;
                unsafe {
                    csMPI.Raw.Allreduce((IntPtr)(&cflhere), (IntPtr)(&cfltotal), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);
                }
                tr.Info("computed CFL timestep: " + cfltotal);
                return cfltotal;

            }
        }
    }

    /// <summary>
    /// used by <see cref="IGridData_Extensions.GetCellNeighbours(IGridData, int, GetCellNeighbours_Mode, out int[], out int[])"/>.
    /// </summary>
    public enum GetCellNeighbours_Mode {

        /// <summary>
        /// in this mode, any cell that shares an edge is considered a neighbor
        /// </summary>
        ViaEdges,

        /// <summary>
        /// in this mode, any cell that shares a vertex is considered a neighbor
        /// </summary>
        ViaVertices
    }

    

}
