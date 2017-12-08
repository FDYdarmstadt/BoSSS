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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using BoSSS.Platform;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Classic {

    partial class GridData {

        /// <summary>
        /// See <see cref="Vertices"/>
        /// </summary>
        private VertexData m_VerticeData;

        /// <summary>
        /// Information about the vertices of the grid elements, see
        /// <see cref="VertexData"/>
        /// </summary>
        public VertexData Vertices {
            get {
                return m_VerticeData;
            }
        }

        /// <summary>
        /// Information about the vertices of the grid elements, see
        /// <see cref="IVertexData"/>.
        /// </summary>
        public IVertexData iVertices {
            get {
                return m_VerticeData;
            }
        }

        /// <summary>
        /// Data about the vertices
        /// </summary>
        public class VertexData : IVertexData {

            /// <summary>
            /// ctor
            /// </summary>
        internal VertexData(GridData owner) {
                m_owner = owner;
            }

            /// <summary>
            /// pointer to owner object
            /// </summary>
            private GridData m_owner;

            /// <summary>
            /// all vertices/nodes of the gird;
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>1st index: vertex index</item>
            ///   <item>2nd index: spatial dimension</item>
            /// </list>
            /// The vertices are sorted in a specific order:
            /// <list type="bullet">
            ///   <item>
            ///     First, all pure local vertices, see
            ///     <see cref="NoOfPurelyLocal"/>
            ///   </item>
            ///   <item>
            ///     Second, all shared vertices owned by this MPI process,
            ///     see <see cref="NoOfRelayed"/>
            ///   </item>
            ///   <item>
            ///     Third, all shared vertices owned by other MPI processes,
            ///     see <see cref="NoOfBorrowed"/>
            ///   </item>
            ///   <item>
            ///     Fourth, all external vertices, see
            ///     <see cref="NoOfExternal"/>
            ///   </item>
            /// </list>
            /// </remarks>
            public MultidimensionalArray Coordinates {
                get;
                internal set;
            }

            /// <summary>
            /// Number of Vertices on local MPI process
            /// </summary>
            public int Count {
                get {
                    return Coordinates.GetLength(0);
                }
            }

            /// <summary>
            /// Number of vertices which are used only by locally updated cells;
            /// </summary>
            public int NoOfPurelyLocal {
                get;
                private set;
            }

            /// <summary>
            /// Number of shared vertices (used by locally updated cells AND
            /// external cells) which are 'owned' by the current MPI process.
            /// </summary>
            public int NoOfRelayed {
                get;
                private set;
            }

            /// <summary>
            /// Number of vertices which are 'owned' by the current MPI process
            /// </summary>
            public int NoOfOwned {
                get {
                    return NoOfPurelyLocal + NoOfRelayed;
                }
            }

            /// <summary>
            /// Number of shared vertices (used by locally updated cells AND
            /// external cells) which are 'owned' by other MPI process.
            /// </summary>
            public int NoOfBorrowed {
                get;
                private set;
            }

            /// <summary>
            /// Number of Vertices which are used only by external cells;
            /// </summary>
            public int NoOfExternal {
                get;
                private set;
            }

            /// <summary>
            /// Number of nodes which are identical to others via the peridicity relations of the grid.
            /// </summary>
            public int NoOfPeriodicElim {
                get;
                private set;
            }

            /// <summary>
            /// Number of nodes that are used by locally updated cells.
            /// </summary>
            public int NoOfNodes4LocallyUpdatedCells {
                get {
                    int R = (this.Count - this.NoOfExternal);
                    Debug.Assert(R == (NoOfPurelyLocal + NoOfBorrowed + NoOfRelayed + NoOfPeriodicElim));
                    return R;
                }
            }

            private Partitioning m_NodePartitioning = null;

            /// <summary>
            /// partitioning of nodes across MPI processes
            /// </summary>
            public Partitioning NodePartitioning {
                get {
                    if (m_NodePartitioning == null) {
                        m_NodePartitioning =
                            new Partitioning(this.NoOfPurelyLocal + this.NoOfRelayed);
                    }
                    return m_NodePartitioning;
                }
            }

            /// <summary>
            /// content: which of the relayed vertices must be send to
            /// other processors
            ///  - 1st index: MPI rank of target processor 'R' 
            ///  - 2nd index: enumeration
            /// </summary>
            public int[][] VertexSendLists {
                get;
                private set;
            }

            /// <summary>
            /// content: where the borrowed vertices received by other
            /// processors must be inserted 
            ///  - 1st index: MPI rank of target processor 'R' 
            ///  - 2nd index: enumeration
            /// </summary>
            public int[][] VertexInsertLists {
                get;
                private set;
            }

            /// <summary>
            /// Computes a unique list of all nodes/vertices of all cells.
            /// </summary>
            /// <param name="vertice">
            /// output: the merged vertices, i.e. each vertex from the
            /// </param>
            /// <param name="cellVertices">
            /// output: indices into the <paramref name="vertice"/>-array.
            /// 1st index: cell index;
            /// 2nd index: vertex index within cell;
            /// </param>
            /// <param name="grdDat"></param>
            /// <param name="h">
            /// cell measure
            /// </param>
            /// <param name="IncludeExt">
            /// true, if also external cells should be cnsidered
            /// </param>
            internal static void CollectVertices(GridData grdDat, out int[][] cellVertices, out MultidimensionalArray vertice, double[] h, bool IncludeExt = false) {
                using (var tr = new FuncTrace()) {
                    // fk: bemerkung: sollte einegermaßen skalieren;
                    // Test am 10sept10: Gitter mit 400x400 dauert cd 8 sec., Gitter mit 800x400 dauert 16 sekunden

                    int D = grdDat.SpatialDimension;
                    int J;
                    J = IncludeExt ? grdDat.Cells.NoOfCells : grdDat.Cells.NoOfLocalUpdatedCells;

                    //int NV = celVtx.GetLength(1);
                    var Krefs = grdDat.Grid.RefElements;

                    // collect vertices
                    // ================

                    int No = 0;
                    {
                        int a, b;
                        for (int i = 0; i < Krefs.Length; i++) {
                            grdDat.GetLocalNoOfCellsRefElement(i, out a, out b);
                            No += (a + (IncludeExt ? b : 0)) * Krefs[i].NoOfVertices;
                        }
                    }
                    MultidimensionalArray Verts = MultidimensionalArray.Create(No, D);

                    // loop over all cells...

                    MultidimensionalArray[] vertsGlobal = (
                        from Kref in Krefs
                        select MultidimensionalArray.Create(1, Kref.NoOfVertices, D)).ToArray();

                    int cnt = 0;
                    for (int j = 0; j < J; j++) {
                        int iKref = grdDat.Cells.GetRefElementIndex(j);
                        var Kref = Krefs[iKref];
                        var Kj = grdDat.Cells.GetCell(j);
                        var vertsGlobal_iKref = vertsGlobal[iKref];

                        //Kref.TransformLocal2Global(Kref.Vertices, vertsGlobal_iKref, 0, Kj.Type, Kj.TransformationParams);
                        grdDat.TransformLocal2Global(Kref.Vertices, j, 1, vertsGlobal_iKref, 0);

                        int NV = Kref.NoOfVertices;

                        // loop over cell vertices...
                        for (int iv = 0; iv < NV; iv++) {
                            for (int d = 0; d < D; d++) {
                                Verts[cnt, d] = vertsGlobal_iKref[0, iv, d];
                            }
                            cnt++;
                        }
                    }

                    // build tree
                    // ==========

                    BoundingBox bb = new BoundingBox(Verts);
                    bb.ExtendByFactor(0.005);
                    int[] Perm = new int[Verts.GetLength(0)];
                    PointLocalization locTree = new PointLocalization(Verts, bb, Perm);
                    Verts = locTree.Points;

                    // eliminate duplicate points
                    // ==========================

                    double[] pt = new double[D];
                    int N = No;
                    List<int> foundPoints = new List<int>();
                    int[] AliasPts = new int[No];
                    AliasPts.SetAll(int.MinValue);
                    int NewVertice = 0;
                    List<int> VerticeTmp = new List<int>();

                    using (new BlockTrace("duplicate Point elimination", tr)) {
                        int n = -1;
                        for (int j = 0; j < J; j++) {
                            int iKref = grdDat.Cells.GetRefElementIndex(j);
                            var Kref = Krefs[iKref];
                            int NV = Kref.NoOfVertices;

                            for (int _n = 0; _n < NV; _n++) {
                                n++;

                                if (AliasPts[n] >= 0)
                                    continue; // point is already assigned.

                                Verts.GetRow(n, pt);  // pt = Verts[n,*]

                                
                                double eps = h[j] * 1.0e-6;

                                locTree.FindNearPoints(foundPoints, eps, pt);
                                if (foundPoints.Count < 1)
                                    throw new ApplicationException("error in algorithm");
                                if (!foundPoints.Contains(n)) {
                                    throw new ApplicationException("error in algorithm");
                                }

                                VerticeTmp.Add(n);

                                for (int k = 0; k < foundPoints.Count; k++) {
                                    AliasPts[foundPoints[k]] = NewVertice;
                                }

                                NewVertice++;

                            }
                        }
                    }

                    // test
                    // ====
                    for (int i = 0; i < AliasPts.Length; i++)
                        if (AliasPts[i] < 0)
                            throw new ApplicationException("error in alg");

                    // store results
                    // =============
                    int[] PermInv = new int[Perm.Length];
                    for (int i = 0; i < PermInv.Length; i++)
                        PermInv[Perm[i]] = i;

                    cellVertices = new int[J][];

                    int m = 0;
                    for (int j = 0; j < J; j++) {
                        int iKref = grdDat.Cells.GetRefElementIndex(j);
                        var Kref = Krefs[iKref];
                        int NV = Kref.NoOfVertices;

                        cellVertices[j] = new int[NV];

                        for (int _n = 0; _n < NV; _n++) {
                            cellVertices[j][_n] = AliasPts[PermInv[m]];
                            m++;
                        }
                    }

                    vertice = MultidimensionalArray.Create(VerticeTmp.Count, D);
                    for (int i = 0; i < NewVertice; i++) {
                        int n = VerticeTmp[i];
                        for (int d = 0; d < D; d++) {
                            vertice[i, d] = Verts[n, d];
                        }
                    }
                }
            }

            /// <summary>
            /// For each vertex, the local indices of the adjacent cells;
            ///  - 1st index: local vertex index
            ///  - 2nd index: collection
            /// </summary>
            public int[][] VerticeToCell {
                get;
                internal set;
            }

            /// <summary>
            /// init code for <see cref="VerticeToCell"/>.
            /// </summary>
            internal void Init_VerticeToCell() {
                Debug.Assert(VerticeToCell == null);

                int K = this.Count;
                int J = this.m_owner.Cells.NoOfCells;
                var CellVertices = this.m_owner.Cells.CellVertices;

                VerticeToCell = new int[K][];

                for (int j = 0; j < J; j++) { // loop over cells
                    var CellVertices_j = CellVertices[j];
                    foreach (int k in CellVertices_j) {
                        ArrayTools.AddToArray(j, ref VerticeToCell[k]);
                    }
                }
            }

            public SortedDictionary<int, int> PeriodicEliminatedPoints {
                get;
                private set;
            }

            List<Tuple<int, int>> PeriodicIdenitiesTmp;

            internal void Init_PeriodicPairings() {
                int NoOfEdges = this.m_owner.Edges.Count;
                int NoOfVtx = this.Count;
                var EdgeTags = this.m_owner.Edges.EdgeTags;
                var E2C = this.m_owner.Edges.CellIndices;
                var E2F = this.m_owner.Edges.FaceIndices;
                var C2V = this.m_owner.Cells.CellVertices;


                var Krefs = this.m_owner.Grid.RefElements;

                // find points which are identic due to periodicity
                // ================================================

                List<int>[] PeriodicPeer = new List<int>[NoOfVtx];

                bool AnythingPeriodic = false;
                for (int iEdg = 0; iEdg < NoOfEdges; iEdg++) {  // loop over edges
                    if (EdgeTags[iEdg] >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                        // found some periodic edge
                        // ++++++++++++++++++++++++


                        if (!(this.m_owner.Edges.IsEdgeConformalWithCell1(iEdg) && this.m_owner.Edges.IsEdgeConformalWithCell2(iEdg))) {
                            throw new NotImplementedException("Periodicity on non-conformal edges not suported yet.");
                        }

                        int jCell1 = E2C[iEdg, 0];
                        int jCell2 = E2C[iEdg, 1];

                        var Cell1 = this.m_owner.Cells.GetCell(jCell1);
                        var Cell2 = this.m_owner.Cells.GetCell(jCell1);

                        int iFace1 = E2F[iEdg, 0];
                        int iFace2 = E2F[iEdg, 1];

                        var Kref1 = m_owner.Cells.GetRefElement(jCell1);
                        var Kref2 = m_owner.Cells.GetRefElement(jCell2);

                        var FVtx1 = Kref1.GetFaceVertices(iFace1); // periodic nodes in local coordinates of first cell
                        var FVtx2 = Kref2.GetFaceVertices(iFace2); // "peers due to periodicity", in local coordinates of second cell

                        var FVtx1G = this.m_owner.GlobalNodes.GetValue_Cell(FVtx1, jCell1, 1).ExtractSubArrayShallow(0, -1, -1);
                        var FVtx2G = this.m_owner.GlobalNodes.GetValue_Cell(FVtx2, jCell2, 1).ExtractSubArrayShallow(0, -1, -1);

                        Debug.Assert(object.ReferenceEquals(this.m_owner.Edges.GetRefElement(iEdg), Kref1.FaceRefElement));
                        Debug.Assert(object.ReferenceEquals(this.m_owner.Edges.GetRefElement(iEdg), Kref2.FaceRefElement));
                        Debug.Assert(FVtx1.NoOfNodes == FVtx2.NoOfNodes);
                        int L = FVtx1.NoOfNodes;

                        int[] GlobIdxFace1 = new int[FVtx1.NoOfNodes];
                        int[] GlobIdxFace2 = new int[FVtx2.NoOfNodes];
                        Debug.Assert(Kref1.FaceToVertexIndices.GetLength(1) == L);
                        Debug.Assert(Kref2.FaceToVertexIndices.GetLength(1) == L);
                        for (int i = 0; i < L; i++) {
                            GlobIdxFace1[i] = C2V[jCell1][Kref1.FaceToVertexIndices[iFace1, i]];
                            GlobIdxFace2[i] = C2V[jCell2][Kref2.FaceToVertexIndices[iFace2, i]];
                        }

                        var Trafo = this.m_owner.Edges.GetPeriodicTrafo(iEdg, true);

                        double Tol = Math.Min(FVtx1G.MindistBetweenRows(), FVtx1G.MindistBetweenRows()) * 1.0e-3;
                        int[] R = NodeCorrespondence(FVtx1G, Trafo, FVtx2G, Tol);

                        // Es sollte gelten:
                        // GlobIdxFace1[R[k]] <==periodic==> GlobIdxFace2[k]

                        for (int k = 0; k < L; k++) {
                            int iVtx1 = GlobIdxFace1[R[k]];
                            int iVtx2 = GlobIdxFace2[k];

                            if (PeriodicPeer[iVtx1] == null)
                                PeriodicPeer[iVtx1] = new List<int>();
                            if (PeriodicPeer[iVtx2] == null)
                                PeriodicPeer[iVtx2] = new List<int>();

                            if (!PeriodicPeer[iVtx1].Contains(iVtx2))
                                PeriodicPeer[iVtx1].Add(iVtx2);
                            if (!PeriodicPeer[iVtx2].Contains(iVtx1))
                                PeriodicPeer[iVtx2].Add(iVtx1);
                        }

                        AnythingPeriodic = true;
                    }
                }

                // eliminate periodic 
                // ==================

                if (AnythingPeriodic) {
                    this.PeriodicIdenitiesTmp = new  List<Tuple<int, int>>();

                    for (int iVtx = 0; iVtx < NoOfVtx; iVtx++) {
                        List<int> pp = PeriodicPeer[iVtx];
                        if (pp != null) {
                            if (!pp.Contains(iVtx))
                                pp.Add(iVtx);

                            for (int i0 = 0; i0 < pp.Count; i0++) {
                                for (int i1 = i0 + 1; i1 < pp.Count; i1++) {
                                    var T = new Tuple<int, int>(pp[i0], pp[i1]);

                                    if(!this.PeriodicIdenitiesTmp.Contains(T, (A,B) => (
                                             (A.Item1 == B.Item1 && A.Item2 == B.Item2) 
                                          || (A.Item2 == B.Item1 && A.Item1 == B.Item2)  ))) {
                                        this.PeriodicIdenitiesTmp.Add(T);
                                    }

                                }
                            }

                        }
                    }
                }
                csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            }


            /// <summary>
            /// Geometrical matching of nodes.
            /// </summary>
            /// <param name="NodesIn">
            /// First set of nodes.
            /// </param>
            /// <param name="ict"></param>
            /// <param name="NodesOut">
            /// second set of nodes.
            /// </param>
            /// <param name="Tol">
            /// Some tolerance
            /// </param>
            /// <returns>
            /// A permutation R, so that for all valid indices k, 
            /// <paramref name="NodesOut"/>[R[k]] == <paramref name="ict"/>(<paramref name="NodesIn"/>[k]).
            /// </returns>
            static int[] NodeCorrespondence(MultidimensionalArray NodesIn, AffineTrafo ict, MultidimensionalArray NodesOut, double Tol) {
                var NodesInTrf = MultidimensionalArray.Create(NodesIn.Lengths);
                ict.Transform(NodesIn, NodesInTrf);

                if (!ArrayTools.ListEquals(NodesInTrf.Lengths, NodesOut.Lengths))
                    throw new ApplicationException("Error in algorithm.");

                int NoOfNodes = NodesIn.GetLength(0);
                var R = new int[NoOfNodes];
                for (int kOut = 0; kOut < NoOfNodes; kOut++) {
                    var vOut = NodesOut.GetRow(kOut);

                    bool bFound = false;

                    for (int kIn = 0; kIn < NoOfNodes; kIn++) {
                        var vIn = NodesInTrf.GetRow(kIn);
                        if (GenericBlas.L2Dist(vOut, vIn) < Tol) {
                            bFound = true;
                            R[kOut] = kIn;
                            break;
                        }
                    }

                    if (bFound == false)
                        throw new ApplicationException("Error in algorithm.");

                }

                return R;
            }

            
            /// <summary>
            /// Work-around for a bug in Mono that has problems
            /// (de)serializing staggered arrays
            /// </summary>
            [Serializable]
            class VertexSharingData {

                /// <summary>
                /// 1st index: cell
                /// 2nd index: cell vertex
                /// 3rd index: collection
                /// </summary>
                public Tuple<int, int>[][][] Data;
            }


            internal void NegogiateSharing() {
                int J = this.m_owner.Cells.NoOfLocalUpdatedCells;
                int[][] SendLists = m_owner.Parallel.SendCommLists;
                int[] InsertIndex = m_owner.Parallel.RcvCommListsInsertIndex;
                int[] RcvNoOfItems = m_owner.Parallel.RcvCommListsNoOfItems;
                int size = this.m_owner.MpiSize;
                int myRank = this.m_owner.MpiRank;
                var Celldat = m_owner.Cells;

                int K = this.Count; // including purly external!

                bool containsPeriodic = this.PeriodicEliminatedPoints != null;

                // PHASE 1: determine vertex indices on other processors
                // =====================================================
                
                // 1st index: local vertex index
                // 2nd index: enum
                // Item1: processor rank 'R'
                // Item2: local vertex index on processor 'R'
                Tuple<int, int>[][] VertexIndicesOnOtherProcessors = new Tuple<int, int>[K][];
                BitArray NonExternalMarker = new BitArray(K); // false on all vertices that are "pure external" (only used by external cells)
                {
                    // mark all vertices on locally updated cells
                    for (int j = 0; j < J; j++) {
                        int[] CV = Celldat.CellVertices[j];
                        foreach (int k in CV)
                            NonExternalMarker[k] = true;
                    }
                    
                    // add the vertex index on this processor
                    for (int k = 0; k < K; k++) {
                        if (NonExternalMarker[k])
                            VertexIndicesOnOtherProcessors[k] = new Tuple<int, int>[] { new Tuple<int, int>(myRank, k) };
                        else
                            VertexIndicesOnOtherProcessors[k] = new Tuple<int, int>[0];
                    }


                    
                    if (this.PeriodicIdenitiesTmp != null) {
                        foreach (var kv in this.PeriodicIdenitiesTmp) {
                            int iVtx1 = kv.Item1;
                            int iVtx2 = kv.Item2;
                            // 'iVtx1' and 'iVtx2' are identical due to periodicity

                            if (NonExternalMarker[iVtx2])
                                (new Tuple<int, int>(myRank, iVtx2)).AddToArray(ref VertexIndicesOnOtherProcessors[iVtx1]);

                            if (NonExternalMarker[iVtx1])
                                (new Tuple<int, int>(myRank, iVtx1)).AddToArray(ref VertexIndicesOnOtherProcessors[iVtx2]);
                        }
                    }


                    // send vertex indices to other processors
                    // =======================================


                    int SomeChange;
                    do {
                        // This has to be done repetitive:
                        // e.g. the center vertex x may be shared by 4 processes, 
                        // 
                        //    *-------*-------*
                        //    | rank0 | rank1 | 
                        //    *-------x-------*
                        //    | rank2 | rank3 | 
                        //    *-------*-------*
                        //
                        // but each process has only two neighbours (with respect to cell-neighbourship).
                        // Thus, in one pass e.g. rank 0 only gets to know that it shares x with rank 1 and rank 2,
                        // while rank 1 gets to know that it shares x with rank 0 and rank 3.
                        // In the secon pass, rank 1 informs rank 0 that x is also shared by rank 3, and vice-versa.

                        SomeChange = 0;

                        Dictionary<int, VertexSharingData> SendData = new Dictionary<int, VertexSharingData>();
                        for (int proc = 0; proc < size; proc++) { // loop over MPI process ranks
                            int[] SendList = SendLists[proc];

                            if (SendList == null)
                                continue;

                            var SndItem = new VertexSharingData();

                            SndItem.Data = new Tuple<int, int>[SendList.Length][][];

                            for (int jj = 0; jj < SendList.Length; jj++) { // loop over all cells that will be send to processor #'proc'
                                int jCell = SendList[jj];
                                //var Kref = Celldat.GetRefElement(jCell);

                                int[] CellVtx = Celldat.CellVertices[jCell];
                                SndItem.Data[jj] = new Tuple<int, int>[CellVtx.Length][];

                                for (int i = 0; i < CellVtx.Length; i++) { // loop over cell vertices
                                    SndItem.Data[jj][i] = VertexIndicesOnOtherProcessors[CellVtx[i]];
                                }
                            }

                            SendData.Add(proc, SndItem);
                        }

                        var RcvData = SerialisationMessenger.ExchangeData(SendData, csMPI.Raw._COMM.WORLD);

                        foreach (var kv in RcvData) {
                            int proc = kv.Key;
                            VertexSharingData RcvItem = kv.Value;
                            int jCell0 = InsertIndex[proc];
                            int L = RcvItem.Data.Length;
                            Debug.Assert(L == RcvNoOfItems[proc]);


                            for (int l = 0; l < L; l++) { // loop over received cells
                                int jCell = jCell0 + l;
                                int[] CellVtx = Celldat.CellVertices[jCell];
                                var VtxIndex_jCell = RcvItem.Data[l];
                                Debug.Assert(VtxIndex_jCell.Length == Celldat.GetRefElement(jCell).NoOfVertices);

                                for (int i = 0; i < VtxIndex_jCell.Length; i++) { // loop over cell vertices
                                    foreach (Tuple<int, int> ProcIvtx in VtxIndex_jCell[i]) {
                                        int k = CellVtx[i];

                                        bool contains = VertexIndicesOnOtherProcessors[k].Contains(ProcIvtx, delegate (Tuple<int, int> a, Tuple<int, int> b) {
                                            bool ret1 = a.Item1 == b.Item1; // equal processor rank
                                            bool ret2 = a.Item2 == b.Item2; // equal vertex index
                                            //Debug.Assert((ret1 == false || a.Item2 == b.Item2)
                                            //    || (containsPeriodic && (this.PeriodicEliminatedPoints.Keys.Contains(a.Item2) || this.PeriodicEliminatedPoints.Keys.Contains(b.Item2))));
                                            return ret1 && ret2;
                                        });

                                        if (!contains) {
                                            SomeChange = 0xFFF;
                                            ProcIvtx.AddToArray(ref VertexIndicesOnOtherProcessors[k]);
                                        }
                                    }
                                }
                            }
                        }

                        if (this.PeriodicIdenitiesTmp != null) {
                            // ++++++++++++++++++
                            //  periodic exchange
                            // ++++++++++++++++++

                            foreach (var kv in this.PeriodicIdenitiesTmp) {
                                int iVtx1 = kv.Item1;
                                int iVtx2 = kv.Item2;
                                // 'iVtx1' and 'iVtx2' are identical due to periodicity

                               
                                // ensure that 'list1' and 'list2' are equal in a set-sense, i.e.
                                // contain the same elements; the sequence does not matter.

                                for (int iii = 0; iii <= 2; iii++) {
                                    
                                    foreach (var T in VertexIndicesOnOtherProcessors[iVtx2]) {

                                        bool contains = VertexIndicesOnOtherProcessors[iVtx1].Contains(T, delegate (Tuple<int, int> a, Tuple<int, int> b) {
                                            bool ret1 = a.Item1 == b.Item1; // equal processor rank
                                            bool ret2 = a.Item2 == b.Item2; // equal vertex index
                                            return ret1 && ret2;
                                        });

                                        if (!contains) {
                                            SomeChange = 0xFFF;
                                            T.AddToArray(ref VertexIndicesOnOtherProcessors[iVtx1]);
                                        }

                                    }

                                    // swap the indices for the second run of the 'iii'-loop
                                    int iv1tmp = iVtx1;
                                    iVtx1 = iVtx2;
                                    iVtx2 = iv1tmp;
                                }
                            }
                        }


                        unsafe {
                            int globSomeChange = 0;
                            csMPI.Raw.Allreduce((IntPtr)(&SomeChange), (IntPtr)(&globSomeChange), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                            SomeChange = globSomeChange;
                        }

                    } while (SomeChange != 0);


                    this.PeriodicIdenitiesTmp = null;
                }

                // PHASE2: negotiate Ownership
                // ===========================

                // VERTEX OWNERSHIP NEGOGIATION:
                // * if a vertex is owned by more than two processors, than the process with the lowest MPI rank will be the owner
                //   (these won't be a lot, so it should produce no significant load imbalance)
                // * vertices shared by exactly two processors are shared between them in zu gleichen Teilen

                // Arten von Vertices
                // (A) purly local
                // (B) shared/gehört mir  == relayed
                // (C) shared/gehört wem anderen  == borrowed
                // (D) external (nur von external cells benutzt);

                // Was will ich wissen?
                // Für jeden (C): Welcher Prozessor besitzt ihn? Welchen Index hat er dort?


                int[] VerticePermuation;
                int NoOfPureLocal;
                int NoOfRelayed;
                int NoOfBorrowed;
                int NoOfExternal;
                int NoOfPeriodicElim;
                int[][] VtxSendLists;
                int[][] VtxInsertLists;
                Tuple<int, int>[] PeriodicElim;

                /*
                using (var diag = new System.IO.StreamWriter("NegOwInput-Proc" + m_owner.MyRank + ".txt")) {
                    for (int k = 0; k < K; k++) {
                        diag.Write(k);
                        diag.Write(" ");
                        diag.Write(this.Coordinates[k, 0] + " " + this.Coordinates[k, 1]);
                        diag.Write(" ");
                        diag.Write(NonExternalMarker[k]);
                        diag.Write(" ");
                        foreach (var t in VertexIndicesOnOtherProcessors[k]) {
                            diag.Write("(" + t.Item1 + "," + t.Item2 + ") "); 
                        }
                        diag.WriteLine();
                    }
                }

                */

                NegogiateOwnership(csMPI.Raw._COMM.WORLD, VertexIndicesOnOtherProcessors, 
                    out VerticePermuation, 
                    out NoOfPureLocal, out NoOfRelayed, out NoOfBorrowed, out NoOfPeriodicElim, out NoOfExternal, 
                    out VtxSendLists, out VtxInsertLists,
                    out PeriodicElim);

                // PHASE 4: apply permutation of vertices
                // ======================================
                {
                    int[] invVerticePermuation = new int[K];
                    for (int k = 0; k < K; k++) {
                        invVerticePermuation[VerticePermuation[k]] = k;
                    }

                    MultidimensionalArray _vertice = this.Coordinates;
                    MultidimensionalArray newVertice = MultidimensionalArray.Create(_vertice.Lengths);
                    Debug.Assert(newVertice.Dimension == 2);

                    int[][] _Vertice2Cell = this.VerticeToCell;
                    int[][] newVertice2Cell = new int[K][];

                    for (int k = 0; k < K; k++) {
                        newVertice2Cell[k] = _Vertice2Cell[VerticePermuation[k]];
                        newVertice.ExtractSubArrayShallow(k, -1).Set(_vertice.ExtractSubArrayShallow(VerticePermuation[k], -1));
                    }
                    
                    this.Coordinates = newVertice;
                    this.VerticeToCell = newVertice2Cell;
                    
                    var CV = this.m_owner.Cells.CellVertices;
                    for (int j = 0; j < CV.Length; j++) {
                        CV[j] = CV[j].Select(k => invVerticePermuation[k]).ToArray();
                    }



                    if (PeriodicElim.Count() > 0) {

                        int i0_Src = NoOfPureLocal + NoOfBorrowed + NoOfRelayed;
                        int iE_Src = i0_Src + NoOfPeriodicElim;

                        int i0_Trg = 0;
                        int iE_Trg = NoOfPureLocal + NoOfBorrowed + NoOfRelayed;

                        //Debugger.Launch();

                        SortedDictionary<int, int> newPeriodicEliminatedPoints = new SortedDictionary<int, int>();
                        foreach (var t in PeriodicElim) {
                            int elimSrc = invVerticePermuation[t.Item1];
                            int elimTrg = invVerticePermuation[t.Item2];
                            Debug.Assert(elimSrc >= i0_Src);
                            Debug.Assert(elimSrc < iE_Src);
                            Debug.Assert(elimTrg >= i0_Trg);
                            Debug.Assert(elimTrg < iE_Trg);

                            if (newPeriodicEliminatedPoints.Keys.Contains(elimSrc))
                                Debugger.Launch();

                            newPeriodicEliminatedPoints.Add(elimSrc, elimTrg);
                        }
                        this.PeriodicEliminatedPoints = newPeriodicEliminatedPoints;

                    } else {
                        PeriodicEliminatedPoints = new SortedDictionary<int, int>();
                    }


                    // transform send lists (based on un-permuded indices)
                    for (int rnk = 0; rnk < size; rnk++) {
                        if (VtxSendLists[rnk] != null) {
                            VtxSendLists[rnk] = VtxSendLists[rnk].Select(k => invVerticePermuation[k]).ToArray();
                        }
                    }

                    // transform insert lists 
                    for (int rnk = 0; rnk < size; rnk++) {
                        if (VtxInsertLists[rnk] != null) {
                            VtxInsertLists[rnk] = VtxInsertLists[rnk].Select(k => invVerticePermuation[k]).ToArray();
                        }
                    }
                }

                // set data
                // ========
                this.NoOfPurelyLocal = NoOfPureLocal;
                this.NoOfRelayed = NoOfRelayed;
                this.NoOfBorrowed = NoOfBorrowed;
                this.NoOfExternal = NoOfExternal;
                this.NoOfPeriodicElim = NoOfPeriodicElim;
                this.VertexSendLists = VtxSendLists;
                this.VertexInsertLists = VtxInsertLists;


                // write data (debugging) 
                // ======================
                /*
                {

                    int D = this.m_owner.SpatialDimension;
                    using (var stw = new System.IO.StreamWriter("Vertex_" + myRank + ".txt")) {
                        stw.WriteLine("Rank: " + myRank);
                        stw.WriteLine("NoOf pure local:     {0}", NoOfPureLocal);
                        stw.WriteLine("NoOf relayed:        {0}", NoOfRelayed);
                        stw.WriteLine("NoOf borrowed:       {0}", NoOfBorrowed);
                        stw.WriteLine("NoOf periodic elim:  {0}", NoOfPeriodicElim);
                        stw.WriteLine("NoOf external:       {0}", NoOfExternal);

                        stw.WriteLine("Vertex List:");
                        
                        for (int k = 0; k < K; k++) {
                            stw.Write(k);
                            stw.Write("\t");

                            for (int d = 0; d < D; d++) {
                                stw.Write(this.Coordinates[k, d]);
                                stw.Write("\t");
                            }

                            int OwnerProcess;
                            if (k < NoOfPureLocal + NoOfRelayed) {
                                OwnerProcess = myRank;
                            } else if (k < NoOfPureLocal + NoOfRelayed + NoOfBorrowed) {
                                OwnerProcess = 22; // BorrowedOwnership[k - (NoOfPureLocal + NoOfRelayed)];
                            } else {
                                OwnerProcess = -1;
                            }
                            stw.Write(OwnerProcess);
                            //stw.Write("\t");

                            stw.WriteLine();
                        }
                        stw.WriteLine("-----------------------");
                        stw.WriteLine("Send Lists:");
                        for (int rnk = 0; rnk < size; rnk++) {
                            if (VtxSendLists[rnk] == null) {
                                continue;
                            }
                            stw.Write(rnk + ":\t");

                            for (int i = 0; i < VtxSendLists[rnk].Length; i++) {
                                stw.Write(VtxSendLists[rnk][i]);
                                if( i < (VtxSendLists[rnk].Length - 1))
                                    stw.Write(",");
                                else
                                    stw.Write(";");
                            }
                            stw.WriteLine();
                        }
                        stw.WriteLine("-----------------------");
                        stw.WriteLine("Insert Lists:");
                        for (int rnk = 0; rnk < size; rnk++) {
                            if (VtxInsertLists[rnk] == null) {
                                continue;
                            }
                            stw.Write(rnk + ":\t");

                            for (int i = 0; i < VtxInsertLists[rnk].Length; i++) {
                                stw.Write(VtxInsertLists[rnk][i]);
                                if (i < (VtxInsertLists[rnk].Length - 1))
                                    stw.Write(",");
                                else
                                    stw.Write(";");
                            }
                            stw.WriteLine();
                        }

                        stw.WriteLine("-----------------------");
                        stw.WriteLine("Local periodic elminiations: {0}", this.PeriodicEliminatedPoints.Count);
                        foreach (var kv in this.PeriodicEliminatedPoints) {
                            stw.WriteLine("{0} -> {1}", kv.Key, kv.Value);
                        }
                    }

                } // */

            }
        }


        /// <summary>
        /// Negogiation of MPI-Ownership for 'items' (e.g. vertices or edges).
        /// </summary>
        private static void NegogiateOwnership(MPI_Comm comm, Tuple<int, int>[][] ItemIndicesOnOtherProcessors, 
            out int[] ItemPermutation, 
            out int NoOfPureLocal, out int NoOfRelayed, out int NoOfBorrowed, out int NoOfPeriodicElim, out int NoOfExternal, 
            out int[][] SendLists, out int[][] InsertLists,
            out Tuple<int,int>[] LocalDuplicates
            ) {

            int size, myRank;
            csMPI.Raw.Comm_Rank(comm, out myRank);
            csMPI.Raw.Comm_Size(comm, out size);

            int K = ItemIndicesOnOtherProcessors.Length; // 
            BitArray ActiveItemsMarker = new BitArray(K);
            for (int k = 0; k < K; k++) {
                if (ItemIndicesOnOtherProcessors[k] == null || ItemIndicesOnOtherProcessors[k].Length <= 0) {
                    ActiveItemsMarker[k] = false;
                } else {
                    var Test = ItemIndicesOnOtherProcessors[k].SingleOrDefault(t => t.Item1 == myRank && t.Item2 == k);
                    ActiveItemsMarker[k] = (Test != null);
                }
            }

            ItemPermutation = new int[K];
            NoOfPureLocal = 0;
            NoOfRelayed = 0;

            // marks all vertices which are owned by this process
            BitArray Owned = new BitArray(K);

            BitArray LocalElim = new BitArray(K); // item 'k' where the 'LocalElim[k]'-flag is set will be moved
            //                                       to the end, i.e. after the borrowed items.
            //                                       It marks representants that are redundant locally.
            List<Tuple<int, int>> LocalIdentities = new List<Tuple<int, int>>();

            // find the local identities
            // =========================
            for (int k = 0; k < K; k++) { // loop over items...
                if (!ActiveItemsMarker[k])
                    // don't bother
                    continue;
                if (LocalElim[k]) {
                    // already eliminated
                    Debug.Assert(LocalIdentities.Where(T => T.Item1 == k).Count() == 1); // the 'LocalIdentities' must already contain item 'k'
                    continue;
                }

                var Equivals = ItemIndicesOnOtherProcessors[k]; // equivalences of item 'k'
                if (Equivals.Length <= 0)
                    throw new ArgumentException();
                if (Equivals.Length == 1) {
                    Debug.Assert(Equivals[0].Item1 == myRank);
                    // no known equivalences
                    continue;
                }

                foreach (var T in Equivals) {
                    if (T.Item1 == myRank && T.Item2 != k) {
                        LocalElim[T.Item2] = true;

                        var lid = new Tuple<int, int>(T.Item2, k);
                        Debug.Assert(lid.Item1 > lid.Item2); // the eliminated representative 'lid.Item2' of some item must be higher than the representative 'k' which is kept
                        Debug.Assert(LocalIdentities.Contains(lid, (A, B) => A.Item1 == B.Item1 && A.Item2 == B.Item2) == false);
                        LocalIdentities.Add(lid); // 'T.Item2' is equivalent to 'k': we select 'k' as a representative of the item.
                    }
                }
            }
            

            // for each MPI process,
            // collect the vertices that are shared with this process
            // ======================================================
            HashSet<int>[] _SharedVertices = new HashSet<int>[size];
            for (int k = 0; k < K; k++) {
                if (!ActiveItemsMarker[k])
                    // don't consider representative k
                    continue;
                if (LocalElim[k])
                    // already eliminated
                    continue;
                
                var VS_k = ItemIndicesOnOtherProcessors[k];
                Debug.Assert(VS_k.Length > 0);
                if (VS_k.Length == 1 || VS_k.Select(T => T.Item1).ToSet().Count == 1) {
                    // Vertex is purly local 
                    Debug.Assert(VS_k[0].Item1 == myRank); // 'k's without 'myRank' should already be excluded by 'ActiveItemsMarker'

                    ItemPermutation[NoOfPureLocal] = k;
                    NoOfPureLocal++;
                    Owned[k] = true;
                } else {
                    foreach (Tuple<int, int> t in VS_k) {
                        int proc = t.Item1;
                        if (proc == myRank)
                            continue;
                        if (_SharedVertices[proc] == null) {
                            _SharedVertices[proc] = new HashSet<int>();
                        }
                        _SharedVertices[proc].Add(k);
                    }
                }
            }

            // sort the shared vertices
            // (i am not sure if that is really necessary)
            List<int>[] SharedVertices = _SharedVertices.Select(hs => hs != null ? (new List<int>(hs)) : default(List<int>)).ToArray();
            foreach (var lst in SharedVertices) {
                if (lst != null)
                    lst.Sort();
            }
            Debug.Assert(SharedVertices[myRank] == null);

            // determine 'RELAYED' and 'BORROWED' items,
            // i.e. determine ownership of shared items
            // ================================================

            BitArray RankDetermined = new BitArray(K);
            List<int> RelayedItems = new List<int>();


            Dictionary<int, int[]> SendData2 = new Dictionary<int, int[]>();
            for (int otherRank = myRank + 1; otherRank < size; otherRank++) { 
                List<int> SV_proc = SharedVertices[otherRank];
                if (SV_proc == null)
                    // sharing no vertex with processor 'otherRank'
                    continue;
                
                // List of items for which this process is responsible to
                // distribute 
                List<int> ItemsToDistribute = new List<int>();


                foreach (int kItem in SV_proc) {  // loop over all items shared with processor 'otherRank'...
                    if (LocalElim[kItem])
                        continue;
                    if (!ActiveItemsMarker[kItem])
                        continue;
                    if (RankDetermined[kItem])
                        // this item is shared by more than two processors, 
                        // but an owner was already assigned.
                        continue;

                    Debug.Assert(Owned[kItem] == false);

                    var OnOther = ItemIndicesOnOtherProcessors[kItem];
                    int minProc = OnOther.Select(T => T.Item1).Min();
                    Debug.Assert(OnOther.Select(T => T.Item1).Max() > minProc); // if not at least two different process ranks, something is worng with the sharing.

                    if (minProc == myRank) {
                        // ++++++++++++++++++++++++++++++++++++++++++++
                        // This process is responsible for determining 
                        // the ownership of the item represented by 'k' 
                        // ++++++++++++++++++++++++++++++++++++++++++++


                        Debug.Assert(!ItemsToDistribute.Contains(kItem));
                        Debug.Assert(LocalElim[kItem] == false);
                        ItemsToDistribute.Add(kItem);
                        RankDetermined[kItem] = true; // the owner rank of this item will be determined immediately (see below)!
                        //                               'kItem' will be owned either by 'myRank' or 'otherRank'
                    } 
                }

                // keep the first half of items for myself... (RELAYED)
                int L2 = ItemsToDistribute.Count / 2;
                int l;
                for (l = 0; l < L2; l++) {

                    int kItem = ItemsToDistribute[l];
                    
                    
                    bool found = false;
                    foreach (var T in ItemIndicesOnOtherProcessors[kItem]) {
                        if (T.Item1 == myRank && ActiveItemsMarker[T.Item2] == true) {
                            if (T.Item2 == kItem) {
                                found = true;
                                Debug.Assert(Owned[kItem] == false);
                                Owned[kItem] = true;
                                Debug.Assert(LocalElim[kItem] == false);
                                ItemPermutation[NoOfPureLocal + NoOfRelayed] = kItem;
                                RelayedItems.Add(kItem);
                                NoOfRelayed++;
                            } else {
                                Debug.Assert(LocalIdentities.Contains(new Tuple<int, int>(T.Item2, kItem), (A, B) => A.Item1 == B.Item1 && A.Item2 == B.Item2) == true);
                                Debug.Assert(LocalElim[T.Item2] == true);
                            }
                        } 
                    }
                    Debug.Assert(found == true);

                }
                // ... and hand the second half of items to the other process (BORROWED)
                int[] OwnedByOther = new int[ItemsToDistribute.Count - L2];
                for (; l < ItemsToDistribute.Count; l++) {
                    int k_myRank = ItemsToDistribute[l];
                    var OnOther = ItemIndicesOnOtherProcessors[k_myRank];

                    RankDetermined[k_myRank] = true;

                    int k_otherRank = OnOther.Where(T => T.Item1 == otherRank).Min(T => T.Item2); // pick the minimum representative 'k_otherRank' on processor 'otherRank' 
                    //                                                                               for the item represented by 'k_myRank' on this processor

                    

                    OwnedByOther[l - L2] = k_otherRank;
                }
                SendData2.Add(otherRank, OwnedByOther);
            }

            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

            // andere prozesse senden uns Indices, 
            // mein schatz,
            // die nur uns gehören.
            var RcvData2 = SerialisationMessenger.ExchangeData(SendData2, comm);
            foreach (var kv in RcvData2) {
                int rcvProc = kv.Key; // data was received from process 'rcvProc'
                int[] OurVertices = kv.Value;

                foreach (int k_mein in OurVertices) { // loop over all items that belong to me!
                    Debug.Assert(ItemIndicesOnOtherProcessors[k_mein].Where(T => T.Item1 == rcvProc).Count() >= 1);
                    Debug.Assert(ItemIndicesOnOtherProcessors[k_mein].Where(T => T.Item1 == myRank && T.Item2 == k_mein).Count() >= 1);
                    Debug.Assert(LocalElim[k_mein] == false);
                    Debug.Assert(ActiveItemsMarker[k_mein] == true);

                    bool found = false;
                    foreach (var T in ItemIndicesOnOtherProcessors[k_mein]) {
                        if (T.Item1 == myRank && ActiveItemsMarker[T.Item2] == true) {
                            Debug.Assert(Owned[T.Item2] == false);

                            if (T.Item2 == k_mein) {
                                found = true;
                                ItemPermutation[NoOfPureLocal + NoOfRelayed] = k_mein;
                                NoOfRelayed++;
                                Debug.Assert(Owned[k_mein] == false);
                                RelayedItems.Add(k_mein);
                                Owned[k_mein] = true;
                                Debug.Assert(LocalElim[k_mein] == false);

                            } else {
                                Debug.Assert(LocalIdentities.Contains(new Tuple<int, int>(T.Item2, k_mein), (A, B) => A.Item1 == B.Item1 && A.Item2 == B.Item2) == true);
                                Debug.Assert(LocalElim[T.Item2] == true);
                            }
                        } 
                    }
                    Debug.Assert(found == true);
                }
            }



            // ++++++++
            // Nun sollten alle Items genau einen Besitzer haben!
            // ++++++++

            // collect borrowed vertices
            // =========================
            {
                NoOfBorrowed = 0;
                for (int k = 0; k < K; k++) {  // loop over vertices...
                    if (Owned[k] || !ActiveItemsMarker[k] || LocalElim[k]) // ignore: owned, external, periodic eliminations
                        continue;

                    // if we reach this point, 'k' is (representing) a borrowed item

                    ItemPermutation[NoOfPureLocal + NoOfRelayed + NoOfBorrowed] = k;// sharedForeign[i].Item1;
                    NoOfBorrowed++;
                }
            }


            // determine send lists
            // ====================
            {
                
                List<int>[] _VtxSendLists = new List<int>[size];
                List<int>[] TargetIndices = new List<int>[size];
                foreach (int k in RelayedItems) {
                    Debug.Assert(Owned[k]);
                    Debug.Assert(ActiveItemsMarker[k]);
                    Debug.Assert(!LocalElim[k]);

                    var OnOther = ItemIndicesOnOtherProcessors[k];
                    ISet<int> Processores = OnOther.Select(T => T.Item1).ToSet();
                    Debug.Assert(Processores.Contains(myRank));
                    Processores.Remove(myRank);

                    foreach (int TargetRank in Processores) {
                        int k_Target = OnOther.Where(T => T.Item1 == TargetRank).Min(T => T.Item2);

                        if (_VtxSendLists[TargetRank] == null) {
                            _VtxSendLists[TargetRank] = new List<int>();
                            TargetIndices[TargetRank] = new List<int>();
                        }
                        _VtxSendLists[TargetRank].Add(k);
                        TargetIndices[TargetRank].Add(k_Target);
                    }
                }

                SendLists = _VtxSendLists.Select(list => list != null ? list.ToArray() : null).ToArray();

                Dictionary<int, int[]> TargetIndicesBla = new Dictionary<int, int[]>();
                for (int p = 0; p < size; p++) {
                    if (TargetIndices[p] != null)
                        TargetIndicesBla.Add(p, TargetIndices[p].ToArray());
                }

                var ReceiveIndices = SerialisationMessenger.ExchangeData(TargetIndicesBla, comm);

                InsertLists = new int[size][];
                foreach (var kv in ReceiveIndices) {
                    int rcvRank = kv.Key;
                    int[] InsList = kv.Value;
                    InsertLists[rcvRank] = InsList;

                    Debug.Assert(InsList.Where(k_Item => Owned[k_Item] == true).Count() <= 0);
                    Debug.Assert(InsList.Where(k_Item => LocalElim[k_Item] == true).Count() <= 0);
                    Debug.Assert(InsList.Where(k_Item => ActiveItemsMarker[k_Item] == false).Count() <= 0);
                }
            }


            // sorting of periodic eliminations
            // ================================
            {
                int ii = 0;
                for (int k = 0; k < K; k++) {
                    if (LocalElim[k] && ActiveItemsMarker[k]) {
                        // vertex k is eliminated due to a periodic identity
                        ItemPermutation[NoOfPureLocal + NoOfRelayed + NoOfBorrowed + ii] = k;
                        ii++;
                    }
                }
                NoOfPeriodicElim = ii;
                LocalDuplicates = LocalIdentities.ToArray();
            }

            // sorting of purly external
            // =========================
            {
                int ii = 0;
                for (int k = 0; k < K; k++) {
                    if (ActiveItemsMarker != null && !ActiveItemsMarker[k]) {
                        // vertex k is used purly by external cells
                        ItemPermutation[NoOfPureLocal + NoOfRelayed + NoOfBorrowed + NoOfPeriodicElim + ii] = k;
                        ii++;
                    }
                }
                NoOfExternal = ii;
            }


        }
    }
}
