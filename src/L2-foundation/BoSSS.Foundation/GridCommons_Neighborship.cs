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
using System.Diagnostics;
using System.Linq;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;

namespace BoSSS.Foundation.Grid.Classic {

    public partial class GridCommons {
        /// <summary>
        /// return values of <see cref="GetCellNeighbourship"/>.
        /// </summary>
        [Serializable]
        public struct Neighbour {

            /// <summary>
            /// global index of neighbor cell.
            /// </summary>
            public long Neighbour_GlobalIndex;

            /// <summary>
            /// if present, a face tag
            /// </summary>
            public CellFaceTag CellFaceTag;

            /// <summary>
            /// true, if the neighbor cell is reached by a periodic boundary
            /// condition
            /// </summary>
            public bool IsPeriodicNeighbour {
                get {
                    return CellFaceTag.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG;
                }
            }

            ///// <summary>
            ///// See <see cref="CellFaceTag."/>.
            ///// </summary>
            //public bool EdgeMayBeEmpty;
        }


        /// <summary>
        /// Helper that wraps a node index with an associated global id of a
        /// cell that uses this particular node.
        /// </summary>
        [Serializable]
        private struct NodeCellIndexPair {
            public long NodeId;
            public long GlobalCellIndex;
        }

        /// <summary>
        /// Helper that wraps a node index with a list of global ids of cells
        /// that share this node
        /// </summary>
        [Serializable]
        private struct NodeCellListPair {
            public long NodeId;
            public long[] CellList;
        }


        class NodeCellIndexPair_ContainerClass {
            public List<NodeCellIndexPair> list = new List<NodeCellIndexPair>();


            /// <summary>
            /// Optimized serialization;
            /// </summary>
            /// <remarks>
            /// Fk, 08sep22: especially in MPI-parallel (64 cores or more) setups, 
            /// the <see cref="SerialisationMessenger"/> produces memory-usage peaks which sometimes even crash the simulation due to out-of-memory.
            /// </remarks>
            public long[] Serialize() {
                var R = new List<long>();
                long L = list.Count;
                R.Add(L);
                for (int l = 0; l < L; l++) {
                    var item_l = list[l];
                    R.Add(item_l.NodeId);
                    R.Add(item_l.GlobalCellIndex);
                }

                return R.ToArray();
            }


            /// <summary>
            /// Optimized de-serialization
            /// </summary>
            public static NodeCellIndexPair_ContainerClass Deserialize(long[] stream) {
                int cnt = 0;
                long L = stream[cnt]; cnt++;
                var R = new NodeCellIndexPair_ContainerClass();
                R.list = new List<NodeCellIndexPair>((int)L);
                for (int l = 0; l < L; l++) {
                    long _NodeID = stream[cnt]; cnt++;
                    long _GlobalCellIndex = stream[cnt]; cnt++;
                    R.list.Add(new NodeCellIndexPair() {
                        NodeId = _NodeID,
                        GlobalCellIndex = _GlobalCellIndex
                    });
                }

                return R;
            }
        }

        [Serializable]
        class NodeCellListPair_ContainerClass {
            public List<NodeCellListPair> list = new List<NodeCellListPair>();

            /// <summary>
            /// Optimized serialization;
            /// </summary>
            /// <remarks>
            /// Fk, 08sep22: especially in MPI-parallel (64 cores or more) setups, 
            /// the <see cref="SerialisationMessenger"/> produces memory-usage peaks which sometimes even crash the simulation due to out-of-memory.
            /// </remarks>
            public long[] Serialize() {
                var R = new List<long>();
                long L = list.Count;
                R.Add(L);
                for (int l = 0; l < L; l++) {
                    var item_l = list[l];
                    R.Add(item_l.NodeId);

                    long LL = item_l.CellList?.Length ?? 0;
                    R.Add(LL);
                    for (int ll = 0; ll < LL; ll++) {
                        R.Add(item_l.CellList[ll]);
                    }
                }

                return R.ToArray();
            }


            /// <summary>
            /// Optimized de-serialization
            /// </summary>
            public static NodeCellListPair_ContainerClass Deserialize(long[] stream) {
                int cnt = 0;
                long L = stream[cnt]; cnt++;
                var R = new NodeCellListPair_ContainerClass();
                R.list = new List<NodeCellListPair>((int)L);
                for (int l = 0; l < L; l++) {
                    long _NodeID = stream[cnt]; cnt++;
                    long LL = stream[cnt]; cnt++;

                    var _CellList = new long[LL];
                    for (int ll = 0; ll < LL; ll++) {
                        _CellList[ll] = stream[cnt]; cnt++;
                    }

                    R.list.Add(new NodeCellListPair() {
                        NodeId = _NodeID,
                        CellList = _CellList
                    });
                }

                return R;
            }
        }

        /// <summary>
        /// Computes the neighbor cells globally (i.e. over all MPI processors) for each local cell.
        /// </summary>
        /// <param name="IncludeBcCells">
        /// If true, also the boundary condition cells (<see cref="BcCells"/>) will be included in the output array.
        /// </param>
        /// <param name="FilterPeriodicDuplicities">
        /// Relevant in case of periodic boundary conditions with one or two cells in periodic direction;
        /// In such cases, the following can occur:
        /// - for one cell in periodic direction: the cell is its own neighbor
        /// - for two cells in periodic direction: the cell has two edges with the same neigbor, i.e. a normal one and a periodic one
        /// Both cases violate the definition of a undirected graph, i.e. there can only be one or zero edge between to cells (aka. vertexes in the sense of a ) and any edge must be between two cells.
        /// 
        /// If set to true, these suckers will be taken out.
        /// </param>
        /// <returns>
        /// Cell-wise neighborship information:
        /// - 1st index: local cell index <em>j</em>, i.e. correlates with <see cref="Cells"/>; if <paramref name="IncludeBcCells"/> is true,
        ///   the information for boundary cells is added after the information for cells.
        /// - 2nd index: enumeration; for the index <em>j</em> the set of neighbor cells. If the global index (<see cref="Neighbour.Neighbour_GlobalIndex"/>)
        ///   is greater or equal than the global number of cells (<see cref="NumberOfCells"/>) the neighbor is a boundary condition cell,
        ///   (<see cref="BcCells"/>).
        /// </returns>
        public Neighbour[][] GetCellNeighbourship(bool IncludeBcCells, bool FilterPeriodicDuplicities) {
            ilPSP.MPICollectiveWatchDog.Watch();
            using (var tr = new FuncTrace()) {
                //tr.InfoToConsole = true;

                var ftNeigh = GetFaceTagsNeigbourIndices(IncludeBcCells);

                var NPart = this.NodePartitioning;
                int K = NPart.LocalLength;
                long k0 = NPart.i0;
                int J = this.NoOfUpdateCells;
                int J_BC = IncludeBcCells ? this.NoOfBcCells : 0;
                long j0 = this.CellPartitioning.i0;
                long Jglob = this.CellPartitioning.TotalLength;
                long j0Bc = this.BcCellPartitioning.i0;
                int mpiRank = this.MyRank;


                Element GetCell(int j) {
                    Element Cell_j;
                    //RefElement Kref;
                    //long jCell_glob;
                    if (j < J) {
                        Cell_j = this.Cells[j];
                        //    Kref = this.m_RefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                        //    jCell_glob = j + j0;
                    } else {
                        Cell_j = this.BcCells[j - J];
                        //    Kref = this.m_EdgeRefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                        //    jCell_glob = (j - J) + j0Bc + Jglob;
                    }
                    return Cell_j;
                }

                RefElement GetRefElement(int j) {
                    RefElement Kref;
                    var Cell_j = GetCell(j);
                    if (j < J) {
                        Kref = this.m_RefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                    } else {
                        Kref = this.m_EdgeRefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                    }
                    return Kref;
                }

                long GetGlobalCellIdx(int j) {
                    long jCell_glob;
                    if (j < J) {
                        jCell_glob = j + j0;
                    } else {
                        jCell_glob = (j - J) + j0Bc + Jglob;
                    }
                    return jCell_glob;
                }

                void AssertNeighborUniqueness(int j, long neighGlIdx) {
                    if (GetCellNeighbours(j).Where(neighEntry => neighEntry.Neighbour_GlobalIndex == neighGlIdx).Count() > 0) {
                        throw new ApplicationException("Neighbor already added");
                    }
                }


                // PART 1: local work
                // ==================
                //
                // In order to reduce MPI communication, we first establish all cell-neighborship which is locally on this processor.
                // This first requires a compression of the global nodes into a local range.


                // compress global node id's to a local range
                // ------------------------------------------
                int[][] Cells_LocalNodeIndices = new int[J + J_BC][]; // 1st index: local cell index;
                                                                      // 2nd index: correlates with respective `Element.NodeIndices`;
                                                                      // content: node index, translated to local index range

                int NoOfLocalNodes;
                {
                    var Global2LocalNode = new SortedDictionary<long, int>(); // key: global node; value: local node idx.; using a SordedDict should lead to an n*log(n) overall runtime behavior
                    int LocalNodeCounter = 0;

                    for (int j = 0; j < (J + J_BC); j++) {
                        Element Cell_j = GetCell(j);
                        var CellNodes = Cell_j.NodeIndices;
                        int L = CellNodes.Length;
                        var CellsLocalNodes_j = new int[L];
                        Cells_LocalNodeIndices[j] = CellsLocalNodes_j;


                        for (int l = 0; l < L; l++) {
                            long GlobalNodeId = CellNodes[l];
                            if (!Global2LocalNode.TryGetValue(GlobalNodeId, out int locIdx)) {
                                Global2LocalNode.Add(GlobalNodeId, LocalNodeCounter);
                                CellsLocalNodes_j[l] = LocalNodeCounter;
                                LocalNodeCounter++;
                            } else {
                                CellsLocalNodes_j[l] = locIdx;
                            }
                        }
                    }

                    NoOfLocalNodes = LocalNodeCounter;
                }

                // build a mapping from local nodes to local cells
                // -----------------------------------------------

                var LocalNodes2Cells = new List<int>[NoOfLocalNodes];
                {
                    for (int j = 0; j < (J + J_BC); j++) {

                        //var CellNodes = Cell_j.NodeIndices;
                        var LocalCellNodes = Cells_LocalNodeIndices[j];
                        foreach (int LocalNodeIdx in LocalCellNodes) {

                            var Cells4Node = LocalNodes2Cells[LocalNodeIdx];
                            if (Cells4Node == null) {
                                Cells4Node = new List<int>();
                                LocalNodes2Cells[LocalNodeIdx] = Cells4Node;
                            }

                            Cells4Node.Add(j);
                        }
                    }
                }

                List<Neighbour>[] CellNeighbours; // index: local cell index
                List<Neighbour> GetCellNeighbours(int j) {
                    List<Neighbour> Cell_j_Neighs = CellNeighbours[j];
                    if (Cell_j_Neighs == null) {
                        Cell_j_Neighs = new List<Neighbour>();
                        CellNeighbours[j] = Cell_j_Neighs;
                    }
                    return Cell_j_Neighs;
                }

                bool[,] FaceIsDone; // 1st index: local cell index; 2nd index: face index
                {
                    CellNeighbours = new List<Neighbour>[J + J_BC];
                    FaceIsDone = new bool[J + J_BC, RefElements.Max(_kref => _kref.NoOfFaces)];

                    List<int> GetNeighborsForCell(int j, int _iface) {
                        if (j >= J)
                            throw new ArgumentOutOfRangeException();
                        Element Cell_j = GetCell(j);
                        RefElement Kref = GetRefElement(j);


                        var faceVtxS = Kref.FaceToVertexIndices;
                        //long[][] B = new long[faceVtx.GetLength(1)][];

                        int NoOfVtxPerFace = faceVtxS.GetLength(1);

                        // the following for-loop intersects the 
                        // cells for each node of the face:
                        int[] LocalCellNodes = Cells_LocalNodeIndices[j];

                        var PossiblePeers = new List<int>();
                        for (int iVtx = 0; iVtx < NoOfVtxPerFace; iVtx++) {
                            int LocalNodeIdx = LocalCellNodes[faceVtxS[_iface, iVtx]];
                            var Cells4Vtx = LocalNodes2Cells[LocalNodeIdx];

                            // since we are dealing with low loop lengths here (maybe 10), the complexity
                            // of the following algorithms should be ok

                            if (iVtx == 0) {
                                PossiblePeers.AddRange(Cells4Vtx);
                                PossiblePeers.Remove(j); // don't consider cell 'j' itself
                            } else {
                                for (int kk = 0; kk < PossiblePeers.Count; kk++) {
                                    int jNeigh = PossiblePeers[kk];
                                    if (!Cells4Vtx.Contains(jNeigh)) {
                                        PossiblePeers.RemoveAt(kk);
                                        kk--;
                                    }
                                }

                            }

                        }

                        return PossiblePeers;
                    }

                    List<int> GetNeighborsForBcCell(int j) {
                        if (j < J || j >= J + J_BC)
                            throw new ArgumentOutOfRangeException();

                        Element Cell_j = GetCell(j);
                        RefElement Kref = GetRefElement(j);

                        int[] LocalCellNodes = Cells_LocalNodeIndices[j];

                        var PossiblePeers = new List<int>();
                        int NoOfVertices = Kref.NoOfVertices;
                        for (int iVtx = 0; iVtx < NoOfVertices; iVtx++) {
                            int LocalNodeIdx = LocalCellNodes[iVtx];
                            var Cells4Vtx = LocalNodes2Cells[LocalNodeIdx];

                            // since we are dealing with low loop lengths here (maybe 10), the complexity
                            // of the following algorithms should be ok

                            if (iVtx == 0) {
                                PossiblePeers.AddRange(Cells4Vtx);
                                PossiblePeers.Remove(j); // don't consider cell 'j' itself
                            } else {
                                for (int kk = 0; kk < PossiblePeers.Count; kk++) {
                                    int jNeigh = PossiblePeers[kk];
                                    if (!Cells4Vtx.Contains(jNeigh)) {
                                        PossiblePeers.RemoveAt(kk);
                                        kk--;
                                    }
                                }

                            }
                        }
                        return PossiblePeers;
                    }

                    bool ConfirmNeigbor(int j, int jNeigh, out int face_neigh) {

                        if (jNeigh < J) {
                            var KrefNeigh = GetRefElement(jNeigh);
                            int NoOfFaces = KrefNeigh.NoOfFaces;

                            for (int iFace = 0; iFace < NoOfFaces; iFace++) {
                                var PossNeigh = GetNeighborsForCell(jNeigh, iFace);
                                if (PossNeigh.Contains(j)) {
                                    face_neigh = iFace;
                                    return true;
                                }

                            }

                            face_neigh = -1;
                            return false;

                        } else {
                            var PossNeigh = GetNeighborsForBcCell(jNeigh);
                            face_neigh = 0;
                            return PossNeigh.Contains(j);
                        }


                    }



                    void AddNeighbor(int j, int jNeig, int _iface) {
                        bool isOk = ConfirmNeigbor(j, jNeig, out int face_neigh);
                        if (isOk) {



                            {
                                Neighbour nCN = default(Neighbour);
                                nCN.Neighbour_GlobalIndex = GetGlobalCellIdx(jNeig);
                                nCN.CellFaceTag.FaceIndex = _iface;
                                nCN.CellFaceTag.ConformalNeighborship = true;

                                AssertNeighborUniqueness(j, nCN.Neighbour_GlobalIndex);


                                GetCellNeighbours(j).Add(nCN);
                                FaceIsDone[j, _iface] = true;
                            }

                            {
                                Neighbour nCN = default(Neighbour);
                                nCN.Neighbour_GlobalIndex = GetGlobalCellIdx(j);
                                nCN.CellFaceTag.FaceIndex = face_neigh;
                                nCN.CellFaceTag.ConformalNeighborship = true;

                                AssertNeighborUniqueness(jNeig, nCN.Neighbour_GlobalIndex);

                                GetCellNeighbours(jNeig).Add(nCN);
                                FaceIsDone[jNeig, face_neigh] = true;
                            }

                        } else {
                            string ErrString;
                            try {
                                ErrString = $"error in mesh; cell {GetCell(j)} is neighbor to cell {GetCell(jNeig)}, but not the other way around.";
                            } catch (Exception e) {
                                ErrString = $"error in mesh; further error in formatting error message for connection between cells {j} and {jNeig}; {e.GetType()} {e.Message}, ";
                            }
                            tr.Error(ErrString);
                            throw new ArgumentException(ErrString);
                        }
                    }


                    // find neighbor cells connected via grid nodes
                    // - - - - - - - - - - - - - - - - - - - - - - - 
                    for (int j = 0; j < J + J_BC; j++) { // loop over cells
                        Element Cell_j = GetCell(j);
                        RefElement Kref = GetRefElement(j);

                        if (j < J) {
                            //normal cells: match faces

                            var faceVtxS = Kref.FaceToVertexIndices;
                            //long[][] B = new long[faceVtx.GetLength(1)][];

                            int NoOfVtxPerFace = faceVtxS.GetLength(1);

                            for (int _iface = 0; _iface < Kref.NoOfFaces; _iface++) { // loop over faces of cell 'j' (local index) resp. 'jCellGlob' (global index)
                                if (FaceIsDone[j, _iface])
                                    continue;

                                var PossiblePeers = GetNeighborsForCell(j, _iface);

                                foreach (int jNeig in PossiblePeers) {
                                    AddNeighbor(j, jNeig, _iface);
                                }
                            }
                        } else {
                            // boundary-condition cell: match the whole element

                            if (FaceIsDone[j, 0])
                                continue;

                            var PossiblePeers = GetNeighborsForBcCell(j);

                            foreach (int jNeig in PossiblePeers) {
                                AddNeighbor(j, jNeig, 0);
                            }
                        }
                    }

                    // find neighbor cells connected via CellFaceTag's
                    // - - - - - - - - - - - - - - - - - - - - - - - - 
                    for (int j = 0; j < J + J_BC; j++) { // loop over cells
                        Element Cell_j = GetCell(j);

                        var otherNeighbours = ftNeigh[j]; // ftNeigh is the result of CellFaceTag-based connectivity
                        if (j < J) {
                            var _Cell_j = (Cell)Cell_j;
                            var Cell_j_Neighs = GetCellNeighbours(j);
                            Debug.Assert(((otherNeighbours == null ? 0 : otherNeighbours.Length) == ((_Cell_j.CellFaceTags == null) ? 0 : _Cell_j.CellFaceTags.Length)));
                            if (otherNeighbours != null) {
                                for (int w = 0; w < otherNeighbours.Length; w++) {
                                    Debug.Assert(_Cell_j.CellFaceTags[w].NeighCell_GlobalID < 0 == otherNeighbours[w] < 0);
                                    if (_Cell_j.CellFaceTags[w].NeighCell_GlobalID >= 0) {
                                        // a connection to 


                                        bool CompareNeighbors_FilterSelfPeriodiodic(Neighbour neigh) {
                                            return neigh.Neighbour_GlobalIndex == otherNeighbours[w];
                                        }

                                        bool CompareNeighbors_IncludeSelfPeriodic(Neighbour neigh) {
                                            if (neigh.Neighbour_GlobalIndex != otherNeighbours[w])
                                                return false;

                                            if (neigh.IsPeriodicNeighbour != _Cell_j.CellFaceTags[w].IsPeriodicNeighbour)
                                                return false;

                                            if (neigh.IsPeriodicNeighbour && _Cell_j.CellFaceTags[w].IsPeriodicNeighbour) {
                                                if (neigh.CellFaceTag.PeriodicInverse != _Cell_j.CellFaceTags[w].PeriodicInverse)
                                                    return false;
                                            }

                                            return true;
                                        }

                                        Func<Neighbour, bool> CompareNeighbors;
                                        if (FilterPeriodicDuplicities)
                                            CompareNeighbors = CompareNeighbors_FilterSelfPeriodiodic;
                                        else
                                            CompareNeighbors = CompareNeighbors_IncludeSelfPeriodic;


                                        if (Cell_j_Neighs.Where(CompareNeighbors).Count() <= 0) { // filter duplicates
                                            var nCN = new Neighbour() {
                                                Neighbour_GlobalIndex = otherNeighbours[w],
                                                CellFaceTag = _Cell_j.CellFaceTags[w],
                                            };

                                            if (FilterPeriodicDuplicities)
                                                AssertNeighborUniqueness(j, nCN.Neighbour_GlobalIndex);
                                            Cell_j_Neighs.Add(nCN);


                                        } else {
                                            //Console.WriteLine("some duplicate found"); artifact from debugging
                                        }
                                    }
                                }
                            }
                        } else {
                            var BcCell_j = (BCElement)Cell_j;
                            var Cell_j_Neighs = GetCellNeighbours(j);
                            Debug.Assert(((otherNeighbours == null ? 0 : otherNeighbours.Length) == ((BcCell_j.NeighCell_GlobalIDs == null) ? 0 : BcCell_j.NeighCell_GlobalIDs.Length)));

                            if (otherNeighbours != null) {
                                for (int w = 0; w < otherNeighbours.Length; w++) {
                                    Cell_j_Neighs.Add(new Neighbour() {
                                        Neighbour_GlobalIndex = otherNeighbours[w],
                                        CellFaceTag = new CellFaceTag() {
                                            EdgeTag = BcCell_j.EdgeTag,
                                            FaceIndex = int.MinValue,
                                            NeighCell_GlobalID = BcCell_j.NeighCell_GlobalIDs[w],
                                            ConformalNeighborship = BcCell_j.Conformal
                                        }
                                    });
                                }
                            }
                        }
                    }
                }

                // PART 2: MPI-global work
                // =======================

                // Which cells make use of a particular node?
                //-------------------------------------------

                // Index: local Node index
                // Entry: Enumeration of global indices of cells that use this particular node
                List<long>[] Nodes2Cells = new List<long>[K];
                {
                    for (int k = 0; k < K; k++) {
                        Nodes2Cells[k] = new List<long>();
                    }

                    // key: MPI processor rank
                    // value: information packet
                    Dictionary<int, NodeCellIndexPair_ContainerClass> Y = new Dictionary<int, NodeCellIndexPair_ContainerClass>();

                    long NodeIdMin = long.MaxValue, NodeIdMax = long.MinValue;

                    for (int j = 0; j < (J + J_BC); j++) {
                        RefElement Kref = GetRefElement(j);
                        int NoOfFaces = Kref.NoOfFaces;

                        if (j < J) {
                            bool Alldone = true;
                            for (int iF = 0; iF < NoOfFaces; iF++) {
                                Alldone = Alldone && FaceIsDone[j, iF];
                            }
                            if (Alldone)
                                continue;
                        } else {
                            if (FaceIsDone[j, 0])
                                continue;
                        }


                        Element Cell_j = GetCell(j);
                        long jCell_glob = GetGlobalCellIdx(j);

                        int NoOfVtx = Kref.NoOfVertices;
                        var CellNodes = Cell_j.NodeIndices;
                        //var CellNodesLocal = Cells_LocalNodeIndices[j];
                        if (CellNodes.Length != NoOfVtx)
                            throw new ApplicationException();
                        //if (CellNodesLocal.Length != NoOfVtx)
                        //    throw new ApplicationException();

                        if (CellNodes.Length != Kref.NoOfVertices) {
                            throw new ApplicationException("error in data structure.");
                        }

                        bool[] VtxToInclude = new bool[NoOfVtx];
                        int[,] Face2Vertex = Kref.FaceToVertexIndices;
                        int NoOfFaceVtx = Face2Vertex.GetLength(1);

                        if (j < J) {
                            for (int iF = 0; iF < NoOfFaces; iF++) {
                                if (!FaceIsDone[j, iF]) {
                                    for (int iFaceVtx = 0; iFaceVtx < NoOfFaceVtx; iFaceVtx++) {
                                        VtxToInclude[Face2Vertex[iF, iFaceVtx]] = true;
                                    }
                                }
                            }
                        } else {
                            VtxToInclude.SetAll(true);
                        }



                        for (int k = 0; k < NoOfVtx; k++) {
                            if (!VtxToInclude[k])
                                continue;

                            //Cells_LocalNodeIndices[j] ;
                            long NodeId = CellNodes[k];
                            //int NodeIdLocal = CellNodesLocal[k];

                            //foreach (long NodeId in CellNodes) {
                            NodeIdMin = Math.Min(NodeIdMin, NodeId);
                            NodeIdMax = Math.Max(NodeIdMax, NodeId);

                            int target_prozi = NPart.FindProcess(NodeId);
                            if (target_prozi == mpiRank) {
                                Nodes2Cells[NodeId - k0].Add(jCell_glob);
                            } else {
                                NodeCellIndexPair Packet;
                                Packet.NodeId = NodeId;
                                Packet.GlobalCellIndex = jCell_glob;

                                if (!Y.TryGetValue(target_prozi, out var Z)) {
                                    Z = new NodeCellIndexPair_ContainerClass();
                                    Y.Add(target_prozi, Z);
                                }

                                Z.list.Add(Packet);
                            }
                        }

                    }

                    tr.Info($"r{mpiRank}: min Node {NodeIdMin} -- max Node {NodeIdMax} // Global No Of Nodes: {NodePartitioning.TotalLength}");


                    foreach (int targProc in Y.Keys) {
                        var item = Y[targProc];
                        tr.Info($"{mpiRank}to{targProc}: nodecellpair {item.list.Count} ");
                    }

                    var W = ArrayMessenger<long>.ExchangeData(Y.Keys, rank => Y[rank].Serialize(), csMPI.Raw._COMM.WORLD);
                    foreach (var wp in W.Values) {
                        var container = NodeCellIndexPair_ContainerClass.Deserialize(wp);

                        foreach (NodeCellIndexPair Packet in container.list) {
                            Nodes2Cells[Packet.NodeId - k0].Add(Packet.GlobalCellIndex);
                        }
                    }


                }

                // For every cell, for every vertex in this cell:
                // Which other cells do also use this node?
                //-----------------------------------------------

                // 1st index: Local cell index
                // 2nd index: Cell vertex index
                // 3rd index: enumeration of 'peer' cells
                // content: global cell index
                long[][][] NodePeers = new long[J + J_BC][][];
                {
                    for (int j = 0; j < J + J_BC; j++) {
                        Element Cell_j = GetCell(j);
                        NodePeers[j] = new long[Cell_j.NodeIndices.Length][];
                    }

                    Dictionary<int, NodeCellListPair_ContainerClass> Y = new Dictionary<int, NodeCellListPair_ContainerClass>();

                    var CPart = this.CellPartitioning;
                    var BcPart = this.BcCellPartitioning;
                    for (int k = 0; k < K; k++) { // loop over locally assigned nodes
                        long k_node = k + k0;

                        if (Nodes2Cells[k] != null) {
                            var cell_list = Nodes2Cells[k].ToArray();
                            foreach (int jCell in cell_list) { // loop over all cells that use node 'k'
                                int cell_proc;
                                long local_offset;
                                if (jCell < Jglob) {
                                    // normal cell
                                    cell_proc = CPart.FindProcess(jCell);
                                    local_offset = j0;
                                } else {
                                    // boundary condition cell
                                    cell_proc = BcPart.FindProcess(jCell - Jglob);
                                    local_offset = Jglob + j0Bc;
                                }


                                if (cell_proc == mpiRank) {
                                    int jCell_loc = checked((int)(jCell - local_offset));
                                    int kC;
                                    bool bfound = false;

                                    Element Cell_j;
                                    int oo;
                                    if (jCell < Jglob) {
                                        // normal cell
                                        Cell_j = this.Cells[jCell_loc];
                                        oo = 0;
                                    } else {
                                        // boundary condition cell
                                        Cell_j = this.BcCells[jCell_loc];
                                        oo = J;
                                    }

                                    for (kC = 0; kC < Cell_j.NodeIndices.Length; kC++) {
                                        if (Cell_j.NodeIndices[kC] == k_node) {
                                            bfound = true;
                                            break;
                                        }
                                    }

                                    if (!bfound)
                                        throw new ApplicationException("error in algorithm.");

                                    NodePeers[jCell_loc + oo][kC] = cell_list;
                                } else {
                                    NodeCellListPair A;
                                    A.NodeId = k_node;
                                    A.CellList = cell_list;


                                    if (!Y.TryGetValue(cell_proc, out var Z)) {
                                        Z = new NodeCellListPair_ContainerClass();
                                        Y.Add(cell_proc, Z);
                                    }

                                    Z.list.Add(A);
                                }
                            }
                        }
                    }

                    //var Yexc = new Dictionary<int, NodeCellListPair_ContainerClass>();
                    //{
                    //    foreach (var kv in Y) {
                    //        Yexc.Add(kv.Key,
                    //            new NodeCellListPair_ContainerClass() { list = kv.Value.ToArray() });
                    //    }
                    //    Y = null;
                    //}


                    //tr.Info($"{mpiRank}: Glob No Of nodes: {NodePartitioning.TotalLength}; local no: {NodePartitioning.LocalLength}");
                    //foreach (int targProc in Yexc.Keys) {
                    //    NodeCellListPair_ContainerClass item = Yexc[targProc];
                    //    int NoOfEntries = item.list.Sum(entry => entry.CellList.Length);
                    //    tr.Info($"{mpiRank}to{targProc}: node-cell-list {item.list.Length} items with {NoOfEntries} entries");
                    //}

                    var W = ArrayMessenger<long>.ExchangeData(Y.Keys, rank => Y[rank].Serialize(), csMPI.Raw._COMM.WORLD);
                    foreach (var wp in W.Values) {
                        var container = NodeCellListPair_ContainerClass.Deserialize(wp);
                        foreach (var P in container.list) {

                            long k_node = P.NodeId;
                            long[] cell_list = P.CellList;

                            foreach (long jCell in cell_list) {
                                int cell_proc;
                                long local_offset;
                                if (jCell < Jglob) {
                                    // normal cell
                                    cell_proc = CPart.FindProcess(jCell);
                                    local_offset = j0;
                                } else {
                                    // boundary condition cell
                                    cell_proc = BcPart.FindProcess(jCell - Jglob);
                                    local_offset = Jglob + j0Bc;
                                }

                                if (cell_proc == mpiRank) {

                                    int jCell_loc = checked((int)(jCell - local_offset));

                                    Element Cell_j;
                                    int oo;
                                    if (jCell < Jglob) {
                                        // normal cell
                                        Cell_j = this.Cells[jCell_loc];
                                        oo = 0;
                                    } else {
                                        // boundary condition cell
                                        Cell_j = this.BcCells[jCell_loc];
                                        oo = J;
                                    }

                                    int kC;
                                    bool bfound = false;
                                    for (kC = 0; kC < Cell_j.NodeIndices.Length; kC++) {
                                        if (Cell_j.NodeIndices[kC] == k_node) {
                                            bfound = true;
                                            break;
                                        }
                                    }

                                    if (!bfound)
                                        throw new ApplicationException("error in algorithm.");

                                    NodePeers[jCell_loc + oo][kC] = cell_list;
                                }
                            }
                        }
                    }
                }

                // Assemble final result
                // ---------------------

                //IEnumerable<Neighbour>[] CellNeighbours;
                {
                    //CellNeighbours = new IEnumerable<Neighbour>[J + J_BC];
                    for (int j = 0; j < J + J_BC; j++) { // loop over cells
                        Element Cell_j = GetCell(j);
                        RefElement Kref = GetRefElement(j);
                        long jCellGlob = GetGlobalCellIdx(j);

                        List<Neighbour> Cell_j_Neighs = GetCellNeighbours(j);

                        // find neighbor cells connected via grid nodes
                        // - - - - - - - - - - - - - - - - - - - - - - - 

                        if (j < J) {
                            //normal cells: match faces


                            var faceVtx = Kref.FaceToVertexIndices;
                            long[][] B = new long[faceVtx.GetLength(1)][]; // 1st index: face vertex; 2nd index: enumeration
                            for (int _iface = 0; _iface < Kref.NoOfFaces; _iface++) { // loop over faces of cell 'j' (local index) resp. 'jCellGlob' (global index)
                                if (FaceIsDone[j, _iface])
                                    continue;
                                for (int iv = 0; iv < B.Length; iv++)
                                    B[iv] = NodePeers[j][faceVtx[_iface, iv]];

                                long NeighIdx = Intersect(B, jCellGlob);
                                if (NeighIdx >= 0) {
                                    int iFound = Cell_j_Neighs.FirstIndexWhere(CN => CN.Neighbour_GlobalIndex == NeighIdx);

                                    if (iFound >= 0) {
                                        Neighbour nCN = Cell_j_Neighs[iFound];
                                        if (nCN.CellFaceTag.FaceIndex != _iface)
                                            throw new ApplicationException("Found two connections between cells with different face indices");

                                        nCN.CellFaceTag.ConformalNeighborship = true;
                                        Cell_j_Neighs[iFound] = nCN;

                                    } else {

                                        Neighbour nCN = default(Neighbour);
                                        nCN.Neighbour_GlobalIndex = NeighIdx;
                                        nCN.CellFaceTag.FaceIndex = _iface;
                                        nCN.CellFaceTag.ConformalNeighborship = true;

                                        Cell_j_Neighs.Add(nCN);
                                    }
                                }
                            }
                        } else {
                            // boundary-condition cell: match the whole element
                            if (FaceIsDone[j, 0])
                                continue;

                            long[][] B = new long[Kref.NoOfVertices][];
                            for (int iv = 0; iv < B.Length; iv++)
                                B[iv] = NodePeers[j][iv];

                            long NeighIdx = Intersect(B, jCellGlob);
                            if (NeighIdx >= 0) {
                                int iFound = Cell_j_Neighs.FirstIndexWhere(CN => CN.Neighbour_GlobalIndex == NeighIdx);

                                if (iFound >= 0) {
                                    Neighbour nCN = Cell_j_Neighs[iFound];

                                    nCN.CellFaceTag.ConformalNeighborship = true;
                                    Cell_j_Neighs[iFound] = nCN;

                                } else {
                                    Neighbour nCN = default(Neighbour);
                                    nCN.Neighbour_GlobalIndex = NeighIdx;
                                    nCN.CellFaceTag.FaceIndex = -1;
                                    nCN.CellFaceTag.ConformalNeighborship = true;

                                    AssertNeighborUniqueness(j, nCN.Neighbour_GlobalIndex);

                                    Cell_j_Neighs.Add(nCN);
                                }
                            }

                        }

                    }
                }

                // PART 3: Return
                // ==============

                {
                    int rL = CellNeighbours.Length;
                    var R = new Neighbour[CellNeighbours.Length][];
                    for (int i = 0; i < rL; i++) {
                        R[i] = CellNeighbours[i]?.ToArray() ?? new Neighbour[0];
                        CellNeighbours[i] = null;

                    }
                    return R;
                }
            }
        }

        /// <summary>
        /// Finds conformal cell neighborings;
        /// </summary>
        /// <param name="B">
        /// For each vertex of some cell face, the global indices of all cells
        /// which use that vertex
        /// - 1st index: face vertex index
        /// - 2nd index: enumeration
        /// - content: global cell index
        /// </param>
        /// <param name="j_cell_myself">
        /// global cell index
        /// </param>
        /// <returns>
        /// Either the global index of a neighbor cell
        /// (this is the set-intersection of all <paramref name="B"/>[i] - sets, excluding <paramref name="j_cell_myself"/> and should have at most 1 element), or a negative number if
        /// there is no neighbor.
        /// </returns>
        static long Intersect(long[][] B, long j_cell_myself) {
            long R;
            long[] B0 = B[0];
            long ret = long.MinValue;
            for (int l = 0; l < B0.Length; l++) { // loop over all 'candidate cell '
                R = B0[l]; // global cell index of 'candidate cell'
                if (R == j_cell_myself)
                    continue;

                bool AllPassed = true;
                for (int j = 1; j < B.Length; j++) { // test if 'candidate cell' is also a neighbor for all other vertices
                    bool bfound = false;
                    long[] Bj = B[j];
                    for (int k = Bj.Length - 1; k >= 0; k--) {
                        if (Bj[k] == R) {
                            bfound = true;
                            break;
                        }
                    }

                    if (!bfound) {
                        AllPassed = false;
                        break;
                    }
                }

                if (AllPassed) {
                    if (ret >= 0)
                        throw new ApplicationException("Found two conformal neighbors - this can't be.");
                    ret = R;
                }
            }

            return ret;
        }

        /// <summary>
        /// converts the GlobalID-based neighborship information (GlobalId's
        /// of the face tags, see<see cref="Cell.CellFaceTags"/> resp.
        /// <see cref="CellFaceTag.NeighCell_GlobalID"/>)
        /// into global cell indices
        /// </summary>
        /// <returns>
        /// content: global cell indices;
        ///  - 1st index: local cell index
        ///  - 2nd index: correlates with the ordering of <see cref="Cell.CellFaceTags"/>.
        /// </returns>
        int[][] GetFaceTagsNeigbourIndices(bool IncludeBcCells) {
            using (new FuncTrace()) {
                var GidInv = this.GetInverseGlobalIDPermutation(IncludeBcCells);
                long GlNoOfCells = this.NumberOfCells_l;
                List<long> gids = new List<long>();
                int J = this.NoOfUpdateCells;
                int J_BC = IncludeBcCells ? this.NoOfBcCells : 0;
                for (int j = 0; j < J; j++) {
                    CellFaceTag[] cfts_j = this.Cells[j].CellFaceTags;

                    if (cfts_j != null && cfts_j.Length > 0) {
                        int K = cfts_j.Length;
                        for (int k = 0; k < K; k++) {
                            long neig_gid = cfts_j[k].NeighCell_GlobalID;
                            if (neig_gid >= 0)
                                gids.Add(neig_gid);
                        }
                    }
                }
                for (int j = 0; j < J_BC; j++) {
                    long[] NeighGids = this.BcCells[j].NeighCell_GlobalIDs;

                    if (NeighGids != null && NeighGids.Length > 0) {
                        int K = NeighGids.Length;

                        for (int k = 0; k < K; k++) {
                            long neig_gid = NeighGids[k];
                            if (neig_gid < 0 || neig_gid >= GlNoOfCells)
                                throw new ApplicationException("illegal GlobalId for boundary-condition cell.");
                            gids.Add(neig_gid);
                        }
                    }
                }


                long[] Gidx = new long[gids.Count];
                GidInv.EvaluatePermutation(gids, Gidx);

                //List<long> work = new List<long>();
                var R = new int[J + J_BC][];
                int cnt = 0;
                for (int j = 0; j < J; j++) {
                    var cfts_j = this.Cells[j].CellFaceTags;

                    if (cfts_j != null && cfts_j.Length > 0) {
                        int K = cfts_j.Length;
                        //work.Clear(); 
                        var Rj = new int[K];
                        R[j] = Rj;
                        for (int k = 0; k < K; k++) {
                            long neig_gid = cfts_j[k].NeighCell_GlobalID;
                            if (neig_gid >= 0) {
                                Rj[k] = (int)Gidx[cnt];
                                cnt++;
                            } else {
                                Rj[k] = int.MinValue;
                            }
                        }

                        //if (work.Count > 0)
                        //    R[j] = work.ToArray();
                        //else
                        //    R[j] = null;

                    }
                }
                for (int j = 0; j < J_BC; j++) {
                    long[] NeighGids = this.BcCells[j].NeighCell_GlobalIDs;

                    if (NeighGids != null && NeighGids.Length > 0) {
                        int K = NeighGids.Length;
                        var Rj = new int[K];
                        R[j + J] = Rj;
                        for (int k = 0; k < K; k++) {
                            Rj[k] = (int)Gidx[cnt];
                            cnt++;
                        }
                    }
                }

                return R;
            }
        }

    }
}
