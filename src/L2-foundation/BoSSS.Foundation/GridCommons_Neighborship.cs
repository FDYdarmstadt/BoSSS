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

namespace BoSSS.Foundation.Grid.Classic {

    public partial class GridCommons {

        /// <summary>
        /// Helper that wraps a node index with an associated global id of a
        /// cell that uses this particular node.
        /// </summary>
        [Serializable]
        private struct NodeCellIndexPair {
            public int NodeId;
            public int GlobalCellIndex;
        }

        /// <summary>
        /// Helper that wraps a node index with a list of global ids of cells
        /// that share this node
        /// </summary>
        [Serializable]
        private struct NodeCellListPair {
            public int NodeId;
            public int[] CellList;
        }

        /// <summary>
        /// return values of <see cref="GetCellNeighbourship"/>.
        /// </summary>
        [Serializable]
        public struct Neighbour {

            /// <summary>
            /// global index of neighbor cell.
            /// </summary>
            public int Neighbour_GlobalIndex;

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

            /// <summary>
            /// See <see cref="CellFaceTag.EdgeMayBeEmpty"/>.
            /// </summary>
            public bool EdgeMayBeEmpty;
        }

        /// <summary>
        /// Computes the neighbor cells globally (i.e. over all MPI processors) for each local cell.
        /// </summary>
        /// <returns>
        /// <param name="IncludeBcCells">
        /// If true, also the boundary condition cells (<see cref="BcCells"/>) will be included in the output array.
        /// </param>
        /// Cell-wise neighborship information:
        /// - index: local cell index <em>j</em>, i.e. correlates with <see cref="Cells"/>; if <paramref name="IncludeBcCells"/> is true,
        ///   the information for boundary cells is added after the information for cells.
        /// - content: for the index <em>j</em> the set of neighbor cells. If the global index (<see cref="Neighbour.Neighbour_GlobalIndex"/>)
        ///   is greater or equal than the global number of cells (<see cref="NumberOfCells"/>) the neighbor is a boundary condition cell,
        ///   (<see cref="BcCells"/>).
        /// </returns>
        public IEnumerable<Neighbour>[] GetCellNeighbourship(bool IncludeBcCells) {
            ilPSP.MPICollectiveWatchDog.Watch();
            using (new FuncTrace()) {

                var ftNeigh = GetFaceTagsNeigbourIndices(IncludeBcCells);

                var NPart = this.NodePartitioning;
                int K = NPart.LocalLength;
                int k0 = NPart.i0;
                int J = this.NoOfUpdateCells;
                int J_BC = IncludeBcCells ? this.NoOfBcCells : 0;
                int j0 = this.CellPartitioning.i0;
                int Jglob = this.CellPartitioning.TotalLength;
                int j0Bc = this.BcCellPartitioning.i0;

                // Which cells make use of a particular node?
                //-------------------------------------------

                // Index: Node index
                // Entry: Enumeration of global indices of cells that use this
                // particular node
                List<int>[] Nodes2Cells = new List<int>[K];
                {
                    for (int k = 0; k < K; k++) {
                        Nodes2Cells[k] = new List<int>();
                    }

                    // key: MPI processor rank
                    // value: information packet
                    Dictionary<int, List<NodeCellIndexPair>> Y =
                        new Dictionary<int, List<NodeCellIndexPair>>();

                    for (int j = 0; j < (J + J_BC); j++) {
                        Element Cell_j;
                        RefElement Kref;
                        int jCell_glob;
                        if (j < J) {
                            Cell_j = this.Cells[j];
                            Kref = this.m_RefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                            jCell_glob = j + j0;
                        } else {
                            Cell_j = this.BcCells[j - J];
                            Kref = this.m_EdgeRefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                            jCell_glob = (j - J) + j0Bc + Jglob;
                        }
                        var CellNodes = Cell_j.NodeIndices;

                        if (CellNodes.Length != Kref.NoOfVertices) {
                            throw new ApplicationException("error in data structure.");
                        }

                        foreach (int NodeId in CellNodes) {
                            int target_prozi = NPart.FindProcess(NodeId);
                            if (target_prozi == MyRank) {
                                Nodes2Cells[NodeId - k0].Add(jCell_glob);
                            } else {
                                NodeCellIndexPair Packet;
                                Packet.NodeId = NodeId;
                                Packet.GlobalCellIndex = jCell_glob;

                                List<NodeCellIndexPair> Z;
                                if (!Y.TryGetValue(target_prozi, out Z)) {
                                    Z = new List<NodeCellIndexPair>();
                                    Y.Add(target_prozi, Z);
                                }

                                Z.Add(Packet);
                            }
                        }
                    }

                    var W = SerialisationMessenger.ExchangeData(Y, csMPI.Raw._COMM.WORLD);
                    foreach (var wp in W.Values) {
                        foreach (NodeCellIndexPair Packet in wp) {
                            Nodes2Cells[Packet.NodeId - k0].Add(Packet.GlobalCellIndex);
                        }
                    }
                }

                // For every cell, for every vertex in this cell:
                // Which other cells die also use this node?
                //-----------------------------------------------

                // 1st index: Local cell index
                // 2nd index: Cell vertex index
                // 3rd index: Collection of 'peer' cells
                int[][][] NodePeers = new int[J + J_BC][][];
                {
                    for (int j = 0; j < J + J_BC; j++) {
                        Element Cell_j;
                        if (j < J) {
                            Cell_j = this.Cells[j];
                        } else {
                            Cell_j = this.BcCells[j - J];
                        }
                        NodePeers[j] = new int[Cell_j.NodeIndices.Length][];
                    }

                    Dictionary<int, List<NodeCellListPair>> Y = new Dictionary<int, List<NodeCellListPair>>();

                    var CPart = this.CellPartitioning;
                    var BcPart = this.BcCellPartitioning;
                    for (int k = 0; k < K; k++) { // loop over locally assigned nodes
                        int k_node = k + k0;

                        var cell_list = Nodes2Cells[k].ToArray();
                        foreach (int jCell in cell_list) { // loop over all cells that use node 'k'
                            int cell_proc;
                            int local_offset;
                            if (jCell < Jglob) {
                                // normal cell
                                cell_proc = CPart.FindProcess(jCell);
                                local_offset = j0;
                            } else {
                                // boundary condition cell
                                cell_proc = BcPart.FindProcess(jCell - Jglob);
                                local_offset = Jglob + j0Bc;
                            }


                            if (cell_proc == MyRank) {
                                int jCell_loc = jCell - local_offset;
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

                                List<NodeCellListPair> Z;
                                if (!Y.TryGetValue(cell_proc, out Z)) {
                                    Z = new List<NodeCellListPair>();
                                    Y.Add(cell_proc, Z);
                                }

                                Z.Add(A);
                            }
                        }
                    }

                    var W = SerialisationMessenger.ExchangeData(Y, csMPI.Raw._COMM.WORLD);
                    foreach (var wp in W.Values) {
                        foreach (var P in wp) {
                            int k_node = P.NodeId;
                            int[] cell_list = P.CellList;

                            foreach (int jCell in cell_list) {
                                int cell_proc;
                                int local_offset;
                                if (jCell < Jglob) {
                                    // normal cell
                                    cell_proc = CPart.FindProcess(jCell);
                                    local_offset = j0;
                                } else {
                                    // boundary condition cell
                                    cell_proc = BcPart.FindProcess(jCell - Jglob);
                                    local_offset = Jglob + j0Bc;
                                }

                                if (cell_proc == MyRank) {

                                    int jCell_loc = jCell - local_offset;

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

                IEnumerable<Neighbour>[] CellNeighbours;
                {
                    CellNeighbours = new IEnumerable<Neighbour>[J + J_BC];
                    for (int j = 0; j < J + J_BC; j++) { // loop over cells
                        //var Cell_j = this.Cells[j];
                        //int jCellGlob = j + j0;
                        //var Kref = this.m_GridSimplices.Single(KK => KK.SupportedTypes.Contains(Cell_j.Type));

                        Element Cell_j;
                        RefElement Kref;
                        int jCellGlob;
                        if (j < J) {
                            Cell_j = this.Cells[j];
                            Kref = this.m_RefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                            jCellGlob = j + j0;
                        } else {
                            Cell_j = this.BcCells[j - J];
                            Kref = this.m_EdgeRefElements.Single(KK => KK.SupportedCellTypes.Contains(Cell_j.Type));
                            jCellGlob = (j - J) + j0Bc + Jglob;
                        }

                        var Cell_j_Neighs = new List<Neighbour>();
                        CellNeighbours[j] = Cell_j_Neighs;

                        // find neighbor cells connected via grid nodes
                        // --------------------------------------------

                        if (j < J) {
                            //normal cells: match faces


                            var faceVtx = Kref.FaceToVertexIndices;
                            int[][] B = new int[faceVtx.GetLength(1)][];
                            for (int _iface = 0; _iface < Kref.NoOfFaces; _iface++) { // loop over faces of cell 'j' (local index) resp. 'jCellGlob' (global index)
                                for (int iv = 0; iv < B.Length; iv++)
                                    B[iv] = NodePeers[j][faceVtx[_iface, iv]];

                                int NeighIdx = Intersect(B, jCellGlob);
                                if (NeighIdx >= 0) {
                                    Neighbour nCN = default(Neighbour);
                                    nCN.Neighbour_GlobalIndex = NeighIdx;
                                    nCN.CellFaceTag.FaceIndex = _iface;
                                    nCN.CellFaceTag.ConformalNeighborship = true;

                                    Cell_j_Neighs.Add(nCN);
                                }
                            }
                        } else {
                            // boundary-condition cell: match the whole element

                            int[][] B = new int[Kref.NoOfVertices][];
                            for (int iv = 0; iv < B.Length; iv++)
                                B[iv] = NodePeers[j][iv];

                            int NeighIdx = Intersect(B, jCellGlob);
                            if (NeighIdx >= 0) {
                                Neighbour nCN = default(Neighbour);
                                nCN.Neighbour_GlobalIndex = NeighIdx;
                                nCN.CellFaceTag.FaceIndex = -1;
                                nCN.CellFaceTag.ConformalNeighborship = true;

                                Cell_j_Neighs.Add(nCN);
                            }

                        }

                        // find neighbor cells connected via CellFaceTag's
                        // -----------------------------------------------

                        var otherNeighbours = ftNeigh[j]; // ftNeigh is the result of CellFaceTag-based connectivity
                        if (j < J) {
                            var _Cell_j = (Cell)Cell_j;
                            Debug.Assert(((otherNeighbours == null ? 0 : otherNeighbours.Length) == ((_Cell_j.CellFaceTags == null) ? 0 : _Cell_j.CellFaceTags.Length)));
                            if (otherNeighbours != null) {
                                for (int w = 0; w < otherNeighbours.Length; w++) {
                                    Debug.Assert(_Cell_j.CellFaceTags[w].NeighCell_GlobalID < 0 == otherNeighbours[w] < 0);
                                    if (_Cell_j.CellFaceTags[w].NeighCell_GlobalID >= 0) {
                                        if (Cell_j_Neighs.Where(neigh => neigh.Neighbour_GlobalIndex == otherNeighbours[w]).Count() <= 0) { // filter duplicates
                                            Cell_j_Neighs.Add(new Neighbour() {
                                                Neighbour_GlobalIndex = otherNeighbours[w],
                                                CellFaceTag = _Cell_j.CellFaceTags[w],
                                            });
                                        }
                                    }
                                }
                            }
                        } else {
                            var BcCell_j = (BCElement)Cell_j;
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

                return CellNeighbours;
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
        /// Either the global index of a neighbor cell, or a negative number if
        /// there is no neighbor.
        /// </returns>
        static int Intersect(int[][] B, int j_cell_myself) {
            int R;
            int[] B0 = B[0];
            int ret = int.MinValue;
            for (int l = 0; l < B0.Length; l++) {
                R = B0[l];
                if (R == j_cell_myself)
                    continue;

                bool AllPassed = true;
                for (int j = 1; j < B.Length; j++) {
                    bool bfound = false;
                    int[] Bj = B[j];
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
