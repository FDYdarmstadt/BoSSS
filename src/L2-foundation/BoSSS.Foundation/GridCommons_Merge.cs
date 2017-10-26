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
using System.Threading.Tasks;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform.LinAlg;
using MPI.Wrappers;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {
    partial class GridCommons {


        /// <summary>
        /// Concatenates grids into one data structure an makes use of <see cref="MergeLogically(GridCommons, GridCommons)"/>
        /// </summary>
        /// <param name="grids">list of different grids</param>
        public static GridCommons MergeLogically(params GridCommons[] grids) {
            GridCommons ret = grids.First();
            for (int i = 1; i < grids.Length; i++) {
                ret = MergeLogically(ret, grids[i]);
            }
            return ret;
        }

        /// <summary>
        /// Concatenates grids <paramref name="A"/> and <paramref name="B"/> into one data structure.
        /// </summary>
        public static GridCommons MergeLogically(GridCommons A, GridCommons B) {
            GridCommons R = new GridCommons();
            R.m_GridGuid = Guid.NewGuid();


            //R.BcCellPartitioning; R.m_BcCellPartitioning
            //R.BcCells;
            //R.BcCellsStorageGuid;
            //R.CellPartitioning;
            //R.Cells;
            //R.Description;
            //R.EdgeRefElements;
            //R.EdgeTagNames;
            //R.GridGuid; == R.ID;
            //R.InversePeriodicTrafo on-demand
            //R.m_ClassNameOfEdgeRefElement
            //R.m_RefElements;
            //R.m_EdgeRefElements;

            /////////////////////////
            // merge ref elements
            /////////////////////////
            {
                R.m_RefElements = A.m_RefElements.CloneAs();
                int[] BrefElm = new int[B.m_RefElements.Length];
                for (int iKrefB = 0; iKrefB < BrefElm.Length; iKrefB++) {
                    RefElement Kref = B.m_RefElements[iKrefB];
                    int iKrefBNew = R.m_RefElements.IndexOf(Kref, (a, b) => a.Equals(b));
                    if (iKrefBNew < 0) {
                        ArrayTools.AddToArray(Kref, ref R.m_RefElements);
                        iKrefBNew = B.m_RefElements.Length - 1;
                    }
                    BrefElm[iKrefB] = iKrefBNew;
                }

                R.m_EdgeRefElements = A.m_EdgeRefElements.CloneAs();
                int[] BEdgeRefElm = new int[B.m_EdgeRefElements.Length];
                for (int iEdgeKrefB = 0; iEdgeKrefB < BEdgeRefElm.Length; iEdgeKrefB++) {
                    RefElement Kref = B.m_EdgeRefElements[iEdgeKrefB];
                    int iEdgeKrefBNew = R.m_EdgeRefElements.IndexOf(Kref, (a, b) => a.Equals(b));
                    if (iEdgeKrefBNew < 0) {
                        ArrayTools.AddToArray(Kref, ref R.m_EdgeRefElements);
                        iEdgeKrefBNew = B.m_EdgeRefElements.Length - 1;
                    }
                    BEdgeRefElm[iEdgeKrefB] = iEdgeKrefBNew;
                }
            }

            //////////////////////////////////////////
            // merge edge tags & periodic transformations 
            //////////////////////////////////////////

            byte[] BNewEdgeTag = new byte[0xff]; // mapping: old edge tag in B (index) => edge tag in R (content)
            {
                R.m_PeriodicTrafo.AddRange(A.m_PeriodicTrafo.Select(at => at.CloneAs()).ToArray());

                R.m_EdgeTagNames.AddRange(A.m_EdgeTagNames);
                foreach (var kv in B.m_EdgeTagNames) {
                    byte _EdgeTagB = kv.Key, _NewEdgeTagB;
                    string NameB = kv.Value;

                    if (R.m_EdgeTagNames.ContainsValue(NameB)) {
                        _NewEdgeTagB = R.m_EdgeTagNames.First(kv2 => kv2.Value.Equals(NameB)).Key;
                    } else {
                        _NewEdgeTagB = R.m_EdgeTagNames.Keys.Where(et => et < FIRST_PERIODIC_BC_TAG).Max();
                    }
                    BNewEdgeTag[_EdgeTagB] = _NewEdgeTagB;
                }


                AffineTrafo[] BperiodicTrafo = B.m_PeriodicTrafo.Select(at => at.CloneAs()).ToArray();

                for (int i = 0; i < BperiodicTrafo.Length; i++) {
                    var Tr = BperiodicTrafo[i];
                    int _BEdgeTag = i + FIRST_PERIODIC_BC_TAG;
                    int _BNewEdgeTag;
                    int idxTr = R.m_PeriodicTrafo.IndexOf(Tr, (TrA, TrB) => TrA.ApproximateEquals(TrB));
                    if (idxTr < 0) {
                        R.m_PeriodicTrafo.Add(Tr);
                        _BNewEdgeTag = R.m_PeriodicTrafo.Count + FIRST_PERIODIC_BC_TAG - 1;
                        R.m_EdgeTagNames.Add((byte)_BNewEdgeTag, B.m_EdgeTagNames[(byte)_BEdgeTag]);
                    } else {
                        _BNewEdgeTag = idxTr + FIRST_PERIODIC_BC_TAG;
                    }

                    BNewEdgeTag[_BEdgeTag] = (byte)_BNewEdgeTag;
                }

            }

            //////////////////////////
            // merge cells (logically)
            //////////////////////////
            {
                int NodeOffset = A.NodePartitioning.TotalLength;
                long JAglb = A.NumberOfCells_l;
                long JBglb = B.NumberOfCells_l;
                long KAglb = A.NoOfBcCells.MPISum();
                long KBglb = B.NoOfBcCells.MPISum();


                int JA = A.Cells.Length;
                int JB = B.Cells.Length;
                R.Cells = new Cell[JA + JB];
                int KA = A.NoOfBcCells;
                int KB = B.NoOfBcCells;
                R.BcCells = new BCElement[KA + KB];

                // check GlobalIds of cells
                for (int ja = 0; ja < JA; ja++) {
                    Cell CA = A.Cells[ja];
                    if (CA.GlobalID < 0 || CA.GlobalID >= JAglb)
                        throw new ArgumentException("Illegal GlobalId in grid A.", "A");
                }
                for (int jb = 0; jb < JB; jb++) {
                    Cell CB = B.Cells[jb];
                    if (CB.GlobalID < 0 || CB.GlobalID >= JBglb)
                        throw new ArgumentException("Illegal GlobalId in grid B.", "B");
                }
                // check globalIDs of boundary elements
                for (int ka = 0; ka < KA; ka++) {
                    BCElement BcA = A.BcCells[ka];
                    if (BcA.GlobalID < JAglb || BcA.GlobalID >= JAglb + KAglb)
                        throw new ArgumentException("Illegal GlobalId for boundary element in grid A.", "A");
                }
                for (int kb = 0; kb < KB; kb++) {
                    BCElement BcB = B.BcCells[kb];
                    if (BcB.GlobalID < JB || BcB.GlobalID >= JBglb + KBglb)
                        throw new ArgumentException("Illegal GlobalId for boundary element in grid B.", "B");
                }

                // cells of A:
                // -----------
                for (int ja = 0; ja < JA; ja++) {
                    Cell CA = A.Cells[ja];

                    R.Cells[ja] = new Cell() {
                        GlobalID = CA.GlobalID,
                        NodeIndices = CA.NodeIndices == null ? null : CA.NodeIndices.CloneAs(),
                        Type = CA.Type,
                        TransformationParams = CA.TransformationParams.CloneAs(),
                        CellFaceTags = CA.CellFaceTags == null ? null : CA.CellFaceTags.Select(tag => new CellFaceTag() {
                            EdgeTag = tag.EdgeTag,
                            FaceIndex = tag.FaceIndex,
                            ConformalNeighborship = tag.ConformalNeighborship,
                            NeighCell_GlobalID = tag.NeighCell_GlobalID < JAglb ?
                                                    tag.NeighCell_GlobalID : // normal cell
                                                    tag.NeighCell_GlobalID + JBglb, // boundary element 
                            PeriodicInverse = tag.PeriodicInverse
                        }).ToArray()
                    };

                    Debug.Assert(R.Cells[ja].GlobalID >= 0 && R.Cells[ja].GlobalID < JAglb);
                    if (CA.CellFaceTags != null) {
                        for (int i = 0; i < CA.CellFaceTags.Length; i++) {
                            var cft = R.Cells[ja].CellFaceTags[i];
                            var orgCft = CA.CellFaceTags[i];

                            Debug.Assert(
                                   (orgCft.NeighCell_GlobalID < JAglb && cft.NeighCell_GlobalID >= 0 && cft.NeighCell_GlobalID < JAglb) // normall cell
                                || (orgCft.NeighCell_GlobalID >= JAglb && cft.NeighCell_GlobalID >= JAglb + JBglb && cft.NeighCell_GlobalID < JAglb + JBglb + KAglb)); // boundary element
                        }
                    }

                }

                // cells of B:
                // -----------
                for (int jb = 0; jb < JB; jb++) {
                    Cell CB = B.Cells[jb];

                    R.Cells[jb + JA] = new Cell() {
                        GlobalID = CB.GlobalID + JAglb,
                        NodeIndices = CB.NodeIndices == null ? null : CB.NodeIndices.Select(idx => idx + NodeOffset).ToArray(),
                        Type = CB.Type,
                        TransformationParams = CB.TransformationParams.CloneAs(),
                        CellFaceTags = CB.CellFaceTags == null ? null : CB.CellFaceTags.Select(tag => new CellFaceTag() {
                            EdgeTag = BNewEdgeTag[tag.EdgeTag],
                            FaceIndex = tag.FaceIndex,
                            ConformalNeighborship = tag.ConformalNeighborship,
                            NeighCell_GlobalID = tag.NeighCell_GlobalID + JAglb, // correct for normal cell and boundary element
                            PeriodicInverse = tag.PeriodicInverse
                        }).ToArray()
                    };

                    Debug.Assert(R.Cells[jb + JA].GlobalID >= JAglb && R.Cells[jb + JA].GlobalID < JAglb + JBglb);
                    if (CB.CellFaceTags != null) {
                        for (int i = 0; i < CB.CellFaceTags.Length; i++) {
                            var cft = R.Cells[jb + JA].CellFaceTags[i];
                            var orgCft = CB.CellFaceTags[i];

                            Debug.Assert(
                                   (orgCft.NeighCell_GlobalID < JBglb && cft.NeighCell_GlobalID >= JAglb && cft.NeighCell_GlobalID < JAglb + JBglb) // normall cell
                                || (orgCft.NeighCell_GlobalID >= JBglb && cft.NeighCell_GlobalID >= JAglb + JBglb + KAglb && cft.NeighCell_GlobalID < JAglb + JBglb + KAglb + KBglb)); // boundary element
                        }
                    }

                }

                // boundary elements of A:
                // -----------------------

                for (int ka = 0; ka < KA; ka++) {
                    BCElement BcA = A.BcCells[ka];

                    R.BcCells[ka] = new BCElement() {
                        Conformal = BcA.Conformal,
                        GlobalID = BcA.GlobalID + JBglb,
                        NeighCell_GlobalIDs = BcA.NeighCell_GlobalIDs == null ? null : BcA.NeighCell_GlobalIDs.Select(
                            gidN => gidN < JAglb ? gidN : gidN + JBglb).ToArray(),
                        EdgeTag = BcA.EdgeTag,
                        NodeIndices = BcA.NodeIndices == null ? null : BcA.NodeIndices.CloneAs(),
                        TransformationParams = BcA.TransformationParams.CloneAs(),
                        Type = BcA.Type
                    };

                    Debug.Assert(R.BcCells[ka].GlobalID >= JAglb + JBglb && R.BcCells[ka].GlobalID < JAglb + JBglb + KAglb);
                    if (BcA.NeighCell_GlobalIDs != null) {
                        for (int i = 0; i < BcA.NeighCell_GlobalIDs.Length; i++) {
                            long Ngid = R.BcCells[ka].NeighCell_GlobalIDs[i];
                            long orgNgid = BcA.NeighCell_GlobalIDs[i];

                            Debug.Assert((orgNgid < JAglb && Ngid >= 0 && Ngid < JAglb) // normall cell
                                || (orgNgid >= JAglb && Ngid >= JAglb + JBglb && Ngid < JAglb + JBglb + KAglb)); // boundary element
                        }
                    }
                }

                // boundary elements of B:
                // -----------------------

                for (int kb = 0; kb < KB; kb++) {
                    BCElement BcB = A.BcCells[kb];

                    R.BcCells[kb + KA] = new BCElement() {
                        Conformal = BcB.Conformal,
                        GlobalID = BcB.GlobalID + JAglb,
                        NeighCell_GlobalIDs = BcB.NeighCell_GlobalIDs == null ? null : BcB.NeighCell_GlobalIDs.Select(
                            gidN => gidN < JBglb ? gidN + JAglb : gidN + JAglb).ToArray(),
                        EdgeTag = BNewEdgeTag[BcB.EdgeTag],
                        NodeIndices = BcB.NodeIndices == null ? null : BcB.NodeIndices.Select(idx => idx + NodeOffset).ToArray(),
                        TransformationParams = BcB.TransformationParams.CloneAs(),
                        Type = BcB.Type
                    };

                    Debug.Assert(R.BcCells[kb].GlobalID >= JAglb + JBglb + KAglb && R.BcCells[kb].GlobalID < JAglb + JBglb + KAglb + KBglb);
                    if (BcB.NeighCell_GlobalIDs != null) {
                        for (int i = 0; i < BcB.NeighCell_GlobalIDs.Length; i++) {
                            long Ngid = R.BcCells[kb + KA].NeighCell_GlobalIDs[i];
                            long orgNgid = BcB.NeighCell_GlobalIDs[i];

                            Debug.Assert((orgNgid < JBglb && Ngid >= JAglb && Ngid < JAglb + JBglb) // normall cell
                                || (orgNgid >= JBglb && Ngid >= JAglb + JBglb + KAglb && Ngid < JAglb + JBglb + KAglb + KBglb)); // boundary element
                        }
                    }
                }
            }

            return R;
        }


        /// <summary>
        /// Usually used after <see cref="MergeLogically(GridCommons, GridCommons)"/>; this method finds element boundaries
        /// which intersect geometrically, but not logically and inserts a <see cref="CellFaceTag"/> which connects those cells.
        /// </summary>
        /// <param name="g"></param>
        /// <param name="upsampling"></param>
        /// <returns></returns>
        public static GridCommons Seal(GridCommons g, int upsampling = 4) {
            GridCommons R = g.CloneAs();
            g = null;
            GridData gdat = new GridData(R);
            int D = gdat.SpatialDimension;
            int J = gdat.Cells.NoOfLocalUpdatedCells;

            if (R.CellPartitioning.MpiSize > 1)
                throw new NotSupportedException("Not supported in MPI-parallel mode.");

            //NodeSet[] TestNodes = gdat.Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.GetSubdivisionTree(upsampling).GlobalVertice).ToArray();
            NodeSet[] TestNodes = gdat.Edges.EdgeRefElements.Select(KrefEdge => KrefEdge.GetBruteForceQuadRule(upsampling, 1).Nodes).ToArray();
            // Its better to use vertices in the interior of the element; if we use vertices at the corners, we might get
            // intersection of edges that just share one point.


            // Define all edges that will be tested (set to boundary edges)
            // ============================================================
            int[] UnknownEdges = gdat.BoundaryEdges.ItemEnum.ToArray();
            int L = UnknownEdges.Sum(iEdg => TestNodes[gdat.Edges.GetRefElementIndex(iEdg)].NoOfNodes);
            
            // Transform nodes on edges (that should be tested) to global coordinates
            // ======================================================================

            MultidimensionalArray TestNodesGlobal = MultidimensionalArray.Create(L, D);
            MultidimensionalArray NormalsGlobal = MultidimensionalArray.Create(L, D);
            int[] NodeToEdge = new int[L]; // pointer l -> Edge index, where l is the first index into 'TestNodesGlobal' & 'NormalsGlobal'


            int[,] E2C = gdat.Edges.CellIndices;
            int cnt = 0;
            foreach (int iEdg in UnknownEdges) {
                int iKref = gdat.Edges.GetRefElementIndex(iEdg);
                NodeSet Ns = TestNodes[iKref];
                int K = Ns.NoOfNodes;

                int[] I0 = new int[] { cnt, 0 };
                int[] IE = new int[] { cnt + K - 1, D - 1 };

                MultidimensionalArray TN = gdat.GlobalNodes.GetValue_EdgeSV(Ns, iEdg, 1);
                TestNodesGlobal.SetSubArray(TN.ExtractSubArrayShallow(0, -1, -1), I0, IE);

                MultidimensionalArray N1 = gdat.Edges.NormalsCache.GetNormals_Edge(Ns, iEdg, 1);
                NormalsGlobal.SetSubArray(N1.ExtractSubArrayShallow(0, -1, -1), I0, IE);

                for (int i = cnt; i < cnt + K; i++)
                    NodeToEdge[i] = iEdg;

                cnt += K;
            }

            // binary tree to speed up point localization
            int[] pl_Permutation = new int[L];
            PointLocalization pl = new PointLocalization(TestNodesGlobal, 0.01, pl_Permutation);
            Debug.Assert(!object.ReferenceEquals(pl.Points, TestNodesGlobal));

            // compare search edges to all other nodes
            // =======================================

            // mapping: cell --> Neighbour cell index, face index
            // 1st index: Cell index;
            // 2nd index: enumeration
            List<Tuple<int, int>>[] FoundPairings = new List<Tuple<int, int>>[gdat.Cells.NoOfCells];

            int[][] C2E = gdat.Cells.Cells2Edges;
            byte[,] E2F = gdat.Edges.FaceIndices;

            int cnt2 = 0;
            for (int iEdgC = 0; iEdgC < UnknownEdges.Length; iEdgC++) { // loop over edges that may get sealed
                int iEdg = UnknownEdges[iEdgC];
                int iKref = gdat.Edges.GetRefElementIndex(iEdg);
                NodeSet Ns = TestNodes[iKref];
                int K = Ns.NoOfNodes;

                int jCell1 = E2C[iEdg, 0];
                Debug.Assert(E2C[iEdg, 1] < 0);
                int iFace1 = E2F[iEdg, 0];
                Debug.Assert(E2F[iEdg, 1] == byte.MaxValue);


                int[] I0 = new int[] { cnt2, 0 };
                int[] IE = new int[] { cnt2 + K - 1, D - 1 };
                MultidimensionalArray TN = TestNodesGlobal.ExtractSubArrayShallow(I0, IE);
                //MultidimensionalArray N1 = NormalsGlobal.ExtractSubArrayShallow(I0, IE);

                // find bounding box for edge
                BoundingBox bbEdge = new BoundingBox(TN);
                if (bbEdge.h_min / bbEdge.h_max < 1.0e-5) {
                    // very thin bounding box, thicken slightly
                    double delta = bbEdge.h_max * 1.0e-5;
                    for (int d = 0; d < D; d++) {
                        bbEdge.Min[d] -= delta;
                        bbEdge.Max[d] += delta;
                    }

                }
                bbEdge.ExtendByFactor(0.01);

                // determine binary code for bounding box
                int bbEdgeSigBits;
                GeomBinTreeBranchCode bbEdgeCode = GeomBinTreeBranchCode.Combine(
                    GeomBinTreeBranchCode.CreateFormPoint(pl.PointsBB, bbEdge.Min),
                    GeomBinTreeBranchCode.CreateFormPoint(pl.PointsBB, bbEdge.Max),
                    out bbEdgeSigBits);

                // determine all points in bounding box
                int iP0, Len;
                pl.GetPointsInBranch(bbEdgeCode, bbEdgeSigBits, out iP0, out Len);

                // determine all edged which potentially overlap with edge 'iEdg'
                HashSet<int> PotOvrlap = new HashSet<int>(); // a set of edge indices
                for (int n = 0; n < Len; n++) {
                    int l = iP0 + n;
                    int iPt = pl_Permutation[l];
                    Debug.Assert(GenericBlas.L2DistPow2(pl.Points.GetRow(l), TestNodesGlobal.GetRow(iPt)) <= 0);

                    int iOvlpEdge = NodeToEdge[iPt];
                    if (iOvlpEdge != iEdg)
                        PotOvrlap.Add(iOvlpEdge);
                }
                //int[] PotOvrlap = UnknownEdges.CloneAs();


                // determine actually overlapping boundary edges:
                foreach (int iOvrlapEdge in PotOvrlap) {
                    int jCell2 = E2C[iOvrlapEdge, 0];
                    Debug.Assert(E2C[iOvrlapEdge, 1] < 0);
                    if (jCell2 == jCell1)
                        continue;
                    int iFace2 = E2F[iOvrlapEdge, 0];
                    Debug.Assert(E2F[iOvrlapEdge, 1] == byte.MaxValue);


                    int AllreadyFound = FoundPairings[jCell1] == null ?
                        0 : FoundPairings[jCell1].Where(tp => tp.Item1 == jCell2 && tp.Item2 == iFace1).Count();
                    if (AllreadyFound > 1)
                        throw new ApplicationException("Error in algorithmus.");
                    if (AllreadyFound > 0)
                        continue;

                    var Kref_j2 = gdat.Cells.GetRefElement(jCell2);
                    double h = Kref_j2.GetMaxDiameter();

                    MultidimensionalArray LocVtx_j2 = MultidimensionalArray.Create(K, D);
                    bool[] NewtonConvervence = new bool[K];
                    gdat.TransformGlobal2Local(TN, LocVtx_j2, jCell2, NewtonConvervence);

                    for (int k = 0; k < K; k++) { // loop over all transformed points
                        if (!NewtonConvervence[k])
                            continue;

                        double[] pt = LocVtx_j2.GetRow(k);
                        double[] ptClose = new double[D];
                        double dist = Kref_j2.ClosestPoint(pt, ptClose);

                        if (dist > h * 1.0e-8)
                            continue;

                        AffineManifold Face = Kref_j2.GetFacePlane(iFace2);
                        double FaceDist = Face.PointDistance(pt);

                        if (FaceDist.Abs() > 1.0e-8 * h)
                            continue;


                        NodeSet Ns2 = new NodeSet(Kref_j2, pt);
                        MultidimensionalArray Normals2 = MultidimensionalArray.Create(1, D);
                        gdat.Edges.GetNormalsForCell(Ns2, jCell2, iFace2, Normals2);

                        double[] N1d = NormalsGlobal.GetRow(cnt2 + k);
                        double[] N2d = Normals2.GetRow(0);


                        //Check if face normals points exactly in the opposite direction, 2 ways:
                        // 1) calculate angle between both normals -> bad choice because of Math.Acos 
                        // 2) inner product of two opposite vectors is -1 
                        //if (Math.Abs(Math.Abs(Math.Acos(GenericBlas.InnerProd(N1d, N2d))) - Math.PI) > 1.0e-8)
                        if (Math.Abs(GenericBlas.InnerProd(N1d, N2d) + 1.0) > 1.0e-8)
                            continue;



                        // if we reach this point, jCell1 should match with jCell2/iFace2
                        if (FoundPairings[jCell1] == null)
                            FoundPairings[jCell1] = new List<Tuple<int, int>>();
                        if (FoundPairings[jCell2] == null)
                            FoundPairings[jCell2] = new List<Tuple<int, int>>();
                        FoundPairings[jCell1].Add(new Tuple<int, int>(jCell2, iFace1));
                        FoundPairings[jCell2].Add(new Tuple<int, int>(jCell1, iFace2));
                        break; // no need to test jCell1 vs. jCell2 anymore
                    }
                }

                cnt2 += K;
            }

            // add the newly found pairings to the grid
            for (int j = 0; j < J; j++) {
                var fp = FoundPairings[j];
                if (fp != null) {
                    foreach (var t in fp) {
                        ArrayTools.AddToArray(new CellFaceTag() {
                            EdgeTag = 0,
                            FaceIndex = t.Item2,
                            ConformalNeighborship = false,
                            NeighCell_GlobalID = gdat.CurrentGlobalIdPermutation.Values[t.Item1]
                        }, ref R.Cells[j].CellFaceTags);
                    }
                }
            }

            return R;
        }



    }
}
