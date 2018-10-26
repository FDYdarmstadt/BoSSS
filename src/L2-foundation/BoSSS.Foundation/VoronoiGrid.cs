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
using ilPSP.Connectors.Matlab;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform.LinAlg;
using System.Diagnostics;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Caching;

namespace BoSSS.Foundation.Grid.Voronoi {
     public partial class VoronoiGrid : IGridData {

        public MultidimensionalArray DelaunayVertices;
        
        /// <summary>
        /// MPI process rank (within world communicator)
        /// </summary>
        public int MpiRank {
            get {
                return this.CellPartitioning.MpiRank;
            }
        }

        /// <summary>
        /// MPI world communicator size 
        /// </summary>
        public int MpiSize {
            get {
                return this.CellPartitioning.MpiSize;
            }
        }

        public int SpatialDimension {
            get {
                return DelaunayVertices.NoOfCols;
            }
        }

        public BasisData ChefBasis {
            get {
                throw new NotImplementedException();
            }
        }

        public IParallelization iParallel {
            get {
                throw new NotImplementedException();
            }
        }

        public Partitioning CellPartitioning {
            get {
                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// The global ID for each cell
        /// </summary>
        public Comm.Permutation CurrentGlobalIdPermutation {
            get {
                throw new NotImplementedException();
            }
        }

        public CacheLogicImplBy_CNs InverseJacobian {
            get {
                throw new NotImplementedException();
            }
        }

        public CacheLogicImplBy_CNs JacobianDeterminat {
            get {
                throw new NotImplementedException();
            }
        }

        public CacheLogic_CNs GlobalNodes {
            get {
                throw new NotImplementedException();
            }
        }

        public CacheLogicImplBy_CNs Jacobian {
            get {
                throw new NotImplementedException();
            }
        }

        public CacheLogicImplBy_CNs AdjungateJacobian {
            get {
                throw new NotImplementedException();
            }
        }

        public Guid GridID {
            get {
                throw new NotImplementedException();
            }
        }

        public IDictionary<byte, string> EdgeTagNames {
            get {
                throw new NotImplementedException();
            }
        }

        public void CreateWithMatlab() {
            // ================================
            // generate voronoi graph in matlab
            // ================================


            int J = this.DelaunayVertices.NoOfRows;
            int D = this.DelaunayVertices.NoOfCols;

            if (D != 2)
                throw new NotSupportedException("todo");


            int[][] OutputVertexIndex = new int[J*5][];
            MultidimensionalArray VertexCoordinates;
            {
                var Matlab = new BatchmodeConnector();

                Matlab.PutMatrix(this.DelaunayVertices, "X");

                // create mirror points
                Matlab.Cmd("[J, D] = size(X);");
                Matlab.Cmd("Xneg = [-X(:, 1), X(:, 2)];");
                Matlab.Cmd("Yneg = [X(:, 1), -X(:, 2)];");
                Matlab.Cmd("X2 = [ones(J, 1) * 2, zeros(J, 1)];");
                Matlab.Cmd("Y2 = [zeros(J, 1), ones(J, 1) * 2];");

                Matlab.Cmd("Xm = X;");
                Matlab.Cmd("Xm = [Xm; Xneg];    % mirror at x = 0");
                Matlab.Cmd("Xm = [Xm; X2 + Xneg]; % mirror at x = 1");
                Matlab.Cmd("Xm = [Xm; Yneg];    % mirror at x = 0");
                Matlab.Cmd("Xm = [Xm; Y2 + Yneg]; % mirror at x = 1");

                // compute Voronoi diagramm
                Matlab.Cmd("[V, C] = voronoin(Xm);");

                // output (export from matlab)
                Matlab.GetStaggeredIntArray(OutputVertexIndex, "C");
                Matlab.GetMatrix(null, "V");

                // run matlab
                Matlab.Execute(false);

                // import here
                VertexCoordinates = (MultidimensionalArray)(Matlab.OutputObjects["V"]);

                // correct indices (1-based index to 0-based index)
                foreach(int[] cell in OutputVertexIndex) {
                    int K = cell.Length;
                    for (int k = 0; k < K; k++) {
                        cell[k]--;
                    }
                }
            }

            // ===============
            // record internal 
            // ===============

            {
                // define Cell data
                {
                    this.m_CellData = new CellData();
                    this.m_CellData.m_Owner = this;
                    this.m_VertexData = new VertexData();
                    this.m_LogEdges = new LogEdgeData();

                    this.m_VertexData.Coordinates = VertexCoordinates;
                    this.m_CellData.CellVertices = OutputVertexIndex.GetSubVector(0, J);

                    this.m_CellData.InfoFlags = new CellInfo[J];
                    ArrayTools.SetAll(this.m_CellData.InfoFlags, CellInfo.CellIsAffineLinear | CellInfo.IsAggregate);
                }

                // decomposition of Voronoi cells to triangles/tetrahedrons
                {
                    m_CellData.AggregateCellToParts = new int[J][];
                    if (D == 2) {
                        var Tri = RefElements.Triangle.Instance;
                        m_CellData.RefElements = new RefElement[] { Tri };

                        int cnt = 0;
                        for (int j = 0; j < J; j++) {
                            int[] VtxIndices = m_CellData.CellVertices[j];

                            int[] PartIdx = new int[VtxIndices.Length - 2];
                            for (int i = 0; i < PartIdx.Length; i++) {
                                PartIdx[i] = cnt;
                                cnt++;
                            }
                            m_CellData.AggregateCellToParts[j] = PartIdx;
                        }

                        int NoOfParts = cnt;
                        m_CellData.PartTransformation = MultidimensionalArray.Create(NoOfParts, D, D);
                        m_CellData.PartCenter = MultidimensionalArray.Create(NoOfParts, D);

                        MultidimensionalArray TriangleVtx = MultidimensionalArray.Create(3, D);
                        for (int j = 0; j < J; j++) {
                            int[] VtxIndices = m_CellData.CellVertices[j];
                            int[] PartIdx = m_CellData.AggregateCellToParts[j];

                            for (int i = 0; i < PartIdx.Length; i++) {
                                int iV0 = VtxIndices[0];
                                int iV1 = VtxIndices[i + 1];
                                int iV2 = VtxIndices[i + 2];

                                TriangleVtx[0, 0] = m_VertexData.Coordinates[iV0, 0];
                                TriangleVtx[0, 1] = m_VertexData.Coordinates[iV0, 1];
                                TriangleVtx[1, 0] = m_VertexData.Coordinates[iV1, 0];
                                TriangleVtx[1, 1] = m_VertexData.Coordinates[iV1, 1];
                                TriangleVtx[2, 0] = m_VertexData.Coordinates[iV2, 0];
                                TriangleVtx[2, 1] = m_VertexData.Coordinates[iV2, 1];

                                var TR = AffineTrafo.FromPoints(Tri.Vertices, TriangleVtx);

                                m_CellData.PartTransformation.ExtractSubArrayShallow(PartIdx[i], -1, -1).Set(TR.Matrix);
                                m_CellData.PartCenter.ExtractSubArrayShallow(PartIdx[i], -1).SetVector(TR.Affine);
                            }
                        }
                    } else if (D == 3) {
                        throw new NotImplementedException("todo");
                    } else {
                        throw new NotSupportedException("Unknown spatial dimension.");
                    }
                }

                // bounding boxes, transformations
                {
                    BoundingBox BB = new BoundingBox(D);
                    m_CellData.BoundingBoxTransformation = MultidimensionalArray.Create(J, D, D);
                    m_CellData.BoundingBoxCenter = MultidimensionalArray.Create(J, D);
                    for (int j = 0; j < J; j++) {
                        m_CellData.GetCellBoundingBox(j, BB);
                        for (int d = 0; d < D; d++) {
                            double lo = BB.Min[d];
                            double hi = BB.Max[d];
                            m_CellData.BoundingBoxCenter[j, d] = 0.5 * (lo + hi);
                            m_CellData.BoundingBoxTransformation[j, d, d] = 0.5 * (hi - lo);
                        }

                    }
                }

                // mapping: vertex to cell
                {
                    List<int>[] VertexToCell = new List<int>[VertexCoordinates.Length];
                    for (int j = 0; j < J; j++) {
                        foreach (int iVtx in OutputVertexIndex[j]) {
                            if (VertexToCell[iVtx] == null)
                                VertexToCell[iVtx] = new List<int>();

                            if (!VertexToCell[iVtx].Contains(j)) {
                                VertexToCell[iVtx].Add(j);
                            }
                        }
                    }

                    m_VertexData.VerticeToCell = new int[VertexToCell.Length][];
                    for (int i = 0; i < VertexToCell.Length; i++) {
                        if (VertexToCell[i] == null)
                            m_VertexData.VerticeToCell[i] = new int[0];
                        else
                            m_VertexData.VerticeToCell[i] = VertexToCell[i].ToArray();
                    }
                    VertexToCell = null;
                }

                // cell neighbors, edges
                {
                    m_CellData.CellNeighbours = new int[J][];
                    var tmpCells2Edges = new List<int>[J];

                    Dictionary<int, int> ShareCount = new Dictionary<int, int>(); // key: cell index; value: number of vertices shared with this cell
                    var EdgesTemp = new List<EdgeTemp>();
                    List<int> Neighs = new List<int>();
                    List<int> EdgeVtx = new List<int>();
                    List<int> IdedEdgsAtOneCell = new List<int>();
                    for (int jCell = 0; jCell < J; jCell++) { // loop over cells
                        ShareCount.Clear();
                        Neighs.Clear();
                        IdedEdgsAtOneCell.Clear();

                        // determine how many vertices 'jCell' shares with other cells 
                        foreach (int iVtx in m_CellData.CellVertices[jCell]) {
                            foreach (int jOtherCell in m_VertexData.VerticeToCell[iVtx]) {
                                if (jOtherCell != jCell) {
                                    if (!ShareCount.ContainsKey(jOtherCell))
                                        ShareCount.Add(jOtherCell, 1);
                                    else
                                        ShareCount[jOtherCell]++;
                                }
                            }
                        }

                        // find faces
                        int[][] FaceIdx = ConvexHullFaces(m_VertexData.Coordinates, m_CellData.CellVertices[jCell]);

                        // determine cell neighbors and edges 
                        int NoOfFacesFound = 0;
                        foreach (var kv in ShareCount) {
                            int jCellNeigh = kv.Key;
                            int NoOfSharedVtx = kv.Value;
                            if (NoOfSharedVtx >= D) {
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                // cells 'jCell' and 'jCellNeigh' share more than 'D' vertices - this is an edge to another cell,
                                // resp. a face of 'jCell'.
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                Debug.Assert(jCellNeigh != jCell);
                                Debug.Assert(!Neighs.Contains(jCellNeigh));
                                Neighs.Add(jCellNeigh);
                                NoOfFacesFound++;

                                EdgeVtx.Clear();
                                foreach (int iVtx in m_CellData.CellVertices[jCell]) {
                                    if (Array.IndexOf(m_VertexData.VerticeToCell[iVtx], jCellNeigh) >= 0)
                                        EdgeVtx.Add(iVtx);
                                }


                                if (jCell < jCellNeigh) {
                                    // the pairing 'jCell'/'jCellNeigh' will be discovered twice;
                                    // we only want to record once

                                    var Etmp = new EdgeTemp() {
                                        jCell1 = jCell,
                                        jCell2 = jCellNeigh,
                                        Vertices = EdgeVtx.ToArray()
                                    };
                                    EdgesTemp.Add(Etmp);
                                    IdedEdgsAtOneCell.Add(EdgesTemp.Count - 1);
                                    if (tmpCells2Edges[jCell] == null) {
                                        tmpCells2Edges[jCell] = new List<int>();
                                    }
                                    if (tmpCells2Edges[jCellNeigh] == null) {
                                        tmpCells2Edges[jCellNeigh] = new List<int>();
                                    }
                                    tmpCells2Edges[jCell].Add(EdgesTemp.Count); // the funky convention for edges-to-cell: the index is 
                                    tmpCells2Edges[jCellNeigh].Add(-EdgesTemp.Count); // shifted by 1, out-cell is negative
                                } else {
                                    Debug.Assert(jCellNeigh < jCell);

                                    int MatchCount = 0;
                                    foreach(int i in tmpCells2Edges[jCellNeigh]) {
                                        int iEdge = Math.Abs(i) - 1;

                                        if (EdgesTemp[iEdge].jCell1 == jCellNeigh && EdgesTemp[iEdge].jCell2 == jCell) {
                                            MatchCount++;
                                            IdedEdgsAtOneCell.Add(iEdge);
                                        }
                                    }
                                    Debug.Assert(MatchCount == 1);
                                }
#if DEBUG
                                if (D == 2) {
                                    Debug.Assert(EdgeVtx.Count == 2);
                                } else if (D == 3) {
                                    // verify that all vertices of the edge are geometrically in one plane

                                    Debug.Assert(EdgeVtx.Count >= 3);
                                    var FacePlane = AffineManifold.FromPoints(
                                        m_VertexData.Coordinates.GetRow(EdgeVtx[0]),
                                        m_VertexData.Coordinates.GetRow(EdgeVtx[1]),
                                        m_VertexData.Coordinates.GetRow(EdgeVtx[2])
                                        );

                                    BoundingBox BB = new BoundingBox(D);
                                    m_CellData.GetCellBoundingBox(jCell, BB);
                                    double h = BB.Diameter;

                                    foreach (int iVtx in EdgeVtx) {
                                        double dist = Math.Abs(FacePlane.PointDistance(m_VertexData.Coordinates.GetRow(iVtx)));
                                        Debug.Assert(dist < h * 1e-8);
                                    }

                                } else {
                                    throw new NotSupportedException("Unknown spatial dimension.");
                                }
#endif

                            }
                        }
                        m_CellData.CellNeighbours[jCell] = Neighs.ToArray();
                        Debug.Assert(NoOfFacesFound <= FaceIdx.Length);
                        Debug.Assert(NoOfFacesFound == IdedEdgsAtOneCell.Count);

                        // boundary edges
                        if(NoOfFacesFound == FaceIdx.Length) {
                            // nothing to do - all faces/edges identified
#if DEBUG
                            for(int i = 0; i < NoOfFacesFound; i++) {
                                int Matches = 0;
                                for(int l = 0; l < NoOfFacesFound; l++) {
                                    if (FaceIdx[i].SetEquals(EdgesTemp[IdedEdgsAtOneCell[l]].Vertices))
                                        Matches++;
                                }
                                Debug.Assert(Matches == 1);
                            }
#endif
                        } else {
                            // missing boundary 

                            for (int i = 0; i < FaceIdx.Length; i++) {
                                int Matches = 0;
                                for (int l = 0; l < NoOfFacesFound; l++) {
                                    if (FaceIdx[i].SetEquals(EdgesTemp[IdedEdgsAtOneCell[l]].Vertices)) {
                                        Matches++;
                                    }
                                }
                                Debug.Assert(Matches <= 1);

                                if(Matches == 0) {
                                    // boundary edge found
                                    var Etmp = new EdgeTemp() {
                                        jCell1 = jCell,
                                        jCell2 = int.MinValue,
                                        Vertices = EdgeVtx.ToArray()
                                    };
                                    EdgesTemp.Add(Etmp);
                                    tmpCells2Edges[jCell].Add(EdgesTemp.Count); // index shifted by 1
                                }
                            }
                        }
                    }

                    // convert temporary data structures to the final ones
                    m_CellData.Cells2Edges = new int[J][];
                    var C2E = m_CellData.Cells2Edges;
                    for(int j = 0; j < J; j++) {
                        C2E[j] = tmpCells2Edges[j] != null ? tmpCells2Edges[j].ToArray() : new int[0];
                    }


                    m_GeomEdges = new GeomEdgeData();
                    int NoOfEdges = EdgesTemp.Count;
                    m_GeomEdges.Info = new EdgeInfo[NoOfEdges];
                    m_LogEdges.CellIndices = new int[NoOfEdges, 2];
                    m_GeomEdges.VertexIndices = new int[NoOfEdges][];
                    var Evtx = m_GeomEdges.VertexIndices;
                    var E2C = m_LogEdges.CellIndices;
                    var Einf = m_GeomEdges.Info;
                    for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                        var Etmp = EdgesTemp[iEdge];
                        E2C[iEdge, 0] = Etmp.jCell1;
                        E2C[iEdge, 1] = Etmp.jCell2;
                        Einf[iEdge] = EdgeInfo.EdgeIsAffineLinear | EdgeInfo.IsAggregate;
                        if (Etmp.jCell2 < 0)
                            Einf[iEdge] |= EdgeInfo.Boundary;
                        Evtx[iEdge] = Etmp.Vertices;
                    }
                }


                // edge metrics
                {
                    if(D == 2) {
                        m_GeomEdges.EdgeRefElements = new RefElement[] { Line.Instance };
                    } else if(D == 3) {
                        m_GeomEdges.EdgeRefElements = new RefElement[] { Triangle.Instance };
                    } else {
                        throw new NotSupportedException("Unknown spatial dimension.");
                    }

                    int[][] Evtx = m_GeomEdges.VertexIndices;

                    int NoOfParts = 0;
                    int NoOfEdges = m_GeomEdges.Count;
                    for(int iEdge = 0; iEdge < NoOfEdges; iEdge++) {
                        NoOfParts += Evtx[iEdge].Length - D + 1;
                    }
                    Debug.Assert(D != 2 || NoOfParts == NoOfEdges);

                    m_GeomEdges.Edge2CellTrafoIndex = new int[NoOfParts, 2];
                    var tmpEdg2CellTrafo = new Dictionary<int, Tuple<int, AffineTrafo>>();

                    // unsolved problem:
                    // (D-1) -- dimensional tesselation of edge is given by D -- dimensional tesselation of adjacent cells
                    // * case D == 2: trivial
                    // * case D == 3: the problem is that two adjacent cells will induce _two different_ tesselations, which 
                    //                wont match in the general case. 


                    MultidimensionalArray PreImage = m_GeomEdges.EdgeRefElements[0].Vertices;
                    MultidimensionalArray Image = MultidimensionalArray.Create(PreImage.NoOfRows, D);

                    //for(int iEdge )



                }

            }

        }

        class EdgeTemp {
            public int[] Vertices;
            public int jCell1;
            public int jCell2;
        }

        static int[][] ConvexHullFaces(MultidimensionalArray C, int[] PointIdx) {
            int D = C.NoOfCols;
            if(D == 2) {
                int L = PointIdx.Length;
                int[][] R = new int[L][];
                for(int l = 0; l < L; l++) {
                    R[l] = new int[] { PointIdx[l], PointIdx[(l + 1) % L] };
                }
#if DEBUG  
                double orientSoll = 0;
                for (int l = 0; l < L; l++) {
                    double[] O1 = C.GetRow(R[l][0]);
                    double[] V1 = C.GetRow(R[l][1]);
                    double[] O2 = V1.CloneAs(); Debug.Assert(R[l][1] == R[(l + 1) % L][0]);
                    double[] V2 = C.GetRow(R[(l + 1) % L][1]);

                    Vector _V1 = new Vector(V1[0] - O1[0], V1[1] - O1[1], 0.0);
                    Vector _V2 = new Vector(V2[0] - O2[0], V2[1] - O2[1], 0.0);
                    Vector Norm = _V1.CrossProduct(_V2);

                    double orient = Math.Sign(Norm.z);
                    if(l == 0) {
                        orientSoll = orient;
                    } else {
                        //if (orient != orientSoll)
                        //    GnuplotCallback(C, PointIdx);

                        Debug.Assert(orient == orientSoll);
                    }
                }
#endif
                return R;
            } else if(D == 3) {
                throw new NotImplementedException("todo");
            } else {
                throw new NotSupportedException("unknown spatial dimension");
            }
        }

        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, MultidimensionalArray GlobalVerticesOut, int jCell) {
            throw new NotImplementedException();
        }

        public void TransformLocal2Global(MultidimensionalArray LocalVerticesIn, int j0, int Len, MultidimensionalArray GlobalVerticesOut, int OutArrayOffset) {
            throw new NotImplementedException();
        }

        public void TransformLocal2Global(MultidimensionalArray NS, int j0, int Len, MultidimensionalArray Nodesglob) {
            throw new NotImplementedException();
        }

        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int j0, int Len, int OutArrayOffset) {
            throw new NotImplementedException();
        }

        public void TransformGlobal2Local(MultidimensionalArray GlobalVerticesIn, MultidimensionalArray LocalVerticesOut, int jCell, bool[] NewtonConvergence) {
            throw new NotImplementedException();
        }


        //public static Action<MultidimensionalArray, int[]> GnuplotCallback;


    }
}
