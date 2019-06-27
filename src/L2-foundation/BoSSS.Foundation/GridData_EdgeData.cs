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
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.Utils.Geom;

namespace BoSSS.Foundation.Grid.Classic {

    partial class GridData {

        /// <summary>
        /// See <see cref="Edges"/>
        /// </summary>
        private EdgeData m_Edges;

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        public EdgeData Edges {
            get {
                return m_Edges;
            }
        }

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        public IGeometricalEdgeData iGeomEdges {
            get {
                return m_Edges;
            }
        }

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        public ILogicalEdgeData iLogicalEdges {
            get {
                return m_Edges;
            }
        }

        /// <summary>
        /// metrics and operations which are associated to edges
        /// </summary>
        public class EdgeData : IGeometricalEdgeData, ILogicalEdgeData {

            /// <summary>
            /// ctor
            /// </summary>
            internal EdgeData(GridData _owner) {
                m_owner = _owner;
                if (m_owner.SpatialDimension <= 2)
                    if (this.EdgeRefElements.Length != 1)
                        throw new ApplicationException("geometry error.");
            }

            /// <summary>
            /// See <see cref="GetEdges4RefElement(RefElement)"/>
            /// </summary>
            private EdgeMask[] m_Edges4RefElement;
            
            ///// <summary>
            ///// For each (edge) reference element, this method provides a
            ///// mask containing all cells which are mapped from the specific
            ///// reference element.
            ///// </summary>
            ///// <param name="iKrefIndex">
            ///// reference element index: <see cref="EdgeRefElements"/>;
            ///// </param>
            //public EdgeMask GetEdges4RefElement(int iKrefIndex) {
            //    return this.GetEdges4RefElement(this.EdgeRefElements[iKrefIndex]);
            //}

            /// <summary>
            /// For each (edge) reference element, this method provides a
            /// mask containing all cells which are mapped from the specific
            /// reference element.
            /// </summary>
            /// <param name="Kref">
            /// Reference element for edges.
            /// </param>
            public EdgeMask GetEdges4RefElement(RefElement Kref) {
                var RefElm = this.EdgeRefElements;
                int iKrefIndex = RefElm.IndexOf(Kref, (A, B) => object.ReferenceEquals(A, B));
                if(iKrefIndex < 0 || iKrefIndex >= RefElm.Length)
                    throw new ArgumentException();

                if(m_Edges4RefElement == null)
                    m_Edges4RefElement = new EdgeMask[RefElm.Length];

                if(m_Edges4RefElement[iKrefIndex] == null) {
                    int E = this.Count;
                    var _mask = new System.Collections.BitArray(E);

                    for(int e = 0; e < E; e++) {
                        int _iKref_edge = this.GetRefElementIndex(e);

                        if(_iKref_edge == iKrefIndex)
                            _mask[e] = true;
                    }

                    m_Edges4RefElement[iKrefIndex] = new EdgeMask(m_owner, _mask, MaskType.Geometrical);
                }

                return m_Edges4RefElement[iKrefIndex];
            }


            /// <summary>
            /// Reference elements for edges.
            /// </summary>
            public RefElement[] EdgeRefElements {
                get {
#if DEBUG
                    // test that each the face reference element, 
                    //     for each cell/volume reference element, 
                    //         is contained in the list of edge reference elements, 
                    // by reference-equality!
                    if(this.m_owner.Cells != null) {
                        IEnumerable<RefElement> Edge_KrefS = m_owner.Grid.EdgeRefElements;
                        foreach(var Cell_Kref in this.m_owner.Cells.RefElements) {
                            RefElement CellFaceKref = Cell_Kref.FaceRefElement;

                            int FoundCount = 0;
                            for(int i = 0; i < Edge_KrefS.Count(); i++) {
                                if(object.ReferenceEquals(CellFaceKref, Edge_KrefS.ElementAt(i)))
                                    FoundCount++;
                            }
                            Debug.Assert(FoundCount == 1);

                            //Debug.Assert(Edge_KrefS.Where(A => object.ReferenceEquals(Cell_Kref.FaceRefElement, A)).Count() == 1);
                        }
                    }
#endif
                    return m_owner.Grid.EdgeRefElements;
                }
            }

            /// <summary>
            /// pointer to owner object
            /// </summary>
            private GridData m_owner;

            /// <summary>
            /// the minimal Euclidean distance between two vertices for each
            /// edge; (Can be used to compute the CFL number);
            /// 1st index: local edge index;
            /// </summary>
            public MultidimensionalArray h_min_Edge;

            /// <summary>
            /// the maximal Euclidean distance between two vertices for each edge;
            /// (Can be used to compute the CFL number);
            /// 1st index: local edge index;
            /// </summary>
            public MultidimensionalArray h_max_Edge;

            /// <summary>
            /// For each edge that is affine-linear (i.e. <see cref="Info"/>[e]
            /// &amp; <see cref="EdgeInfo.EdgeIsAffineLinear"/> != 0), the
            /// square root of the Gram determinant; NaN for nonlinear edges.<br/>
            /// 1st index: local edge index;
            /// </summary>
            /// <remarks>
            /// Let be 
            /// \f[ 
            ///   \mathbb{R}^{D-1} 
            ///     \ni \vec{\xi} 
            ///       \mapsto
            ///         \gamma(\vec{\xi}) \in
            ///           \mathbb{R}^{D-1}
            /// \f]
            /// the mapping from the edge coordinate system to the physical coordinate system.
            /// Then the integral of \f$  f \f$  over the edge 
            /// \f$ \gamma(K_\textrm{ref}) \f$ 
            /// is given as 
            /// \f[ 
            ///   \int_{\vec{x} \in \gamma(K_\textrm{ref})} f(\vec{x}) \ \textrm{dS}
            ///   =
            ///   \int_{\xi \in K_\textrm{ref}} f(\gamma(\xi)) g(\vec{\xi}) \ \textrm{d} \vec{\xi}
            /// \f]
            /// with the square-root of the Gram determinant
            /// \f[ 
            ///   g(\vec{xi}) = \sqrt{ 
            ///      \textrm{det} ( (\partial \gamma)^T \cdot (\partial \gamma) )  
            ///   }.
            /// \f]
            /// If the transformation 
            /// \f$ \gamma \f$
            /// of the edge to the global coordinate system 
            /// is affine-linear, the Jacobian 
            /// \f$ \partial \gamma \f$
            /// is constant and 
            /// \f$ g \f$
            /// can be precomputed.
            /// (see Analysis 2, Königsberger, Springer-Verlag 2000, pp. 343)
            /// </remarks>
            public MultidimensionalArray SqrtGramian {
                get;
                internal set;
            }

            /// <summary>
            /// true, if edge <paramref name="e"/> is affine-linear, false if
            /// not
            /// </summary>
            public bool IsEdgeAffineLinear(int e) {
                return (Info[e] & EdgeInfo.EdgeIsAffineLinear) != 0;
            }

           

            /// <summary>
            /// returns the area (to be more exact: the (D-1) - dimensional
            /// measure) of the edge <paramref name="e"/>;
            /// </summary>
            /// <param name="e">local edge index</param>
            /// <returns></returns>
            public double GetEdgeArea(int e) {
                if (IsEdgeAffineLinear(e)) {
                    var Kref = this.EdgeRefElements[GetRefElementIndex(e)];
                    double area = Kref.Volume;
                    return area * SqrtGramian[e];
                } else {

                    RefElement KrefEdge = GetRefElement(e);
                    int jCell = this.CellIndices[e, 0];
                    var Kref = m_owner.Cells.GetRefElement(jCell);
                    int D = Kref.SpatialDimension;
                    double deg = m_owner.Cells.GetInterpolationDegree(jCell);
                    var Cj = m_owner.Cells.GetCell(jCell);

                    if(deg > 1)
                        deg -= 1;
                    var qr = KrefEdge.GetQuadratureRule((int)(deg * D));

                    MultidimensionalArray Metrics = this.NormalsCache.GetIntegrationMetric(qr.Nodes, e, 1);


                    //Loop over quadrature nodes
                    double area = 0;
                    for(int n = qr.Weights.GetLength(0) - 1; n >= 0; n--) {
                        area += qr.Weights[n] * Metrics[0, n];
                    }

                    return area;
                }
            }

            /// <summary>
            /// <see cref="NoOfBoundaryEdges"/>;
            /// </summary>
            private int m_NoOfBoundaryEdges;

            /// <summary>
            /// Number of edges which lie on the boundary of the physical
            /// domain i.e. edges which bound to only one cell. This are all
            /// <see cref="CellIndices"/> from index 0 to this value -1.
            /// All edges at higher indices bound to two cells.
            /// </summary>
            public int NoOfBoundaryEdges {
                get {
                    return m_NoOfBoundaryEdges;
                }
            }

            /// <summary>
            /// internal edges are all edges which do not lie on an external
            /// cell; All <see cref="CellIndices"/> from index 0 to this value
            /// are internal, edges at higher indices lie on the inter-process
            /// border.
            /// </summary>
            public int NoOfInternalEdges {
                get {
                    return Count - m_NoOfBoundaryEdges;
                }
            }

            /// <summary>
            /// total number of all edges handled on this processor;
            /// </summary>
            public int Count {
                get {
                    return CellIndices.GetLength(0);
                }
            }

            /// <summary>
            /// helper structure used in <see cref="CollectEdges"/>;
            /// </summary>
            private struct ComputeEdgesHelper {
                public int Cell1;
                public int Cell2;
                public byte FaceIndex1;
                public byte FaceIndex2;
                public byte EdgeTag;
                public EdgeInfo info;
                public int Cell1TrafoIdx;
                public int Cell2TrafoIdx;
                public bool IsPeriodic;
                public int Cell1_PeriodicTrafoIdx;
                public int Cell2_PeriodicTrafoIdx;
                
                public int EdgeKrefIndex {
                    get {
                        return (int)(EdgeInfo.EdgeSimplexIdxMask & info);
                    }
                    set {
                        Debug.Assert(value >= 0);
                        Debug.Assert(value <= (int)(EdgeInfo.EdgeSimplexIdxMask));
                        int _info = (int)info;
                        _info |= value;
                        info = (EdgeInfo)_info;
                    }
                }
            }

            /// <summary>
            /// For edge number <paramref name="e"/>, the index into
            /// <see cref="EdgeRefElements"/>.
            /// </summary>
            public int GetRefElementIndex(int e) {
                int info = (int)this.Info[e];
                int ret = ((int)(EdgeInfo.EdgeSimplexIdxMask)) & info;
                return ret;
            }

            /// <summary>
            /// For edge number <paramref name="e"/>, the respective reference element.
            /// </summary>
            public RefElement GetRefElement(int e) {
                return this.EdgeRefElements[this.GetRefElementIndex(e)];
            }
            
            /// <summary>
            /// Temporary edge data structure during assembly process.
            /// </summary>
            private List<ComputeEdgesHelper> m_EdgesTmp;

            /// <summary>
            /// temporary data structure during assembly process;
            /// 1st index: local cell index j;
            /// 2nd index: collection
            /// content: CellsToEdges[j,0] to CellsToEdges[j,N] are the N+1
            /// edges that bound to cell j.
            /// </summary>
            private List<int>[] m_CellsToEdgesTmp;

            /// <summary>
            /// this method initializes the <see cref="CellIndices"/>,
            /// <see cref="FaceIndices"/>, <see cref="NormalsForAffine"/> and
            /// <see cref="EdgeTags"/> arrays (from simpler/more
            /// straightforward data structures).
            /// </summary>
            internal void CollectEdges() {
                using (new ilPSP.Tracing.FuncTrace()) {
                    int[][] CellNeighbours = m_owner.m_Cells.CellNeighbours;
                    var CellNeighbours_Global = m_owner.m_Cells.CellNeighbours_global_tmp;

                    int Je = m_owner.Cells.Count;
                    int J = m_owner.Cells.NoOfLocalUpdatedCells;
                    int j0 = m_owner.CellPartitioning.i0;
                    long[] GlidxExternal = m_owner.Parallel.GlobalIndicesExternalCells;

                    if (m_EdgesTmp != null)
                        throw new ApplicationException("internal error.");
                    m_EdgesTmp = new List<ComputeEdgesHelper>(J * 2);
                   
#if DEBUG
                    for (int j = 0; j < J; j++) {
                        int[] CellNeigh = CellNeighbours[j];
                        foreach(int jN in CellNeigh) {
                            Debug.Assert(jN >= 0);
                            Debug.Assert(jN != j);
                        }
                    }
#endif

                    int mask;
                    unchecked {
                        mask = (int)0x80000000;
                    }

                    //Debug.Assert(false, "break0" + m_owner.MyRank);

                    // loop over all locally updated cells...
                    for (int j = 0; j < J; j++) {

                        int[] CellNeigh = CellNeighbours[j];
                        int K = CellNeigh.Length;


                        // loop over neighbor...
                        for (byte e = 0; e < K; e++) {

                            if ((CellNeigh[e] & mask) != 0)
                                continue; // this edge is already in the list
                            // we use the most significant bit to mark processed edges

                            int jNeig = CellNeigh[e];
                            int jNeigGlob;
                            if(jNeig < J) {
                                Debug.Assert(jNeig >= 0);
                                jNeigGlob = jNeig + j0;
                            } else {
                                Debug.Assert(jNeig >= J);
                                Debug.Assert(jNeig < Je);
                                jNeigGlob = (int)(GlidxExternal[jNeig - J]);
                            }

                            //Single(dnsjdkvnskj)
                            //GridCommons.Neighbour cn_je = CellNeighbours_Global[j].ElementAt(e);
                            GridCommons.Neighbour cn_je = CellNeighbours_Global[j].Single(Neigh => Neigh.Neighbour_GlobalIndex == jNeigGlob);

                            // if we reach this point, we've found a new edge
                            ComputeEdgesHelper ceh = default(ComputeEdgesHelper);
                            ceh.info = EdgeInfo.Default;
                            ceh.Cell1 = j;
                            ceh.Cell2 = CellNeigh[e];
                            ceh.FaceIndex1 = (byte)cn_je.CellFaceTag.FaceIndex;
                            ceh.FaceIndex2 = byte.MaxValue; // still unknown
                            ceh.IsPeriodic = cn_je.IsPeriodicNeighbour;
                            ceh.EdgeTag = cn_je.CellFaceTag.EdgeTag;
                            if (!cn_je.CellFaceTag.ConformalNeighborship)
                                ceh.info |= EdgeInfo.Cell1_Nonconformal;
                            
                            if (ceh.Cell2 < J) {
                                // ++++++++++++++++++++++++++++++++++++
                                // edge between cells on this processor
                                // ++++++++++++++++++++++++++++++++++++

                                int K2 = CellNeighbours[ceh.Cell2].Length;
                                int found = 0;
                                for (int e2 = 0; e2 < K2; e2++) {
                                    if (CellNeighbours[ceh.Cell2][e2] == j) {
                                        found++;
                                        CellNeighbours[ceh.Cell2][e2] |= mask;
                                    }
                                }

                                if (found <= 0 || found > 1) {
                                    long GlId0 = m_owner.m_Cells.GetCell(ceh.Cell1).GlobalID;
                                    long GlId1 = m_owner.m_Cells.GetCell(ceh.Cell2).GlobalID;

                                    if (found <= 0)
                                        throw new ApplicationException("Fatal error in cell graph: cell " + GlId0
                                            + " bounds to cell " + GlId1 + ", but not the other way around.");
                                    else
                                        throw new ApplicationException("Fatal error in cell graph: edge between " + GlId0
                                            + " and cell " + GlId1 + ", is defined multiple times.");
                                }
                            } else {
                                // +++++++++++++++++++
                                // interprocess - edge
                                // +++++++++++++++++++

                                ceh.info = EdgeInfo.Interprocess;
                            }

                            {
                                GridCommons.Neighbour cn_je2 = CellNeighbours_Global[ceh.Cell2].Single(x => x.Neighbour_GlobalIndex == (j + j0));
                                ceh.FaceIndex2 = (byte)cn_je2.CellFaceTag.FaceIndex;
                                if (!cn_je2.CellFaceTag.ConformalNeighborship)
                                    ceh.info |= EdgeInfo.Cell2_Nonconformal;


                                if (cn_je2.CellFaceTag.EdgeTag != ceh.EdgeTag) {
                                    throw new ApplicationException(string.Format("Inconsistent edge tag."));
                                }

                                if (ceh.EdgeTag > 0 && ceh.EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG)
                                    throw new ApplicationException("Found internal edge with invalid edge tag.");


                                if (cn_je.IsPeriodicNeighbour != cn_je2.IsPeriodicNeighbour) {
                                    throw new ApplicationException("inconsistent specification of periodic boundaries.");
                                }

                                if (ceh.IsPeriodic) {
                                    if (cn_je.CellFaceTag.PeriodicInverse == cn_je2.CellFaceTag.PeriodicInverse)
                                        throw new ApplicationException("inconsistent specification of periodic boundaries.");

                                    ceh.Cell1_PeriodicTrafoIdx = ((int)cn_je.CellFaceTag.EdgeTag) - GridCommons.FIRST_PERIODIC_BC_TAG;
                                    if (cn_je.CellFaceTag.PeriodicInverse) {
                                        ceh.Cell1_PeriodicTrafoIdx++;
                                        ceh.Cell1_PeriodicTrafoIdx *= -1;
                                    }

                                    ceh.Cell2_PeriodicTrafoIdx = ((int)cn_je2.CellFaceTag.EdgeTag) - GridCommons.FIRST_PERIODIC_BC_TAG;
                                    if (cn_je2.CellFaceTag.PeriodicInverse) {
                                        ceh.Cell2_PeriodicTrafoIdx++;
                                        ceh.Cell2_PeriodicTrafoIdx *= -1;
                                    }
                                }

                            }

                            // mark edges that touch ghost cells
                            if (ceh.Cell1 >= J || ceh.Cell2 >= J)
                                ceh.info |= EdgeInfo.Interprocess;

                            // add edge to list
                            m_EdgesTmp.Add(ceh);

                           
                        }
                    }


                    // remove the marking
                    int invmask = ~mask;
                    Je = CellNeighbours.Length;
                    for (int j = 0; j < Je; j++) {

                        int[] CellNeigh = CellNeighbours[j];
                        int K = CellNeigh.Length;

                        // loop over neighbours of a cell...
                        for (byte e = 0; e < K; e++) {
                            //faul
                            CellNeigh[e] &= invmask;
                        }
                    }


#if DEBUG
                    for (int j = 0; j < J; j++) {

                        int[] CellNeigh = CellNeighbours[j];
                        int K = CellNeigh.Length;

                        Debug.Assert(CellNeigh.Where(i => i < 0).Count() == 0);
                    }
#endif
                }
            }

            internal void DetermineEdgeTrafo() {
                using (new FuncTrace()) {

                    m_CellsToEdgesTmp = new List<int>[m_owner.Cells.Count];


                    // preparation: helper vars
                    // ========================

                    var GridSimplices = m_owner.Grid.RefElements;
                    int NS = GridSimplices.Length;
                    int D = this.m_owner.SpatialDimension;

                    List<Tuple<int,AffineTrafo>> e2cTrafo = new List<Tuple<int, AffineTrafo>>();

                    MultidimensionalArray[] VerticesFor_KrefEdge = this.EdgeRefElements.Select(KrefEdge => KrefEdge.Vertices).ToArray();

                    // work
                    // ====

                    // loop over all edges....
                    int skippedEdgesCount = 0;
                    List<int> skippedEdges = new List<int>();
                    for (int e = 0; e < this.m_EdgesTmp.Count; e++) {

                        // cache some vars..
                        // -----------------
                        var Edge = this.m_EdgesTmp[e];
                        int j1 = Edge.Cell1;
                        int j2 = Edge.Cell2;

                        bool conformal1 = ((Edge.info & EdgeInfo.Cell1_Nonconformal) == 0);
                        bool conformal2 = ((Edge.info & EdgeInfo.Cell2_Nonconformal) == 0);
                        
                        var K_j1 = this.m_owner.Cells.GetCell(j1);
                        var K_j2 = this.m_owner.Cells.GetCell(j2);

                        var Kref1 = this.m_owner.Cells.GetRefElement(j1); // GridSimplices[K_j1.MajorCellTypeIndex];
                        var Kref2 = this.m_owner.Cells.GetRefElement(j2); // GridSimplices[K_j2.MajorCellTypeIndex];

                        int L1 = Kref1.FaceRefElement.NoOfVertices;
                        int E1 = Kref1.NoOfFaces;
                        int L2 = Kref2.FaceRefElement.NoOfVertices;
                        int E2 = Kref2.NoOfFaces;

                        int[,] E1vtxIdx = Kref1.FaceToVertexIndices;
                        int[,] E2vtxIdx = Kref2.FaceToVertexIndices;
                        Debug.Assert(E1vtxIdx.GetLength(1) == L1);
                        Debug.Assert(E2vtxIdx.GetLength(1) == L2);

                        int face1 = Edge.FaceIndex1;
                        int face2 = Edge.FaceIndex2;


                        // first face of edge: cell j1, face e1, in local coordinates of cell j1
                        NodeSet V_f1 = Kref1.GetFaceVertices(face1);
                        //for (int l1 = 0; l1 < L1; l1++) {
                        //    for (int d = 0; d < D; d++)
                        //        V_f1[l1, d] = Kref1.Vertices[E1vtxIdx[face1, l1], d];
                        //}

                        // second face of edge: cell j2, face e2 ...
                        var V_f2 = MultidimensionalArray.Create(L2, D); //        ... in local coordinates of cell j1, obtained by transformation
                        NodeSet V_f2_in_K2 = Kref2.GetFaceVertices(face2);  // ... in local coordinates of cell j2, where they are defined
                        {
                            //for (int l2 = 0; l2 < L2; l2++) {
                            //    for (int d = 0; d < D; d++)
                            //        V_f2_in_K2[l2, d] = Kref2.Vertices[E2vtxIdx[face2, l2], d];
                            //}

                            // transform cell2 -> global
                            var V_e2_in_Global = MultidimensionalArray.Create(L2, D);
                            m_owner.TransformLocal2Global(V_f2_in_K2, V_e2_in_Global, j2);

                            // apply periodic trafo, if required
                            MultidimensionalArray _V_e2_in_Global;
                            if (Edge.IsPeriodic) {
                                // if necessary, apply periodic transformation
                                AffineTrafo PerT;
                                if (Edge.Cell2_PeriodicTrafoIdx >= 0) {
                                    PerT = m_owner.Grid.PeriodicTrafo[Edge.Cell2_PeriodicTrafoIdx];
                                } else {
                                    PerT = m_owner.Grid.InversePeriodicTrafo[Edge.Cell2_PeriodicTrafoIdx * (-1) - 1];
                                }

                                _V_e2_in_Global = MultidimensionalArray.Create(L2, D);
                                PerT.Transform(V_e2_in_Global, _V_e2_in_Global);
                            } else {
                                _V_e2_in_Global = V_e2_in_Global;
                            }
                                                        
                            // transform global -> cell1
                            if (conformal1 && conformal2) {
                                // if both edges are known to be conformal, V_e2 is (resp. should be) just a permutation of V_e1.
                                // therefore we can avoid the global-to-local transformation,
                                // which is potentially problematic for high-order curved elements.

                                if (L1 != L2)
                                    throw new ApplicationException("Cell1_Conformal, Cell1_Conformal -- flags in edge tags are fucked up.");

                                var V_e1_in_Global = MultidimensionalArray.Create(L1, D);
                                m_owner.TransformLocal2Global(V_f1, V_e1_in_Global, j1);

                                double Threshold = 1.0e-5 * Math.Min(V_e1_in_Global.MindistBetweenRows(), _V_e2_in_Global.MindistBetweenRows());

                                for (int i = 0; i < L2; i++) {
                                    double[] pt = _V_e2_in_Global.GetRow(i);

                                    double dist;
                                    int k;
                                    V_e1_in_Global.MindistRow(pt, out dist, out k);

                                    if (dist > Threshold)
                                        throw new ArgumentException("Edge is specified to be conformal, but vertex distance is far above acceptable threshold.");

                                    V_f2.SetRow(i, V_f1.GetRow(k));
                                }
                            } else {
                                bool[] Converged = new bool[_V_e2_in_Global.NoOfRows];
                                m_owner.TransformGlobal2Local(_V_e2_in_Global, V_f2, j1, Converged);
                                if (Converged.Any(t => t == false))
                                    throw new ArithmeticException("Newton divergence");
                            }
                        }


                        AffineTrafo newTrafo = null;

                        int Edg_idx;

                        if (conformal1 && conformal2) {

                            Debug.Assert(L1 == L2);

                            newTrafo = Kref1.GetFaceTrafo(face1);
                            Edg_idx = VerticesFor_KrefEdge.IndexWhere(EdgeVtx => EdgeVtx.GetLength(0) == L1); 

                        } else {

                           


                            bool bFoundIntersection = FaceIntersect(V_f1, V_f2,
                                Kref1.GetFaceTrafo(face1), Kref1.GetInverseFaceTrafo(face1),
                                VerticesFor_KrefEdge,
                                out conformal1, out conformal2, out newTrafo, out Edg_idx);

                            if(bFoundIntersection) {
                                // intersection found

                                //Debug.Assert(Edge.FaceIndex1 == byte.MaxValue || Edge.FaceIndex1 == e1);
                                //Debug.Assert(Edge.FaceIndex2 == byte.MaxValue || Edge.FaceIndex2 == e2);

                                Edge.EdgeKrefIndex = Edg_idx;


                                face1 = int.MaxValue - 10;
                                face2 = int.MaxValue - 10;
                            } else {
                                MultidimensionalArray Vtx1 = MultidimensionalArray.Create(Kref1.Vertices.Lengths);
                                m_owner.TransformLocal2Global(Kref1.Vertices, Vtx1, j1);
                                MultidimensionalArray Vtx2 = MultidimensionalArray.Create(Kref2.Vertices.Lengths);
                                m_owner.TransformLocal2Global(Kref2.Vertices, Vtx2, j2);

                                BoundingBox BB1 = new BoundingBox(Vtx1);
                                BoundingBox BB2 = new BoundingBox(Vtx2);
                                BB1.ExtendByFactor(1e-8);
                                BB2.ExtendByFactor(1e-8);
                                bool bbIntersect = BB1.Overlap(BB2);

                                throw new ApplicationException("Fatal error in grid: Cell " + K_j1.GlobalID + " and Cell "
                                + K_j2.GlobalID + " are specified to be neighbors, but the do not seem "
                                + "to have a matching edge geometrically. Bounding box overlap is " + bbIntersect + ".");
                            }
                        }

                        Edge.info &= ~EdgeInfo.Cell1_Nonconformal;
                        Edge.info &= ~EdgeInfo.Cell2_Nonconformal;
                        if (!conformal1)
                            Edge.info |= EdgeInfo.Cell1_Nonconformal;
                        if (!conformal2)
                            Edge.info |= EdgeInfo.Cell2_Nonconformal;

                        // collect transformations
                        // -----------------------
                        var Trafo1 = newTrafo; // transformation from edge to cell 1;
                        AffineTrafo Trafo2; // transformation form edge to cell 2;
                        {
                            int[] st = new int[] { 0, 0 };
                            int[] en = new int[] { D - 1, D - 1 };

                            var Urbild = MultidimensionalArray.Create(D + 1, D);
                            var Bild = MultidimensionalArray.Create(D + 1, D);
                            Urbild.ExtractSubArrayShallow(st, en).Set(V_f2.ExtractSubArrayShallow(st, en));
                            Bild.ExtractSubArrayShallow(st, en).Set(V_f2_in_K2.ExtractSubArrayShallow(st, en));
                            var K_j2ToK_j1 = AffineTrafo.FromPoints(Urbild, Bild);
                            Trafo2 = K_j2ToK_j1 * Trafo1; // transformation form edge to cell 2;
                        }

                        // cells to edges
                        // --------------
                        {
                            var cell1_edges = m_CellsToEdgesTmp[j1];
                            if(cell1_edges == null)
                                m_CellsToEdgesTmp[j1] = cell1_edges = new List<int>();
                            var cell2_edges = m_CellsToEdgesTmp[j2];
                            if(cell2_edges == null)
                                m_CellsToEdgesTmp[j2] = cell2_edges = new List<int>();

                            cell1_edges.Add(e - skippedEdgesCount);
                            cell2_edges.Add(e - skippedEdgesCount);
                        }

                        {
                            var KrefEdge = this.EdgeRefElements[Edge.EdgeKrefIndex];
                            
                            NodeSet V1 = new NodeSet(Kref1, KrefEdge.NoOfVertices, D);
                            Trafo1.Transform(KrefEdge.Vertices, V1);
                            V1.LockForever();
                            
                            NodeSet V2 = new NodeSet(Kref2, KrefEdge.NoOfVertices, D);
                            Trafo2.Transform(KrefEdge.Vertices, V2);
                            V2.LockForever();

                            var V1G = MultidimensionalArray.Create(V1.Lengths);
                            m_owner.TransformLocal2Global(V1, V1G, j1);

                            var V2G = MultidimensionalArray.Create(V1.Lengths);
                            m_owner.TransformLocal2Global(V2, V2G, j2);

                            //var JacDet1 = MultidimensionalArray.Create(1, V1.GetLength(0));
                            //Kref1.JacobianDetTransformation(V1, JacDet1, 0, K_j1.Type, K_j1.TransformationParams);
                            var JacDet1 = m_owner.JacobianDeterminat.GetValue_Cell(V1, j1, 1);
                           
                            //var JacDet2 = MultidimensionalArray.Create(1, V1.GetLength(0));
                            //Kref1.JacobianDetTransformation(V2, JacDet2, 0, K_j2.Type, K_j2.TransformationParams);
                            var JacDet2 = m_owner.JacobianDeterminat.GetValue_Cell(V2, j2, 1);

                            if(JacDet1.Min() <= 0)
                                throw new ArithmeticException("Non-positive Jacobian found in cell " + j1 + ".");
                            if(JacDet2.Min() <= 0)
                                throw new ArithmeticException("Non-positive Jacobian found in cell " + j2 + ".");
                            
                            var RelScale = Math.Max(Math.Max(V1G.MaxdistBetweenRows(), V2G.MaxdistBetweenRows()), Math.Max(JacDet1.Max(), JacDet2.Max()));

                            var Diff = V1G.CloneAs();
                            Diff.Acc(-1.0, V2G);
                            var err = Diff.L2Norm()/RelScale;
                            if (!(err <= 1.0e-8 || Edge.IsPeriodic))
                                throw new ArithmeticException("Edges do not match geometrically.");
                        }

                        int iKref1 = this.m_owner.Cells.GetRefElementIndex(j1); 
                        int iKref2 = this.m_owner.Cells.GetRefElementIndex(j2); 

                        Tuple<int,AffineTrafo> Trafo1Pair = new Tuple<int,AffineTrafo>(iKref1, Trafo1);
                        Tuple<int,AffineTrafo> Trafo2Pair = new Tuple<int,AffineTrafo>(iKref2, Trafo2);

                        Edge.Cell1TrafoIdx = e2cTrafo.IndexOf(Trafo1Pair, (a,b) => (a.Item1 == b.Item1 && a.Item2.ApproximateEquals(b.Item2)));
                        if (Edge.Cell1TrafoIdx < 0) {
                            e2cTrafo.Add(Trafo1Pair);
                            Edge.Cell1TrafoIdx = e2cTrafo.Count - 1;
                        }
                        Edge.Cell2TrafoIdx = e2cTrafo.IndexOf(Trafo2Pair, (a, b) => (a.Item1 == b.Item1 && a.Item2.ApproximateEquals(b.Item2)));
                        if (Edge.Cell2TrafoIdx < 0) {
                            e2cTrafo.Add(Trafo2Pair);
                            Edge.Cell2TrafoIdx = e2cTrafo.Count - 1;
                        }

                        // store
                        // -----
                        this.m_EdgesTmp[e] = Edge;
                    }

                    if(skippedEdgesCount > 0) {
                        int sk = 1;
                        for(int e = skippedEdges[0]; e < m_EdgesTmp.Count - skippedEdges.Count; e++) {
                            while(sk < skippedEdges.Count && e + sk == skippedEdges[sk]) {
                                sk++;
                            }
                            m_EdgesTmp[e] = m_EdgesTmp[e + sk];
                        }

                        m_EdgesTmp.RemoveRange(m_EdgesTmp.Count - skippedEdges.Count, skippedEdges.Count);

                    }
                    
                    lock(padlock) {
                        this.e2C_offet = offset_counter;
                        this.Edge2CellTrafos = e2cTrafo.Select(tt => tt.Item2).ToArray().ToList();
                        this.Edge2CellTrafosRefElementIndices = e2cTrafo.Select(tt => tt.Item1).ToArray().ToList();
                        
                        offset_counter += this.Edge2CellTrafos.Count;
                    }
                }
            }


            private static int offset_counter = 0;
            private static readonly object padlock = new object();

            /// <summary>
            /// Some hack, used by <see cref="NodeSet.GetVolumeNodeSet"/>; 
            /// only effective (un-equal 0), if more than one grid is used in the application.
            /// </summary>
            public int e2C_offet {
                get;
                private set;
            }
            
            /// <summary>
            /// Transformations from edge coordinate system to local cell
            /// coordinate systems.
            /// </summary>
            public IList<AffineTrafo> Edge2CellTrafos {
                get;
                private set;
            }

            /// <summary>
            /// Square-root of the Gramian determinat for each transformation in <see cref="Edge2CellTrafos"/>, i.e.
            /// if \f$ \myMatrix{M} \f$ 
            /// is the matrix of the transformation, this number is 
            /// \f$ \sqrt{ \operatorname{det} ( \myMatrix{M}^T \cdot \myMatrix{M} ) } \f$.
            /// </summary>
            public MultidimensionalArray Edge2CellTrafos_SqrtGramian {
                get {
                    if (m_Edge2CellTrafos_SqrtGramian == null) {
                        this.m_Edge2CellTrafos_SqrtGramian = MultidimensionalArray.Create(Edge2CellTrafos.Count);
                        for (int i = 0; i < Edge2CellTrafos.Count; i++) {
                            var tr = this.Edge2CellTrafos[i];
                            this.Edge2CellTrafos_SqrtGramian[i] = IMatrixExtensions.GEMM(tr.Matrix.Transpose(), tr.Matrix).Determinant().Sqrt();
                        }
                    }
                    return m_Edge2CellTrafos_SqrtGramian;
                }
            }

            MultidimensionalArray m_Edge2CellTrafos_SqrtGramian;

            /// <summary>
            /// For each edge-to-cell transformation, see <see cref="Edge2CellTrafos"/>,
            /// the index of the cell reference element, i.e. an index into <see cref="EdgeRefElements"/>.
            /// </summary>
            public IList<int> Edge2CellTrafosRefElementIndices {
                get;
                private set;
            }

            /// <summary>
            /// Geometric intersection of non-conformal edges.
            /// </summary>
            /// <param name="VtxFace1">
            /// vertices of one face of the first cell (face 1)
            /// </param>
            /// <param name="VtxFace2">
            /// vertices of one face of the second cell (face 2)
            /// </param>
            /// <param name="TrafoEdge">
            /// transformation from face 1 to the reference coordinate system
            /// of cell 1
            /// </param>
            /// <param name="InvTrafoEdge">
            /// inverse of <paramref name="TrafoEdge"/>
            /// </param>
            /// <param name="conformal1">
            /// true, if face 1 is fully enclosed by face  2
            /// </param>
            /// <param name="conformal2">
            /// true, if face 2 is fully enclosed by face 1
            /// </param>
            /// <param name="NewTrafo">
            /// transformation from the intersection of both faces to the 
            /// reference coordinate system (in which
            /// <paramref name="VtxFace1"/> and <paramref name="VtxFace2"/>
            /// are defined).
            /// </param>
            /// <param name="EdgeRefElementIndex">
            /// on exit, the index of the reference element for the edge,
            /// i.e. index into <see cref="EdgeRefElements"/>.
            /// </param>
            /// <returns>
            /// true, if the intersection of face 1 and 2 is non-empty.
            /// </returns>
            /// <param name="VerticesFor_KrefEdge">
            /// How many vertices the i-th edge element have?
            /// </param>
            internal static bool FaceIntersect(
                MultidimensionalArray VtxFace1,
                MultidimensionalArray VtxFace2,
                AffineTrafo TrafoEdge,
                AffineTrafo InvTrafoEdge,
                MultidimensionalArray[] VerticesFor_KrefEdge,
                out bool conformal1,
                out bool conformal2,
                out AffineTrafo NewTrafo,
                out int EdgeRefElementIndex) {

                int D = VtxFace1.GetLength(1);
                Debug.Assert(VtxFace1.GetLength(1) == VtxFace2.GetLength(1), "mismatch in spatial dimension.");

                conformal1 = false;
                conformal2 = false;
                NewTrafo = null;
                EdgeRefElementIndex = -1;

                // test if both edges match
                // ========================
                if (D == 1) {
                    // 1D -- case
                    // ++++++++++

                    double Tol = Math.Max(Math.Abs(VtxFace1[0, 0]), Math.Abs(VtxFace2[0, 0])) * 1.0e-12;

                    // (edges are points)

                    if (Math.Abs(VtxFace1[0, 0] - VtxFace2[0, 0]) < Tol) {

                        // in 1D, edges are always conformal
                        conformal1 = true;
                        conformal2 = true;
                        NewTrafo = TrafoEdge;

                        // in 1D, there is only one possible edge element (the point)
                        EdgeRefElementIndex = 0;
                        return true;
                    } else {
                        return false;
                    }

                } else {
                    // 2D, 3D -- case
                    // --------------

                    // pre-test
                    AffineManifold A1 = AffineManifold.FromPoints(VtxFace1);
                    AffineManifold A2 = AffineManifold.FromPoints(VtxFace2);

                    if (!A1.Equals(A2, 1.0e-8))
                        return false;

                    // main test
                    var VtxEdge1_Loc = InvTrafoEdge.Transform(VtxFace1);
                    var VtxEdge2_Loc = InvTrafoEdge.Transform(VtxFace2);

                    if (D == 2) {
                        // 2D: edges are lines
                        // +++++++++++++++++++

                        // check
                        Debug.Assert(VtxEdge1_Loc.GetLength(0) == 2, "number of vertices check failed, edge 1");
                        Debug.Assert(VtxEdge2_Loc.GetLength(0) == 2, "number of vertices check failed, edge 2");
                        Debug.Assert(VtxEdge1_Loc.GetLength(1) == 1, "spatial dim. check failed, edge 1");
                        Debug.Assert(VtxEdge2_Loc.GetLength(1) == 1, "spatial dim. check failed, edge 2");

                        // in 2D, there is only one possible edge element (the line)
                        EdgeRefElementIndex = 0;

                        // intersection between cells...
                        double l1 = VtxEdge1_Loc[0, 0];
                        double h1 = VtxEdge1_Loc[1, 0];
                        if (l1 > h1) {
                            double buf = l1;
                            l1 = h1;
                            h1 = buf;
                        }

                        double l2 = VtxEdge2_Loc[0, 0];
                        double h2 = VtxEdge2_Loc[1, 0];
                        if (l2 > h2) {
                            double buf = l2;
                            l2 = h2;
                            h2 = buf;
                        }

                        double lc = Math.Max(l1, l2);
                        double hc = Math.Min(h1, h2);

                        // define test tolerance
                        double Tol_Dist = ((h1 - l1) + (h2 - l2)) * 1.0e-10;
                        
                        // test for no intersection... 
                        if (h1 < l2 + Tol_Dist || h2 < l1 + Tol_Dist)
                            return false;

                        if (Math.Abs((h1 - l1) / (h2 - l2)) <= 1.0e-9 || Math.Abs((h2 - l2) / (h1 - l1)) < 1.0e-9)
                            throw new NotSupportedException("Cell size ratios too high for reliable non-conformal matching.");

                        // so there is some intersection ...
                        if (Math.Abs(l1 - l2) < Tol_Dist && Math.Abs(h1 - h2) < Tol_Dist) {
                            // no hanging node
                            // +++++++++++++++
                            conformal1 = true;
                            conformal2 = true;
                            NewTrafo = TrafoEdge;

                        } else {
                            // hanging node/nonconformal edge
                            // ++++++++++++++++++++++++++++++
                            conformal1 = ((Math.Abs(l1 - lc) <= Tol_Dist) && (Math.Abs(h1 - hc) <= Tol_Dist));
                            conformal2 = ((Math.Abs(l2 - lc) <= Tol_Dist) && (Math.Abs(h2 - hc) <= Tol_Dist));

                            var Intersect_edge = new double[,] { { lc }, { hc } };
                            var Intersect_vol = TrafoEdge.Transform(Intersect_edge);

                            Debug.Assert(VerticesFor_KrefEdge.Length == 1);
                            Debug.Assert(VerticesFor_KrefEdge[0].GetLength(0) == 2);
                            Debug.Assert(VerticesFor_KrefEdge[0].GetLength(1) == 1);

                            NewTrafo = AffineTrafo.FromPoints(VerticesFor_KrefEdge[0].To2DArray(), Intersect_vol);
                        }


                        return true;

                    } else if (D == 3) {
                        // 3D: edges are quads or triangles
                        // ++++++++++++++++++++++++++++++++


                        // For each edge V1-V2 in the first polygon,
                        //     Let H := Half-plane tangenting V1-V2, with the remaining
                        //         vertices on the "inside".
                        //     Let C := New empty polygon.
                        //     For each edge V3-V4 in the second polygon,
                        //         Let X := The intersection between V3-V4 and H.
                        //         If V3 inside H, and V4 is outside H then,
                        //             Add V3 to C.
                        //             Add X to C.
                        //         Else if both V3 and V4 lies outside H then,
                        //             Skip.
                        //         Else if V3 outside H, and V4 is inside H then,
                        //             Add X to C.
                        //         Else
                        //             Add V3 to C.
                        //     Replace the second polygon with C.



                        var first_polygon = CheckFace(VtxEdge1_Loc);
                        var second_polygon = CheckFace(VtxEdge2_Loc);
                        double ref_area = Math.Min(SpanArea(first_polygon), SpanArea(second_polygon));

                        var intersect = new List<Vector>(second_polygon);
                        bool _conformal2 = true;
                        for (int i = 0; i < first_polygon.Count; i++) {
                            Vector V1 = first_polygon[i];
                            int ip1 = (i + 1) % first_polygon.Count;
                            Vector V2 = first_polygon[ip1];

                            var H = AffineManifold.FromPoints(V1, V2);
                            var Tol_Dist = Math.Max(V1.Abs(), V2.Abs()) * 1e-10;

                            {
                                int ip2 = (i + 2) % first_polygon.Count;
                                var Pin = first_polygon[ip2];
                                if (H.PointDistance(Pin) < 0)
                                    H.Normal.Scale(-1);
                            }

                            var C = new List<Vector>();
                            for (int j = 0; j < intersect.Count; j++) {
                                var V3 = intersect[j];
                                int jp1 = (j + 1) % intersect.Count;
                                var V4 = intersect[jp1];
                                
                                bool V3_in = H.PointDistance(V3) > (-1*Tol_Dist);
                                bool V4_in = H.PointDistance(V4) > (-1*Tol_Dist);

                                if (V3_in && !V4_in) {
                                    C.Add(V3);
                                    var Hx = AffineManifold.FromPoints(V3, V4);
                                    var X = AffineManifold.Intersect2D(H, Hx);
                                    C.Add(X);
                                    _conformal2 = false;
                                } else if (!V3_in && !V4_in) {
                                    // skip
                                    _conformal2 = false;
                                } else if (!V3_in && V4_in) {
                                    var Hx = AffineManifold.FromPoints(V3, V4);
                                    var X = AffineManifold.Intersect2D(H, Hx);
                                    C.Add(X);
                                    _conformal2 = false;
                                } else {
                                    C.Add(V3);
                                }
                            }

                            intersect = C;
                        }


                        if (second_polygon.Count <= 2) {
                            return false;
                        }

                        double span_area = SpanArea(intersect);
                        if (span_area < ref_area * (1.0e-10)) {
                            return false;
                        }

                        bool _conformal1 = true;
                        {
                            int N2 = second_polygon.Count;
                            int N1 = first_polygon.Count;
                            for (int i = 0; i < N2; i++) {
                                int ip1 = (i + 1) % N2;
                                var H = AffineManifold.FromPoints(second_polygon[i], second_polygon[ip1]);
                                var Tol_Dist = Math.Max(second_polygon[i].L2Norm(), second_polygon[ip1].L2Norm()) * 1e-10;
                                {
                                    int ip2 = (i + 2) % N2;
                                    var Pin = second_polygon[ip2];
                                    if (H.PointDistance(Pin) < 0)
                                        H.Normal.Scale(-1);
                                }


                                for (int j = 0; j < N1; j++) {
                                    if (H.PointDistance(first_polygon[j]) < -1*Tol_Dist) {
                                        _conformal1 = false;
                                        break;
                                    }
                                }
                            }
                        }


                        conformal1 = _conformal1;
                        conformal2 = _conformal2;


                        for (EdgeRefElementIndex = 0; EdgeRefElementIndex < VerticesFor_KrefEdge.Length; EdgeRefElementIndex++)
                            if (VerticesFor_KrefEdge[EdgeRefElementIndex].GetLength(0) == intersect.Count)
                                break;

                        if (_conformal1 && _conformal2) {
                            NewTrafo = TrafoEdge;
                        } else {
                            // compute transformation from edge to computed intersection
                            var Intersect_edge = new double[intersect.Count, 2];
                            for (int i = 0; i < intersect.Count; i++) {
                                Intersect_edge[i, 0] = intersect[i][0];
                                Intersect_edge[i, 1] = intersect[i][1];
                            }
                            var Intersect_vol = TrafoEdge.Transform(Intersect_edge);

                            var VtxKref = CheckFace(VerticesFor_KrefEdge[EdgeRefElementIndex]);
                            double[,] VtxKrefSpan = new double[3, 2]  {
                                { VtxKref[0][0], VtxKref[0][1] },
                                { VtxKref[1][0], VtxKref[1][1] },
                                { VtxKref[2][0], VtxKref[2][1] }
                            };

                            NewTrafo = AffineTrafo.FromPoints(
                                VtxKrefSpan,
                                Intersect_vol.GetSubMatrix(0, 3, 0, 3));
                        }




                        return true;

                    } else {
                        throw new NotSupportedException("spatial dimension " + D + " is not implemented.");
                    }
                }
            }

            private static double SpanArea(List<Vector> second_polygon) {
                double span_area;
                int N2 = second_polygon.Count;
                double sp1x = second_polygon[1][0] - second_polygon[0][0];
                double sp1y = second_polygon[1][1] - second_polygon[0][1];
                double sp2x = second_polygon[N2 - 1][0] - second_polygon[0][0];
                double sp2y = second_polygon[N2 - 1][1] - second_polygon[0][1];
                span_area = sp1x * sp2y - sp1y * sp2x;
                return Math.Abs(span_area);
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="VtxEdge1_Loc"></param>
            /// <returns></returns>
            private static List<Vector> CheckFace(MultidimensionalArray _VtxEdge1_Loc) {

                MultidimensionalArray VtxEdge1_Loc = MultidimensionalArray.Create(_VtxEdge1_Loc.Lengths);
                VtxEdge1_Loc.Set(_VtxEdge1_Loc);
                

                int N1;
                N1 = VtxEdge1_Loc.GetLength(0);
                if (N1 == 4) {
                    // quad: due to the strange vertex numbering in bosss, exchange vertex 2 and 3

                    for (int i = 0; i < 2; i++) {
                        double tmp;
                        tmp = VtxEdge1_Loc[2, i];
                        VtxEdge1_Loc[2, i] = VtxEdge1_Loc[3, i];
                        VtxEdge1_Loc[3, i] = tmp;
                    }
                } else if (N1 == 3) {

                } else {
                    throw new NotSupportedException();
                }


                var ret = new List<Vector>();
                for (int i = 0; i < N1; i++) {
                    ret.Add(VtxEdge1_Loc.GetRowPt(i));
                }
                return ret;
            }

            /// <summary>
            /// sets, for each entry in <see cref="m_EdgesTmp"/>, it sets
            /// <see cref="ComputeEdgesHelper.EdgeTag"/>;
            /// </summary>
            internal void SetEdgeTags() {
                using (new FuncTrace()) {

                    var Cells = m_owner.Grid.Cells;

                    for (int jEdge = 0; jEdge < m_EdgesTmp.Count; jEdge++) {
                        ComputeEdgesHelper Edge_j = m_EdgesTmp[jEdge];
                        byte edgeTag_cell1 = 0;
                        var Cell1 = m_owner.Cells.GetCell(Edge_j.Cell1);


                        if (Edge_j.Cell2 >= 0) {
                            // an internal edge
                            // that can be:
                            //  - 'normal' (edge tag must be 0)
                            //  - periodic

                            /*
                            long globalId_1 = Cell1.GlobalID;
                            var Cell2 = m_owner.Cells.GetCell(Edge_j.Cell2);
                            long globalId_2 = Cell2.GlobalID;

                            var N1 = default(CellFaceTag);
                            if (Cell1.CellFaceTags != null)
                                N1 = Cell1.CellFaceTags.SingleOrDefault(
                                    X => X.NeighCell_GlobalID < Cell2.GlobalID && X.FaceIndex == Edge_j.FaceIndex1);
                            var N2 = default(CellFaceTag);
                            if (Cell2.CellFaceTags != null)
                                N2 = Cell2.CellFaceTags.SingleOrDefault(
                                    X => X.NeighCell_GlobalID < Cell1.GlobalID && X.FaceIndex == Edge_j.FaceIndex2);

                            if (N1.EdgeTag != N2.EdgeTag) {
                                throw new ApplicationException(string.Format("Mismatch in edge tags"));
                            }

                            if (N1.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                                if (!(N1.PeriodicInverse == false && N2.PeriodicInverse == true))
                                    throw new ApplicationException("periodic mapping fucked up.");
                            }
                            */
                        } else {
                            // boundary edge
                            // +++++++++++++

                            if (Cell1.CellFaceTags != null) {
                                try {
                                    var bnd = Cell1.CellFaceTags.SingleOrDefault(X => X.NeighCell_GlobalID < 0 && X.FaceIndex == Edge_j.FaceIndex1);
                                    edgeTag_cell1 = bnd.EdgeTag;
                                } catch (InvalidOperationException) {
                                    throw new ApplicationException(string.Format("Unable to uniquely determine boundary conditions for cell (globalID = {0}, FaceIndex = {1}) -- boundary condition defined multiple times?.", Cell1.GlobalID, Edge_j.FaceIndex1));
                                }
                            }

                        }

                        if (edgeTag_cell1 != 0) {
                            Edge_j.EdgeTag = edgeTag_cell1;
                            m_EdgesTmp[jEdge] = Edge_j; // necessary, because we are dealing with structs
                        }
                    }
                }
            }

           

            /// <summary>
            /// finds all boundary edges
            /// </summary>
            internal void CollectBoundaryEdges() {
                using (new FuncTrace()) {
                    int J = this.m_owner.Cells.NoOfLocalUpdatedCells;
                    int Jglob = this.m_owner.Grid.NumberOfCells;

                    var CellNeighbours_global = m_owner.Cells.CellNeighbours_global_tmp;
                    

                    // loop over cells...
                    for (int j = 0; j < J; j++) {
                        var Kj = this.m_owner.Cells.GetCell(j);
                        var Kref = this.m_owner.Cells.GetRefElement(j);
                        int S = Kref.NoOfFaces;


                        // all edges (indices into 'm_TemporaryEdges')
                        // that are connected with cell 'j'
                        if (m_CellsToEdgesTmp[j] == null) {
                            m_CellsToEdgesTmp[j] = new List<int>();
                        }
                        var allEdges = m_CellsToEdgesTmp[j];

                        //if (Kj.CellFaceTags != null && Kj.CellFaceTags.Length != S) {
                        //    throw new ApplicationException("Fatal error in grid edge tags: "
                        //        + "cell " + Kj.GlobalID + " has " + S + " faces, but only "
                        //        + Kj.CellFaceTags.Length + " edge tags are defined. "
                        //        + "Alternatively, 'EdgeTags' may be set to null.");
                        //}

                        bool allconformal = true;
                        int markCounter = 0;
                        for (int s = 0; s < allEdges.Count; s++) { // loop over all edges that are 
                            //                                        currently attached to cell 'j'
                            var Edge = m_EdgesTmp[allEdges[s]];

                            if (j == Edge.Cell1) {
                                if ((Edge.info & EdgeInfo.Cell1_Nonconformal) != 0)
                                    allconformal = false;
                                else
                                    markCounter++;

                            } else if (j == Edge.Cell2) {
                                if ((Edge.info & EdgeInfo.Cell2_Nonconformal) != 0)
                                    allconformal = false;
                                else
                                    markCounter++;
                            } else {
                                throw new ApplicationException("internal error.");
                            }
                        }

                        Debug.Assert(!(allconformal && markCounter > S));

                        if (allconformal && markCounter >= S) {
                            // inner cell, all edges conformal
                            // -> NO NEW EDGES WILL BE FOUND
                            // ++++++++++++++++++++++++++++++++
                            continue;
                        }


                        var CellNeighbours_global_j = CellNeighbours_global[j];


                        // loop over cell faces ...
                        for (int s = 0; s < S; s++) {

                            // select all edges that bound to face 's'
                            var neigs = allEdges.Where(x => (
                                   (m_EdgesTmp[x].FaceIndex1 == s && m_EdgesTmp[x].Cell1 == j)
                                || (m_EdgesTmp[x].FaceIndex2 == s && m_EdgesTmp[x].Cell2 == j)));

                            byte EdgeTag;
                            {
                                IEnumerable<CellFaceTag> _Q = (Kj.CellFaceTags == null) ? new CellFaceTag[0] : (Kj.CellFaceTags.Where(x => x.FaceIndex == s && x.EdgeTag != 0));

                                IEnumerable<GridCommons.Neighbour> _R = CellNeighbours_global_j.Where(delegate(GridCommons.Neighbour cn) {
                                    return ((cn.Neighbour_GlobalIndex >= Jglob) && (cn.CellFaceTag.FaceIndex == s));
                                });

                                if (_Q.Count() > 0 && _R.Count() > 0) {
                                    throw new NotSupportedException("Boundary conditions must be specified either by edge tags or by boundary-condition elements. Both options simultaneously is not supported.");
                                }

                                if (_Q.Count() > 1)
                                    throw new NotSupportedException("Found more than one EdgeTag for a boundary condition, for some face.");

                                if (_Q.Count() > 0) {
                                    EdgeTag = _Q.First().EdgeTag;
                                } else if (_R.Count() > 0) {

                                    if (_R.Count() > 1) {
                                        throw new NotImplementedException("currently, there is no support for non-conformal boundary-condition elements.");
                                    }

                                    var cn = _R.First();
                                    if (!cn.CellFaceTag.ConformalNeighborship)
                                        throw new NotImplementedException("currently, there is no support for non-conformal boundary-condition elements.");


                                    BCElement bcc = m_owner.m_BcCells_tmp[cn.Neighbour_GlobalIndex];

                                    EdgeTag = bcc.EdgeTag;

                                } else {
                                    // no def for boundary edge tag
                                    EdgeTag = 0;
                                }
                            }


                            if (neigs.Count() <= 0) {
                                // no edge found that is attached to face 's' -->
                                // face 's' of cell 'j' will form a conformal boundary edge
                                // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                var Edge = default(ComputeEdgesHelper);

                                Edge.Cell1 = j;
                                Edge.FaceIndex1 = (byte)s;

                                Edge.Cell2 = -12345;
                                Edge.FaceIndex2 = byte.MaxValue;


                                var trafo = Kref.GetFaceTrafo(s);
                                int iKref = this.m_owner.Cells.GetRefElementIndex(j);

                                Edge.Cell1TrafoIdx = -1;
                                for(int iTrafo = 0; iTrafo < this.Edge2CellTrafos.Count; iTrafo++) {
                                    if(this.Edge2CellTrafosRefElementIndices[iTrafo] != iKref)
                                        continue;

                                    if(!trafo.ApproximateEquals(this.Edge2CellTrafos[iTrafo]))
                                        continue;

                                    Debug.Assert(Edge.Cell1TrafoIdx < 0, "multiple matches");
                                    Edge.Cell1TrafoIdx = iTrafo;
                                }
                                                                
                                if (Edge.Cell1TrafoIdx < 0) {
                                    // not in list: must add
                                    // +++++++++++++++++++++
                                    this.Edge2CellTrafos.Add(trafo);
                                    Edge.Cell1TrafoIdx = this.Edge2CellTrafos.Count - 1;
                                    this.Edge2CellTrafosRefElementIndices.Add(iKref);
                                                                        
                                    lock(padlock) {
                                        offset_counter += 1;
                                    }
                                }

                                Edge.Cell2TrafoIdx = -456778;

                                Edge.info |= EdgeInfo.Boundary;

                                Edge.EdgeTag = EdgeTag;

                                m_EdgesTmp.Add(Edge);
                                allEdges.Add(m_EdgesTmp.Count - 1);
                                this.m_NoOfBoundaryEdges++;

                            } else {
                                var Edge = m_EdgesTmp[neigs.First()];

                                if ((j == Edge.Cell1 && ((Edge.info & EdgeInfo.Cell1_Nonconformal) == 0))
                                    || (j == Edge.Cell2 && ((Edge.info & EdgeInfo.Cell2_Nonconformal) == 0))) {
                                    // face 's' of cell 'j' is completely covered by 'Edge' -> no boundary edge
                                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                    Debug.Assert(neigs.Count() == 1); // in case of conformal edge, there can be only one neighbour

                                    if (Kj.CellFaceTags != null && (EdgeTag > 0 && EdgeTag < GridCommons.FIRST_PERIODIC_BC_TAG))
                                        throw new ApplicationException("Fatal error in grid edge tags: "
                                            + "cell " + Kj.GlobalID + " has edge tag " + Kj.CellFaceTags[s]
                                            + " for face " + s + ", but this face is completely internal. "
                                            + "(EdgeTag is expected to be 0.)");

                                    continue;
                                } else {
                                    // Difficult part: nonconformal edges:
                                    // Determine whether face 's' of cell 'j' is completely covered
                                    // by edges, or if there is some 'leftover'!
                                    // The leftover forms one or more boundary edges.
                                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                    continue;
                                    //throw new NotImplementedException("todo");
                                    // deactivate if you're sure that there are absolutely no hanging nodes on the boundary,
                                    // i.e. noting like:
                                    //
                                    //   x-------xx-----x
                                    //   |       ||     |
                                    //   |       ||     |
                                    //   |       |x-----x
                                    //   |       |
                                    //   |       |
                                    //   |       |
                                    //   x-------x

                                }

                            }
                        }
                    }
                }
            }

            /// <summary>
            /// sharing of edges between processors
            /// </summary>
            internal void NegogiateNeighbourship() {
                using (new FuncTrace()) {
                    int Size = m_owner.MpiSize;
                    int myRank = m_owner.MpiRank;
                    int E = m_EdgesTmp.Count;
                    int Jupdt = this.m_owner.Cells.NoOfLocalUpdatedCells;
                    int[][] SendLists = m_owner.Parallel.SendCommLists;
                    int[] InsertIndex = m_owner.Parallel.RcvCommListsInsertIndex;
                    int[] RcvNoOfItems = m_owner.Parallel.RcvCommListsNoOfItems;
                    long[] GidxExt = m_owner.m_Parallel.GlobalIndicesExternalCells;
                    var CellPart = m_owner.CellPartitioning;
                    int D = m_owner.SpatialDimension;

                    // put boundary edges to the beginning
                    // ===================================
                    {
                        int Ebnd = this.NoOfBoundaryEdges;
                        int Eint = E - Ebnd;

#if DEBUG
                        for (int e = 0; e < E; e++) {
                            if (e < Eint)
                                Debug.Assert(m_EdgesTmp[e].Cell2 >= 0);
                            else
                                Debug.Assert(m_EdgesTmp[e].Cell2 < 0);
                        }
#endif

                        var BndEdges = m_EdgesTmp.GetSubVector(Eint, Ebnd);
                        m_EdgesTmp.RemoveRange(Eint, Ebnd);
                        m_EdgesTmp.InsertRange(0, BndEdges);

                        var C2E = this.m_CellsToEdgesTmp;
                        for (int j = 0; j < C2E.Length; j++) {
                            var C2E_j = C2E[j];
                            int K = C2E_j.Count;
                            for (int k = 0; k < K; k++) {
                                int iEdge = C2E_j[k];
                                Debug.Assert(iEdge >= 0 && iEdge < E);

                                if (iEdge >= Eint) {
                                    // correct re-ordering of boundary-edges
                                    iEdge -= Eint;
                                    Debug.Assert(iEdge >= 0 && iEdge < Ebnd);
                                } else {
                                    // correct re-ordering of internal-edges
                                    iEdge += Ebnd;
                                    Debug.Assert(iEdge >= Ebnd && iEdge < E);
                                }

                                C2E_j[k] = iEdge;
                            }
                        }
                    }


                    // find local indices of edges on other processors
                    // (and synchronize Edge-To-Cell transformations on different processors)
                    // ======================================================================
                    
                    // 1st index: edge index
                    // 2nd index: enumeration
                    // Item1: process rank R
                    // Item2: index of edge on rank R
                    Tuple<int, int>[][] EdgeIndicesOnOtherProcessors = new Tuple<int, int>[E][];
                    {
                        for (int e = 0; e < E; e++) {
                            EdgeIndicesOnOtherProcessors[e] = new Tuple<int, int>[] { new Tuple<int, int>(myRank, e) };
                        }


                        //var SendData = new Dictionary<int, Tuple<int, long, long, double[,], double[,]>[]>();
                        //var Data = new List<Tuple<int, long, long, double[,], double[,]>>();
                        var SendData = new Dictionary<int, Tuple<int, long, long, AffineTrafo, AffineTrafo>[]>();
                        var Data = new List<Tuple<int, long, long, AffineTrafo, AffineTrafo>>();
                        var C2E = this.m_CellsToEdgesTmp;
                        var EdgesTmp = this.m_EdgesTmp;

                        double[,] OriginAndNframe = new double[Math.Max(D, 1), Math.Max(D - 1, 1)];
                        for (int d = 0; d < (D - 1); d++) {
                            OriginAndNframe[d + 1, d] = 1.0;
                        }

                        for (int targRank = 0; targRank < Size; targRank++) {
                            var SndList = SendLists[targRank];

                            if (SndList == null)
                                continue;

                            Data.Clear();

                            foreach (int jCell in SndList) {
                                var C2E_j = C2E[jCell];
                                foreach (int iEdge in C2E_j) {
                                    var Edge = EdgesTmp[iEdge];

                                    if (Edge.Cell1 >= Jupdt || Edge.Cell2 >= Jupdt) {
                                        // found some edge 

                                        Debug.Assert(Edge.Cell1 < Jupdt || Edge.Cell2 < Jupdt);
                                        Debug.Assert(Edge.Cell1 == jCell || Edge.Cell2 == jCell);


                                        long GIdx1, GIdx2;
                                        int targRankE;
                                        if (Edge.Cell1 >= Jupdt) {
                                            // Cell1 is external
                                            // +++++++++++++++++
                                            Debug.Assert(Edge.Cell2 < Jupdt);

                                            GIdx1 = GidxExt[Edge.Cell1 - Jupdt];
                                            GIdx2 = CellPart.i0 + Edge.Cell2;
                                            targRankE = CellPart.FindProcess(GIdx1);
                                            Debug.Assert(targRankE != myRank);
                                            Debug.Assert(CellPart.FindProcess(GIdx2) == myRank);
                                        } else {
                                            // Cell2 is external
                                            // +++++++++++++++++
                                            Debug.Assert(Edge.Cell2 >= Jupdt);

                                            GIdx1 = CellPart.i0 + Edge.Cell1;
                                            GIdx2 = GidxExt[Edge.Cell2 - Jupdt];
                                            targRankE = CellPart.FindProcess(GIdx2);
                                            Debug.Assert(targRankE != myRank);
                                            Debug.Assert(CellPart.FindProcess(GIdx1) == myRank);
                                        }

                                        if (targRankE != targRank)
                                            continue;

                                        AffineTrafo T1 = this.Edge2CellTrafos[Edge.Cell1TrafoIdx];
                                        AffineTrafo T2 = this.Edge2CellTrafos[Edge.Cell2TrafoIdx];

                                        //double[,] T1_OriginAndNframe = T1.Transform(OriginAndNframe);
                                        //double[,] T2_OriginAndNframe = T2.Transform(OriginAndNframe);

                                        Data.Add(new Tuple<int, long, long, AffineTrafo, AffineTrafo>(iEdge, GIdx1, GIdx2, T1, T2));
                                    }
                                }

                            }

                            SendData.Add(targRank, Data.ToArray());
                        }

                        var RcvData = SerialisationMessenger.ExchangeData(SendData, csMPI.Raw._COMM.WORLD);
                        foreach (var kv in RcvData) {
                            int originRank = kv.Key;
                            var RData = kv.Value;

                            foreach (var t in RData) {
                                int iEdge_foreign = t.Item1;
                                long Gidx1 = t.Item2;
                                long Gidx2 = t.Item3;
                                //double[,] T1_OriginAndNframe = t.Item4;
                                //double[,] T2_OriginAndNframe = t.Item5;
                                AffineTrafo T1 = t.Item4;
                                AffineTrafo T2 = t.Item5;

                                int jCell_loc;
                                if (CellPart.IsInLocalRange(Gidx1)) {
                                    Debug.Assert(CellPart.FindProcess(Gidx2) == originRank);
                                    jCell_loc = CellPart.TransformIndexToLocal((int)Gidx1);
                                } else {
                                    Debug.Assert(CellPart.IsInLocalRange(Gidx2));
                                    Debug.Assert(CellPart.FindProcess(Gidx1) == originRank);
                                    jCell_loc = CellPart.TransformIndexToLocal((int)Gidx2);
                                }


                                var C2E_jCell_loc = C2E[jCell_loc];
                                int matchCount = 0;
                                int iEdge_local = int.MinValue;
                                foreach (int iEdge in C2E_jCell_loc) {
                                    var Edge = EdgesTmp[iEdge];
                                    bool bEdgeChanged = false;
                                    long _Gidx1 = Edge.Cell1 < Jupdt ? Edge.Cell1 + CellPart.i0 : GidxExt[Edge.Cell1 - Jupdt];
                                    long _Gidx2 = Edge.Cell2 < Jupdt ? Edge.Cell2 + CellPart.i0 : GidxExt[Edge.Cell2 - Jupdt];

                                    bool evenMatch = (_Gidx1 == Gidx1 && _Gidx2 == Gidx2);
                                    bool oddMatch = (_Gidx1 == Gidx2 && _Gidx2 == Gidx1);

                                    if (!(evenMatch || oddMatch))
                                        continue;

                                    // in certain situations, which can occur e.g. with curved elements (e.g. a ring consisting of two curved quads)
                                    // there may actually be more than one edge between the same pair of cells
                                    // therefore, we have to check also the cell centers

                                    AffineTrafo _T1 = null, _T2 = null;
                                    if (evenMatch) {
                                        _T1 = this.Edge2CellTrafos[Edge.Cell1TrafoIdx];
                                        _T2 = this.Edge2CellTrafos[Edge.Cell2TrafoIdx];
                                    } else {
                                        Debug.Assert(oddMatch == true);
                                        _T1 = this.Edge2CellTrafos[Edge.Cell2TrafoIdx];
                                        _T2 = this.Edge2CellTrafos[Edge.Cell1TrafoIdx];
                                    }

                                    double[,] _T1_OriginAndNframe = _T1.Transform(OriginAndNframe);
                                    //double[,] _T2_OriginAndNframe = _T2.Transform(OriginAndNframe);
                                    double[,] T1_OriginAndNframe = T1.Transform(OriginAndNframe);

                                    double h = 1.0;
                                    if (D > 1) {
                                        h = 0;
                                        for (int d = D - 1; d >= 0; d--) {
                                            double dx = _T1_OriginAndNframe[0, d] - _T1_OriginAndNframe[1, d];
                                            h += dx * dx;
                                        }
                                    }
                                    double dist = 0;
                                    for (int d = D - 1; d >= 0; d--) {
                                        double dx = _T1_OriginAndNframe[0, d] - T1_OriginAndNframe[0, d];
                                        dist += dx * dx;
                                    }
                                    double normDist = dist / h;

                                    if (normDist > 1.0e-12)
                                        // center of received edge differs from local edge => must be different
                                        continue;

                                    // if we reached this point, we have a match
                                    matchCount++;

                                    if (originRank < myRank) {
                                        // on the higher rank, we ensure that the transformations match the lower rank

                                        int iT1 = this.Edge2CellTrafos.IndexOf(T1);
                                        if (iT1 < 0) {
                                            this.Edge2CellTrafos.Add(T1);
                                            iT1 = this.Edge2CellTrafos.Count - 1;
                                            this.Edge2CellTrafosRefElementIndices.Add(0);
                                            Debug.Assert(this.m_owner.Grid.RefElements.Length == 1);
                                        }
                                        int iT2 = this.Edge2CellTrafos.IndexOf(T2);
                                        if (iT2 < 0) {
                                            this.Edge2CellTrafos.Add(T2);
                                            iT2 = this.Edge2CellTrafos.Count - 1;
                                            this.Edge2CellTrafosRefElementIndices.Add(0);
                                            Debug.Assert(this.m_owner.Grid.RefElements.Length == 1);
                                        }

                                        if (evenMatch) {
                                            if (Edge.Cell1TrafoIdx != iT1) {
                                                bEdgeChanged = true;
                                                Edge.Cell1TrafoIdx = iT1;
                                            }
                                            if (Edge.Cell2TrafoIdx != iT2) {
                                                bEdgeChanged = true;
                                                Edge.Cell2TrafoIdx = iT2;
                                            }
                                        } else {
                                            Debug.Assert(oddMatch == true);
                                            if (Edge.Cell1TrafoIdx != iT2) {
                                                bEdgeChanged = true;
                                                Edge.Cell1TrafoIdx = iT2;
                                            }
                                            if (Edge.Cell2TrafoIdx != iT1) {
                                                bEdgeChanged = true;
                                                Edge.Cell2TrafoIdx = iT1;
                                            }
                                        }
                                    }

                                    if (bEdgeChanged) {
                                        EdgesTmp[iEdge] = Edge;
                                    }

                                    iEdge_local = iEdge;
                                }

                                if (matchCount != 1)
                                    throw new ApplicationException("error in algorithm");

                                (new Tuple<int, int>(originRank, iEdge_foreign)).AddToArray(ref EdgeIndicesOnOtherProcessors[iEdge_local]);
                            }
                        }
                    }

                    // Ownership negogiation
                    // ======================================================================
                                        
                    int[] EdgePermuation;
                    int NoOfPureLocal;
                    int NoOfShOwned;
                    int NoOfShForeign;
                    int NoOfExternal;
                    int NoOfPeriodicElim;
                    int[][] EdgeSendLists;
                    int[][] EdgeInsertLists;
                    Tuple<int, int>[] LocalId;
                    NegogiateOwnership(csMPI.Raw._COMM.WORLD, EdgeIndicesOnOtherProcessors,
                        out EdgePermuation, out NoOfPureLocal, out NoOfShOwned, out NoOfShForeign, out NoOfPeriodicElim, out NoOfExternal,
                        out EdgeSendLists, out EdgeInsertLists,
                        out LocalId);
                    Debug.Assert(NoOfPeriodicElim == 0); // there is no duplicate representation of edges due to periodicity
                    Debug.Assert(LocalId.Length == 0); // there is no duplicate representation of edges due to periodicity

                    // apply permutation
                    // ======================================================================
                    {
                        int[] invEdgePermuation = new int[E];
                        for (int e = 0; e < E; e++) {
                            invEdgePermuation[EdgePermuation[e]] = e;
                        }

                        // permute edges
                        var newEdgesTmp = new List<ComputeEdgesHelper>(E);
                        for (int e = 0; e < E; e++) {
                            newEdgesTmp.Add(this.m_EdgesTmp[EdgePermuation[e]]);
                        }
                        this.m_EdgesTmp = newEdgesTmp;

                        // permute cell-to-edge transformation
                        var C2E = this.m_CellsToEdgesTmp;
                        for (int j = 0; j < C2E.Length; j++) {
                            var C2E_j = C2E[j];
                            int K = C2E_j.Count;
                            for (int k = 0; k < K; k++) {
                                C2E_j[k] = invEdgePermuation[C2E_j[k]];
                            }
                        }

                        // transform send lists (were based on old indices)
                        for (int rnk = 0; rnk < Size; rnk++) {
                            if (EdgeSendLists[rnk] != null) {
                                EdgeSendLists[rnk] = EdgeSendLists[rnk].Select(k => invEdgePermuation[k]).ToArray();
                            }
                        }
                        
                        // transform insert lists (were based on old indices)
                        for (int rnk = 0; rnk < Size; rnk++) {
                            if (EdgeInsertLists[rnk] != null) {
                                EdgeInsertLists[rnk] = EdgeInsertLists[rnk].Select(k => invEdgePermuation[k]).ToArray();
                            }
                        }
                    }

                    // set data
                    // ========
                    this.NoOfPurelyLocal = NoOfPureLocal;
                    this.NoOfRelayed = NoOfShOwned;
                    this.NoOfBorrowed = NoOfShForeign;
                    //this.NoOfExternal = NoOfExternal;
                    this.EdgeSendLists = EdgeSendLists;
                    this.EdgeInsertLists = EdgeInsertLists;
                }
            }

            /// <summary>
            /// Number of edges which are used only by locally updated cells
            /// and do not bound to any cell which is exchanged over MPI;
            /// </summary>
            public int NoOfPurelyLocal {
                get;
                private set;
            }

            /// <summary>
            /// Number of shared edges (between by locally updated cells and
            /// external cells) which are 'owned' by the current MPI process.
            /// </summary>
            public int NoOfRelayed {
                get;
                private set;
            }

            /// <summary>
            /// Number of edges 
            /// which are 'owned' by the current MPI process.
            /// </summary>
            public int NoOfOwned {
                get {
                    return NoOfPurelyLocal + NoOfRelayed;
                }
            }

            /// <summary>
            /// Number of shared edges (between by locally updated cells and
            /// external cells) which are 'owned' by other MPI process.
            /// </summary>
            public int NoOfBorrowed {
                get;
                private set;
            }

            /// <summary>
            /// content: which of the shared/foreign edges must be send to
            /// other processors
            ///  - 1st index: MPI rank of target processor 'R' 
            ///  - 2nd index: enumeration
            /// </summary>
            public int[][] EdgeSendLists {
                get;
                private set;
            }


            /// <summary>
            /// content: where the shared/foreign edges received by other
            /// processors must be inserted 
            /// 1st index: MPI rank of target processor 'R' <br/>
            /// 2nd index: enumeration
            /// </summary>
            public int[][] EdgeInsertLists {
                get;
                private set;
            }

            /// <summary>
            /// converts temporary data structures in permanent ones
            /// </summary>
            internal void FinalizeAssembly() {
                using (new FuncTrace()) {
                    int E = m_EdgesTmp.Count;

                    this.CellIndices = new int[E, 2];
                    this.FaceIndices = new byte[E, 2];
                    this.Edge2CellTrafoIndex = new int[E, 2];
                    this.EdgeTags = new byte[E];
                    this.Info = new EdgeInfo[E];

                    int D = m_owner.SpatialDimension;

                    this.Edge2CellTrafosRefElementIndices = this.Edge2CellTrafosRefElementIndices.ToList().AsReadOnly();
                    this.Edge2CellTrafos = this.Edge2CellTrafos.ToList().AsReadOnly();

                    for (int e = 0; e < E; e++) {
                        ComputeEdgesHelper edg = m_EdgesTmp[e];

                        this.CellIndices[e, 0] = edg.Cell1;
                        this.CellIndices[e, 1] = edg.Cell2;

                        this.Edge2CellTrafoIndex[e, 0] = edg.Cell1TrafoIdx;
                        this.Edge2CellTrafoIndex[e, 1] = edg.Cell2TrafoIdx;

                        this.FaceIndices[e, 0] = edg.FaceIndex1;
                        this.FaceIndices[e, 1] = edg.FaceIndex2;

                        this.EdgeTags[e] = edg.EdgeTag;
                        this.Info[e] = edg.info;


                        var K_j1 = m_owner.Cells.GetCell(edg.Cell1);
                        if (K_j1.Type.IsLinear() || D == 1) {
                            this.Info[e] |= EdgeInfo.EdgeIsAffineLinear;
                        }

                        if (edg.Cell2 >= 0) {
                            var K_j2 = m_owner.Cells.GetCell(edg.Cell2);
                            if (K_j2.Type.IsLinear() || D == 1) {
                                this.Info[e] |= EdgeInfo.EdgeIsAffineLinear;
                            }
                        }

                        //if ((this.Info[e] & EdgeInfo.Boundary) != 0)
                        //    this.m_NoOfBoundaryEdges++;

                        Debug.Assert(((this.Info[e] & EdgeInfo.Boundary) != 0) == (edg.Cell2 < 0));
                    }

                    m_EdgesTmp = null;
                }
            }

            /// <summary>
            /// initializes <see cref="CellData.Cells2Edges"/>;
            /// </summary>
            internal void InitCells2Edges() {
                using (new FuncTrace()) {
                    var c2e_tmp = this.m_CellsToEdgesTmp;
                    int J = m_owner.Cells.NoOfLocalUpdatedCells;
                    Debug.Assert(c2e_tmp.Length >= J);

                    var EdgTmp = m_owner.Edges.m_EdgesTmp;

                    var Cells2Edges = new int[J][];
                    m_owner.Cells.Cells2Edges = Cells2Edges;
                    for (int j = 0; j < J; j++) {
                        Cells2Edges[j] = c2e_tmp[j].ToArray();
                        var Cells2Edges_j = Cells2Edges[j];

                        for (int k = Cells2Edges_j.Length - 1; k >= 0; k--) {
                            int iEdge = Cells2Edges_j[k];
                            var Edge = EdgTmp[iEdge];

                            Debug.Assert((Edge.Cell1 == j) != (Edge.Cell2 == j));

                            Cells2Edges_j[k]++;
                            if (Edge.Cell2 == j) {
                                Cells2Edges_j[k] *= -1;
                            }
                        }
                    }

                    this.m_CellsToEdgesTmp = null;
                }
            }

            /// <summary>
            /// Normals for all affine-linear edges
            /// - 1st index: edge index
            /// - 2nd index: spatial direction
            /// </summary>
            public MultidimensionalArray NormalsForAffine {
                get;
                private set;
            }


            Caching.EdgeNormalsCacheLogic_CNsFace m_NormalsCache;

            /// <summary>
            /// Cached normals at nodes.
            /// </summary>
            public Caching.EdgeNormalsCacheLogic_CNsFace NormalsCache {
                get {
                    if(m_NormalsCache == null) {
                        m_NormalsCache = new Caching.EdgeNormalsCacheLogic_CNsFace(this.m_owner);
                    }
                    return m_NormalsCache;
                }
            }


            /// <summary>
            /// Computes the normals on face <paramref name="iFace"/> in the
            /// volume coordinate system of a given cell
            /// <paramref name="jCell"/> at the given <paramref name="Nodes"/>
            /// and writes the result into <paramref name="NormalsOut"/>
            /// </summary>
            /// <param name="jCell">Cell index</param>
            /// <param name="iFace">Face index</param>
            /// <param name="Nodes">
            /// Evaluation nodes
            /// <list type="bullet">
            ///   <item>1st index: Node index</item>
            ///   <item>2nd index: Spatial dimension</item>
            /// </list>
            /// </param>
            /// <param name="NormalsOut">
            /// <list type="bullet">
            ///     <item>1st index: cell index</item>
            ///     <item>2nd index: Node index</item>
            ///     <item>3rd index: Spatial dimension</item>
            /// </list>
            /// </param>
            /// <param name="QuadMetric">
            /// A by-product: the integral transformation metric.
            /// </param>
            /// <param name="Offset">
            /// An offset into the first entry of <paramref name="NormalsOut"/> and <paramref name="QuadMetric"/>.
            /// </param>
            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut, MultidimensionalArray QuadMetric, int Offset) {


                var Cj = m_owner.Cells.GetCell(jCell); //0 for 1st neighbour of edge e
                var Kref = m_owner.Cells.GetRefElement(jCell);
                int D = m_owner.SpatialDimension;
                double metric = Kref.FaceTrafoGramianSqrt[iFace];
                int NoOfNodes = Nodes.GetLength(0);

                Debug.Assert(object.ReferenceEquals(Nodes.RefElement, Kref));
                Debug.Assert(NormalsOut.Dimension == 3);
                Debug.Assert(NormalsOut.GetLength(2) == D);
                Debug.Assert(NormalsOut.GetLength(1) == NoOfNodes);

                if(QuadMetric != null) {
                    Debug.Assert(QuadMetric.Dimension == 2);
                    Debug.Assert(NormalsOut.GetLength(1) == NoOfNodes);
                }


                MultidimensionalArray RefNormals = Kref.FaceNormals;
                MultidimensionalArray AdjJacOut = m_owner.AdjungateJacobian.GetValue_Cell(Nodes, jCell, 1);

                double[] Acc = new double[D];

                for (int q = 0; q < NoOfNodes; q++) {// loop over the quadrature nodes
                    double normLen = 0.0;

                    // compute AdjJacOut[0,q,:,:]^T * RefNormal
                    for (int d1 = 0; d1 < D; d1++) {
                        double nX = 0;

                        for (int dS = 0; dS < D; dS++) {
                            nX += AdjJacOut[0, q, dS, d1] * RefNormals[iFace, dS];

                        }
                        normLen += nX * nX;

                        Acc[d1] = nX;
                    }

                    // normalize:
                    normLen = Math.Sqrt(normLen);
                    double sc = 1.0 / normLen;
                    for (int i = 0; i < D; i++) {
                        NormalsOut[Offset, q, i] = Acc[i]*sc;
                    }

                    if(QuadMetric != null)
                        QuadMetric[Offset, q] = normLen * metric;
                }
            }

            public void GetNormalsForCell(NodeSet Nodes, int jCell, int iFace, MultidimensionalArray NormalsOut) {
                this.GetNormalsForCell(Nodes, jCell, iFace, NormalsOut.ResizeShallow(1, Nodes.NoOfNodes, Nodes.SpatialDimension), null, 0);
            }


            /// <summary>
            /// computes <see cref="SqrtGramian"/>
            /// </summary>
            internal void ComputeSqrtGramian() {
                using (new FuncTrace()) {
                    var E = this.Count;
                    this.SqrtGramian = MultidimensionalArray.Create(E);
                    int D = this.m_owner.SpatialDimension;

                    if (D == 1) {
                        this.SqrtGramian.SetAll(1.0);
                    } else {

                        //var simplices = this.m_owner.Grid.RefElements;
                        
                        MultidimensionalArray Jac = MultidimensionalArray.Create(1, 1, D, D);
                        MultidimensionalArray JacSh = Jac.ExtractSubArrayShallow(0, 0, -1, -1);

                        int DMinusOne = Math.Max(D - 1, 1);
                        MultidimensionalArray JacTedge = MultidimensionalArray.Create(D, DMinusOne);
                        MultidimensionalArray JacTj = MultidimensionalArray.Create(D, D);
                        MultidimensionalArray JacFull = MultidimensionalArray.Create(D, DMinusOne);
                        MultidimensionalArray JacFullTranp = MultidimensionalArray.Create(DMinusOne, D);
                        MultidimensionalArray Gramian = MultidimensionalArray.Create(DMinusOne, DMinusOne);

                        for (int e = 0; e < E; e++) {
                            if ((this.Info[e] & EdgeInfo.EdgeIsAffineLinear) != 0) {
                                JacTedge.Clear();
                                //JacTedgeTransp.Clear();
                                JacTj.Clear();
                                //JacTedgeTransp.Clear();

                                Jac.Clear();

                                RefElement KrefEdge = this.GetRefElement(e);

                                var trafo = this.Edge2CellTrafos[this.Edge2CellTrafoIndex[e, 0]];
                                JacTedge.Acc(1.0, trafo.Matrix);
                                //JacTedge.Transpose(JacTedgeTransp);

                                int jCell1 = this.CellIndices[e, 0];
                                var Cell1 = m_owner.Cells.GetCell(jCell1);
                                var Kref = m_owner.Cells.GetRefElement(jCell1);
                                int iKref = m_owner.Cells.GetRefElementIndex(jCell1);
                                //Kref.JacobianOfTransformation(
                                //    EdgeCenters[iKref][this.FaceIndices[e, 0]],
                                //    Jac,
                                //    0, Cell1.Type, Cell1.TransformationParams);
                                m_owner.EvaluateJacobian(KrefEdge.Center.GetVolumeNodeSet(this.m_owner, this.Edge2CellTrafoIndex[e, 0]), jCell1, 1, Jac);
                                JacTj.Acc(1.0, JacSh);


                                JacFull.GEMM(1.0, JacTj, JacTedge, 0.0);
                                JacFull.TransposeTo(JacFullTranp);

                                Gramian.GEMM(1.0, JacFullTranp, JacFull, 0.0);
                                this.SqrtGramian[e] = Math.Sqrt(Gramian.Determinant());
                                if(double.IsInfinity(this.SqrtGramian[e]) || double.IsNaN(this.SqrtGramian[e]) || this.SqrtGramian[e] == 0)
                                    throw new ArithmeticException(string.Format("Illegal Gramian determint for some edge at cell {0}; value is {1}.", jCell1, this.SqrtGramian[e]));

                            } else {
                                this.SqrtGramian[e] = double.NaN;
                            }
                        }
                    }
                }
            }


            public double GetSqrtGramianForNonConformEdge(int iEdge, int _inOut) {

                double scaling = this.SqrtGramian[iEdge];

                RefElement KrefEdge = this.GetRefElement(iEdge);
                NodeSet KRefVert = KrefEdge.Vertices.GetVolumeNodeSet(this.m_owner, this.Edge2CellTrafoIndex[iEdge, _inOut]);

                int D = this.m_owner.SpatialDimension;
                double len = 0.0;
                for (int d = 0; d < D; d++) {
                    len += (KRefVert[1, d] - KRefVert[0, d]).Pow2();
                }
                len = len.Sqrt();
                if (len < 0 || len > 2)
                    throw new ArithmeticException();

                return scaling * (2.0 / len);

            }



            /// <summary>
            /// Edge-to-Cell - transformation index, i.e. index into <see cref="Edge2CellTrafos"/>;
            /// - 1st index: local edge index;
            /// - 2nd index: 0,1 first and second neighbor;
            /// </summary>
            public int[,] Edge2CellTrafoIndex {
                get;
                private set;
            }

            /// <summary>
            /// local cell indices of cells that belong to an edge;
            /// - 1st index: local edge index
            /// - 2nd index: 0,1 first and second neighbor;
            /// </summary>
            /// <remarks>
            /// Example: Let be <see cref="CellIndices"/>[i,0] = 123 and
            /// <see cref="CellIndices"/>[i,1] = 321; Then edge i is located
            /// on the intersection of (the closure of) cell 123 and cell 321;
            /// A negative cell index indicates that an edge is only subset of
            /// one cell (cells on the border of an domain). The negative cell
            /// index is always stored at the 2nd entry.
            /// </remarks>
            public int[,] CellIndices {
                get;
                private set;
            }

            /// <summary>
            /// Equal to <see cref="CellIndices"/>.
            /// </summary>
            public int[,] LogicalCellIndices {
                get {
                    return CellIndices;
                }
            }

            /// <summary>
            /// Face index, where the numbering of faces is defined by the reference element, see e.g. <see cref="RefElement.FaceToVertexIndices"/>.
            /// - 1st index: local edge index; 
            /// - 2nd index: 0 and 1 for first and second neighbor;
            /// </summary>
            /// <remarks>
            /// Example: let be <see cref="FaceIndices"/>[i,0] = 1 and
            /// <see cref="FaceIndices"/>[i,1] = 1 and
            /// <see cref="CellIndices"/>[i,0] = 123 and
            /// <see cref="CellIndices"/>[i,1] = 321; Then edge i is on the 1st
            /// face of cell 123 and also on the 1st edge of cell 321; If edge
            /// i lies on a border entry [i,1] is negative;
            /// </remarks>
            public byte[,] FaceIndices {
                get;
                private set;
            }

            /// <summary>
            /// Edge Tags
            /// index: local edge index;
            /// </summary>
            public byte[] EdgeTags {
                get;
                private set;
            }

            /// <summary>
            /// additional edge information
            /// </summary>
            public EdgeInfo[] Info {
                get;
                private set;
            }

            /// <summary>
            /// Not required for the classic grid; therefore, null.
            /// </summary>
            public int[][] EdgeToParts {
                get {
                    return null;
                }
            }

            /// <summary>
            /// fills <see cref="h_min_Edge"/>;
            /// </summary>
            internal void Initialize_h_Edge() {
                using (new FuncTrace()) {

                    int D = m_owner.SpatialDimension;
                    int E = Count;
                    double[] vk = new double[D];
                    double[] vi = new double[D];
                    h_min_Edge = MultidimensionalArray.Create(E);
                    h_max_Edge = MultidimensionalArray.Create(E);

                    
                    MultidimensionalArray verticesGlob = MultidimensionalArray.Create(1, 1, D);



                    for (int e = 0; e < E; e++) {

                        int jCell;
                        jCell = CellIndices[e, 0];
                        
                        RefElement Kref_edge = this.EdgeRefElements[GetRefElementIndex(e)];
                        NodeSet verticesLoc = Kref_edge.Vertices.GetVolumeNodeSet(this.m_owner, this.Edge2CellTrafoIndex[e, 0]);

                        int N = verticesLoc.NoOfNodes;
                        if(verticesGlob.GetLength(1) != N) {
                            verticesGlob = MultidimensionalArray.Create(1, N, D);
                        }
                        
                        m_owner.TransformLocal2Global(verticesLoc, jCell, 1, verticesGlob, 0);

                        double dist_min = double.MaxValue;
                        double dist_max = 0;
                        // für punkte vi der kante ....
                        for (int i = 0; i < N; i++) {
                            for (int d = 0; d < D; d++)
                                vi[d] = verticesGlob[0, i, d];

                            // vergleiche mit allen anderen punkten vk der kante ...
                            for (int k = i + 1; k < N; k++) {
                                for (int d = 0; d < D; d++)
                                    vk[d] = verticesGlob[0, k, d];

                                double dist_ik = 0;
                                for (int d = 0; d < D; d++) {
                                    double hh = vi[d] - vk[d];
                                    dist_ik += hh * hh;
                                }
                                if (dist_ik < dist_min)
                                    dist_min = dist_ik;
                                if (dist_ik > dist_max)
                                    dist_max = dist_ik;
                            }
                        }

                        h_min_Edge[e] = Math.Sqrt(dist_min);
                        h_max_Edge[e] = Math.Sqrt(dist_max);
                    }
                }

            }

            /// <summary>
            /// sets <see cref="NormalsForAffine"/>.
            /// </summary>
            internal void InitNormals() {
                int E = this.Count;
                int D = m_owner.SpatialDimension;
                var __Normals = MultidimensionalArray.Create(E, D);

                var Krefs = this.m_owner.Grid.RefElements;

                MultidimensionalArray[,] FaceCenters = new MultidimensionalArray[Krefs.Length, Krefs.Max(Kref => Kref.NoOfFaces)];

                for (int e = 0; e < E; e++) {
                    var Normal_e = __Normals.ExtractSubArrayShallow(new int[] { e, 0 }, new int[] { e, D - 1 });

                    if (this.IsEdgeAffineLinear(e)) {

                        int jCell = this.CellIndices[e, 0];
                        //var Kref = this.m_owner.Cells.GetRefElement(jCell);
                        int iKref = this.m_owner.Cells.GetRefElementIndex(jCell);
                        int iFace = this.FaceIndices[e, 0];

                        //this.GetNormals(e, 1, Kref.FaceCenters.ExtractSubArrayShallow(new int[] { iFace, 0 }, new int[] { iFace, D - 1 }), Normal_e.ResizeShallow(1, 1, D));
                        //todo

                        //if(FaceCenters[iKref, iFace] == null)
                        //    FaceCenters[iKref, iFace] = Krefs[iKref].FaceCenters.ExtractSubArrayShallow(new int[] { iFace, 0 }, new int[] { iFace, D - 1 });


                        this.GetNormalsForCell(Krefs[iKref].GetFaceCenter(iFace), jCell, iFace, Normal_e);


                    } else {
                        // 
                        Normal_e.SetAll(double.NaN);
                    }
                }

                this.NormalsForAffine = __Normals;
            }

            /// <summary>
            /// Returns the periodic transformation for edge <paramref name="iEdge"/>.
            /// </summary>
            /// <param name="iEdge"></param>
            /// <param name="InToOut">
            /// If true, the transformation from the in- to the out-cell, if false the other way around.
            /// </param>
            /// <returns></returns>
            public AffineTrafo GetPeriodicTrafo(int iEdge, bool InToOut) {
                if (EdgeTags[iEdge] < GridCommons.FIRST_PERIODIC_BC_TAG)
                    throw new ArgumentException("Edge is not periodic.");

                int iii = InToOut ? 0 : 1;
                int jCell = CellIndices[iEdge, iii];
                var Cell = this.m_owner.Cells.GetCell(jCell);
                int iFace = FaceIndices[iEdge, iii];
                var FaceTag = Cell.CellFaceTags.Single(cft => cft.FaceIndex == iFace);
                int iTrafo = ((int)EdgeTags[iEdge]) - GridCommons.FIRST_PERIODIC_BC_TAG;
                Debug.Assert(iTrafo == ((int)FaceTag.EdgeTag) - GridCommons.FIRST_PERIODIC_BC_TAG);

                AffineTrafo Trafo;
                if (FaceTag.PeriodicInverse)
                    Trafo = this.m_owner.Grid.InversePeriodicTrafo[iTrafo];
                else
                    Trafo = this.m_owner.Grid.PeriodicTrafo[iTrafo];

                return Trafo;
            }


        }
    }
}
