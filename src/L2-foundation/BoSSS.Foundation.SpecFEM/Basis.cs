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
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform.LinAlg;
using MPI.Wrappers;
using ilPSP.Tracing;
using System.Collections;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.SpecFEM {
    
    /// <summary>
    /// A Spectral Element (high order continuous FEM) basis
    /// </summary>
    public class SpecFemBasis {


        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="p">
        /// The number of nodes per edge - this, in consequence, determines the SpecFEM-space, resp. the polynomial degree of the interpolation.
        /// </param>
        /// <param name="grdDat">
        /// The grid, upon which the SpecFEM-Basis is build.
        /// </param>
        /// <param name="__MassSolver">
        /// optional specification of the mass matrix solver.
        /// </param>
        public SpecFemBasis(Grid.Classic.GridData grdDat, int p, ilPSP.LinSolvers.ISparseSolver __MassSolver = null) {
            ConstructorCreateNodalPolynomials(grdDat, p);
            ConstructorCommon(p, __MassSolver);
        }

        private void ConstructorCreateNodalPolynomials(Grid.Classic.GridData grdDat, int p) {
            this.GridDat = grdDat;

            var Krefs = grdDat.Grid.RefElements;
            m_CellNodes = new NodeSet[Krefs.Length];
            m_NodalBasis = new PolynomialList[Krefs.Length];
            m_Type = new int[Krefs.Length][];
            m_EntityIndex = new int[Krefs.Length][];
            m_Nodal2Modal = new MultidimensionalArray[Krefs.Length];
            m_Modal2Nodal = new MultidimensionalArray[Krefs.Length];
            for (int iKref = 0; iKref < Krefs.Length; iKref++) {
                Krefs[iKref].SelectNodalPolynomials(p, out m_CellNodes[iKref], out m_NodalBasis[iKref], out m_Type[iKref], out m_EntityIndex[iKref], out m_Nodal2Modal[iKref], out m_Modal2Nodal[iKref]);
            }
        }

        private void ConstructorCommon(int p, ilPSP.LinSolvers.ISparseSolver __MassSolver) {
            using(new FuncTrace()) {
                Grid.Classic.GridData grdDat = this.GridDat;
                var Krefs = grdDat.Grid.RefElements;

                int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
                CellNode_To_Node = new int[J, m_CellNodes.Max(x => x.GetLength(0))];
                CellNode_To_Node.SetAll(int.MinValue);

                int D = this.GridDat.SpatialDimension;
                this.MaxFullDegree = (p - 1);

                this.UNodes = new MultidimensionalArray[Krefs.Length][][];

                for(int iKref = 0; iKref < Krefs.Length; iKref++) { // for all reference elements in the grid...

                    SortNodes(D, CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], out UNodes[iKref],
                        out m_NodesPerTypePerUnit // Bem.: noch nicht angepasst für mehrere Ref-Elemente im Gitter
                        );

                    if(this.GridDat.Grid.RefElements.Count() != 1)
                        throw new NotImplementedException();
                    if(this.GridDat.Edges.EdgeRefElements.Count() != 1)
                        throw new NotImplementedException();

                    {
                        if(D == 3) {
                            this.IdentifyCell(CellNodes[iKref], m_Type[iKref]);
                            this.IdentifyEdge(CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], iKref);
                            this.IdentifyVertice(CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], iKref);
                            // we still miss the co-edges (edges of edges)
                            throw new NotImplementedException("3D SpecFEM may come with GridOfTomorrow.");
                        } else if(D == 2) {
                            this.IdentifyCell(CellNodes[iKref], m_Type[iKref]);
                            this.IdentifyEdge(CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], iKref);
                            this.IdentifyVertice(CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], iKref);
                        } else if(D == 1) {
                            this.IdentifyCell(CellNodes[iKref], m_Type[iKref]);
                            this.IdentifyVertice(CellNodes[iKref], m_Type[iKref], m_EntityIndex[iKref], iKref);
                        } else {
                            throw new NotSupportedException("Unknown spatial dimension.");
                        }
                    }
                }

#if DEBUG
                {
                    // test if all local nodes are used.
                    
                    int K = this.NoOfLocalNodes;
                    int[] isNodeUsed = new int[K];
                    int B = CellNode_To_Node.GetLength(1);
                    int[] Ks = CellNodes.Select(NS => NS.NoOfNodes).ToArray();


                    for (int j = 0; j < J; j++) {
                        int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                        int Kj = Ks[iKref];

                        for (int b = 0; b < B; b++) {
                            int iNode = CellNode_To_Node[j, b];

                            if (b < Kj) {
                                isNodeUsed[iNode]++;
                            } else {
                                if (iNode >= 0)
                                    throw new ApplicationException();
                            }
                        }
                    }

                    for (int k = 0; k < K; k++) {
                        if (isNodeUsed[k] <= 0)
                            throw new ApplicationException("");
                    }

                }
#endif


                InitGlobalNodes();

                List<int>[] SendList = new List<int>[grdDat.MpiSize];
                List<int>[] InsList = new List<int>[grdDat.MpiSize];
                if(D == 3) {
                    throw new NotImplementedException("3D SpecFEM may come with GridOfTomorrow.");
                } else if(D == 2) {
                    this.BuildCommLists_forEdgeNodes(1, SendList, InsList);
                    this.BuildCommLists_forVertexNodes(2, SendList, InsList);
                } else if(D == 1) {
                    this.BuildCommLists_forVertexNodes(1, SendList, InsList);
                } else {
                    throw new NotSupportedException("Unknown spatial dimension.");
                }
                this.MPI_SendLists = SendList.Select(ll => ll != null ? ll.ToArray() : null).ToArray();
                this.MPI_InsertLists = InsList.Select(ll => ll != null ? ll.ToArray() : null).ToArray();


                ContainingDGBasis = new Foundation.Basis(grdDat, m_NodalBasis.Max(poly_list => poly_list.Max(poly => poly.AbsoluteDegree)));

                this.m_MassSolver = __MassSolver;
                if(m_MassSolver != null) {
                    m_MassSolver.DefineMatrix(this.MassMatrix);
                }

                InitMultiplicity();
                ForeignNodeMapping();
            }
        }


        void InitMultiplicity() {
            int Kg = this.NoOfLocalNodes;
            m_NodeMultiplicity = MultidimensionalArray.Create(Kg);
            int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
            var C2N = this.CellNode_To_Node;

            int[] _K = this.NodesPerCell;
            

            for (int j = 0; j < J; j++) {
                int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                int K = _K[iKref];
                for (int k  = 0; k < K; k++) {
                    m_NodeMultiplicity[C2N[j, k]] += 1.0;
                }
            }

            using(var trx = new Transceiver(this)) {
                trx.AccumulateGather(m_NodeMultiplicity);
            }
#if DEBUG
            // assert that the multiplicity of all cell nodes is 1
            int k0C = this.GetLocalOwnedNodesOffset(0);
            int k0E = this.GetNoOfLocalNodes(0);
            for(int k = 0; k < k0E; k++) {
                Debug.Assert(m_NodeMultiplicity[k + k0C] == 1.0);
            }

            int k1C = this.GetBorrowedNodesOffset(0);
            int k1E = this.GetNoOfBorrowedNodes(0);
            Debug.Assert(k1E == 0); // there should not be any cell node which is shared with other processors
#endif
        }

        MultidimensionalArray m_NodeMultiplicity;

        /// <summary>
        /// 
        /// </summary>
        public MultidimensionalArray NodeMultiplicity {
            get {
                return m_NodeMultiplicity;
            }
        }


        NodeSet[] m_CellNodes;

        /// <summary>
        /// Nodes for the reference cell in local coordinates<br/>
        ///  - 1st index: reference element <br/>
        ///  - 2nd index: node index <br/>
        ///  - 3rd index: spatial dimension
        /// </summary>
        public NodeSet[] CellNodes {
            get {
                return m_CellNodes;
            }
        }

        /// <summary>
        ///  - 1st index: reference element.<br/>
        ///  - 2nd index: type <br/>
        ///  - 3rd index: unit (face index, vertex index, ...)<br/>
        /// </summary>
        MultidimensionalArray[][][] UNodes;



        /// <summary>
        /// Number of Nodes in one cell;<br/>
        ///  - index: reference element <br/>
        /// </summary>
        public int[] NodesPerCell {
            get {
                return CellNodes.Select(nodes => nodes.GetLength(0)).ToArray();
            }
        }
        
        PolynomialList[] m_NodalBasis;

        /// <summary>
        /// the nodal basis polynomials.<br/>
        ///  - 1st index: reference element index<br/>
        ///  - 2nd index: node index.
        /// </summary>
        public PolynomialList[] NodalBasis {
            get {
                return m_NodalBasis.CloneAs();
            }
        }

        /// <summary>
        /// The Node type for each node of a cell, i.e. whether the node is a cell- (1D,2D,3D), face- (only 2D and 3D), co-face- (only 3D), or vertex node (1D,2D,3D);<br/>
        ///  - 1st index: reference element;<br/>
        ///  - 2nd index: cell node;
        /// </summary>
        int[][] m_Type;

        /// <summary>
        /// The entity index, see also <see cref="RefElement.GetNodeSet(int, out NodeSet, out int[], out int[], int[])"/>.
        /// If the k-th node is 
        /// * a cell node, the EntityIndex[k] is 0
        /// * a face/edge node, the EntityIndex[k] is the index of the face within the reference element
        /// * a co-face/co-edge node, the EntityIndex[k] is the index of the co-face within the reference element
        /// * a vertex node, the EntityIndex[k] is the index of the vertex within the reference element
        ///  - 1st index: reference element;<br/>
        ///  - 2nd index: cell node;
        /// </summary>
        int[][] m_EntityIndex;
        


        internal MultidimensionalArray[] m_Nodal2Modal;


        internal MultidimensionalArray[] m_Modal2Nodal;


        /// <summary>
        /// Returns the matrix to transform from
        /// a coordinate vector with respect to DG basis <see cref="Basis"/>
        /// int the nodal DG basis.
        /// </summary>
        public MsrMatrix GetModal2NodalOperator(Basis ModalBasis) {
            if(this.GridDat.Grid.RefElements.Length != 1)
                throw new NotImplementedException();
            
            
            int J = this.GridDat.Cells.NoOfCells;
            int N = this.NodesPerCell[0];
            int M = ModalBasis.Length;

            MsrMatrix R = new MsrMatrix(new Partitioning(J * N), new Partitioning(J * M));

            MultidimensionalArray M2N = this.m_Nodal2Modal[0];
            Debug.Assert(M2N.NoOfRows == N);
            Debug.Assert(M2N.NoOfCols <= M);
            M = Math.Min(M, M2N.NoOfCols);
            var Trafo = this.GridDat.ChefBasis.Scaling;
            

            for(int j = 0; j < J; j++) {
                int i0 = j * N + R.RowPartitioning.i0; // row offset into R for cell 'j'
                int j0 = j * ModalBasis.Length + R.ColPartition.i0; // column offset into R for cell 'j'
                double tr = Trafo[j];

                for(int n = 0; n < N; n++) { // loop over rows
                    for(int m = 0; m < M; m++) { // lop over columns
                        R[i0 + n, j0 + m] = tr * M2N[n, m];
                    }
                }
            }
            return R;
        }

        /// <summary>
        /// Returns the matrix to transform from
        /// a coordinate vector with respect to the nodal DG basis
        /// into the DG basis <see cref="Basis"/>.
        /// </summary>
        public MsrMatrix GetNodal2ModalOperator(Basis ModalBasis) {
            if(this.GridDat.Grid.RefElements.Length != 1)
                throw new NotImplementedException();
            
            int J = this.GridDat.Cells.NoOfCells;
            int N = this.NodesPerCell[0];
            int M = ModalBasis.Length;

            MsrMatrix R = new MsrMatrix(new Partitioning(J * M), new Partitioning(J * N));

            MultidimensionalArray N2M = this.m_Modal2Nodal[0];
            Debug.Assert(N2M.NoOfCols == N);
            Debug.Assert(N2M.NoOfRows <= M);
            M = Math.Min(M, N2M.NoOfRows);
            var Trafo = this.GridDat.ChefBasis.Scaling;


            for(int j = 0; j < J; j++) {
                int i0 = j * N + R.RowPartitioning.i0; // column offset into R for cell 'j'
                int j0 = j * ModalBasis.Length + R.ColPartition.i0; // row offset into R for cell 'j'
                double tr = Trafo[j];

                for(int m = 0; m < M; m++) { // loop over rows
                    for(int n = 0; n < N; n++) { // lop over columns
                        R[j0 + m, i0 + n] = (1.0/tr) * N2M[m, n];
                    }
                }
            }
            return R;
        }

        /// <summary>
        /// Scattering operator from global nodes to cell-wise nodes
        /// </summary>
        public MsrMatrix GetNodeScatterMatrix() {
            if(this.GridDat.Grid.RefElements.Length != 1)
                throw new NotImplementedException();

            int J = this.GridDat.Cells.NoOfCells;
            int N = this.NodesPerCell[0];
            int[] _K = this.NodesPerCell;
            var C2N = this.CellNode_To_Node;
                        
            MsrMatrix R = new MsrMatrix(new Partitioning(J * N), this.NodePartition);

            var CellData = this.GridDat.Cells;

            int i0Row = 0;
            for(int j = 0; j < J; j++) { // loop over cells...
                int iKref = CellData.GetRefElementIndex(j);
                
                //double[] NodalCoordinates = _NodalCoordinates[iKref];
                int K = _K[iKref];

                // collect coordinates for cell 'j':
                for(int k = 0; k < K; k++) {
                    int _c2n = C2N[j, k];
                    R[i0Row + k, _c2n] = 1.0; // NodalCoordinates[k] = m_Coordinates[_c2n];
                }


                i0Row += K;
            }

            return R;
        }

        /// <summary>
        /// An implicit description of the SEM space;
        /// </summary>
        public MsrMatrix GetNullSpaceMatrix() {
            int[][] dummy;
            return GetNullSpaceMatrix(out dummy);
        }
        

        /// <summary>
        /// An implicit description of the SEM space;
        /// </summary>
        public MsrMatrix GetNullSpaceMatrix(out int[][] Cells2Rows) {
            int[,] C2N = this.CellNode_To_Node;
            int K = this.NodePartition.LocalLength;
            int J = this.GridDat.Cells.NoOfCells;
            int N = this.NodesPerCell[0];
            Cells2Rows = J.ForLoop(j => new int[0]);

            
            // invert the cell_node-to-node mapping
            // ====================================
            
            List<Tuple<int,int>>[] N2C = K.ForLoop(k => new List<Tuple<int, int>>());

            Debug.Assert(J == C2N.GetLength(0));
            Debug.Assert(N == C2N.GetLength(1));
            for(int j = 0; j < J; j++) { // loop over cells...
                for(int n = 0; n < N; n++) { // loop over cell-nodes...
                    N2C[C2N[j, n]].Add(new Tuple<int, int>(j, n));
                }
            }

            // filter all global nodes which impose a restriction
            // ==================================================

            List<Tuple<int,int>>[] Restrictions = N2C.Where(list => list.Count > 1).ToArray();
            int NoOfRestrictions = Restrictions.Sum(list => list.Count - 1);

            // write Null-space matrix
            // =======================

            MsrMatrix MtxNS = new MsrMatrix(new Partitioning(NoOfRestrictions), new Partitioning(J * N));

            int irest = 0;
            foreach(var list in Restrictions) {
                for(int i = 0; i < list.Count - 1; i++) {
                    int j1 = list[i].Item1;  // first cell
                    int n1 = list[i].Item2;  // node in first cell
                    int j2 = list[i + 1].Item1; // second cell
                    int n2 = list[i + 1].Item2; // node in second cell

                    int jj1 = j1 * N + n1;
                    int jj2 = j2 * N + n2;

                    irest.AddToArray(ref Cells2Rows[j1]);
                    irest.AddToArray(ref Cells2Rows[j2]);

                    MtxNS[irest, jj1] = +1;
                    MtxNS[irest, jj2] = -1;
                    
                    irest++;
                }
            }

            return MtxNS;
        }


        /// <summary>
        /// the minimal DG Basis which contains this SpecFEM Basis
        /// </summary>
        public BoSSS.Foundation.Basis ContainingDGBasis {
            get;
            private set;
        }

        /// <summary>
        /// Nodes in global coordinates.
        /// </summary>
        public MultidimensionalArray GlobalNodes {
            get;
            private set;
        }

        static void SortNodes(int D, MultidimensionalArray Nodes, int[] Type, int[] Unit, out MultidimensionalArray[][] _UNodes, out int[] NodesPerTypePerUnit) {
            //int D = this.GridDat.SpatialDimension;
            int NoOfNodes = Nodes.GetLength(0);

          

            NodesPerTypePerUnit = new int[D + 1];
            for (int k = 0; k < NoOfNodes; k++) {
                if (Unit[k] == 0)
                    NodesPerTypePerUnit[Type[k]]++;
            }

            _UNodes = new MultidimensionalArray[D + 1][];
            var lUNodes = new List<MultidimensionalArray>[D + 1];

            int TYPE = -1;
            int U = -1; 
            int k0 = int.MinValue;
            for(int k = 0; k <= NoOfNodes; k++) {
                if (k >= NoOfNodes || Type[k] != TYPE || Unit[k] != U) {
                    if (k0 >= 0) {
                        int kE = k - 1;

                        var __uNodes = Nodes.ExtractSubArrayShallow(new int[] { k0, 0 }, new int[] { kE, D - 1 });
                        var _uNodes = MultidimensionalArray.Create(__uNodes.Lengths);
                        _uNodes.Set(__uNodes);

                        int _t = Type[kE];

                        if (lUNodes[_t] == null)
                            lUNodes[_t] = new List<MultidimensionalArray>();
                        lUNodes[_t].Add(_uNodes);
                    }
                    if (k < NoOfNodes) {
                        k0 = k;
                        TYPE = Type[k];
                        U = Unit[k];
                    }
                }
            }

            for (int h = 0; h < lUNodes.Length; h++) {
                if (lUNodes[h] != null)
                    _UNodes[h] = lUNodes[h].ToArray();
                else
                    _UNodes[h] = new MultidimensionalArray[0];
            }
        }

               

        void IdentifyCell(MultidimensionalArray CellNodes, int[] NodeType) {
            int NoOfNodes = GetNoOfLocalNodes(0);
            int Offset = 0;

            int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
            var C2N = CellNode_To_Node;

            int cnt = Offset;
            int NoOfNodesPerCell = NodeType.Length;
            for (int j = 0; j < J; j++) { // loop over cells

                for (int k = 0; k < NoOfNodesPerCell; k++) {
                    if (NodeType[k] == 0) { // loop over all cell nodes of cell j ...
                        C2N[j, k] = cnt; // asociate global node index 'cnt' with cell node 'k' of cell 'j'
                        cnt++;
                    }
                }
            }

            if (cnt != NoOfNodes - Offset)
                throw new ApplicationException("Error in Algorithm.");
        }

        int[] m_NodesPerTypePerUnit;
        
        void IdentifyEdge(MultidimensionalArray CellNodes, int[] NodeType, int[] EntityIndex, int iKrefFilter) {
            int TYPE = this.GridDat.SpatialDimension - 1; // type of edge nodes.
            int Offset_Owned = GetLocalOwnedNodesOffset(TYPE);
            int Offset_Foreign = GetBorrowedNodesOffset(TYPE);
            int D = this.GridDat.SpatialDimension;
            int L = m_NodesPerTypePerUnit[TYPE]; // Nodes per edge

            int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
            var C2N = CellNode_To_Node;

            var Cell2Edge = this.GridDat.Cells.Cells2Edges;
            var FaceIndices = this.GridDat.Edges.FaceIndices;
            var GidxExtCell = this.GridDat.Parallel.GlobalIndicesExternalCells;
            var Edges2Cell = this.GridDat.Edges.CellIndices;
            var EdgeTags = this.GridDat.Edges.EdgeTags;
            int j0 = this.GridDat.CellPartitioning.i0;
            

            
            int Y = this.GridDat.Edges.Edge2CellTrafos.Count;
            AffineTrafo[,] InterCellTrafos = new AffineTrafo[Y, Y];

            int NoOfOwnedEdges = this.GridDat.Edges.NoOfOwned;
            
            //int cnt = Offset;
            int NoOfNodesPerCell = NodeType.Length;
            for (int j = 0; j < J; j++) { // loop over cells...
                int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                if (iKref != iKrefFilter)
                    continue;

                var Cell2Edge_j = Cell2Edge[j];

                for (int k = 0; k < NoOfNodesPerCell; k++) { // loop over cell-nodes...
                    if (NodeType[k] == TYPE) {
                        // select EDGE-nodes

                        // find face indices, edge indices, in- and out-cell

                        int iFace = EntityIndex[k]; // face index
                        int iEdge = -1;
                        int iPeriodic = -1;
                        bool isInCell;
                        int jDefCell, jCell1, jCell2;
                        {
                            int InOrOut = -1;
                            isInCell = false;
                            for (int kk = Cell2Edge_j.Length - 1; kk >= 0; kk--) {
                                int q = Cell2Edge_j[kk];
                                int _iEdge = Math.Abs(q) - 1;
                                InOrOut = q > 0 ? 0 : 1;
                                isInCell = q > 0;
                                int _iFace = FaceIndices[_iEdge, InOrOut];
                                if (_iFace == iFace) {
                                    iEdge = _iEdge;
                                    break;
                                }
                            }
                            if (iEdge < 0)
                                throw new ApplicationException("error in alg.");
                            if (!this.GridDat.Edges.IsEdgeConformal(iEdge, InOrOut)) {
                                throw new NotSupportedException("currently, no hanging nodes for SpecFEM");
                            }

                            if (EdgeTags[iEdge] >= GridCommons.FIRST_PERIODIC_BC_TAG) {
                                iPeriodic = EdgeTags[iEdge] - GridCommons.FIRST_PERIODIC_BC_TAG;
                            }

                            jCell1 = Edges2Cell[iEdge,0];
                            jCell2 = Edges2Cell[iEdge,1];

                            // identify the defining cell:
                            // ======================================
                            if (jCell2 >= 0) {
                                // inner edge:
                                // the defining cell should be the one with smaller GlobalId
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                                long GidxCell1 = jCell1 < J ? j0 + jCell1 : GidxExtCell[jCell1 - J];
                                long GidxCell2 = jCell2 < J ? j0 + jCell2 : GidxExtCell[jCell2 - J];

                                if (GidxCell1 < GidxCell2) {
                                    jDefCell = jCell1;
                                } else {
                                    Debug.Assert(GidxCell2 < GidxCell1);
                                    jDefCell = jCell2;
                                }
                            } else {
                                // boundary edge (no out-cell):
                                // the defining cell is obviously the inner cell
                                // +++++++++++++++++++++++++++++++++++++++++++++
                                jDefCell = jCell1;
                            }
                        }

                        int i0;
                        if (iEdge < NoOfOwnedEdges) {
                            i0 = iEdge*L + Offset_Owned;
                        } else {
                            i0 = (iEdge - NoOfOwnedEdges)*L + Offset_Foreign;
                        }


                        if (j == jDefCell) {
                            // ++++++++++++++++++++++++++++++++++++++++++++++++
                            // cell j is responsible for defining the nodes
                            // ++++++++++++++++++++++++++++++++++++++++++++++++

                            //if (iPeriodic >= 0) {
                            //    Console.WriteLine("periodic node " + i0 + "\t(MPI rank " + this.GridDat.Partitioning.Rank  + ")");
                            //}

                            for (int i = 0; i < L; i++) {
                                C2N[j, k] = i0 + i;
                                Debug.Assert(iFace == EntityIndex[k]);
                                Debug.Assert(TYPE == NodeType[k]);
                                k++;
                            }
                            k--;
                        } else {
                            //var ict = GridDat.InterCellTransformations[GridDat.TrafoTo2ndCell[iEdge]];

                            int _0, _1;
                            if (isInCell) {
                                _0 = 1;
                                _1 = 0;
                            } else {
                                _0 = 0;
                                _1 = 1;
                            }
                            
                            int j_IN = Edges2Cell[iEdge, _0];
                            int iFace_IN = FaceIndices[iEdge, _0];

                            int j_OUT = Edges2Cell[iEdge, _1];
                            Debug.Assert(j_OUT == j);
                            int iFace_OT = FaceIndices[iEdge, _1];
                            Debug.Assert(iFace == iFace_OT);
                            
                            var NodesIn = UNodes[iKref][TYPE][iFace_IN];
                            var NodesOut = UNodes[iKref][TYPE][iFace_OT];

                            var ict = InterCellTrafos[iFace_IN, iFace_OT];
                            if(ict == null) {
                                var trafo_IN = this.GridDat.Edges.Edge2CellTrafos[this.GridDat.Edges.Edge2CellTrafoIndex[iEdge, _0]];
                                var trafo_OT = this.GridDat.Edges.Edge2CellTrafos[this.GridDat.Edges.Edge2CellTrafoIndex[iEdge, _1]];

                                var N1 = MultidimensionalArray.Create(D + 1, D);
                                var N2 = MultidimensionalArray.Create(D + 1, D);

                                var NE = MultidimensionalArray.Create(D, Math.Min(D - 1, 1));
                                if(D == 1) {
                                    // nop
                                } else if(D == 2) {
                                    NE[1, 0] = 1;
                                } else if(D == 3) {
                                    NE[1, 0] = 1;
                                    NE[2, 1] = 1;
                                } else
                                    throw new NotSupportedException("unknown spatial dimension");
                                trafo_IN.Transform(NE, N1.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { D - 1, D - 1 }));
                                trafo_OT.Transform(NE, N2.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { D - 1, D - 1 }));
                                ict = AffineTrafo.FromPoints(N1, N2);

#if DEBUG
                                {
                                    var Kref1 = this.GridDat.Cells.GetRefElement(jCell1);
                                    var Kref2 = this.GridDat.Cells.GetRefElement(jCell2);
                                    var FacePlane1 = Kref1.GetFacePlane(iFace_IN);
                                    var FacePlane2 = Kref1.GetFacePlane(iFace_OT);

                                    for(int k2 = 0; k2 < (D - 1); k2++) {
                                        double[] pt1 = N1.GetRow(k2);
                                        double[] pt2 = N2.GetRow(k2);
                                        double pt1Dist = FacePlane1.PointDistance(pt1).Abs();
                                        double pt2Dist = FacePlane2.PointDistance(pt2).Abs();
                                        Debug.Assert(pt1Dist < 1.0e-6);
                                        Debug.Assert(pt2Dist < 1.0e-6);
                                    }
                                }
#endif

                                //ict = trafo_Out * (trafo_IN.Invert());
                                InterCellTrafos[iFace_IN, iFace_OT] = ict;
                            }

                            

                            var R = KnotenZuordnung(NodesIn, ict, NodesOut);

                            //if (iPeriodic >= 0) {
                            //    Console.WriteLine("periodic node " + i0 + "\t(MPI rank " + this.GridDat.Partitioning.Rank + ") ext");
                            //}


                            for (int i = 0; i < L; i++) {
                                C2N[j_OUT, k] = i0 + R[i];
                                Debug.Assert(iFace == EntityIndex[k]);
                                Debug.Assert(TYPE == NodeType[k]);
                                k++;
                            }
                            k--;
                        }
                    }
                }
            }
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
        /// <returns>
        /// A permutation R, so that for all valid indices k, 
        /// <paramref name="NodesOut"/>[R[k]] == <paramref name="ict"/>(<paramref name="NodesIn"/>[k]).
        /// </returns>
        int[] KnotenZuordnung(MultidimensionalArray NodesIn, AffineTrafo ict, MultidimensionalArray NodesOut) {
            var NodesInTrf = MultidimensionalArray.Create(NodesIn.Lengths);
            ict.Transform(NodesIn, NodesInTrf);

            if(!ArrayTools.ListEquals(NodesInTrf.Lengths, NodesOut.Lengths))
                throw new ApplicationException("Error in algorithm.");

            int NoOfNodes = NodesIn.GetLength(0);
            var R = new int[NoOfNodes];
            for (int kOut = 0; kOut < NoOfNodes; kOut++) {
                var vOut = NodesOut.GetRow(kOut);

                bool bFound = false;

                for (int kIn = 0; kIn < NoOfNodes; kIn++) {
                    var vIn = NodesInTrf.GetRow(kIn);
                    if (GenericBlas.L2Dist(vOut, vIn) < 1.0e-10) {
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
        /// builds the cell-to-node (<see cref="CellNode_To_Node"/>) mapping for Vertex nodes
        /// </summary>
        void IdentifyVertice(MultidimensionalArray CellNodes, int[] NodeType, int[] EntityIndex, int iKrefFilter) {
            int TYPE = this.GridDat.SpatialDimension;
            int Offset_Owned = GetLocalOwnedNodesOffset(TYPE);
            int Offset_Foreign = GetBorrowedNodesOffset(TYPE);
            int Kmax = this.NoOfLocalNodes;
            int NoOfNodes = this.GetNoOfLocalNodes(TYPE);

            int L = m_NodesPerTypePerUnit[TYPE];

            int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
            var C2N = CellNode_To_Node;

            var Cell2Vertex = this.GridDat.Cells.CellVertices;
            int NoOfOwned = this.GridDat.Vertices.NoOfOwned;
            int Periodic0 = NoOfOwned + this.GridDat.Vertices.NoOfBorrowed;
            int PeriodicE = Periodic0 + this.GridDat.Vertices.NoOfPeriodicElim;


            /*
            int[] PeriodicIndexTrafo;
            {
                BitArray PeriodicFlag = null; //this.GridDat.Vertices.PeriodicFlag;
                if (this.GridDat.Vertices.PeriodicEliminatedPoints != null) {
                    throw new NotImplementedException("todo");
                    var PeriodicElim = this.GridDat.Vertices.PeriodicEliminatedPoints;
                    PeriodicIndexTrafo = new int[PeriodicFlag.Length];

                    int cnt = 0;
                    for (int k = 0; k < PeriodicIndexTrafo.Length; k++) {
                        if (PeriodicFlag[k]) {
                            PeriodicIndexTrafo[k] = PeriodicElim[k];
                        } else {
                            PeriodicIndexTrafo[k] = cnt;
                            cnt++;
                        }
                    }
                } else {
                    PeriodicIndexTrafo = null;
                }
            }
            */

            
            int NoOfNodesPerCell = NodeType.Length;
            for (int j = 0; j < J; j++) { // loop over locally updated cells
                int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                if (iKref != iKrefFilter)
                    continue;
                
                for (int k = 0; k < NoOfNodesPerCell; k++) {
                    if (NodeType[k] == TYPE) {
                        int iVtx = EntityIndex[k]; // index within cell
                        int k_Vertex = Cell2Vertex[j][iVtx]; // global index

                        if(k_Vertex >= Periodic0) {
                            Debug.Assert(k_Vertex < PeriodicE);
                            k_Vertex = GridDat.Vertices.PeriodicEliminatedPoints[k_Vertex];
                        }

                        //if (PeriodicIndexTrafo != null) {
                        //    k_Vertex = PeriodicIndexTrafo[k_Vertex];
                        //}
                        Debug.Assert(k_Vertex < NoOfNodes);

                        if (k_Vertex < NoOfOwned) {
                            // vertex is owned by this MPI process
                            C2N[j, k] = k_Vertex + Offset_Owned;
                        } else {
                            // some borrowed vertex
                            Debug.Assert(this.GridDat.MpiSize > 1, "");
                            C2N[j, k] = (k_Vertex - NoOfOwned) + Offset_Foreign;
                        }

                        if (C2N[j, k] >= Kmax)
                            throw new ArgumentException();
                    }
                }
            }
        }

        void BuildCommLists_forVertexNodes(int TYPE_vertexNodes, List<int>[] SendLists, List<int>[] InsertLists) {
            int myMpiRk = this.GridDat.MpiRank;
            int mpiSz = this.GridDat.MpiSize;
            
            var VtxSendLists = this.GridDat.Vertices.VertexSendLists;
            var VtxInsertLists = this.GridDat.Vertices.VertexInsertLists;
            int NoOfOwnedVtx =  this.GridDat.Vertices.NoOfOwned;
            int NoOfForeignVtx = this.GridDat.Vertices.NoOfBorrowed;
            int NoOfPureLocalVtx = this.GridDat.Vertices.NoOfPurelyLocal;

            int SEM_ForeignNodeOffset = this.GetBorrowedNodesOffset(TYPE_vertexNodes);
            int SEM_NoOfForeignNodes = this.GetNoOfBorrowedNodes(TYPE_vertexNodes);
            int SEM_OwnedNodeOffset = this.GetLocalOwnedNodesOffset(TYPE_vertexNodes);
            int SEM_NoOfOwnedNodes = this.GetNoOfOwnedNodes(TYPE_vertexNodes);

            Debug.Assert(SEM_NoOfOwnedNodes == NoOfOwnedVtx);

            // build send lists
            // ----------------

            for (int p = 0; p < mpiSz; p++) { // loop over MPI processes
                int[] VtxSL = VtxSendLists[p];
                if (VtxSL == null)
                    continue;
                if (SendLists[p] == null)
                    SendLists[p] = new List<int>();

                int L = VtxSL.Length;
                for (int l = 0; l < L; l++) { // loop over all vertices that should be sent to process 'p'...
                    int k_send = VtxSL[l];
                    Debug.Assert(k_send >= NoOfPureLocalVtx);
                    Debug.Assert(k_send < NoOfOwnedVtx);

                    int k_send_SEM = k_send + SEM_OwnedNodeOffset;
                    Debug.Assert(k_send_SEM >= SEM_OwnedNodeOffset);
                    Debug.Assert(k_send_SEM < SEM_OwnedNodeOffset + SEM_NoOfOwnedNodes);
                                        
                    SendLists[p].Add(k_send_SEM);
                }
            }

            // build receive lists
            // -------------------

            for (int p = 0; p < mpiSz; p++) {
                int[] VtxIL = VtxInsertLists[p];
                if (VtxIL == null)
                    continue;
                if (InsertLists[p] == null)
                    InsertLists[p] = new List<int>();

                int L = VtxIL.Length;
                for (int l = 0; l < L; l++) {
                    int k_ins = VtxIL[l];
                    Debug.Assert(k_ins >= NoOfOwnedVtx);
                    Debug.Assert(k_ins < NoOfOwnedVtx + NoOfForeignVtx);

                    int k_ins_SEM = (k_ins - NoOfOwnedVtx) + SEM_ForeignNodeOffset;
                    Debug.Assert(k_ins_SEM >= SEM_ForeignNodeOffset);
                    Debug.Assert(k_ins_SEM < SEM_ForeignNodeOffset + SEM_NoOfForeignNodes);
                    
                    InsertLists[p].Add(k_ins_SEM);
                }
            }
        }

        void BuildCommLists_forEdgeNodes(int TYPE_edgeNodes, List<int>[] SendLists, List<int>[] InsertLists) {
            int myMpiRk = this.GridDat.MpiRank;
            int mpiSz = this.GridDat.MpiSize;

            var EdgeSendLists = this.GridDat.Edges.EdgeSendLists;
            var EdgeInsertLists = this.GridDat.Edges.EdgeInsertLists;
            int NoOfOwnedEdges =  this.GridDat.Edges.NoOfOwned;
            int NoOfBorrowedEdges = this.GridDat.Edges.NoOfBorrowed;
            int NoOfPureLocaledges = this.GridDat.Edges.NoOfPurelyLocal;

            int SEM_ForeignNodeOffset = this.GetBorrowedNodesOffset(TYPE_edgeNodes);
            int SEM_NoOfForeignNodes = this.GetNoOfBorrowedNodes(TYPE_edgeNodes);
            int SEM_OwnedNodeOffset = this.GetLocalOwnedNodesOffset(TYPE_edgeNodes);
            int SEM_NoOfOwnedNodes = this.GetNoOfOwnedNodes(TYPE_edgeNodes);

            int KU = this.m_NodesPerTypePerUnit[TYPE_edgeNodes];

            // build send lists
            // ----------------

            for (int p = 0; p < mpiSz; p++) {
                int[] EdgeSL = EdgeSendLists[p];
                if (EdgeSL == null)
                    continue;
                if (SendLists[p] == null)
                    SendLists[p] = new List<int>();
                var SendList = SendLists[p];

                int L = EdgeSL.Length;
                for (int l = 0; l < L; l++) {
                    int k_send = EdgeSL[l];
                    Debug.Assert(k_send >= NoOfPureLocaledges);
                    Debug.Assert(k_send < NoOfOwnedEdges);


                    for (int ku = 0; ku < KU; ku++) {
                        int k_send_SEM = k_send * KU + SEM_OwnedNodeOffset + ku;
                        Debug.Assert(k_send_SEM >= SEM_OwnedNodeOffset);
                        Debug.Assert(k_send_SEM < SEM_OwnedNodeOffset + SEM_NoOfOwnedNodes);

                        SendList.Add(k_send_SEM);
                    }
                }
            }

            // build receive lists
            // -------------------

            for (int p = 0; p < mpiSz; p++) {
                int[] EdgeIL = EdgeInsertLists[p];
                if (EdgeIL == null)
                    continue;
                if (InsertLists[p] == null)
                    InsertLists[p] = new List<int>();
                var InsertList = InsertLists[p];

                int L = EdgeIL.Length;
                for (int l = 0; l < L; l++) {
                    int k_ins = EdgeIL[l];
                    Debug.Assert(k_ins >= NoOfOwnedEdges);
                    Debug.Assert(k_ins < NoOfOwnedEdges + NoOfBorrowedEdges);
                    
                    for (int ku = 0; ku < KU; ku++) {
                        int k_ins_SEM = (k_ins - NoOfOwnedEdges) * KU + SEM_ForeignNodeOffset + ku;
                        Debug.Assert(k_ins_SEM >= SEM_ForeignNodeOffset);
                        Debug.Assert(k_ins_SEM < SEM_ForeignNodeOffset + SEM_NoOfForeignNodes);
                        InsertLists[p].Add(k_ins_SEM);
                    }
                }
            }
        }

        /// <summary>
        /// compute global indices and processor ranks (of respective owner processors) for foreign nodes
        /// </summary>
        void ForeignNodeMapping() {
            int Kloc = this.NoOfLocalNodes;
            int Kown = this.NoOfLocalOwnedNodes;
            int myrank = this.GridDat.MpiRank;
            var AbusedArray = MultidimensionalArray.Create(Kloc);
            int k0 = this.NodePartition.i0;

            for (int k = 0; k < Kown; k++) { // loop over 'owned' nodes ...
                double pseudo_struct = 0;
                Debug.Assert(sizeof(double) == sizeof(int)*2);
                unsafe {
                    int* p = (int*)(&pseudo_struct);
                    p[0] = k + k0; // store the global index of the node on the first 4 bytes ...
                    p[1] = myrank; // ... and the rank on the last 4 bytes of the double
                }

                AbusedArray[k] = pseudo_struct;
            }

            using (var trx = new Transceiver(this)) {
                trx.Scatter(AbusedArray);
            }

            ForeignNodes_GlobalIndex = new int[Kloc - Kown];
            ForeignNodes_OwnerRank = new int[Kloc - Kown];

            for (int k = Kown; k < Kloc; k++) { // loop over 'foreign' nodes ...
                double pseudo_struct = AbusedArray[k];
                Debug.Assert(sizeof(double) == sizeof(int)*2);
                int kGlob, ownerRank;
                unsafe {
                    int* p = (int*)(&pseudo_struct);
                    kGlob = p[0];
                    ownerRank = p[1];
                }
                Debug.Assert(ownerRank == this.NodePartition.FindProcess(kGlob));
                ForeignNodes_GlobalIndex[k - Kown] = kGlob;
                ForeignNodes_OwnerRank[k - Kown] = ownerRank;
            }
        }

        /// <summary>
        /// foreach foreign node, the global node index
        /// </summary>
        public int[] ForeignNodes_GlobalIndex;

        /// <summary>
        /// foreach foreign node, the MPI-rank of the owner process
        /// </summary>
        public int[] ForeignNodes_OwnerRank;


        /// <summary>
        /// the associated grid
        /// </summary>
        public Grid.Classic.GridData GridDat {
            get;
            private set;
        }

        public int GetLocalOwnedNodesOffset(int type) {
            int N = 0;
            for (int i = 0; i < type; i++)
                N += GetNoOfOwnedNodes(i);
            return N;
        }


        public int GetBorrowedNodesOffset(int type) {
            int N = this.NoOfLocalOwnedNodes;
            for (int i = 0; i < type; i++)
                N += GetNoOfBorrowedNodes(i);
            return N;
        }


        /// <summary>
        /// For a global node index <paramref name="k"/>, this returns either the respective cell-, edge- or vertex-index
        /// which which the node is associated.
        /// </summary>
        /// <param name="k">
        /// a global node index
        /// </param>
        /// <param name="EntityIndex">
        /// On exit, 
        /// either a cell index, if <paramref name="k"/> is in the range of volume/cell nodes,
        /// or an edge index, if <paramref name="k"/> is in ther range of edge nodes,
        /// or a vertex index, if <paramref name="k"/> is in the range of vertex nodes.
        /// </param>
        /// <param name="_type">
        /// the node type
        /// </param>
        public void GetEntityIndex(int k, out int EntityIndex, out int _type) {
            if(k < 0)
                throw new IndexOutOfRangeException("Global node index out of range.");

            int MaxType = GridDat.SpatialDimension;
            for(int type = 0; type <= MaxType; type++) {
                int Offset = this.GetLocalOwnedNodesOffset(type);
                int NoOf = this.GetNoOfLocalNodes(type);

                int R = k - Offset;
                if(R >= 0 && R < NoOf) {
                    EntityIndex = R;
                    _type = type;
                    return;
                }
            }

            throw new IndexOutOfRangeException("Global node index out of range.");
        }


        /// <summary>
        /// Number of local nodes (on current MPI process), for different node types (<paramref name="type"/>);
        /// includes shared/foreign/owned but NOT periodic/external.
        /// </summary>
        /// <param name="type">
        /// node type index:
        /// <list type="bullet">
        /// <item>in 1D: 0 for volume/cell nodes, 1 for vertex nodes/edges (in 1D, the edges of the grid are 0-dimensional nodes)</item>
        /// <item>in 2D: 0 for volume/cell nodes, 1 for edge nodes and 2 for vertex nodes. </item>
        /// <item>in 3D: 0 for volume/cell nodes, 1 for edge nodes (located on faces), 2 for co-edges (edges of faces), 3 for vertex nodes.</item>
        /// </list>
        /// </param>
        /// <returns></returns>
        public int GetNoOfLocalNodes(int type) {
            int D = GridDat.SpatialDimension;
            if (type < 0 || type > D)
                throw new ArgumentOutOfRangeException("Unknown node type index;");

            if (type == 0) {
                // cell nodes
                // ++++++++++
                return this.GridDat.Cells.NoOfLocalUpdatedCells*m_NodesPerTypePerUnit[type];
            }

            if (type == D) {
                // vertex nodes
                // ++++++++++++

                // This assertion is not true, when running 1D and parallel
                //if ((D == 1) && (this.GridDat.Edges.Count != this.GridDat.Vertices.Count))
                //    throw new ApplicationException("Error in grid data structure.");

                Debug.Assert(m_NodesPerTypePerUnit[type] == 1); // vertex nodes should always have only one node per Vertex.
                int Periodic = this.GridDat.Vertices.NoOfPeriodicElim;
                
                return (this.GridDat.Vertices.NoOfNodes4LocallyUpdatedCells - Periodic)*m_NodesPerTypePerUnit[type];
            }

            if (type == 1) {
                // edge nodes
                // ++++++++++
                return this.GridDat.Edges.Count*m_NodesPerTypePerUnit[type];
            }

            if (D == 3 && type == 2)
                throw new NotImplementedException("3D SpecFEM may come somewhen...");

            throw new ApplicationException("Error in algorithm.");
        }

        /// <summary>
        /// local Number Of Nodes on this MPI process
        /// (includes shared/foreign/owned but NOT periodic).
        /// </summary>
        public int NoOfLocalNodes {
            get {
                int N = 0, D = this.GridDat.SpatialDimension;
                for (int d = 0; d <= D; d++) {
                    N += GetNoOfLocalNodes(d);
                }
                return N;
            }
        }

        /// <summary>
        /// local Number Of Nodes owned by MPI process
        /// </summary>
        public int NoOfLocalOwnedNodes {
            get {
                int N = 0, D = this.GridDat.SpatialDimension;
                for (int d = 0; d <= D; d++) { // loop over vertex types
                    N += GetNoOfOwnedNodes(d);
                }
                return N;
            }
        }

        public int GetNoOfOwnedNodes(int type) {
            int D = GridDat.SpatialDimension;
            if (type < 0 || type > D)
                throw new ArgumentOutOfRangeException("Unknown node type index;");

            if (type == 0)
                // all cell nodes are locally owned
                return this.GetNoOfLocalNodes(type);

            if (type == D) {
                // vertices

                Debug.Assert(m_NodesPerTypePerUnit[type] == 1); // vertex nodes should always have only one node per Vertex.
                
                return this.GridDat.Vertices.NoOfOwned*m_NodesPerTypePerUnit[type];
            }
            
            if( type == 1) {
                // edges
                return this.GridDat.Edges.NoOfOwned*this.m_NodesPerTypePerUnit[type];
            }

            if (D == 3 && type == 2)
                // co-edges
                throw new NotImplementedException("3D SpecFEM may come with GridOfTomorrow.");

            throw new ApplicationException("Error in algorithm.");
        }

        public int GetNoOfBorrowedNodes(int type) {
            int foreign = (GetNoOfLocalNodes(type) - GetNoOfOwnedNodes(type));

            

            Debug.Assert(this.GridDat.MpiSize > 1 || foreign == 0); // when running with 1 process, there should not be 
            //                                                      any foreign node
            return foreign;
        }

        Partitioning m_NodePartition = null;

        /// <summary>
        /// partition of nodes among MPI processes
        /// </summary>
        public Partitioning NodePartition {
            get {
                if (m_NodePartition == null) {
                    m_NodePartition = new Partitioning(this.NoOfLocalOwnedNodes);
                }
                return m_NodePartition;
            }
        }
        
        /// <summary>
        /// Global Number of Nodes (over all MPI processes), i.e. dimension of the SEM vector space.
        /// </summary>
        public int GlobalNoOfNodes {
            get {
                return NodePartition.TotalLength;
            }
        }



        /// <summary>
        /// Transformation from cell-node index to local node index:<br/>
        ///  - 1st index: cell index<br/>
        ///  - 2nd index: node index within cell
        /// </summary>
        public int[,] CellNode_To_Node;

        /// <summary>
        /// For each MPI processor rank <em>p</em>, a list of indices which determines the DOFs that must be send
        /// to rank <em>p</em> during the scatter-operation <see cref="Transceiver.Scatter(MultidimensionalArray)"/>. <br/>
        ///  - 1st index: MPI rank <em>p</em>, if an entry is null, no data is send to the corresponding processor.<br/>
        ///  - 2nd index: enumeration.<br/>
        /// content: indices of locally owned DOFs.
        /// </summary>
        public int[][] MPI_SendLists;

        /// <summary>
        /// For each MPI processor rank <em>p</em>, a list of indices which determines 
        /// where the DOFs which this rank receices from rank  <em>p</em> should be inserted
        /// during the scatter-operation <see cref="Transceiver.Scatter(MultidimensionalArray)"/>. <br/>
        ///  - 1st index: MPI rank <em>p</em>, if an entry is null, no data is received from processor <em>p</em>.<br/>
        ///  - 2nd index: enumeration.<br/>
        /// content: indices of borrowed DOFs.
        /// </summary>
        public int[][] MPI_InsertLists;


        /// <summary>
        /// Computes the coordinates of the global nodes in physical coordinates.
        /// </summary>
        void InitGlobalNodes() {
            using(new FuncTrace()) {
                int D = GridDat.SpatialDimension;
                this.GlobalNodes = MultidimensionalArray.Create(this.NoOfLocalNodes, D);

                int vecsize = -100;
                int j = 0;
                int J = GridDat.Cells.NoOfLocalUpdatedCells;
                int _K = this.CellNodes[0].GetLength(0);

                MultidimensionalArray GlobalNodesPerCell = MultidimensionalArray.Create(100, _K, D);

                while(j < J) {  // loop over cells

                    // find appropriate vector size
                    int maxSize = Math.Min(100, J - j);
                    Debug.Assert(maxSize > 0);
                    vecsize = GridDat.Cells.GetNoOfSimilarConsecutiveCells(CellInfo.RefElementIndex_Mask, j, maxSize);


                    int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                    int K = this.CellNodes[iKref].GetLength(0); // number of nodes per cell

                    // transform cell nodes to global coordinates
                    if(GlobalNodesPerCell.GetLength(0) < vecsize || GlobalNodesPerCell.GetLength(1) != K)
                        GlobalNodesPerCell.Allocate(vecsize, K, D);
                    GridDat.TransformLocal2Global(this.CellNodes[iKref], j, vecsize, GlobalNodesPerCell, 0);
                    
                    for(int i = 0; i < vecsize; i++) {
                        int jCell = i + j;

                        for(int k = 0; k < K; k++) {
                            int iTarg = CellNode_To_Node[jCell, k];

                            for(int d = 0; d < D; d++)
                                this.GlobalNodes[iTarg, d] = GlobalNodesPerCell[i, k, d];
                        }
                    }
                    j += vecsize;
                }
            }
        }

        MsrMatrix m_MassMatrix;


        /// <summary>
        /// the SpecFEM mass matrix;
        /// </summary>
        public MsrMatrix MassMatrix {
            get {
                if (m_MassMatrix == null)
                    m_MassMatrix = ComputeMassMatrix();
                return m_MassMatrix;
            }
        }


        MultidimensionalArray[] ComputeCellLocalMassMatrix() {
            var ret = new MultidimensionalArray[this.GridDat.Grid.RefElements.Length];

            for (int iKref = 0; iKref < ret.Length; iKref++) {
                int K = NodesPerCell[iKref];
                var Mass = MultidimensionalArray.Create(K, K);

                var MR = m_Modal2Nodal[iKref];
                var ML = MR.Transpose();

                Mass.GEMM(1.0, ML, MR, 0.0);

                ret[iKref] = Mass;
            }

            return ret;
        }


        /// <summary>
        /// computes the Mass-Matrix of this SpecFEM-Basis.
        /// </summary>
        /// <param name="cm">
        /// optional restriction onto a subgrid
        /// </param>
        /// <returns></returns>
        internal MsrMatrix ComputeMassMatrix(CellMask cm = null) {
            var MassMatrix = new MsrMatrix(this.NodePartition);
            MassMatrix.AssumeSymmetric = true;

            if (cm == null) {
                cm = CellMask.GetFullMask(this.GridDat);
            }

            int[] _K = NodesPerCell;
            int J = GridDat.Cells.NoOfLocalUpdatedCells;
            MultidimensionalArray[] _CellLocalMass = ComputeCellLocalMassMatrix();
            var scalings = GridDat.Cells.JacobiDet;

            var C2N = this.CellNode_To_Node;

            int k0 = this.NodePartition.i0;
            int Kown = this.NoOfLocalOwnedNodes;


            var SendData = new Dictionary<int, List<Tuple<int, int, double>>>();


            foreach (Chunk cnk in cm) {
                int j0 = cnk.i0;
                int jE = cnk.JE;
                for (int j = j0; j < jE; j++) { // loop over cells

                    int iKref = this.GridDat.Cells.GetRefElementIndex(j);
                    var CellLocalMass = _CellLocalMass[iKref];
                    int K = _K[iKref];
                    Debug.Assert(CellLocalMass.NoOfCols == K);
                    Debug.Assert(CellLocalMass.NoOfRows == K);

                    double tr = scalings[j];

                    if (!this.GridDat.Cells.IsCellAffineLinear(j))
                        throw new NotImplementedException("Curved elements currently not supported.");
                    
                    //
                    
                    for (int l = 0; l < K; l++) {
                        int rowLoc = C2N[j, l];
                        int rowGlob;
                        List<Tuple<int, int, double>> SendData_p = null;
                        if (rowLoc < Kown) {
                            rowGlob = rowLoc + k0;
                        } else {
                            int ownRank = this.ForeignNodes_OwnerRank[rowLoc - Kown];
                            if (!SendData.TryGetValue(ownRank, out SendData_p)) {
                                SendData_p = new List<Tuple<int, int, double>>();
                                SendData.Add(ownRank, SendData_p);
                            }
                            rowGlob = this.ForeignNodes_GlobalIndex[rowLoc - Kown];
                        }

                        for (int k = 0; k < K; k++) {
                            int colLoc = C2N[j, k];
                            int colGlob;
                            if (colLoc < Kown) {
                                colGlob = colLoc + k0;
                            } else {
                                colGlob = this.ForeignNodes_GlobalIndex[colLoc - Kown];
                            }

                            if (SendData_p == null)
                                MassMatrix[rowGlob, colGlob] += CellLocalMass[l, k]*tr;
                            else
                                SendData_p.Add(new Tuple<int, int, double>(rowGlob, colGlob, CellLocalMass[l, k]*tr));
                        }
                    }
                }
            }

            var ReceiveData = SerialisationMessenger.ExchangeData(SendData, csMPI.Raw._COMM.WORLD);
            foreach (var kv in ReceiveData) {
                foreach (var ee in kv.Value) {
                    MassMatrix[ee.Item1, ee.Item2] += ee.Item3;
                }
            }

            return MassMatrix;
        }

        ilPSP.LinSolvers.ISparseSolver m_MassSolver;

        /// <summary>
        /// a linear solver for the omnipresent problem <br/>
        /// <see cref="MassMatrix"/>*x = b. <br/>
        /// </summary>
        public ilPSP.LinSolvers.ISparseSolver MassSolver {
            get {
                if (m_MassSolver == null) {
                    var solver = new ilPSP.LinSolvers.monkey.CG();
                    solver.MatrixType = ilPSP.LinSolvers.monkey.MatrixType.CCBCSR;
                    solver.DevType = ilPSP.LinSolvers.monkey.DeviceType.CPU;
                    solver.Tolerance = 1.0e-12;
                    solver.DefineMatrix(this.MassMatrix);
                    m_MassSolver = solver;
                }

                return m_MassSolver;
            }
        }

        /// <summary>
        /// test code for mass matrix computation
        /// </summary>
        public void ComputeMassMatrixChk() {
            var _polys = this.m_NodalBasis;
            int K = _polys.Length;
            int J = GridDat.Cells.NoOfCells;

            MultidimensionalArray[] CellMass = new MultidimensionalArray[J];
            MsrMatrix MassChk = new MsrMatrix(NoOfLocalNodes, 1);


            CellQuadrature.GetQuadrature(new int[] { K, K },
                GridDat,
                (new CellQuadratureScheme()).Compile(GridDat, this.ContainingDGBasis.Degree*2),
                delegate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    NodeSet Nodes = QR.Nodes;
                    int _NoOfNodes = Nodes.NoOfNodes;
                    var PolyAtNode = MultidimensionalArray.Create(K, _NoOfNodes);

                    int iKref = this.GridDat.Cells.GetRefElementIndex(i0);
                    var polys = _polys[iKref];
                    for(int k = 0; k < K; k++) {
                        polys[k].Evaluate(PolyAtNode.ExtractSubArrayShallow(k, -1), Nodes);
                    }


                    EvalResult.Multiply(1.0, PolyAtNode, PolyAtNode, 0.0, "jnkl", "kn", "ln");

                    double errSum = 0;
                    for (int i = 0; i < Length; i++) {
                        Debug.Assert(this.GridDat.Cells.GetRefElementIndex(i0 + i) == iKref);
                        for (int n = 0; n < _NoOfNodes; n++) {
                            for (int k = 0; k < K; k++) {
                                for (int l = 0; l < K; l++) {
                                    double soll = PolyAtNode[k, n]*PolyAtNode[l, n];
                                    errSum += Math.Abs(soll - EvalResult[i, n, k, l]);
                                }
                            }
                        }
                    }
                    Console.WriteLine("errsum = " + errSum);
                },
                delegate(int i0, int Length, MultidimensionalArray ResultsOfIntegration) { // SaveIntegrationResults
                    for (int i = 0; i < Length; i++) {
                        int jCell = i + i0;
                        CellMass[jCell] = MultidimensionalArray.Create(K, K);
                        CellMass[jCell].Set(ResultsOfIntegration.ExtractSubArrayShallow(i, -1, -1));
                    }
                }).Execute();

            var Andere = this.ComputeCellLocalMassMatrix();

            var C2N = this.CellNode_To_Node;
            for (int j = 0; j < J; j++) {
                var M = CellMass[j];

                //
                for (int l = 0; l < K; l++) {
                    int rowGlob = C2N[j, l];
                    for (int k = 0; k < K; k++) {
                        int colGlob = C2N[j, k];
                        MassChk[rowGlob, colGlob] += M[l, k];
                    }
                }
            }

            MassChk.Acc(-1.0, this.MassMatrix);
            double test = MassChk.InfNorm();
            if (test > 1.0e-10)
                throw new Exception();
        }


        /// <summary>
        /// maximum polynomial degree of the nodal basis
        /// </summary>
        public int MaxDegree {
            get {
                return this.m_NodalBasis.Max(pList => pList.Max(p => p.AbsoluteDegree));
            }
        }

        /// <summary>
        /// the maximum degree that is available in all directions
        /// </summary>
        public int MaxFullDegree {
            get;
            private set;
        }

        Basis m_MaxFullDGBasis;

        /// <summary>
        /// a DG basis of degree <see cref="MaxFullDegree"/>.
        /// </summary>
        public Basis MaxFullDGBasis {
            get {
                if (m_MaxFullDGBasis == null)
                    m_MaxFullDGBasis = new Basis(this.GridDat, this.MaxFullDegree);
                return m_MaxFullDGBasis;
            }
        }
        
    }
}
