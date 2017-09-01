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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// Specialized version of <see cref="Quadrature"/> for cell boundary
    /// integrals. That is, it builds volume quad rules based on the edge rules
    /// of all edges of a given cell.
    /// </summary>
    /// <remarks>
    /// This functionality has been separated from the <see cref="Quadrature"/>
    /// because having different numbers of quadrature points per edge was hard
    /// to realize using the existing structures.
    /// </remarks>
    public abstract class CellBoundaryQuadrature<TQuadRule> : Quadrature<TQuadRule, CellMask>
        where TQuadRule : CellBoundaryQuadRule {

        /// <summary>
        /// constructor, for a quadrature rule which is already compiled (<paramref name="domNrule"/>)
        /// </summary>
        /// <param name="noOfIntegralsPerCell">tensor dimension of the integrand</param>
        /// <param name="context"></param>
        /// <param name="domNrule">quadrature rule and domain</param>
        /// <param name="cs">Physical or reference coordinate system?</param>
        public CellBoundaryQuadrature(
            int[] noOfIntegralsPerCell, Grid.Classic.GridData context, ICompositeQuadRule<TQuadRule> domNrule, CoordinateSystem cs = Quadrature.CoordinateSystem.Physical)
            : base(noOfIntegralsPerCell, context, domNrule, cs) //
        {
            foreach(IChunkRulePair<QuadRule> crp in domNrule) {
                NodeCoordinateSystem ncs = crp.Rule.Nodes.GetNodeCoordinateSystem(context);
                if(ncs != NodeCoordinateSystem.CellCoord) {
                    throw new ArgumentException("Illegal node set for cell boundary quadrature. Found some node set defined for: " + ncs.ToString() + ".");
                }
            }
        }


        /// <summary>
        /// Sweeps whether cell <paramref name="i0"/> is linear/nonlinear and how many cells of the same type are going to come after it.
        /// </summary>
        protected override void NextPart(out bool Linear, out int NoOfElm, int i0, int Len) {
            var cells = base.gridData.iGeomCells;

            bool bLinear = cells.IsCellAffineLinear(i0);
            Linear = bLinear;
            int L;
            for (L = 1; L < Len; L++) {
                if (cells.IsCellAffineLinear(i0 + L) != bLinear) {
                    NoOfElm = L;
                    return;
                }
            }

            NoOfElm = Len;
        }

        /// <summary>
        /// Some DEBUG checks on the current quadrature chunk.
        /// </summary>
        protected override void CheckQuadratureChunk(int j0, int Len, int iKRef) {
            // check1: Ref element index:
            int JE = Len + j0;
            for (int j = 0; j < JE; j++) {
                int _iKref = gridData.iGeomCells.GetRefElementIndex(j);
                if (iKRef != _iKref)
                    throw new Exception("Internal error: mismatch between reference element index for given quadrature rule and specific cell that should be integrated.");

            }

        }

        /// <summary>
        /// Index of the current rule's reference element into
        /// <see cref="BoSSS.Foundation.Grid.GridCommons.RefElements"/>.
        /// </summary>
        public override int CurrentRuleRefElementIndex {
            get {
                var Kref = base.CurrentRule.RefElement;
                Debug.Assert(Kref != null);
                int ret = Array.IndexOf(gridData.iGeomCells.RefElements, Kref);
                Debug.Assert(ret >= 0, "unable to identify reference element");


                if (ret < 0) {
                    throw new Exception();
                }


                return ret;
            }
        }

        /// <summary>
        /// integration metrics for linear mapped cells
        /// </summary>
        /// <returns>
        ///  - 1st index: cell 
        ///  - 2nd index: face 
        /// </returns>
        protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {

        
            int NoOFFaces = this.CurrentRule.RefElement.NoOfFaces;
            var ret = MultidimensionalArray.Create(L, NoOFFaces);
            var Cell2Edge = this.GridDat.iLogicalCells.Cells2Edges;
            var FaceIndices = this.GridDat.iGeomEdges.FaceIndices;
            var Scalings = this.GridDat.iGeomEdges.SqrtGramian;
            var faceGrm = this.GridDat.iGeomCells.GetRefElement(i0).FaceTrafoGramianSqrt;
                       
            for (int i = 0; i < L; i++) {
                int jCell = i + i0;

                int[] Cell2edge_j = Cell2Edge[jCell];
                Debug.Assert(this.GridDat.iGeomCells.IsCellAffineLinear(jCell) == true);

                for (int iF = 0; iF < NoOFFaces; iF++) {
                    ret[i, iF] = 0;
                }

                for (int e = Cell2edge_j.Length - 1; e >= 0; e--) {
                    int iEdge = Math.Abs(Cell2edge_j[e]) - 1;
                    int iFace = FaceIndices[iEdge, Cell2edge_j[e] > 0 ? 0 : 1];

                    ret[i, iFace] += Scalings[iEdge]/faceGrm[iFace];
                }
            }

#if DEBUG
            MultidimensionalArray NonLinScaling = GetScalingsForNonlinElements(i0, L);
            int[] NodesPerFace = CurrentRule.NumbersOfNodesPerFace;
            for (int l = 0; l < L; l++) {
                int k = 0;
                for (int iFace = 0; iFace < NoOFFaces; iFace++) {
                    for (int n = 0; n < NodesPerFace[iFace]; n++) {
                        Debug.Assert(Math.Abs(NonLinScaling[l, k] - ret[l, iFace]) <= 1.0e-8);
                        k++;
                    }
                }
            }
#endif

            return ret;
        }

        /// <summary>
        /// integration metrics for non-linear mapped cells
        /// </summary>
        /// <returns>
        ///  - 1st index: cell
        ///  - 2nd index: quadrature node
        /// </returns>
        protected override MultidimensionalArray GetScalingsForNonlinElements(int i0, int L) {

            // Nodes in the edge coordinate system
            NodeSet Nodes = CurrentRule.Nodes;
            int[] NodesPerFace = CurrentRule.NumbersOfNodesPerFace;

            // Cell-Reference element:
            // we can assume that cells {i0, ... , i0 + Len - 1}
            // map to the same reference element
            var grd = gridData;
            var RefElmCell = grd.iGeomCells.GetRefElement(i0);

            // bla
            int N = CurrentRule.NoOfNodes;
            int D = gridData.SpatialDimension;
            var Trafos = GridDat.iGeomEdges.Edge2CellTrafos;
            int NoOfFaces = RefElmCell.NoOfFaces;
            Debug.Assert(NoOfFaces == NodesPerFace.Length);

            // normals in reference space
            MultidimensionalArray RefNormals = RefElmCell.FaceNormals;

            // return value
            MultidimensionalArray _R = MultidimensionalArray.Create(L, N);



            NodeSet[] NodesForFace = new NodeSet[NoOfFaces];
            {
                int n = 0;
                for (int iFace = 0; iFace < NoOfFaces; iFace++) {
                    int Ne = NodesPerFace[iFace];
                    if(Ne > 0)
                        NodesForFace[iFace] = new NodeSet(RefElmCell, Nodes.ExtractSubArrayShallow(new int[] { n, 0 }, new int[] { n + Ne - 1, D - 1 }));
                    n += Ne;
                }
                Debug.Assert(n == N);
            }

            MultidimensionalArray Jac = grd.Jacobian.GetValue_Cell(Nodes, i0, L);
            MultidimensionalArray DetJ = grd.JacobianDeterminat.GetValue_Cell(Nodes, i0, L);
            
            // loop over cells:
            for (int j = 0; j < L; j++) {
                int jCell = j + i0;  // cell index

               
                // return value for Edge # jEdge
                var R = _R.ExtractSubArrayShallow(new int[] { j, 0 }, new int[] { j - 1, N - 1 });

                //// transform nodes to cell coordinates
                //int TrafoIdx = grd.Edges.TrafoIndex[jEdge, 0];
                //var Trafo = grd.Edges.EdgeToCellTrafos[TrafoIdx];
                //Trafo.Transform(NodesEdge, Nodes);
                ////RefElmCell.


                // the Jacobian and it's determinat
                //RefElmCell.JacobianDetTransformation(
                //    Nodes, DetJ, 0, grd.Grid.Cells[jCell].Type, grd.Grid.Cells[jCell].TransformationParams);
                //RefElmCell.JacobianOfTransformation(
                //    Nodes, Jac, 0, grd.Grid.Cells[jCell].Type, grd.Grid.Cells[jCell].TransformationParams);
                

                int n0 = 0;
                for (int RefFace = 0; RefFace < NoOfFaces; RefFace++) {
                    int Ne = NodesPerFace[RefFace];
                    if (Ne <= 0)
                        continue;

                    // compute normals in physical space
                    MultidimensionalArray CurvedNormals = MultidimensionalArray.Create(Ne, D);
                    grd.iGeomEdges.GetNormalsForCell(NodesForFace[RefFace], jCell, RefFace,  CurvedNormals);


                    for (int _n = 0; _n < Ne; _n++) { // loop over quadrature nodes
                        int n = n0 + _n;
                        double r = 0.0;
                        double sc = 0.0;
                        for (int d = 0; d < D; d++) {//calculate the scalar product between the normals in ref. and phyic. domain (pro quad. node)
                            double acc = 0.0;
                            for (int s = 0; s < D; s++) {
                                acc += RefNormals[RefFace, s] * Jac[j, n, d, s];
                            }
                            sc += acc * CurvedNormals[_n, d];
                        }

                        r += DetJ[j, n] * (1 / (sc));

                        R[n] = r;
                    }

                    n0 += Ne;
                }

            }

            return _R;
        }

        /// <summary>
        /// performs the quadrature
        /// </summary>
        protected override void DoQuadrature(TQuadRule quadRule, int j0, int _Bulksize) {
            var currentRuleWeights = quadRule.Weights;
            int[] NumbersOfNodesPerEdge = quadRule.NumbersOfNodesPerFace;
            int NoOfFaces = base.CurrentRule.RefElement.NoOfFaces;

            switch (CoordinateSystem) {

                case CoordinateSystem.Physical: {

                        int JE = j0 + _Bulksize;
                        int j = j0;
                        while (j < JE) {

                            // determine sub-chunk
                            bool Linear;
                            int Bulksize;
                            NextPart(out Linear, out Bulksize, j, _Bulksize - (j - j0));


                            if (Linear) {
                                // codepath for integration of linear elements
                                // +++++++++++++++++++++++++++++++++++++++++++



                                var scalings = GetScalingsForLinearElements(j0, Bulksize);

                                Debug.Assert(scalings.Dimension == 2);
                                Debug.Assert(scalings.IsContinious);
                                Debug.Assert(quadRule.Weights.IsContinious);
                                Debug.Assert(quadRule.Weights.Dimension == 1);


                                for (int jj = 0; jj < Bulksize; jj++) {
                                    int j_cell = jj + j;

                                    // loop over integral components:
                                    for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {
                                        int n = 0; // <--- node counter

                                        // loop over the edges of a cell
                                        for (int e = 0; e < NoOfFaces; e++) {
                                            //int iEdge = Math.Abs(jCell2Edge[j_cell, e]) - 1;
                                            double re = 0.0;

                                            // loop over the nodes of an edge
                                            for (int ne = 0; ne < NumbersOfNodesPerEdge[e]; ne++) {
                                                re += m_EvalResultsCollapsed[jj, n, m] * quadRule.Weights[n];
                                                n++;
                                            }

                                            re *= scalings[jj, e];

                                            //Debug.Assert(m_QuadResultsCollapsed[jj, e, m] == 0.0);
                                            m_QuadResultsCollapsed[jj, e, m] = re;
                                        }
                                    }
                                }


                            } else {
                                // codepath for integration of nonlinear/curved elements
                                // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                                var scalings = GetScalingsForNonlinElements(j0, Bulksize);

                                Debug.Assert(scalings.Dimension == 2);
                                Debug.Assert(scalings.IsContinious);
                                Debug.Assert(quadRule.Weights.IsContinious);
                                Debug.Assert(quadRule.Weights.Dimension == 1);


                                for (int jj = 0; jj < Bulksize; jj++) {
                                    int j_cell = jj + j;

                                    // loop over integral components:
                                    for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {
                                        int n = 0; // <--- node counter

                                        // loop over the edges of a cell
                                        for (int e = 0; e < NoOfFaces; e++) {
                                            //int iEdge = Math.Abs(jCell2Edge[j_cell, e]) - 1;
                                            double re = 0.0;

                                            // loop over the nodes of an edge
                                            for (int ne = 0; ne < NumbersOfNodesPerEdge[e]; ne++) {
                                                re += m_EvalResultsCollapsed[jj, n, m] * quadRule.Weights[n] * scalings[jj, n];
                                                n++;
                                            }

                                            //Debug.Assert(m_QuadResultsCollapsed[jj, e, m] == 0.0);
                                            m_QuadResultsCollapsed[jj, e, m] = re;
                                        }
                                    }
                                }

                            }


                            j += Bulksize;
                        }


                        break;
                    }

                case CoordinateSystem.Reference: {
                        Debug.Assert(quadRule.Weights.IsContinious);
                        Debug.Assert(quadRule.Weights.Dimension == 1);

                        for (int jj = 0; jj < _Bulksize; jj++) {

                            //int j_cell = jj + j;

                            // loop over integral components:
                            for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {
                                int n = 0; // <--- node counter

                                // loop over the edges of a cell
                                for (int e = 0; e < NoOfFaces; e++) {
                                    //int iEdge = Math.Abs(jCell2Edge[j_cell, e]) - 1;
                                    double re = 0.0;

                                    // loop over the nodes of an edge
                                    for (int ne = 0; ne < NumbersOfNodesPerEdge[e]; ne++) {
                                        re += m_EvalResultsCollapsed[jj, n, m] * quadRule.Weights[n];
                                        n++;
                                    }

                                    //Debug.Assert(m_QuadResultsCollapsed[jj, e, m] == 0.0);
                                    m_QuadResultsCollapsed[jj, e, m] = re;
                                }
                            }
                        }

                        break;
                    }

                default:
                    throw new NotImplementedException("not known: " + CoordinateSystem);
            }
        }

        /// <summary>
        /// 2nd phase of quadrature: allocation of memory for 
        /// the <see cref="Quadrature{S, T}.Evaluate"/>-method;
        /// Called whenever the node set or the number of cells per evaluation is changed;
        /// </summary>
        /// <param name="NoOfItems">number of edges or cells to integrate</param>
        /// <param name="rule">The quadrature rule applied to each item</param>
        protected override void AllocateBuffersInternal(int NoOfItems, NodeSet rule) {
            int NoOfNodes = rule.GetLength(0);
            if (m_EvalResults == null)
                m_EvalResults = new MultidimensionalArray(2 + IntegralCompDim.Length);
            m_EvalResults.Allocate(((new int[] { NoOfItems, NoOfNodes }).Concat(IntegralCompDim)).ToArray());

            if (m_QuadResults == null)
                m_QuadResults = new MultidimensionalArray(2 + IntegralCompDim.Length);
            m_QuadResults.Allocate(((new int[] { NoOfItems, base.CurrentRule.RefElement.NoOfFaces }).Concat(IntegralCompDim)).ToArray());

            m_EvalResultsCollapsed = m_EvalResults.ResizeShallow(
                (m_EvalResults.Lengths.Take(2).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());
            m_QuadResultsCollapsed = m_QuadResults.ResizeShallow(
                (m_QuadResults.Lengths.Take(2).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());


            if (m_AllocateBuffers != null)
                m_AllocateBuffers(NoOfItems, rule);
        }

        /// <summary>
        /// creates a cell quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="CellQuadrature"/>.
        /// </summary>
        static public CellBoundaryQuadrature<TQuadRule> GetQuadrature(
            int[] noOfIntegralsPerCell,
            Grid.Classic.GridData context,
            ICompositeQuadRule<TQuadRule> domNrule,
            Del_Evaluate _Evaluate,
            Del_SaveIntegrationResults _SaveIntegrationResults,
            Del_AllocateBuffers _AllocateBuffers = null,
            Del_QuadNodesChanged _PostLockNodes = null,
            CoordinateSystem cs = CoordinateSystem.Physical) {

            var ret = new CellBoundaryQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_Evaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }


        private class CellBoundaryQuadratureImpl : CellBoundaryQuadrature<TQuadRule> {

            public CellBoundaryQuadratureImpl(
                int[] noOfIntegralsPerCell, Grid.Classic.GridData context, ICompositeQuadRule<TQuadRule> domNrule, CoordinateSystem cs)
                : base(noOfIntegralsPerCell, context, domNrule, cs) {
            }

        }

    }







}
