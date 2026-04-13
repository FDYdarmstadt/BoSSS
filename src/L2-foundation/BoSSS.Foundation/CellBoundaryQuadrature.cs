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
            int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<TQuadRule> domNrule, CoordinateSystem cs = Quadrature.CoordinateSystem.Physical)
            : base(noOfIntegralsPerCell, context, domNrule, cs) //
        {
            foreach(IChunkRulePair<QuadRule> crp in domNrule) {
                NodeCoordinateSystem ncs = crp.Rule.Nodes.GetNodeCoordinateSystem(context);
                if(ncs != NodeCoordinateSystem.CellCoord) {
                    throw new ArgumentException("Illegal node set for cell boundary quadrature. Found some node set defined for: " + ncs.ToString() + ".");
                }
            }
        }

        protected override IIntegrationMetric GetDefaultIntegrationMetric() {
            return new CellBoundaryIntegrationMetric();
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
        /// <see cref="BoSSS.Foundation.Grid.IGeometricalCellsData.RefElements"/>.
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



                                var scalings = IntegrationMetric.GetScalingsForLinearElements(this.gridData, quadRule, j0, Bulksize);

                                Debug.Assert(scalings.Dimension == 2);
                                Debug.Assert(scalings.IsContinuous);
                                Debug.Assert(quadRule.Weights.IsContinuous);
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

                                var scalings = IntegrationMetric.GetScalingsForNonlinElements(this.gridData, quadRule, j0, Bulksize);

                                Debug.Assert(scalings.Dimension == 2);
                                Debug.Assert(scalings.IsContinuous);
                                Debug.Assert(quadRule.Weights.IsContinuous);
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
                        Debug.Assert(quadRule.Weights.IsContinuous);
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
        /// <param name="ThreadRank"></param>
        /// <param name="iThread"></param>
        protected override void AllocateBuffersInternal(int NoOfItems, NodeSet rule, int iThread, int ThreadRank) {
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
                m_AllocateBuffers(NoOfItems, rule, iThread, ThreadRank);
        }

        /// <summary>
        /// creates a cell quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="CellQuadrature"/>.
        /// </summary>
        static public CellBoundaryQuadrature<TQuadRule> GetQuadrature(
            int[] noOfIntegralsPerCell,
            IGridData context,
            ICompositeQuadRule<TQuadRule> domNrule,
            Del_Evaluate _Evaluate,
            Del_SaveIntegrationResults _SaveIntegrationResults,
            //Del_AllocateBuffers _AllocateBuffers = null,
            //Del_QuadNodesChanged _PostLockNodes = null,
            CoordinateSystem cs = CoordinateSystem.Physical) {

            var ret = new CellBoundaryQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_Evaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                //m_AllocateBuffers = _AllocateBuffers,
                //m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }


        private class CellBoundaryQuadratureImpl : CellBoundaryQuadrature<TQuadRule> {

            public CellBoundaryQuadratureImpl(
                int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<TQuadRule> domNrule, CoordinateSystem cs)
                : base(noOfIntegralsPerCell, context, domNrule, cs) {
            }

            public override Quadrature<TQuadRule, CellMask> CloneForThreadParallelization(int iThread, int NumThreads) {
                return new CellBoundaryQuadratureImpl(
                    this.IntegralCompDim, this.GridDat, this.m_compositeRule, this.CoordinateSystem) {
                    m_AllocateBuffers = this.m_AllocateBuffers,
                    m_SaveIntegrationResults = this.m_SaveIntegrationResults,
                    m_quadNodesChanged = this.m_quadNodesChanged,
                    m_Evaluate = this.m_Evaluate,
                    m_ExEvaluate = this.m_ExEvaluate,
                };
            }
        
        }

    }







}
