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
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using System.Diagnostics;
using ilPSP;

namespace BoSSS.Foundation.Quadrature {


    /// <summary>
    /// A variant of the edge quadrature that uses different 
    /// quadrature nodes and weights on both sides of an edge (different quadrature for 'in'- and 'out'-cell).  
    /// </summary>
    public abstract class DoubleEdgeQuadrature : Quadrature<DoubleEdgeQuadRule,EdgeMask> {

        /// <summary>
        /// constructor, for a quadrature rule which is already compiled (<paramref name="domNrule"/>)
        /// </summary>
        /// <param name="noOfIntegralsPerCell">tensor dimension of the integrand</param>
        /// <param name="context"></param>
        /// <param name="domNrule">quadrature rule and domain</param>
        /// <param name="cs">Physical or reference coordinate system?</param>
        public DoubleEdgeQuadrature(
            int[] noOfIntegralsPerCell, Grid.Classic.GridData context, ICompositeQuadRule<DoubleEdgeQuadRule> domNrule, CoordinateSystem cs = Quadrature.CoordinateSystem.Physical)
            : base(noOfIntegralsPerCell, context, domNrule, cs) //
        {
            foreach(IChunkRulePair<QuadRule> crp in domNrule) {
                NodeCoordinateSystem ncs = crp.Rule.Nodes.GetNodeCoordinateSystem(context);
                if(ncs != NodeCoordinateSystem.EdgeCoord) {
                    throw new ArgumentException("Illegal node set for edge quadrature. Found some node set defined for: " + ncs.ToString() + ".");
                }
            }
        }

       

        /// <summary>
        /// not used/suitable for cell boundary quadrature
        /// </summary>
        protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {
            return gridData.iGeomEdges.SqrtGramian.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + L - 1 });
        }

        /// <summary>
        /// not used/suitable for cell boundary quadrature
        /// </summary>
        protected override MultidimensionalArray GetScalingsForNonlinElements(int i0, int L) {
            throw new NotImplementedException("todo.");
        }

        /// <summary>
        /// see <see cref="Quadrature{A,B}.CheckQuadratureChunk"/>
        /// </summary>
        protected override void CheckQuadratureChunk(int j0, int Len, int iKRef) {
            var edges = gridData.iGeomEdges;
            for (int e = j0; e < Len; e++) {
                int _iKref = edges.GetRefElementIndex(e);
                if (iKRef != _iKref)
                    throw new Exception("Internal error: mismatch between reference element index for given quadrature rule and specific edge that should be integrated.");
            }
            
        }

        /// <summary>
        /// see <see cref="Quadrature{A,B}.NextPart"/>
        /// </summary>
        protected override void NextPart(out bool Linear, out int NoOfElm, int i0, int Len) {
            var edges = base.gridData.iGeomEdges;

            bool bLinear = edges.IsEdgeAffineLinear(i0);
            Linear = bLinear;
            int L;
            for (L = 1; L < Len; L++) {
                if (edges.IsEdgeAffineLinear(i0+L) != bLinear) {
                    NoOfElm = L;
                    return;
                }
            }

            NoOfElm = Len;
        }

        /// <summary>
        /// not implemented
        /// </summary>
        public override int CurrentRuleRefElementIndex {
            get { throw new NotImplementedException(); }
        }

        /// <summary>
        /// performs the quadrature
        /// </summary>
        protected override void DoQuadrature(DoubleEdgeQuadRule quadRule, int j0, int _Bulksize) {
            var currentRuleWeights = quadRule.Weights;
            

            int NoOfNodes = quadRule.NoOfNodes;
            int Median = quadRule.Median;

            switch (CoordinateSystem) {
                case CoordinateSystem.Physical: {
                    
                    int JE = j0 + _Bulksize;
                    int j = j0;
                    while (j < JE) {

                        // determine sub-chunk
                        bool Linear; int Bulksize;
                        NextPart(out Linear, out Bulksize, j, _Bulksize - (j - j0));

                        if (Linear) {


                            MultidimensionalArray scalings = GetScalingsForLinearElements(j0, Bulksize);
                            Debug.Assert(scalings.Dimension == 1);
                            Debug.Assert(scalings.Lengths.Length == 1);
                            Debug.Assert(scalings.IsContinious);
                            Debug.Assert(quadRule.Weights.IsContinious);
                            Debug.Assert(quadRule.Weights.Dimension == 1);

                            for (int jj = 0; jj < Bulksize; jj++) {
                                int iEdge = jj + j;
                                double _scaling = scalings[jj];

                                // loop over integral components:
                                for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {

                                    double reIn = 0.0;
                                    double reOt = 0.0;

                                    // loop over the nodes of an edge
                                    int ne;
                                    for (ne = 0; ne < Median; ne++) {
                                        reIn += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                                    }
                                    for (; ne < NoOfNodes; ne++) {
                                        reOt += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                                    }
                                    //Debug.Assert(reIn == reOt);


                                    m_QuadResultsCollapsed[jj, 0, m] = reIn*_scaling;
                                    m_QuadResultsCollapsed[jj, 1, m] = reOt*_scaling;

                                }
                            }
                        } else {
                            var scalings = GetScalingsForNonlinElements(j0, Bulksize);
                            Debug.Assert(scalings.Dimension == 2);
                            Debug.Assert(scalings.GetLength(0) == Bulksize);
                            Debug.Assert(scalings.GetLength(1) == currentRuleWeights.GetLength(0));
                            Debug.Assert(quadRule.Weights.IsContinious);
                            Debug.Assert(quadRule.Weights.Dimension == 1);

                            for (int jj = 0; jj < Bulksize; jj++) {
                                int iEdge = jj + j;
                                double _scaling = scalings[jj];

                                // loop over integral components:
                                for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {

                                    double reIn = 0.0;
                                    double reOt = 0.0;

                                    // loop over the nodes of an edge
                                    int ne;
                                    for (ne = 0; ne < Median; ne++) {
                                        reIn += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                                    }
                                    for (; ne < NoOfNodes; ne++) {
                                        reOt += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                                    }
                                    //Debug.Assert(reIn == reOt);


                                    m_QuadResultsCollapsed[jj, 0, m] = reIn*_scaling;
                                    m_QuadResultsCollapsed[jj, 1, m] = reOt*_scaling;

                                }
                            }
                        }

                        j += Bulksize; 
                    }

                    break;
                }

                case CoordinateSystem.Reference:
                Debug.Assert(quadRule.Weights.IsContinious);
                Debug.Assert(quadRule.Weights.Dimension == 1);

                for (int jj = 0; jj < _Bulksize; jj++) {

                    // loop over integral components:
                    for (int m = 0; m < m_TotalNoOfIntegralsPerItem; m++) {

                        double reIn = 0.0;
                        double reOt = 0.0;

                        // loop over the nodes of an edge
                        int ne;
                        for (ne = 0; ne < Median; ne++) {
                            reIn += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                        }
                        for (; ne < NoOfNodes; ne++) {
                            reOt += m_EvalResultsCollapsed[jj, ne, m] * quadRule.Weights[ne];
                        }

                        m_QuadResultsCollapsed[jj, 0, m] = reIn;
                        m_QuadResultsCollapsed[jj, 1, m] = reOt;
                    }

                }

                break;

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
        protected override void AllocateBuffers(int NoOfItems, NodeSet rule) {
            int NoOfNodes = rule.GetLength(0);
            if (m_EvalResults == null)
                m_EvalResults = new MultidimensionalArray(2 + base.IntegralCompDim.Length);
            m_EvalResults.Allocate(((new int[] { NoOfItems, NoOfNodes }).Concat(IntegralCompDim)).ToArray());

            if (m_QuadResults == null)
                m_QuadResults = new MultidimensionalArray(2 + IntegralCompDim.Length);
            m_QuadResults.Allocate(((new int[] { NoOfItems, 2 }).Concat(IntegralCompDim)).ToArray());

            m_EvalResultsCollapsed = m_EvalResults.ResizeShallow(
                (m_EvalResults.Lengths.Take(2).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());
            m_QuadResultsCollapsed = m_QuadResults.ResizeShallow(
                (m_QuadResults.Lengths.Take(2).Concat(new int[] { m_TotalNoOfIntegralsPerItem })).ToArray());


            if (m_AllocateBuffers != null)
                m_AllocateBuffers(NoOfItems, rule);
        }


        /// <summary>
        /// creates a double-edge quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="CellQuadrature"/>.
        /// </summary>
        static public DoubleEdgeQuadrature GetQuadrature(int[] noOfIntegralsPerCell,
                              Grid.Classic.GridData context,
                              ICompositeQuadRule<DoubleEdgeQuadRule> domNrule,
                              Del_Evaluate _Evaluate,
                              Del_SaveIntegrationResults _SaveIntegrationResults,
                              Del_AllocateBuffers _AllocateBuffers = null,
                              Del_QuadNodesChanged _PostLockNodes = null,
                              CoordinateSystem cs = CoordinateSystem.Physical) {
            var ret = new DoubleEdgeQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_Evaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }


        private class DoubleEdgeQuadratureImpl : DoubleEdgeQuadrature {

            public DoubleEdgeQuadratureImpl(
                int[] noOfIntegralsPerCell, Grid.Classic.GridData context, ICompositeQuadRule<DoubleEdgeQuadRule> domNrule, CoordinateSystem cs)
                : base(noOfIntegralsPerCell, context, domNrule, cs) {
            }

        }

    }
}
