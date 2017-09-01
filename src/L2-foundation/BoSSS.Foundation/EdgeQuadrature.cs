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
    /// quadrature over edges.
    /// </summary>
    public abstract class EdgeQuadrature : Quadrature<QuadRule, EdgeMask> {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="noOfIntegralsPerCell"></param>
        /// <param name="context"></param>
        /// <param name="rule"></param>
        /// <param name="cs">integrate in physical or reference coordinates?</param>
        public EdgeQuadrature(int[] noOfIntegralsPerCell, IGridData context, ICompositeQuadRule<QuadRule> rule, CoordinateSystem cs = Quadrature.CoordinateSystem.Physical)
            : base(noOfIntegralsPerCell, context, rule, cs) //
        {
            foreach(IChunkRulePair<QuadRule> crp in rule) {
                NodeCoordinateSystem ncs = crp.Rule.Nodes.GetNodeCoordinateSystem(context);
                if(ncs != NodeCoordinateSystem.EdgeCoord) {
                    throw new ArgumentException("Illegal node set for edge quadrature. Found some node set defined for: " + ncs.ToString() + ".");
                }
            }
        }

        

        /// <summary>
        /// the square root of the Gramian determinant
        /// </summary>
        protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {
            var R = gridData.iGeomEdges.SqrtGramian.ExtractSubArrayShallow(new int[] { i0 }, new int[] { i0 + L - 1 });
            return R;
        }

        /// <summary>
        /// Sweeps whether edge <paramref name="i0"/> is linear/nonlinear and how many cells of the same type are going to come after it.
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
        /// the transformation metric for nonlinear edges.
        /// </summary>
        protected override MultidimensionalArray GetScalingsForNonlinElements(int i0, int L) {
            return base.GridDat.iGeomEdges.NormalsCache.GetIntegrationMetric(base.CurrentRule.Nodes, i0, L);
        }

        /// <summary>
        /// creates an edge quadrature, where integrand evaluation (<paramref name="_Evaluate"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="EdgeQuadrature"/>.
        /// </summary>
        static public EdgeQuadrature GetQuadrature(int[] noOfIntegralsPerCell,
                              Grid.Classic.GridData context,
                              ICompositeQuadRule<QuadRule> domNrule,
                              Del_Evaluate _Evaluate,
                              Del_SaveIntegrationResults _SaveIntegrationResults,
                              Del_AllocateBuffers _AllocateBuffers = null,
                              Del_QuadNodesChanged _PostLockNodes = null,
                              CoordinateSystem cs = CoordinateSystem.Physical) {
            var ret = new EdgeQuadratureImpl(noOfIntegralsPerCell, context, domNrule, cs) {
                m_Evaluate = _Evaluate,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }

        /// <summary>
        /// creates an edge quadrature, where integrand evaluation AND quadrature (!!!) (<paramref name="_EvaluateEx"/>) and other methods
        /// can be passed as delegates. Use this, if you do not want to derive from <see cref="EdgeQuadrature"/>.
        /// </summary>
        static public EdgeQuadrature GetQuadrature2(int[] noOfIntegralsPerEdge,
                              IGridData context,
                              ICompositeQuadRule<QuadRule> domNrule,
                              Del_Evaluate _EvaluateEx,
                              Del_SaveIntegrationResults _SaveIntegrationResults,
                              Del_AllocateBuffers _AllocateBuffers = null,
                              Del_QuadNodesChanged _PostLockNodes = null,
                              CoordinateSystem cs = CoordinateSystem.Physical) {
            var ret = new EdgeQuadratureImpl(noOfIntegralsPerEdge, context, domNrule, cs) {
                m_ExEvaluate = _EvaluateEx,
                m_SaveIntegrationResults = _SaveIntegrationResults,
                m_AllocateBuffers = _AllocateBuffers,
                m_quadNodesChanged = _PostLockNodes,

            };
            return ret;
        }


        private class EdgeQuadratureImpl : EdgeQuadrature {

            public EdgeQuadratureImpl(int[] noOfIntegralsPerCell,
                              IGridData context,
                              ICompositeQuadRule<QuadRule> domNrule,
                              CoordinateSystem cs)
                : base(noOfIntegralsPerCell, context, domNrule, cs) {
            }

        }


        /// <summary>
        /// Index of the current rule's reference element into
        /// <see cref="BoSSS.Foundation.Grid.GridData.EdgeData.EdgeRefElements"/>.
        /// </summary>
        public override int CurrentRuleRefElementIndex {
            get {
                var Kref = base.CurrentRule.RefElement;
                int ret = Array.IndexOf(gridData.iGeomEdges.EdgeRefElements, Kref);
                Debug.Assert(ret >= 0, "unable to identify reference element");
                return ret;
            }
        }

        /// <summary>
        /// Some DEBUG checks on the current quadrature chunk.
        /// </summary>
        protected override void CheckQuadratureChunk(int j0, int Len, int iKRef) {
            var edges = gridData.iGeomEdges;
            for (int e = j0; e < Len; e++) {
                int _iKref = edges.GetRefElementIndex(e);
                if (iKRef != _iKref)
                    throw new Exception("Internal error: mismatch between reference element index for given quadrature rule and specific edge that should be integrated.");
            }
        }
    }
}
