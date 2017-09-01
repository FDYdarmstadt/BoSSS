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
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Quadrature {

    /// <summary>
    /// a double-edge-quadrature <see cref="DoubleEdgeQuadrature"/> which uses the same rule for all edges.
    /// </summary>
    public class StandardDoubleEdgeRuleFactory : IQuadRuleFactory<DoubleEdgeQuadRule>  {

        /// <summary>
        /// ctor.
        /// </summary>
        public StandardDoubleEdgeRuleFactory(IGridData g, RefElement KrefEdge) {
            if (!g.iGeomEdges.EdgeRefElements.Contains(KrefEdge, (a, b) => object.ReferenceEquals(a, b)))
                throw new ArgumentException("The reference element provided is not an edge reference element of the given grid.", "KrefEdge");
            this.RefElement = KrefEdge;
        }

        /// <summary>
        /// the edge simplex of the grid
        /// </summary>
        public RefElement RefElement {
            get;
            private set;
        }

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        /// <summary>
        /// construction of the quadrature rule.
        /// </summary>
        public IEnumerable<IChunkRulePair<DoubleEdgeQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            var orgQr = this.RefElement.GetQuadratureRule(order);
            var dblQr = new DoubleEdgeQuadRule();
            dblQr.OrderOfPrecision = orgQr.OrderOfPrecision;
            int L = orgQr.NoOfNodes;
            int D = orgQr.Nodes.GetLength(1);
            dblQr.Median = L;
            dblQr.Nodes = new NodeSet(this.RefElement, L * 2, D);
            dblQr.Weights = MultidimensionalArray.Create(L*2);
            dblQr.Nodes.SetSubArray(orgQr.Nodes, new int[] { 0, 0 }, new int[] { L-1, D-1 });
            dblQr.Nodes.SetSubArray(orgQr.Nodes, new int[] { L, 0 }, new int[] { 2*L-1, D-1 });
            dblQr.Weights.SetSubArray(orgQr.Weights, new int[] { 0 }, new int[] { L-1 });
            dblQr.Weights.SetSubArray(orgQr.Weights, new int[] { L }, new int[] { 2*L-1 });
            dblQr.Nodes.LockForever();

            return mask.Select(chunk => new ChunkRulePair<DoubleEdgeQuadRule>(chunk,dblQr));
        }
    }
}
