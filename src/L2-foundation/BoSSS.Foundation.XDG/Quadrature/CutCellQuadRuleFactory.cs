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

using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.Subdivision {

    /// <summary>
    /// A quad rule factory that has special support for cut cells. Uses a
    /// subdivision strategy (see <see cref="ISubdivisionStrategy"/>) to
    /// subdivide elements (i.e., cell or edges) and builds a quad rule that is
    /// (hopefully) suitable for the integration in cut cells.
    /// </summary>
    public class CutCellQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        /// <summary>
        /// Quadrature rule to be used in cut leaves (see
        /// <see cref="SubdivisionNode.IsCut"/>).
        /// </summary>
        private QuadRule cutNodeRule;

        /// <summary>
        /// Constructs a quad rule factory using the given subdivision
        /// strategy.
        /// </summary>
        /// <param name="subdivisionStrategy">
        /// <see cref="SubdivisionStrategy"/>
        /// </param>
        /// <param name="leafDivisions">
        /// The additional number of subdivisions for the quadrature rule to be
        /// used in cut leaves (see <see cref="SubdivisionNode.IsCut"/>).
        /// </param>
        public CutCellQuadRuleFactory(ISubdivisionStrategy subdivisionStrategy, int leafDivisions) {
            SubdivisionStrategy = subdivisionStrategy;
            if (leafDivisions >= 0) {
                cutNodeRule = RefElement.GetBruteForceQuadRule(leafDivisions, 1);
            }
        }

        /// <summary>
        /// The strategy that should be used for the construction of the
        /// cut-cell quadrature rules.
        /// </summary>
        public ISubdivisionStrategy SubdivisionStrategy {
            get;
            private set;
        }

        #region IQuadRuleFactory<QuadRule> Members

        /// <summary>
        /// The simplex to be integrated over.
        /// </summary>
        public RefElement RefElement {
            get {
                return SubdivisionStrategy.RefElement;
            }
        }

        /// <summary>
        /// Uses <see cref="SubdivisionStrategy"/> to construct the subdivision
        /// nodes of the given simplex. Then, each node is supplied with a
        /// standard quadrature rule of the desired order.
        /// </summary>
        /// <param name="mask">
        /// <see cref="CutCellQuadRuleFactory.GetQuadRuleSet"/>
        /// </param>
        /// <param name="order">
        /// <see cref="CutCellQuadRuleFactory.GetQuadRuleSet"/>
        /// </param>
        /// <returns>
        /// A quadrature rule that is <paramref name="order"/>th order accurate
        /// in every node constructed by <see cref="SubdivisionStrategy"/> that
        /// is not cut by the interface.
        /// </returns>
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            using (new FuncTrace()) {
                QuadRule baseRule = RefElement.GetQuadratureRule(order);

                var result = new List<ChunkRulePair<QuadRule>>();
                List<SubdivisionNode> lastNodeList = null;
                QuadRule lastRule = null;
                foreach (var chunkNodePair in SubdivisionStrategy.GetSubdivisionNodes(mask)) {
                    List<SubdivisionNode> nodeList = chunkNodePair.Value.ToList();

                    // Make sure consecutive chunks share the same quad rule
                    // instance if quad rules are equal.
                    if (lastNodeList != null && nodeList.SequenceEqual(lastNodeList)) {
                        result.Add(new ChunkRulePair<QuadRule>(
                            chunkNodePair.Key, lastRule));
                        continue;
                    }

                    QuadRule[] quadRules = new QuadRule[nodeList.Count];
                    int noOfNodes = 0;
                    for (int i = 0; i < nodeList.Count; i++) {
                        QuadRule leafRule = baseRule;
                        if (nodeList[i].IsCut) {
                            leafRule = cutNodeRule ?? baseRule;
                        }

                        quadRules[i] = leafRule;
                        noOfNodes += leafRule.NoOfNodes;
                    }

                    QuadRule rule = QuadRule.CreateEmpty(RefElement, noOfNodes, RefElement.SpatialDimension);
                    int offset = 0;
                    for (int k = 0; k < nodeList.Count; k++) {
                        QuadRule leaveRule = quadRules[k];
                        SubdivisionNode node = nodeList[k];
                        double det = node.Transformation.Matrix.Determinant();

                        for (int i = 0; i < leaveRule.NoOfNodes; i++) {
                            double[] vertex = node.Transformation.Transform(
                                leaveRule.Nodes.ExtractSubArrayShallow(i, -1).To1DArray());
                            rule.Nodes.SetSubVector(vertex, offset + i, -1);
                            rule.Weights[offset + i] = det * leaveRule.Weights[i];
                        }

                        offset += leaveRule.NoOfNodes;
                    }

                    result.Add(new ChunkRulePair<QuadRule>(
                        chunkNodePair.Key, rule));

                    lastRule = rule;
                    rule.Nodes.LockForever();
                    lastNodeList = nodeList;
                }

                return result;
            }
        }


        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        #endregion
    }
}
