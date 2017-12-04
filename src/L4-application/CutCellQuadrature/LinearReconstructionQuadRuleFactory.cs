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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Linq;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Grid.RefElements;

namespace CutCellQuadrature {

    class LinearReconstructionQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        private LevelSetTracker tracker;

        private IQuadRuleFactory<CellBoundaryQuadRule> rootFactory;

        public LinearReconstructionQuadRuleFactory(LevelSetTracker tracker, LineAndPointQuadratureFactory lineAndPointFactory) {
            this.tracker = tracker;
            this.rootFactory = lineAndPointFactory.GetPointFactory();
        }

        public RefElement RefElement {
            get {
                return tracker.GridDat.Grid.RefElements[0];
            }
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            QuadRule baseRule = Line.Instance.GetQuadratureRule(order);
            var pointRules = rootFactory.GetQuadRuleSet(mask, order);

            var result = new List<ChunkRulePair<QuadRule>>();
            foreach (Chunk chunk in mask) {
                foreach (int cell in chunk.Elements) {
                    var pointRule = pointRules.Single(pair => pair.Chunk.i0 <= cell && pair.Chunk.JE > cell).Rule;

                    QuadRule reconstructedRule;
                    switch (pointRule.NoOfNodes) {
                        case 0:
                        case 1:
                            // Cell not really cut
                            reconstructedRule = QuadRule.CreateEmpty(RefElement, 1, 2);
                            break;

                        case 2: {
                                LineSegment linearReconstruction = new LineSegment(
                                    RefElement.SpatialDimension,
                                    RefElement,
                                    pointRule.Nodes.GetRow(0),
                                    pointRule.Nodes.GetRow(1));

                                reconstructedRule = QuadRule.CreateEmpty(
                                    RefElement, baseRule.NoOfNodes, RefElement.SpatialDimension);
                                for (int i = 0; i < baseRule.NoOfNodes; i++) {
                                    reconstructedRule.Weights[i] =
                                        baseRule.Weights[i] * linearReconstruction.Length / Line.Instance.Volume;

                                    double[] point = linearReconstruction.GetPointOnSegment(baseRule.Nodes[i, 0]);
                                    for (int d = 0; d < RefElement.SpatialDimension; d++) {
                                        reconstructedRule.Nodes[i, d] = point[d];
                                    }
                                }
                            }
                            break;

                        case 4: {
                                // Assume only single interface
                                double[] xCoords = pointRule.Nodes.GetColumn(0);
                                double maxXDist = xCoords.Max() - xCoords.Min();

                                double[] yCoords = pointRule.Nodes.GetColumn(1);
                                double maxYDist = yCoords.Max() - yCoords.Min();

                                int orderingDirection;
                                if (maxXDist > maxYDist) {
                                    orderingDirection = 0;
                                } else {
                                    orderingDirection = 1;
                                }

                                double[][] roots = new double[pointRule.NoOfNodes][];
                                for (int i = 0; i < roots.Length; i++) {
                                    roots[i] = new double[] { pointRule.Nodes[i, 0], pointRule.Nodes[i, 1] };
                                }
                                var orderedRoots = roots.OrderBy(t => t[orderingDirection]).ToArray();

                                // Now, construct two separate rules
                                LineSegment linearReconstruction1 = new LineSegment(
                                    RefElement.SpatialDimension,
                                    RefElement,
                                    orderedRoots[0],
                                    orderedRoots[1]);
                                LineSegment linearReconstruction2 = new LineSegment(
                                    RefElement.SpatialDimension,
                                    RefElement,
                                    orderedRoots[2],
                                    orderedRoots[3]);

                                reconstructedRule = QuadRule.CreateEmpty(
                                    RefElement, 2 * baseRule.NoOfNodes, RefElement.SpatialDimension);
                                for (int i = 0; i < baseRule.NoOfNodes; i++) {
                                    reconstructedRule.Weights[i] =
                                        baseRule.Weights[i] * linearReconstruction1.Length / Line.Instance.Volume;

                                    double[] point = linearReconstruction1.GetPointOnSegment(baseRule.Nodes[i, 0]);
                                    for (int d = 0; d < RefElement.SpatialDimension; d++) {
                                        reconstructedRule.Nodes[i, d] = point[d];
                                    }

                                    int offset = baseRule.NoOfNodes;
                                    reconstructedRule.Weights[offset + i] =
                                        baseRule.Weights[i] * linearReconstruction2.Length / Line.Instance.Volume;

                                    point = linearReconstruction2.GetPointOnSegment(baseRule.Nodes[i, 0]);
                                    for (int d = 0; d < RefElement.SpatialDimension; d++) {
                                        reconstructedRule.Nodes[offset + i, d] = point[d];
                                    }

                                }
                            }
                            break;

                        default:
                            throw new NotImplementedException();
                    }

                    reconstructedRule.Nodes.LockForever();

                    MultidimensionalArray metrics =
                        tracker.DataHistories[0].Current.GetLevelSetNormalReferenceToPhysicalMetrics(reconstructedRule.Nodes, cell, 1);
                    for (int i = 0; i < reconstructedRule.NoOfNodes; i++) {
                        reconstructedRule.Weights[i] /= metrics[0, i];
                    }

                    result.Add(new ChunkRulePair<QuadRule>(
                        Chunk.GetSingleElementChunk(cell), reconstructedRule));
                }
            }

            return result;
        }
    }
}
