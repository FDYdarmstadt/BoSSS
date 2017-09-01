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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP;
using System.Linq;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// Quadrature for the boundary of the boundary (also known as the co-faces) of a three-dimensional
    /// cell. That is, a quadrature rule for the line elements bounding the
    /// three-dimensional reference element.
    /// </summary>
    public class CellEdgeBoundaryQuadRule : CellBoundaryQuadRule {

        /// <summary>
        /// Creates an empty rule
        /// </summary>
        /// <param name="noOfNodes">
        /// The desired number of nodes
        /// </param>
        /// <param name="element">
        /// The selected reference element
        /// </param>
        /// <returns>
        /// A quadrature rule with <paramref name="noOfNodes"/> nodes where all
        /// relevant entries are zero.
        /// </returns>
        public static CellEdgeBoundaryQuadRule CreateEmpty(int noOfNodes, RefElement element) {
            if (element.SpatialDimension < 3) {
                throw new ArgumentException("Only makes sense for 3D objects");
            }

            CellEdgeBoundaryQuadRule rule = new CellEdgeBoundaryQuadRule() {
                Nodes = new NodeSet(element, noOfNodes, element.SpatialDimension),
                Weights = MultidimensionalArray.Create(noOfNodes),
                OrderOfPrecision = 0,
                NumbersOfNodesPerFace = new int[element.NoOfFaces],
                NumbersOfNodesPerFaceOfFace =
                    new int[element.NoOfFaces, element.FaceRefElement.NoOfFaces]
            };
            return rule;
        }

        /// <summary>
        /// 1st index: edge index, which is a (D-1)-dimensional manifold.
        /// 2nd index: edge of edge, which is a (D-2)-dimensional manifold.
        /// </summary>
        public int[,] NumbersOfNodesPerFaceOfFace {
            get;
            set;
        }
    }

    /// <summary>
    /// A factory for quadrature rule that creates Gaussian quadrature rules on
    /// sub-sections of the lines bounding a reference element (cf.
    /// <see cref="CellEdgeBoundaryQuadRule"/>) in the style of
    /// <see cref="CutLineQuadRuleFactory"/>
    /// </summary>
    class CutLineOnEdgeQuadRuleFactory : IQuadRuleFactory<CellEdgeBoundaryQuadRule>, IObserver<LevelSetTracker.LevelSetRegionsInfo> {

        /// <summary>
        /// Minimal distance between two points
        /// </summary>
        private const double EPSILON = 1e-14;

        /// <summary>
        /// Static line element
        /// </summary>
        private static readonly Line lineSimplex = Line.Instance;

        /// <summary>
        /// Tracks the level set location
        /// </summary>
        private LevelSetTracker tracker;

        /// <summary>
        /// The line segments of the reference element.
        /// </summary>
        private LineSegment[,] referenceLineSegments;

        /// <summary>
        /// The considered jump type, i.e. whether to integrate over the
        /// positive region, the negative region, or over both
        /// </summary>
        private JumpTypes jumpType;

        /// <summary>
        /// The last integration order used in <see cref="GetQuadRuleSet"/>.
        /// </summary>
        private int lastOrder = -1;

        /// <summary>
        /// The index of the level set to be considered
        /// </summary>
        private int levelSetIndex;

        /// <summary>
        /// Cache for the quadrature rules. The key denotes the cell index
        /// </summary>
        /// <remarks>
        /// Is cleared every time the level set degree (cf.
        /// <see cref="lastOrder"/>) is changed level set moves (see
        /// <see cref="OnNext"/>)
        /// </remarks>
        private Dictionary<int, CellEdgeBoundaryQuadRule> cache =
            new Dictionary<int, CellEdgeBoundaryQuadRule>();

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return cache.Keys.ToArray();
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="tracker"></param>
        /// <param name="levelSetIndex"></param>
        /// <param name="rootFindingAlgorithm"></param>
        /// <param name="jumpType"></param>
        public CutLineOnEdgeQuadRuleFactory(LevelSetTracker tracker, int levelSetIndex, LineSegment.IRootFindingAlgorithm rootFindingAlgorithm = null, JumpTypes jumpType = JumpTypes.Heaviside) {
            if (tracker.GridDat.SpatialDimension < 3) {
                throw new ArgumentException("Only applicable in 3d", "tracker");
            }

            this.tracker = tracker;
            this.RootFindingAlgorithm = rootFindingAlgorithm ?? LineSegment.DefaultRootFindingAlgorithm;
            this.jumpType = jumpType;
            if (levelSetIndex >= tracker.LevelSets.Count) {
                throw new ArgumentOutOfRangeException("Please specify a valid index for the level set function");
            }
            this.levelSetIndex = levelSetIndex;
            this.referenceLineSegments = GetReferenceLineSegments();

            tracker.Subscribe(this);
        }

        /// <summary>
        /// Algorithm used for the root search on the individual line segments
        /// </summary>
        public LineSegment.IRootFindingAlgorithm RootFindingAlgorithm {
            get;
            private set;
        }

        #region IQuadRuleFactory<CellBndQuadRule> Members

        /// <summary>
        /// Returns a set of <see cref="CellEdgeBoundaryQuadRule"/>s that
        /// enables the integration over sub-segments of the edges of the edges
        /// of a (three-dimensional) domain. This is obviously only useful if
        /// the integrand has a discontinuity that is aligned with the zero
        /// iso-contour of the level set function.
        /// </summary>
        /// <param name="mask"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        public IEnumerable<IChunkRulePair<CellEdgeBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask == null) {
                mask = CellMask.GetFullMask(tracker.GridDat);
            }

            if (mask is CellMask == false) {
                throw new ArgumentException("Edge mask required", "mask");
            }

            if (lastOrder != order) {
                cache.Clear();
            }

            QuadRule baseRule = lineSimplex.GetQuadratureRule(order);
            int D = tracker.GridDat.SpatialDimension;
            int noOfEdges = tracker.GridDat.Grid.RefElements[0].NoOfFaces;
            int noOfEdgesOfEdge = RefElement.FaceRefElement.NoOfFaces;

            var result = new List<ChunkRulePair<CellEdgeBoundaryQuadRule>>(mask.NoOfItemsLocally);
            foreach (Chunk chunk in mask) {
                for (int i = 0; i < chunk.Len; i++) {
                    int cell = i + chunk.i0;

                    if (cache.ContainsKey(cell)) {
                        result.Add(new ChunkRulePair<CellEdgeBoundaryQuadRule>(
                            Chunk.GetSingleElementChunk(cell),
                            cache[cell]));
                        continue;
                    }

                    List<double[]> nodes = new List<double[]>();
                    List<double> weights = new List<double>();

                    if (tracker.GridDat.Cells.Cells2Edges[cell].Length != noOfEdges) {
                        throw new NotImplementedException("Not implemented for hanging nodes");
                    }

                    int[] noOfNodesPerEdge = new int[noOfEdges];
                    int[,] noOfNodesPerEdgeOfEdge = new int[noOfEdges, noOfEdgesOfEdge];
                    for (int e = 0; e < noOfEdges; e++) {
                        int edge = Math.Abs(tracker.GridDat.Cells.Cells2Edges[cell][e]) - 1;
                        double edgeDet = tracker.GridDat.Edges.SqrtGramian[edge];

                        for (int ee = 0; ee < noOfEdgesOfEdge; ee++) {
                            LineSegment refSegment = referenceLineSegments[e, ee];
                            double edgeOfEdgeDet = RefElement.FaceRefElement.FaceTrafoGramianSqrt[ee];

                            double[] roots = refSegment.GetRoots(tracker.LevelSets[levelSetIndex], cell, 0);
                            LineSegment[] subSegments = refSegment.Split(roots);

                            for (int k = 0; k < subSegments.Length; k++) {
                                // Evaluate sub segment at center to determine sign
                                NodeSet _point = new NodeSet(this.RefElement, subSegments[k].GetPointOnSegment(0.0));
                                
                                double scaling = edgeOfEdgeDet * subSegments[k].Length / refSegment.Length;

                                if(jumpType != JumpTypes.Implicit) {
                                    //using (tracker.GridDat.NSC.CreateLock(
                                    //    MultidimensionalArray.CreateWrapper(point, 1, D), 0, -1.0)) {
                                    MultidimensionalArray levelSetValue = tracker.GetLevSetValues(levelSetIndex, _point, cell, 1);

                                    switch(jumpType) {
                                        case JumpTypes.Heaviside:
                                        if(levelSetValue[0, 0] <= 0.0) {
                                            continue;
                                        }
                                        break;

                                        case JumpTypes.OneMinusHeaviside:
                                        if(levelSetValue[0, 0] > 0.0) {
                                            continue;
                                        }
                                        break;

                                        case JumpTypes.Sign:
                                        scaling *= levelSetValue[0, 0].Sign();
                                        break;

                                        default:
                                        throw new NotImplementedException();
                                    }

                                }

                                for (int m = 0; m < baseRule.NoOfNodes; m++) {
                                    // Base rule _always_ is a line rule, thus Nodes[*, _0_]
                                    double[] point = subSegments[k].GetPointOnSegment(baseRule.Nodes[m, 0]);

                                    weights.Add(baseRule.Weights[m] * scaling);
                                    nodes.Add(point);

                                    noOfNodesPerEdge[e]++;
                                    noOfNodesPerEdgeOfEdge[e, ee]++;
                                }
                            }
                        }
                    }

                    if (weights.Count == 0) {
                        CellEdgeBoundaryQuadRule emptyRule =
                            CellEdgeBoundaryQuadRule.CreateEmpty(1, RefElement);
                        emptyRule.Nodes.LockForever();
                        cache.Add(cell, emptyRule);
                        result.Add(new ChunkRulePair<CellEdgeBoundaryQuadRule>(
                            Chunk.GetSingleElementChunk(cell), emptyRule));
                        continue;
                    }

                    NodeSet localNodes = new NodeSet(this.RefElement, nodes.Count, D);
                    for (int j = 0; j < nodes.Count; j++) {
                        for (int d = 0; d < D; d++) {
                            localNodes[j, d] = nodes[j][d];
                        }
                    }
                    localNodes.LockForever();

                    CellEdgeBoundaryQuadRule subdividedRule = new CellEdgeBoundaryQuadRule() {
                        OrderOfPrecision = order,
                        Weights = MultidimensionalArray.Create(weights.Count),
                        Nodes = localNodes,
                        NumbersOfNodesPerFace = noOfNodesPerEdge,
                        NumbersOfNodesPerFaceOfFace = noOfNodesPerEdgeOfEdge
                    };
                    subdividedRule.Weights.SetSubVector(weights, -1);

                    cache.Add(cell, subdividedRule);

                    result.Add(new ChunkRulePair<CellEdgeBoundaryQuadRule>(
                        Chunk.GetSingleElementChunk(cell), subdividedRule));
                }
            }

            return result;
        }

        /// <summary>
        /// The reference element
        /// </summary>
        public RefElement RefElement {
            get {
                return tracker.GridDat.Grid.RefElements[0];
            }
        }

        #endregion

        /// <summary>
        /// Gathers the line segments bounding <see cref="RefElement"/> in
        /// reference coordinates.
        /// </summary>
        /// <returns>
        /// An array of line segments
        /// <list type="bullet">
        ///     <item>1st index: Edge index</item>
        ///     <item>2nd index: Edge of edge index</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Each line segment will appear twice in the result, since it is
        /// stored for each edge separately. However, this method ensures that
        /// corresponding segments are represented by the same object (reference
        /// equality) which ensures that this does not affect performance.
        /// </remarks>
        private LineSegment[,] GetReferenceLineSegments() {
            int D = tracker.GridDat.SpatialDimension;
            int noOfEdges = tracker.GridDat.Grid.RefElements[0].NoOfFaces;
            int noOfEdgesOfEdge = RefElement.FaceRefElement.NoOfFaces;

            MultidimensionalArray edgeOfEdgeVertices = MultidimensionalArray.Create(2, 1);
            edgeOfEdgeVertices[0, 0] = -1.0;
            edgeOfEdgeVertices[1, 0] = 1.0;

            MultidimensionalArray edgeVertices = MultidimensionalArray.Create(
                edgeOfEdgeVertices.GetLength(0), RefElement.FaceRefElement.SpatialDimension);
            MultidimensionalArray volumeVertices = MultidimensionalArray.Create(
                edgeVertices.GetLength(0), tracker.GridDat.SpatialDimension);

            // Remember encountered segments so that $lineSegments contains no
            // duplicates and roots can be cached efficiently
            Dictionary<LineSegment, LineSegment> seenSegments = new Dictionary<LineSegment, LineSegment>(
                noOfEdges * noOfEdgesOfEdge / 2);

            LineSegment[,] lineSegments = new LineSegment[noOfEdges, noOfEdgesOfEdge];
            LevelSet levelSetField = tracker.LevelSets[levelSetIndex] as LevelSet;
            for (int ee = 0; ee < noOfEdgesOfEdge; ee++) {
                RefElement.FaceRefElement.TransformFaceCoordinates(ee, edgeOfEdgeVertices, edgeVertices);

                for (int e = 0; e < noOfEdges; e++) {
                    tracker.GridDat.Grid.RefElements[0].TransformFaceCoordinates(
                        e, edgeVertices, volumeVertices);

                    double[] start = new double[D];
                    double[] end = new double[D];
                    for (int d = 0; d < D; d++) {
                        start[d] = volumeVertices[0, d];
                        end[d] = volumeVertices[1, d];
                    }
                    LineSegment newSegment = new LineSegment(D, this.RefElement, start, end, rootFindingAlgorithm: RootFindingAlgorithm);

                    // Assert that the segment does not already exist
                    LineSegment segment;
                    if (!seenSegments.TryGetValue(newSegment, out segment)) {
                        segment = newSegment;
                        seenSegments.Add(segment, segment);

                        if (levelSetField != null) {
                            segment.ProjectBasisPolynomials(levelSetField.Basis);
                        }

                        tracker.Subscribe(segment);
                    }

                    lineSegments[e, ee] = segment;
                }
            }

            return lineSegments;
        }

        #region IObserver<LevelSetData> Members

        /// <summary>
        /// Empty.
        /// </summary>
        public void OnCompleted() {
        }

        /// <summary>
        /// Empty.
        /// </summary>
        /// <param name="error"></param>
        public void OnError(Exception error) {
        }

        /// <summary>
        /// Clears cached quadrature rules.
        /// </summary>
        /// <param name="value"></param>
        public void OnNext(LevelSetTracker.LevelSetRegionsInfo value) {
            cache.Clear();
        }

        #endregion
    }
}
