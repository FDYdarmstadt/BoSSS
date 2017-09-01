using System;
using System.Collections.Generic;
using System.Diagnostics;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    public class MakeshiftCutLineQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        private const double EPSILON = 1e-6;

        private const int MAX_SIMPLEX_LEVELS = 3;

        private static Line lineSimplex = new Line();

        private LevelSetTracker tracker;

        private int levSetIndex;

        private LineSegment[] referenceLineSegments;

        public MakeshiftCutLineQuadRuleFactory(LevelSetTracker tracker, int levSetIndex)
            : this(tracker, levSetIndex, new LineSegment.SafeGuardedNewtonMethod(1e-14)) {
        }

        public MakeshiftCutLineQuadRuleFactory(LevelSetTracker tracker, int levSetIndex, LineSegment.IRootFindingAlgorithm rootFindingAlgorithm) {
            this.tracker = tracker;
            if (tracker.LevelSets.Count <= levSetIndex)
                throw new ArgumentOutOfRangeException("Please provide a valid index for the level set.");
            this.levSetIndex = levSetIndex;
            this.RootFindingAlgorithm = rootFindingAlgorithm;
            this.referenceLineSegments = GetReferenceLineSegments();
        }

        public LineSegment.IRootFindingAlgorithm RootFindingAlgorithm {
            get;
            private set;
        }

        #region IQuadRuleFactory<QuadRule> Members

        public Simplex Simplex {
            get {
                return tracker.Ctx.Grid.GridSimplex;
            }
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask == null) {
                mask = EdgeMask.GetFullMask(tracker.Ctx.GridDat);
            }

            if (mask is EdgeMask == false) {
                throw new Exception();
            }

            QuadRule baseRule = lineSimplex.GetQuadratureRule(order);
            int D = tracker.Ctx.Grid.GridSimplex.SpatialDimension;

            var result = new List<ChunkRulePair<QuadRule>>(mask.NoOfItemsLocally);
            foreach (Chunk chunk in mask) {
                for (int i = 0; i < chunk.Len; i++) {
                    List<double[]> nodes = new List<double[]>();
                    List<double> weights = new List<double>();
                    int[] noOfNodesPerEdge = new int[referenceLineSegments.Length];

                    int edge = i + chunk.i0;

                    // Always choose 'left' edge
                    int cell = tracker.Ctx.GridDat.Edges[edge, 0];
                    int localEdge = tracker.Ctx.GridDat.EdgeIndices[edge, 0];

                    LineSegment referenceSegment = referenceLineSegments[localEdge];
                    double[] roots = referenceSegment.GetRoots(tracker.LevelSets[levSetIndex], cell);

                    LineSegment[] subSegments = referenceSegment.Split(roots);

                    for (int k = 0; k < subSegments.Length; k++) {
                        for (int m = 0; m < baseRule.NoOfNodes; m++) {
                            // Base rule _always_ is a line rule, thus Nodes[*, _0_]
                            double[] point = subSegments[k].GetPointOnSegment(baseRule.Nodes[m, 0]);

                            uint lh = tracker.Ctx.NSC.LockNodeSetFamily(
                                tracker.Ctx.NSC.CreateContainer(
                                    MultidimensionalArray.CreateWrapper(point, 1, D), -1.0));
                            MultidimensionalArray levelSetValue = tracker.GetLevSetValues(levSetIndex, 0, cell, 1);
                            tracker.Ctx.NSC.UnlockNodeSetFamily(lh);

                            // Only positive volume
                            if (levelSetValue[0, 0] <= 0.0) {
                                continue;
                            }

                            weights.Add(baseRule.Weights[m] * subSegments[k].Length / referenceSegment.Length);
                            nodes.Add(point);
                        }
                    }

                    if (weights.Count == 0) {
                        continue;
                    }

                    MultidimensionalArray localNodes = MultidimensionalArray.Create(nodes.Count, D);
                    for (int j = 0; j < nodes.Count; j++) {
                        for (int d = 0; d < D; d++) {
                            localNodes[j, d] = nodes[j][d];
                        }
                    }

                    MultidimensionalArray localEdgeNodes = MultidimensionalArray.Create(nodes.Count, 1);
                    tracker.Ctx.Grid.GridSimplex.VolumeToEdgeCoordinates(localEdge, localNodes, localEdgeNodes);

                    QuadRule subdividedRule = new QuadRule() {
                        OrderOfPrecision = order,
                        Weights = MultidimensionalArray.Create(weights.Count),
                        Nodes = localEdgeNodes.CloneAs()
                    };
                    subdividedRule.Weights.SetV(weights, -1);

                    result.Add(new ChunkRulePair<QuadRule>(
                        Chunk.GetSingleElementChunk(edge), subdividedRule));
                }
            }

            return result;
        }

        #endregion

        public LineSegment[] GetReferenceLineSegments() {
            Stack<Simplex> simplexHierarchy = new Stack<Simplex>();
            Simplex currentSimplex = Simplex;
            int spatialDimension = Simplex.SpatialDimension;

            int n = 0;
            while (currentSimplex.GetType() != lineSimplex.GetType()) {
                if (n > MAX_SIMPLEX_LEVELS) {
                    throw new ApplicationException("Something went terribly wrong. Please contact Björn");
                }

                simplexHierarchy.Push(currentSimplex);

                currentSimplex = currentSimplex.EdgeSimplex;
                n++;
            }

            MultidimensionalArray vertexCoordinates = MultidimensionalArray.Create(2, 1);
            vertexCoordinates[0, 0] = -1.0;
            vertexCoordinates[1, 0] = 1.0;

            while (simplexHierarchy.Count > 0) {
                currentSimplex = simplexHierarchy.Pop();

                int noOfVertices = vertexCoordinates.GetLength(0);
                int D = currentSimplex.SpatialDimension;
                MultidimensionalArray volumeCoordinates = MultidimensionalArray.Create(
                    noOfVertices * currentSimplex.NoOfEdges, currentSimplex.SpatialDimension);

                for (int e = 0; e < currentSimplex.NoOfEdges; e++) {
                    MultidimensionalArray coordinates = MultidimensionalArray.Create(noOfVertices, D);
                    currentSimplex.EdgeToVolumeCoordinates(e, vertexCoordinates, coordinates);

                    for (int i = 0; i < noOfVertices; i++) {
                        for (int d = 0; d < D; d++) {
                            volumeCoordinates[e * noOfVertices + i, d] = coordinates[i, d];
                        }
                    }
                }

                vertexCoordinates = volumeCoordinates;
            }

            Debug.Assert(
                vertexCoordinates.GetLength(0) % 2 == 0,
                "Even number of vertices expected");
            int initialNumberOfLineSegments = vertexCoordinates.GetLength(0) / 2;

            List<LineSegment> lineSegments = new List<LineSegment>(initialNumberOfLineSegments);
            for (int i = 0; i < initialNumberOfLineSegments; i++) {
                LineSegment newSegment = new LineSegment(spatialDimension, RootFindingAlgorithm);
                for (int d = 0; d < spatialDimension; d++) {
                    newSegment.Start[d] = vertexCoordinates[2 * i + 0, d];
                    newSegment.End[d] = vertexCoordinates[2 * i + 1, d];
                }

                if (!lineSegments.Contains(newSegment)) {
                    lineSegments.Add(newSegment);
                    tracker.Subscribe(newSegment);
                }
            }

            foreach (LineSegment segment in lineSegments) {
                LevelSet levelSetField = tracker.LevelSets[levSetIndex] as LevelSet;
                if (levelSetField != null) {
                    segment.ProjectBasisPolynomials(levelSetField.Basis);
                }
            }

            return lineSegments.ToArray();
        }
    }
}
