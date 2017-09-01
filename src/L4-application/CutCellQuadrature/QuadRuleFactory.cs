using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.LevelSet;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Utils;

namespace RegularizedQuadrature {

    public class QuadRuleFactory : IQuadRuleFactoryOld {

        private Context context;

        private LevelSetTracker tracker;

        private int maxLevels;

        private int minLevels;

        private NestedVertexSet baseVertexSet;

        private SimplexSubdivisionTree subdivisionTree;

        public QuadRuleFactory(Context context, LevelSetTracker tracker, int maxLevels, int minLevels) {
            this.context = context;
            this.tracker = tracker;
            this.maxLevels = maxLevels;
            this.minLevels = minLevels;

            baseVertexSet = new NestedVertexSet(context.Grid.SpatialDimension);

            double[,] vertices = context.Grid.GridSimplex.Vertices;
            int verticesPerCell = vertices.GetLength(0);
            int[] polyhedronVertices = new int[verticesPerCell];
            for (int i = 0; i < verticesPerCell; i++) {
                double[] vertex = ArrayTools.GetRow(vertices, i);
                polyhedronVertices[i] = baseVertexSet.RegisterVertex(vertex);
            }

            subdivisionTree = new SimplexSubdivisionTree(context.Grid.GridSimplex, baseVertexSet, polyhedronVertices);

            for (int i = 0; i < minLevels; i++) {
                subdivisionTree.Subdivide(baseVertexSet, false);
            }
            subdivisionTree.SetSavePoint();
        }

        private static double[,] NodesToRectangularArray(List<double[]> nodes) {
            int N = nodes.Count;
            int D = nodes.First().Length;

            double[,] result = new double[N, D];
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < D; j++) {
                    result[i, j] = nodes[i][j];
                }
            }

            return result;
        }

        #region IQuadRuleFactory Members

        public QuadRule GetStandardQuadRule(int order) {
            return context.Grid.GridSimplex.GetQuadratureRule(order);
        }

        public QuadRule GetCutCellQuadRule(LevelSet levelSet, int baseOrder, int cell) {
            // Build tree
            subdivisionTree.ResetToSavePoint();

            NestedVertexSet currentSet = baseVertexSet;
            for (int i = minLevels; i < maxLevels; i++) {
                if (currentSet.LocalNumberOfVertices == 0) {
                    // No new vertices were added during last subdivision
                    break;
                }

                NodeSetController.NodeSetContainer nsc = context.NSC.CreateContainer(currentSet.Vertices, -1.0);
                uint lh = context.NSC.LockNodeSetFamily(nsc);
                MultidimensionalArray levelSetValues = MultidimensionalArray.Create(1, nsc.NodeSet.GetLength(0));
                levelSet.Evaluate(cell, 1, 0, levelSetValues);
                context.NSC.UnlockNodeSetFamily(lh);

                subdivisionTree.ReadLevelSetValues(levelSetValues.ExtractSubArrayShallow(0, -1));

                currentSet = new NestedVertexSet(currentSet);
                subdivisionTree.Subdivide(currentSet, true);
            }

            // Read level set values of leaves (only if IsCut is used!)
            if (currentSet.Vertices != null) {
                NodeSetController.NodeSetContainer nsc2 = context.NSC.CreateContainer(currentSet.Vertices, -1.0);
                uint lh2 = context.NSC.LockNodeSetFamily(nsc2);
                MultidimensionalArray levelSetValues2 = MultidimensionalArray.Create(1, nsc2.NodeSet.GetLength(0));
                levelSet.Evaluate(cell, 1, 0, levelSetValues2);
                context.NSC.UnlockNodeSetFamily(lh2);
                subdivisionTree.ReadLevelSetValues(levelSetValues2.ExtractSubArrayShallow(0, -1));
            }

            // Construct rule
            List<double[]> nodes = new List<double[]>();
            List<double> weights = new List<double>();
            foreach (SimplexSubdivisionTree.Node leave in subdivisionTree.Leaves) {
                double det = leave.TransformationFromRoot.Matrix.Determinat();
                QuadRule rule;
                if (leave.IsCut) {
                    rule = GetStandardQuadRule(baseOrder);
                } else {
                    rule = GetStandardQuadRule(1);
                }

                for (int i = 0; i < rule.NoOfNodes; i++) {
                    double[] vertex = ArrayTools.GetRow(rule.Nodes, i);
                    nodes.Add(leave.TransformationFromRoot.Transform(vertex));
                    weights.Add(det * rule.Weights[i]);
                }
            }

            QuadRule result = new QuadRule();
            result.Nodes = NodesToRectangularArray(nodes);
            result.Weights = weights.ToArray();

            return result;
        }

        #endregion
    }
}
