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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CutCellQuadrature {

    class EquivalentPolynomialQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        private IQuadRuleFactory<QuadRule> baseFactory;

        private LevelSetTracker tracker;

        private IQuadRuleFactory<CellBoundaryQuadRule> rootFactory;

        public EquivalentPolynomialQuadRuleFactory(IQuadRuleFactory<QuadRule> baseFactory, LevelSetTracker tracker, LineAndPointQuadratureFactory lineAndPointFactory) {
            this.baseFactory = baseFactory;
            this.tracker = tracker;
            this.rootFactory = lineAndPointFactory.GetPointFactory();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            if (mask.MaskType != MaskType.Geometrical)
                throw new ArgumentException("Expecting a geometrical mask.");
            
            var standardVolumeRules = baseFactory.GetQuadRuleSet(mask, order);
            var pointRules = rootFactory.GetQuadRuleSet(mask, order);

            var result = new List<ChunkRulePair<QuadRule>>();
            foreach (var chunkRulePair in standardVolumeRules) {
                foreach (int cell in chunkRulePair.Chunk.Elements) {
                    QuadRule standardRule = chunkRulePair.Rule;

                    var pointRule = pointRules.Single(pair => pair.Chunk.i0 <= cell && pair.Chunk.JE > cell).Rule;
                    MultidimensionalArray gradientsAtRoots = tracker.DataHistories[0].Current.GetLevelSetGradients(pointRule.Nodes, cell, 1);

                    QuadRule modifiedRule;
                    switch (pointRule.NoOfNodes) {
                        case 0:
                        case 1: {
                                // Cell not really cut
                                MultidimensionalArray levelSetValue = tracker.DataHistories[0].Current.GetLevSetValues(
                                    new NodeSet(RefElement, new double[2]), cell, 1);
                                if (Math.Sign(levelSetValue.Storage[0]) < 0) {
                                    // Cell is completely void
                                    QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, 2);
                                    emptyRule.Nodes.LockForever();
                                    modifiedRule = emptyRule;
                                } else {
                                    // Cell is completely non-void
                                    modifiedRule = standardRule;
                                }
                            }
                            break;

                        case 2: {
                                double a, b, c;
                                GetLevelSetApproximationCofficients(cell, pointRule.Nodes.GetRow(0), pointRule.Nodes.GetRow(1), out a, out b, out c);
                                double[] eqcv = Heqpol_coefficients(a, b, c);

                                modifiedRule = standardRule.CloneAs();
                                for (int i = 0; i < modifiedRule.NoOfNodes; i++) {
                                    double x = modifiedRule.Nodes[i, 0];
                                    double y = modifiedRule.Nodes[i, 1];
                                    // double[] v = new double[] {
                                    //     1.0, x, x * x, y, x * y, y * y };
                                    double v = eqcv[0] + eqcv[1] * x + eqcv[2] * x * x
                                        + eqcv[3] * y + eqcv[4] * x * y + eqcv[5] * y * y;

                                    modifiedRule.Weights[i] *= v;
                                }

                                modifiedRule.Nodes.LockForever();
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

                                // Now, do the equivalent polynomial thing
                                double a, b, c;
                                GetLevelSetApproximationCofficients(cell, orderedRoots[0], orderedRoots[1], out a, out b, out c);
                                double[] eqcv1 = Heqpol_coefficients(a, b, c);

                                GetLevelSetApproximationCofficients(cell, orderedRoots[2], orderedRoots[3], out a, out b, out c);
                                double[] eqcv2 = Heqpol_coefficients(a, b, c);

                                modifiedRule = standardRule.CloneAs();
                                for (int i = 0; i < modifiedRule.NoOfNodes; i++) {
                                    double x = modifiedRule.Nodes[i, 0];
                                    double y = modifiedRule.Nodes[i, 1];
                                    // double[] v = new double[] {
                                    //     1.0, x, x * x, y, x * y, y * y };
                                    double v1 = eqcv1[0] + eqcv1[1] * x + eqcv1[2] * x * x
                                        + eqcv1[3] * y + eqcv1[4] * x * y + eqcv1[5] * y * y;
                                    double v2 = eqcv2[0] + eqcv2[1] * x + eqcv2[2] * x * x
                                        + eqcv2[3] * y + eqcv2[4] * x * y + eqcv2[5] * y * y;

                                    modifiedRule.Weights[i] *= v1 * v2;
                                }
                            }
                            break;

                        default:
                            throw new NotImplementedException();
                    }

                    result.Add(new ChunkRulePair<QuadRule>(
                        Chunk.GetSingleElementChunk(cell), modifiedRule));
                }
            }

            return result;
        }

        private void GetLevelSetApproximationCofficients(int cell, double[] firstNode, double[] secondNode, out double a, out double b, out double c) {
            Vector2D p = new Vector2D(firstNode[0], firstNode[1]);
            Vector2D q = new Vector2D(secondNode[0], secondNode[1]);
            Vector2D n = new Vector2D(-(q[1] - p[1]), q[0] - p[0]);

            // Use gradient of first root to determine correct direction of normal vector
            NodeSet bla = new NodeSet(RefElement, firstNode);
            MultidimensionalArray gradient = tracker.DataHistories[0].Current.GetLevelSetGradients(bla, cell, 1);
            Vector2D gradientVector = new Vector2D(
                gradient[0, 0, 0], gradient[0, 0, 1]);
            n.Scale(Math.Sign(gradientVector * n));

            a = n[0];
            b = n[1];
            c = p * n * -1.0;
        }

        public RefElement RefElement {
            get {
                return baseFactory.RefElement;
            }
        }

        public int[] GetCachedRuleOrders() {
            return new int[0];
        }


        private static readonly double TOL = 1e-4;

        private static readonly double[,] InvA21 = new double[,] {
             { 7.0/8.0,    0.0,     -15.0/16.0, 0.0,     0.0,     -15.0/16.0 },
             { 0.0,        3.0/4.0, 0.0,        0.0,     0.0,     0.0 },
             { -15.0/16.0, 0.0,     45.0/16.0,  0.0,     0.0,     0.0 },
             { 0.0,        0.0,     0.0,        3.0/4.0, 0.0,     0.0 },
             { 0.0,        0.0,     0.0,        0.0,     9.0/4.0, 0.0 },
             { -15.0/16.0, 0.0,     0.0,        0.0,     0.0,     45.0/16.0 }
        };

        private double[] Heqpol_coefficients(double a, double b, double c) {
            double t = Math.Sqrt(a * a + b * b);
            a = a / t;
            b = b / t;
            c = c / t;
            double am = Math.Abs(a);
            double bm = Math.Abs(b);

            double a2 = a * a;
            double b2 = b * b;
            double c2 = c * c;
            double a3 = a * a * a;
            double b3 = b * b * b;
            double c3 = c * c * c;

            double[] BV = new double[6];
            if (am < bm * TOL) {
                // a = 0, b <> 0
                double abs01 = Math.Abs(b - c);
                double abs02 = Math.Abs(b + c);

                // include 'quad_ph_6_b.f90';
                BV[0] = (2.0 * b - abs01 + abs02) / b;
                BV[1] = 0.0;
                BV[2] = (2.0 * b - abs01 + abs02) / (3.0 * b);
                BV[3] = ((b + c) * abs01 + (b - c) * abs02) / (2.0 * b2);
                BV[4] = 0.0;
                BV[5] = (2.0 * b3 - (b2 + b * c + c2) * abs01 + (b2 - b * c + c2) * abs02) / (3.0 * b3);
            } else if (bm < am * TOL) {
                // a <> 0, b = 0
                double abs01 = Math.Abs(a - c);
                double abs02 = Math.Abs(a + c);

                //  include 'quad_ph_6_a.f90';
                BV[0] = (2.0 * a - abs01 + abs02) / a;
                BV[1] = ((a + c) * abs01 + (a - c) * abs02) / (2.0 * a2);
                BV[2] = (2.0 * a3 - (a2 + a * c + c2) * abs01 + (a2 - a * c + c2) * abs02) / (3.0 * a3);
                BV[3] = 0.0;
                BV[4] = 0.0;
                BV[5] = (2.0 * a - abs01 + abs02) / (3.0 * a);
            } else {
                // a <> 0, b <> 0
                double abs01 = Math.Abs(a - b - c);
                double abs02 = Math.Abs(a + b - c);
                double abs03 = Math.Abs(a - b + c);
                double abs04 = Math.Abs(a + b + c);

                // include 'quad_ph_6_ab.f90';
                BV[0] = 1.0 / (4.0 * a * b) * (
                    8.0 * a * b
                    + (a - b - c) * abs01
                    - (a + b - c) * abs02
                    - a * abs03
                    + b * abs03
                    - c * abs03
                    + a * abs04
                    + b * abs04
                    + c * abs04);
                BV[1] = 1.0 / (12.0 * a2 * b) * (
                    (-2.0 * a2 + a * (b + c) + (b + c) * (b + c)) * abs01
                    + (2.0 * a2 + a * (b - c) - (b - c) * (b - c)) * abs02
                    - 2.0 * a2 * abs03
                    + a * b * abs03
                    + b2 * abs03
                    - a * c * abs03
                    - 2.0 * b * c * abs03
                    + c2 * abs03
                    + 2.0 * a2 * abs04
                    + a * b * abs04
                    - b2 * abs04
                    + a * c * abs04
                    - 2.0 * b * c * abs04
                    - c2 * abs04);
                BV[2] = 1.0 / (24.0 * a3 * b) * (
                    16.0 * a3 * b
                    + (3.0 * a3 - a2 * (b + c) - a * (b + c) * (b + c) - (b + c) * (b + c) * (b + c)) * abs01
                    + (-3.0 * a3 + a * (b - c) * (b - c) - (b - c) * (b - c) * (b - c) + a2 * (-b + c)) * abs02
                    - 3.0 * a3 * abs03
                    + a2 * b * abs03
                    + a * b2 * abs03
                    + b3 * abs03
                    - a2 * c * abs03
                    - 2.0 * a * b * c * abs03
                    - 3.0 * b2 * c * abs03
                    + a * c2 * abs03
                    + 3.0 * b * c2 * abs03
                    - c3 * abs03
                    + 3.0 * a3 * abs04
                    + a2 * b * abs04
                    - a * b2 * abs04
                    + b3 * abs04
                    + a2 * c * abs04
                    - 2.0 * a * b * c * abs04
                    + 3.0 * b2 * c * abs04
                    - a * c2 * abs04
                    + 3.0 * b * c2 * abs04
                    + c3 * abs04);
                BV[3] = 1.0 / (12.0 * a * b2) * (
                    (a2 - 2.0 * b2 + a * (b - 2.0 * c) - b * c + c2) * abs01
                    - (a2 - 2.0 * b2 + b * c + c2 - a * (b + 2.0 * c)) * abs02
                    + a2 * abs03
                    + a * b * abs03
                    - 2.0 * b2 * abs03
                    + 2.0 * a * c * abs03
                    + b * c * abs03
                    + c2 * abs03
                    - a2 * abs04
                    + a * b * abs04
                    + 2.0 * b2 * abs04
                    - 2.0 * a * c * abs04
                    + b * c * abs04
                    - c2 * abs04);
                BV[4] = 1.0 / (48.0 * a2 * b2) * (
                    (-3.0 * a3 + (3.0 * b - c) * (b + c) * (b + c) + a2 * (-3.0 * b + 5.0 * c) + a * (3.0 * b2 + 2.0 * b * c - c2)) * abs01
                    + (3.0 * a3 + (b - c) * (b - c) * (3.0 * b + c) - a2 * (3.0 * b + 5.0 * c) + a * (-3.0 * b2 + 2.0 * b * c + c2)) * abs02
                    + 3.0 * a3 * abs03
                    + 3.0 * a2 * b * abs03
                    - 3.0 * a * b2 * abs03
                    - 3.0 * b3 * abs03
                    + 5.0 * a2 * c * abs03
                    + 2.0 * a * b * c * abs03
                    + 5.0 * b2 * c * abs03
                    + a * c2 * abs03
                    - b * c2 * abs03
                    - c3 * abs03
                    - 3.0 * a3 * abs04
                    + 3.0 * a2 * b * abs04
                    + 3.0 * a * b2 * abs04
                    - 3.0 * b3 * abs04
                    - 5.0 * a2 * c * abs04
                    + 2.0 * a * b * c * abs04
                    - 5.0 * b2 * c * abs04
                    - a * c2 * abs04
                    - b * c2 * abs04
                    + c3 * abs04);
                BV[5] = -1.0 / (24.0 * a * b3) * (
                    -16.0 * a * b3
                    - (a3 - 3.0 * b3 + a2 * (b - 3.0 * c) - b2 * c + b * c2 - c3 + a * (b2 - 2.0 * b * c + 3.0 * c2)) * abs01
                    - (-a3 - 3.0 * b3 + b2 * c + b * c2 + c3 + a2 * (b + 3.0 * c) - a * (b2 + 2.0 * b * c + 3.0 * c2)) * abs02
                    + a3 * abs03
                    + a2 * b * abs03
                    + a * b2 * abs03
                    - 3.0 * b3 * abs03
                    + 3.0 * a2 * c * abs03
                    + 2.0 * a * b * c * abs03
                    + b2 * c * abs03
                    + 3.0 * a * c2 * abs03
                    + b * c2 * abs03
                    + c3 * abs03
                    - a3 * abs04
                    + a2 * b * abs04
                    - a * b2 * abs04
                    - 3.0 * b3 * abs04
                    - 3.0 * a2 * c * abs04
                    + 2.0 * a * b * c * abs04
                    - b2 * c * abs04
                    - 3.0 * a * c2 * abs04
                    + b * c2 * abs04
                    - c3 * abs04);
            }

            //eqcv = matmul(InvA21, BV);
            double[] coefficients = new double[InvA21.GetLength(0)];
            for (int i = 0; i < InvA21.GetLength(0); i++) {
                for (int j = 0; j < InvA21.GetLength(1); j++) {
                    coefficients[i] += InvA21[i, j] * BV[j];
                }
            }
            return coefficients;
        }
    }
}
