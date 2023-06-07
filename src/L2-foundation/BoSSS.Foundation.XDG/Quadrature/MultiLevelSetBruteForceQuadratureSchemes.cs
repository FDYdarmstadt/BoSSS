using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System;
using System.Collections;

namespace BoSSS.Foundation.XDG.Quadrature {

    internal interface IScheme {
        RefElement ReferenceElement { get; }

        QuadRule GetQuadRule(int j);

        void Initialize(int resolution);
    }

    internal class TensorGrid {
        private double[] xTics;

        private double[] yTics;

        public NodeSet VolumeNodes => volumeNodes;

        private NodeSet volumeNodes;

        public NodeSet HorizontalEdgeNodes => horizontalEdgeNodes;

        private NodeSet horizontalEdgeNodes;

        public NodeSet VerticalEdgeNodes => verticalEdgeNodes;

        private NodeSet verticalEdgeNodes;

        public TensorGrid(double[] xTics, double[] yTics) {
            this.xTics = xTics;
            this.yTics = yTics;
            CreateVolumeNodes();
            CreateInnerEdgeNodes();
        }

        private void CreateVolumeNodes() {
            volumeNodes = new NodeSet(Square.Instance, xTics.Length * yTics.Length, 2, true);
            for (int i = 0; i < xTics.Length; ++i) {
                for (int j = 0; j < yTics.Length; ++j) {
                    volumeNodes[VolumeNodeIndice(i, j), 0] = xTics[i];
                    volumeNodes[VolumeNodeIndice(i, j), 1] = yTics[j];
                }
            }
            volumeNodes.LockForever();
        }

        private int VolumeNodeIndice(int i, int j) {
            return i * yTics.Length + j;
        }

        public (int x1, int x2, int y1, int y2) NeighborNodes(int node) {
            int i = node / yTics.Length;
            int j = node % yTics.Length;
            int x1 = VolumeNodeIndice(i - 1, j);
            int x2 = VolumeNodeIndice(i + 1, j);
            int y1 = VolumeNodeIndice(i, j - 1);
            int y2 = VolumeNodeIndice(i, j + 1);

            x1 = CheckBounds(x1);
            x2 = CheckBounds(x2);
            y1 = CheckBounds(y1);
            y2 = CheckBounds(y2);

            return (x1, x2, y1, y2);

            int CheckBounds(int k) {
                if (k > VolumeNodes.NoOfNodes - 1) {
                    return -1;
                } else {
                    return k;
                }
            }
        }

        public (int x1, int x2, int y1, int y2, int x1y1, int x1y2, int x2y1, int x2y2) SixNeighborNodes(int node) {
            int i = node / yTics.Length;
            int j = node % yTics.Length;
            int x1 = VolumeNodeIndice(i - 1, j);
            int x2 = VolumeNodeIndice(i + 1, j);
            int y1 = VolumeNodeIndice(i, j - 1);
            int y2 = VolumeNodeIndice(i, j + 1);

            int x1y1 = VolumeNodeIndice(i - 1, j -1);
            int x1y2 = VolumeNodeIndice(i - 1, j + 1);
            int x2y1 = VolumeNodeIndice(i + 1, j - 1);
            int x2y2 = VolumeNodeIndice(i + 1, j + 1);

            x1 = CheckBounds(x1);
            x2 = CheckBounds(x2);
            y1 = CheckBounds(y1);
            y2 = CheckBounds(y2);
            x1y1 = CheckBounds(x1y1);
            x1y2 = CheckBounds(x1y2);
            x2y1 = CheckBounds(x2y1);
            x2y2 = CheckBounds(x2y2);

            return (x1, x2, y1, y2, x1y1, x1y2, x2y1, x2y2);

            int CheckBounds(int k) {
                if (k > VolumeNodes.NoOfNodes - 1) {
                    return -1;
                } else {
                    return k;
                }
            }
        }

        private void CreateInnerEdgeNodes() {
            verticalEdgeNodes = new NodeSet(Square.Instance, xTics.Length * (yTics.Length - 1), 2, false);
            for (int i = 0; i < xTics.Length - 1; ++i) {
                for (int j = 0; j < yTics.Length; ++j) {
                    verticalEdgeNodes[i * yTics.Length + j, 0] = (xTics[i] + xTics[i + 1]) / 2;
                    verticalEdgeNodes[i * yTics.Length + j, 1] = yTics[j];
                }
            }
            verticalEdgeNodes.LockForever();
            horizontalEdgeNodes = new NodeSet(Square.Instance, (xTics.Length - 1) * yTics.Length, 2, false);
            for (int i = 0; i < xTics.Length; ++i) {
                for (int j = 0; j < yTics.Length - 1; ++j) {
                    horizontalEdgeNodes[i * (yTics.Length - 1) + j, 0] = xTics[i];
                    horizontalEdgeNodes[i * (yTics.Length - 1) + j, 1] = (yTics[j] + yTics[j + 1]) / 2.0;
                }
            }
            horizontalEdgeNodes.LockForever();
        }

        public (int nodeA, int nodeB) NodesOfHorizontalEdge(int edge) {
            int nodeA, nodeB;
            int i = edge / (yTics.Length - 1);
            int j = edge % (yTics.Length - 1);
            nodeA = VolumeNodeIndice(i, j);
            nodeB = VolumeNodeIndice(i, j + 1);
            return (nodeA, nodeB);
        }

        public (int nodeA, int nodeB) NodesOfVerticalEdge(int edge) {
            int nodeA, nodeB;
            int i = edge / yTics.Length;
            int j = edge % yTics.Length;
            nodeA = VolumeNodeIndice(i, j);
            nodeB = VolumeNodeIndice(i + 1, j);
            return (nodeA, nodeB);
        }
    }

    internal class BruteForceVolumeScheme : IScheme {
        private TensorGrid grid;

        private double volumeWeight;

        private Action<int, NodeSet, MultidimensionalArray> phi;

        public RefElement ReferenceElement => Square.Instance;

        public BruteForceVolumeScheme(Action<int, NodeSet, MultidimensionalArray> phi) {
            this.phi = phi;
        }

        public void Initialize(int resolution) {
            grid = new TensorGrid(CenteredLinSpace(resolution), CenteredLinSpace(resolution));
            volumeWeight = 4.0 / (double)(resolution * resolution);
            phiValues = MultidimensionalArray.Create(1, resolution * resolution);
        }

        private static double[] CenteredLinSpace(int resolution) {
            double[] nodes = new double[resolution];
            double increment = 2.0 / resolution;
            double first = increment / 2 - 1.0;
            for (int i = 0; i < resolution; ++i) {
                nodes[i] = increment * i + first;
            }
            return nodes;
        }

        private MultidimensionalArray phiValues;

        public QuadRule GetQuadRule(int cell) {
            phi(cell, grid.VolumeNodes, phiValues);
            BitArray nodeMap = new BitArray(grid.VolumeNodes.NoOfNodes);
            int numberOfVolumeNodes = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (phiValues[0, i] > 0) {
                    nodeMap[i] = true;
                    ++numberOfVolumeNodes;
                }
            }
            QuadRule rule;
            if (numberOfVolumeNodes == 0) {
                rule = QuadRule.CreateEmpty(Square.Instance, 1, 2);
                rule.Nodes.LockForever();
            } else {
                rule = ExtractQuadRule(nodeMap, numberOfVolumeNodes);
            }
            return rule;
        }

        private QuadRule ExtractQuadRule(BitArray nodeMap, int count) {
            QuadRule rule = QuadRule.CreateEmpty(Square.Instance, count, 2, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = grid.VolumeNodes[i, 0];
                    rule.Nodes[j, 1] = grid.VolumeNodes[i, 1];
                    rule.Weights[j] = volumeWeight;
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }
    }

    internal class BruteForceSurfaceScheme : IScheme {
        private double surfaceWeight;

        private TensorGrid grid;

        private FiniteDifferenceOperator differentialOp;

        private Action<int, NodeSet, MultidimensionalArray> phi0;
        private Action<int, NodeSet, MultidimensionalArray> phi1;
        private Func<int, Vector> jacobian;
        private Action<int, NodeSet, MultidimensionalArray> gradient;

        public RefElement ReferenceElement => Square.Instance;

        public BruteForceSurfaceScheme(
            Action<int, NodeSet, MultidimensionalArray> phi0,
            Action<int, NodeSet, MultidimensionalArray> phi1,
            Func<int, Vector> jacobian,
            Action<int, NodeSet, MultidimensionalArray> gradient) {
            this.phi0 = phi0;
            this.phi1 = phi1;
            this.jacobian = jacobian;
            this.gradient = gradient;
        }

        public void Initialize(int resolution) {
            grid = new TensorGrid(CenteredLinSpace(resolution), CenteredLinSpace(resolution));
            surfaceWeight = 2.0 / (double)(resolution);
            phiValues0 = MultidimensionalArray.Create(1, resolution * resolution);
            phiValues1 = MultidimensionalArray.Create(1, resolution * resolution);
            differentialOp = new FiniteDifferenceOperator(phiValues0, grid);
        }

        private static double[] CenteredLinSpace(int resolution) {
            double[] nodes = new double[resolution];
            double increment = 2.0 / resolution;
            double first = increment / 2 - 1.0;
            for (int i = 0; i < resolution; ++i) {
                nodes[i] = increment * i + first;
            }
            return nodes;
        }

        private MultidimensionalArray phiValues0;
        private MultidimensionalArray phiValues1;
        private Vector scales;
        private MultidimensionalArray grad;

        public QuadRule GetQuadRule(int cell) {
            phi0(cell, grid.VolumeNodes, phiValues0);
            phi1(cell, grid.VolumeNodes, phiValues1);
            (BitArray horizontalEdgeMap, int horizontalEdgeCount) = FindHorizontalEdges();
            (BitArray verticalEdgeMap, int verticalEdgeCount) = FindVerticalEdges();
            scales = jacobian(cell);
            QuadRule rule;
            if (horizontalEdgeCount == 0 && verticalEdgeCount == 0) {
                rule = QuadRule.CreateEmpty(Square.Instance, 1, 2);
                rule.Nodes.LockForever();
            } else if (horizontalEdgeCount > verticalEdgeCount) {
                rule = ExtractHorizontalQuadRule(horizontalEdgeMap, horizontalEdgeCount, grid.HorizontalEdgeNodes);
                grad = MultidimensionalArray.Create(1, rule.NoOfNodes * rule.NoOfNodes, 2);
                gradient(cell, rule.Nodes, grad);
                for (int i = 0; i < rule.NoOfNodes; ++i) {
                    Vector gradient = new Vector(grad[0, i, 0], grad[0, i, 1]);
                    rule.Weights[i] *= gradient.Abs() / Math.Abs(gradient[1]);
                }
            } else {
                rule = ExtractVerticalQuadRule(verticalEdgeMap, verticalEdgeCount, grid.VerticalEdgeNodes);
                grad = MultidimensionalArray.Create(1, rule.NoOfNodes * rule.NoOfNodes, 2);
                gradient(cell, rule.Nodes, grad);
                for (int i = 0; i < rule.NoOfNodes; ++i) {
                    Vector gradient = new Vector(grad[0, i, 0], grad[0, i, 1]);
                    rule.Weights[i] *= gradient.Abs() / Math.Abs(gradient[0]);
                }
            }            
            return rule;
        }

        private (BitArray edgeMap, int edgeCount) FindHorizontalEdges() {
            int edgeCount = 0;
            BitArray edgeMap = new BitArray(grid.HorizontalEdgeNodes.NoOfNodes);
            for (int i = 0; i < grid.HorizontalEdgeNodes.NoOfNodes; ++i) {
                (int adjacentA, int adjacentB) = grid.NodesOfHorizontalEdge(i);
                if (phiValues0[0, adjacentA] * phiValues0[0, adjacentB] < 0 && phiValues1[0, adjacentA] * phiValues1[0, adjacentB] < 0) {
                    edgeMap[i] = true;
                    ++edgeCount;
                }
            }
            return (edgeMap, edgeCount);
        }

        private (BitArray edgeMap, int edgeCount) FindVerticalEdges() {
            int edgeCount = 0;
            BitArray edgeMap = new BitArray(grid.VerticalEdgeNodes.NoOfNodes);
            for (int i = 0; i < grid.VerticalEdgeNodes.NoOfNodes; ++i) {
                (int adjacentA, int adjacentB) = grid.NodesOfVerticalEdge(i);
                if (phiValues0[0, adjacentA] * phiValues0[0, adjacentB] < 0 && phiValues1[0, adjacentA] * phiValues1[0, adjacentB] < 0) {
                    edgeMap[i] = true;
                    ++edgeCount;
                }
            }
            return (edgeMap, edgeCount);
        }

        private QuadRule ExtractVerticalQuadRule(BitArray nodeMap, int count, NodeSet nodes) {
            QuadRule rule = QuadRule.CreateEmpty(Square.Instance, count, 2, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = nodes[i, 0];
                    rule.Nodes[j, 1] = nodes[i, 1];
                    rule.Weights[j] = surfaceWeight * scales[0];
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }

        private double VerticalWeight(int edgeNode) {
            (int adjacentA, int adjacentB) = grid.NodesOfVerticalEdge(edgeNode);
            Vector gradiendA = differentialOp.Gradient(adjacentA);
            Vector gradiendB = differentialOp.Gradient(adjacentB);
            Vector grad = (gradiendA + gradiendB) / 2;
            grad[0] *= scales[0];
            grad[1] *= scales[1];
            double weight = grad.Abs() / Math.Abs(grad[0]);
            return weight;
        }

        private QuadRule ExtractHorizontalQuadRule(BitArray nodeMap, int count, NodeSet nodes) {
            QuadRule rule = QuadRule.CreateEmpty(Square.Instance, count, 2, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = nodes[i, 0];
                    rule.Nodes[j, 1] = nodes[i, 1];
                    rule.Weights[j] = surfaceWeight * scales[1];
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }

        private double HorizontalWeight(int edgeNode) {
            (int adjacentA, int adjacentB) = grid.NodesOfHorizontalEdge(edgeNode);
            Vector gradiendA = differentialOp.Gradient(adjacentA);
            Vector gradiendB = differentialOp.Gradient(adjacentB);
            Vector grad = (gradiendA + gradiendB) / 2;
            grad[0] *= scales[0];
            grad[1] *= scales[1];
            double weight = grad.Abs() / Math.Abs(grad[1]);
            return weight;
        }
    }

    internal class FiniteDifferenceOperator {
        private MultidimensionalArray values;

        private TensorGrid grid;

        public FiniteDifferenceOperator(MultidimensionalArray values, TensorGrid grid) {
            this.grid = grid;
            this.values = values;
        }

        public Vector Gradient(int nodeA) {
            (int x1, int x2, int y1, int y2) = grid.NeighborNodes(nodeA);
            double dx = 0;
            int dxCounter = 0;
            if (x1 > -1) {
                dx += FiniteDiff(x1, nodeA, 0);
                ++dxCounter;
            }
            if (x2 > -1) {
                dx += FiniteDiff(nodeA, x2, 0);
                ++dxCounter;
            }

            double dy = 0;
            int dyCounter = 0;
            if (y1 > -1) {
                dy += FiniteDiff(y1, nodeA, 1);
                ++dyCounter;
            }
            if (y2 > -1) {
                dy += FiniteDiff(nodeA, y2, 1);
                ++dyCounter;
            }

            return new Vector(dx / dxCounter, dy / dyCounter);

            double FiniteDiff(int i, int j, int dim) {
                double diff = values[0, j] - values[0, i];
                diff /= (grid.VolumeNodes[j, dim] - grid.VolumeNodes[i, dim]);
                return diff;
            }
        }
    }

    internal class BruteForceEdgeScheme : IScheme {
        public RefElement ReferenceElement => Line.Instance;

        private NodeSet edgeNodes;

        private Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi;

        private double weight;

        public BruteForceEdgeScheme(Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi) {
            this.phi = phi;
        }

        public void Initialize(int resolution) {
            weight = 2.0 / (resolution);
            phiValuesIn = MultidimensionalArray.Create(1, resolution);
            phiValuesOut = MultidimensionalArray.Create(1, resolution);
            CreateEdgeNodes(resolution);
        }

        private void CreateEdgeNodes(int resolution) {
            edgeNodes = new NodeSet(Line.Instance, resolution, 1, true);
            double increment = 2.0 / resolution;
            double first = increment / 2 - 1.0;
            for (int i = 0; i < resolution; ++i) {
                edgeNodes[i, 0] = increment * i + first;
            }
            edgeNodes.LockForever();
        }

        private MultidimensionalArray phiValuesIn;

        private MultidimensionalArray phiValuesOut;

        public QuadRule GetQuadRule(int edge) {
            phi(edge, edgeNodes, phiValuesIn, phiValuesOut);

            BitArray nodeMap = new BitArray(edgeNodes.NoOfNodes);
            int numberOfEdgeNodes = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (phiValuesIn[0, i] > 0) {
                    nodeMap[i] = true;
                    ++numberOfEdgeNodes;
                }
            }
            QuadRule rule;
            if (numberOfEdgeNodes == 0) {
                rule = QuadRule.CreateEmpty(Line.Instance, 1, 1);
                rule.Nodes.LockForever();
            } else {
                rule = ExtractQuadRule(nodeMap, numberOfEdgeNodes);
            }
            return rule;
        }

        private QuadRule ExtractQuadRule(BitArray nodeMap, int count) {
            QuadRule rule = QuadRule.CreateEmpty(ReferenceElement, count, 1, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = edgeNodes[i, 0];
                    rule.Weights[j] = weight;
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }
    }

    internal class BruteForceEdgePointScheme : IScheme {
        public RefElement ReferenceElement => Line.Instance;

        private NodeSet edgeNodes;

        private Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi;

        private Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi1;

        Func<int, double> gram;

        public BruteForceEdgePointScheme(Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi,
            Action<int, NodeSet, MultidimensionalArray, MultidimensionalArray> phi1,
            Func<int, double> gram) {
            this.phi = phi;
            this.phi1 = phi1;
            this.gram = gram;
        }

        public void Initialize(int resolution) {
            phiValuesIn = MultidimensionalArray.Create(1, resolution);
            phiValuesOut = MultidimensionalArray.Create(1, resolution);
            phiValuesIn1 = MultidimensionalArray.Create(1, 1);
            phiValuesOut1 = MultidimensionalArray.Create(1, 1);
            CreateEdgeNodes(resolution);
        }

        private void CreateEdgeNodes(int resolution) {
            edgeNodes = new NodeSet(Line.Instance, resolution, 1, true);
            double increment = 2.0 / (resolution - 1);
            double first = - 1.0;
            for (int i = 0; i < resolution; ++i) {
                edgeNodes[i, 0] = increment * i + first;
            }
            edgeNodes.LockForever();
        }

        private MultidimensionalArray phiValuesIn;

        private MultidimensionalArray phiValuesOut;

        private MultidimensionalArray phiValuesIn1;

        private MultidimensionalArray phiValuesOut1;

        public QuadRule GetQuadRule(int edge) {
            phi(edge, edgeNodes, phiValuesIn, phiValuesOut);

            BitArray nodeMap = new BitArray(edgeNodes.NoOfNodes);
            int numberOfEdgeNodes = 0;
            for (int i = 0; i < nodeMap.Length - 1; ++i) {
                if (phiValuesIn[0, i] * phiValuesIn[0, i+1] < 0) {
                    NodeSet edgeNodes1 = new NodeSet(Line.Instance, 1, 1, false);
                    edgeNodes1[0] = edgeNodes[i];
                    edgeNodes1.LockForever();
                    phi1(edge, edgeNodes1, phiValuesIn1, phiValuesOut1);
                    if(phiValuesIn1[0] > 0) {
                        nodeMap[i] = true;
                        ++numberOfEdgeNodes;
                    }
                }
            }
            QuadRule rule;
            if (numberOfEdgeNodes == 0) {
                rule = QuadRule.CreateEmpty(Line.Instance, 1, 1);
                rule.Nodes.LockForever();
            } else {
                double weight = gram(edge);
                rule = ExtractQuadRule(nodeMap, numberOfEdgeNodes, weight);
            }
            return rule;
        }

        private QuadRule ExtractQuadRule(BitArray nodeMap, int count, double weight) {
            QuadRule rule = QuadRule.CreateEmpty(ReferenceElement, count, 1, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = edgeNodes[i, 0];
                    rule.Weights[j] = weight;
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }
    }

    internal class BruteForceZeroScheme : IScheme {

        private TensorGrid grid;
        private Action<int, NodeSet, MultidimensionalArray> phi0;
        private Func<int, double> gram;

        public RefElement ReferenceElement => Square.Instance;

        ///Zeros must be also only minima of phi
        public BruteForceZeroScheme(
            Action<int, NodeSet, MultidimensionalArray> phi0,
            Func<int, double> gram) {
            this.phi0 = phi0;
            this.gram = gram;
        }

        public void Initialize(int resolution) {
            grid = new TensorGrid(LinSpace(resolution), LinSpace(resolution));
            phiValues0 = MultidimensionalArray.Create(1, resolution * resolution);
        }

        private static double[] LinSpace(int resolution) {
            double[] nodes = new double[resolution ];
            double increment = 2.0 / (resolution - 1);
            double first = - 1.0;
            for (int i = 0; i < resolution; ++i) {
                nodes[i] = increment * i + first;
            }
            return nodes;
        }

        private MultidimensionalArray phiValues0;

        public QuadRule GetQuadRule(int cell) {
            phi0(cell, grid.VolumeNodes, phiValues0);

            (BitArray nodeMap, int numberOfVolumeNodes) = FindZeroNodes();
            
            QuadRule rule;
            if (numberOfVolumeNodes == 0) {
                rule = QuadRule.CreateEmpty(Square.Instance, 1, 2);
                rule.Nodes.LockForever();
            } else {
                if(numberOfVolumeNodes != 1) {
                    Console.WriteLine($"Contact Line Operator Warning: More than one levelSetintersection in Cell{cell}");
                }
                double weight = gram(cell);
                rule = ExtractQuadRule(nodeMap, numberOfVolumeNodes, weight);
            }
            return rule;
        }
        
        private (BitArray, int) FindZeroNodes() {
            BitArray nodeMap = new BitArray(grid.VolumeNodes.NoOfNodes);
            int numberOfVolumeNodes = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                (int a, int b, int c, int d, int e, int f, int g, int h) = grid.SixNeighborNodes(i);
                if(a > -1 && b > -1 && c > -1 && d > -1 && e > -1 && f > -1 && g > -1 && h > -1) {
                    bool isMinimum = true;
                    if (phiValues0[0, i] > phiValues0[0, a]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, b]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, c]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, d]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, e]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, f]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, g]) {
                        isMinimum = false;
                    }
                    if (phiValues0[0, i] > phiValues0[0, h]) {
                        isMinimum = false;
                    }
                    if (isMinimum) {
                        if (Math.Abs(phiValues0[0, i]) < 1e-4) {
                            nodeMap[i] = true;
                            ++numberOfVolumeNodes;
                        }
                    }
                }
            }
            return (nodeMap, numberOfVolumeNodes);
        }

        private QuadRule ExtractQuadRule(BitArray nodeMap, int count, double weight) {
            QuadRule rule = QuadRule.CreateEmpty(Square.Instance, count, 2, true);
            int j = 0;
            for (int i = 0; i < nodeMap.Length; ++i) {
                if (nodeMap[i]) {
                    rule.Nodes[j, 0] = grid.VolumeNodes[i, 0];
                    rule.Nodes[j, 1] = grid.VolumeNodes[i, 1];
                    rule.Weights[j] = weight;
                    ++j;
                }
            }
            rule.Nodes.LockForever();
            return rule;
        }
    }
}