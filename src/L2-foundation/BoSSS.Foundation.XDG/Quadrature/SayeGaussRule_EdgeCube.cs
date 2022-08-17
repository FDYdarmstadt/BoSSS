using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;
using ilPSP;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.XDG.Quadrature {
    class SayeGaussRule_EdgeCube : SayeGaussRule_Cube, ISayeGaussEdgeRule {

        QuadratureMode mode;

        IGridData grid;

        public SayeGaussRule_EdgeCube(LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm rootFinder,
            QuadratureMode mode) : base(_lsData, rootFinder, mode) {
            this.mode = mode;
            grid = _lsData.GridDat;
        }

        LinearSayeSpace<Cube> activeSpace;

        public CellBoundaryQuadRule EvaluateEdges(int Cell) {
            LinearSayeSpace<Cube>[] edgeSubspace = CreateEdgeSubspaces(mode == QuadratureMode.Surface);
            QuadRule[] edgeRules = new QuadRule[edgeSubspace.Length];
            for (int i = 0; i < edgeSubspace.Length; ++i) {
                edgeSubspace[i].Reset();
                activeSpace = edgeSubspace[i];
                QuadRule edgeRule = Evaluate(Cell, edgeSubspace[i]);
                edgeRules[i] = edgeRule;
            }
            CellBoundaryQuadRule combinedRules = Combine(edgeRules);
            return combinedRules;
        }

        //Map from bosss-cube indices to saye edges (subspaces)
        //RowNumber: index of edge in bosss-cube 
        //Each Row: first entry: height direction, second entry: offset from center
        static int[,] subspace2CubeMap = new int[6,2]{ 
            { 0, -1 },
            { 0,  1 },
            { 1,  1 },
            { 1, -1},
            { 2, 1 },
            { 2, -1 }};

        LinearSayeSpace<Cube>[] CreateEdgeSubspaces(bool isSurface) {
            //6 Faces in a cube
            int domainSign = 1;
            if (mode == QuadratureMode.NegativeVolume) {
                domainSign = -1;
            }
            LinearSayeSpace<Cube>[] edges = new LinearSayeSpace<Cube>[6];
            LinearPSI<Cube> psi = new LinearPSI<Cube>(Cube.Instance);
            for (int i = 0; i < 6; ++i) {
                LinearSayeSpace<Cube> fullCube = new LinearSayeSpace<Cube>(Cube.Instance, isSurface);
                fullCube.RemoveDimension(subspace2CubeMap[i, 0]);
                LinearPSI<Cube> edge = psi.ReduceDim(subspace2CubeMap[i,0], subspace2CubeMap[i, 1]);
                fullCube.PsiAndS.Add(new Tuple<LinearPSI<Cube>, int>(edge, domainSign));
                edges[i] = fullCube;
            }
            return edges;
        }

        CellBoundaryQuadRule Combine(QuadRule[] rules) {
            int numberOfNodes = 0;
            foreach (QuadRule rule in rules) {
                numberOfNodes += rule.NoOfNodes;
            }
            CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(Cube.Instance, numberOfNodes, 3, 6, true);
            int subarrayPointer = 0;
            for (int i = 0; i < rules.Length; ++i) {
                int subNumberOfNodes = rules[i].NoOfNodes;
                combinedRule.NumbersOfNodesPerFace[i] = subNumberOfNodes;
                for(int j = 0; j < subNumberOfNodes; ++j) {
                    for(int d = 0; d < 3; ++d) {
                        combinedRule.Nodes[subarrayPointer + j, d] = rules[i].Nodes[j, d];
                    }
                    combinedRule.Weights[subarrayPointer + j] = rules[i].Weights[j];
                }
                subarrayPointer += subNumberOfNodes;
            }
            combinedRule.Nodes.LockForever();
            return combinedRule;
        }

        void RestrictToActiveSpace(MultidimensionalArray gradient) {
            for(int i = 0; i < 3; ++i) {
                if (!activeSpace.DimActive(i)) {
                    gradient[i] = 0;
                }
            }
        }

        protected override SayeQuadRule BuildSurfaceQuadRule(MultidimensionalArray X, double X_weight, int heightDirection, int cell) {
            double weight = X_weight;

            NodeSet node = new NodeSet(RefElement, X.To2DArray(), true);
            MultidimensionalArray gradient = ReferenceGradient(node, cell);
            RestrictToActiveSpace(gradient);

            weight *= gradient.L2Norm() / Math.Abs(gradient[heightDirection]);

            MultidimensionalArray jacobian = grid.Jacobian.GetValue_Cell(node, cell, 1).ExtractSubArrayShallow(0, 0, -1, -1);
            //Scale weight
            if (IsScalingMatrix(jacobian)) {
                weight /= jacobian[heightDirection, heightDirection];
            } else {
                throw new NotImplementedException("To do");
            }

            MultidimensionalArray weightArr = new MultidimensionalArray(1);
            weightArr.Allocate(1);
            weightArr[0] = weight;
            return new SayeQuadRule(node, weightArr, RefElement);
        }
    }
}
