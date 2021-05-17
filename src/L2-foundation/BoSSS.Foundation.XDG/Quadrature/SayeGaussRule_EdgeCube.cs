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

        public SayeGaussRule_EdgeCube(LevelSetTracker.LevelSetData _lsData,
            IRootFindingAlgorithm rootFinder,
            QuadratureMode mode) : base(_lsData, rootFinder, mode) {
            this.mode = mode;
        }

        public CellBoundaryQuadRule EvaluateEdges(int Cell) {
            LinearSayeSpace<Cube>[] edgeSubspace = CreateEdgeSubspaces(false);
            QuadRule[] edgeRules = new QuadRule[edgeSubspace.Length];
            for (int i = 0; i < edgeSubspace.Length; ++i) {
                edgeSubspace[i].Reset();
                QuadRule edgeRule = Evaluate(Cell, edgeSubspace[i]);
                edgeRules[i] = edgeRule;
            }
            CellBoundaryQuadRule combinedRules = Combine(edgeRules);
            return combinedRules;
        }

        LinearSayeSpace<Cube>[] CreateEdgeSubspaces(bool surface) {
            //6 Faces in a cube
            int domainSign = 1;
            if (mode == QuadratureMode.NegativeVolume) {
                domainSign = -1;
            }
            LinearSayeSpace<Cube>[] edges = new LinearSayeSpace<Cube>[6];
            LinearPSI<Cube> psi = new LinearPSI<Cube>(Cube.Instance);
            for (int i = 0; i < 3; ++i) {
                for(int j = 0; j < 2; ++j) {
                    LinearSayeSpace<Cube> fullCube = new LinearSayeSpace<Cube>(Cube.Instance, surface);
                    fullCube.RemoveDimension(i);
                    LinearPSI<Cube> edge = psi.ReduceDim(i, 1 - 2 * j);
                    fullCube.PsiAndS.Add(new Tuple<LinearPSI<Cube>, int>(edge, domainSign));
                    edges[2 * i + j] = fullCube;
                }
            }
            return edges;
        }

        CellBoundaryQuadRule Combine(QuadRule[] rules) {
            int numberOfNodes = 0;
            foreach (QuadRule rule in rules) {
                numberOfNodes += rule.NoOfNodes;
            }
            CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(Cube.Instance, numberOfNodes, 3, 6);
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
    }
}
