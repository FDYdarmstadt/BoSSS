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
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using ilPSP.Tracing;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    class LevelSetEdgeSurfaceQuadRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {

        private LevelSetTracker tracker;

        private IQuadRuleFactory<CellEdgeBoundaryQuadRule> edgeRuleFactory;

        private CellBoundaryQuadRule baseRule;

        private int currentOrder = -1;

        private DivergenceFreeFaceBasis[] phiBasis;

        private JumpTypes jumpType;

        private Dictionary<int, CellBoundaryQuadRule> cache =
            new Dictionary<int, CellBoundaryQuadRule>();

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return cache.Keys.ToArray();
        }

        private int levelSetIndex;

        LevelSetTracker.LevelSetData LevelSetData;

        /// <summary>
        /// Ctor
        /// </summary>
        public LevelSetEdgeSurfaceQuadRuleFactory(LevelSetTracker.LevelSetData lsData,  IQuadRuleFactory<CellEdgeBoundaryQuadRule> edgeRuleFactory, JumpTypes jumpType) {
            this.tracker = lsData.Tracker;
            this.levelSetIndex = lsData.LevelSetIndex;
            this.edgeRuleFactory = edgeRuleFactory;
            this.jumpType = jumpType;
            this.LevelSetData = lsData;
        }

        #region IQuadRuleFactory<CellBoundaryQuadRule> Members

        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            using (var tr = new FuncTrace()) {
                if (mask != null && mask is CellMask == false) {
                    throw new ArgumentException("Cell mask required", "mask");
                }

                Stopwatch totalTimer = new Stopwatch();
                Stopwatch projectionTimer = new Stopwatch();
                Stopwatch optimizationTimer = new Stopwatch();

                totalTimer.Start();
                if (order != currentOrder) {
                    cache.Clear();
                    SwitchOrder(order);
                }

                var result = new List<ChunkRulePair<CellBoundaryQuadRule>>(mask.NoOfItemsLocally);
                foreach (Chunk chunk in mask) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int cell = chunk.i0 + i;

                        if (cache.ContainsKey(cell)) {
                            result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                Chunk.GetSingleElementChunk(cell),
                                cache[cell]));
                            continue;
                        }

                        optimizationTimer.Start();
                        CellBoundaryQuadRule optimizedRule = GetOptimizedRule(chunk.i0 + i, order);
                        optimizationTimer.Stop();

                        cache.Add(cell, optimizedRule);
                        result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                            Chunk.GetSingleElementChunk(i + chunk.i0), optimizedRule));
                    }
                }
                totalTimer.Stop();

                double totalTicks = (double)totalTimer.ElapsedTicks;
                double percentageProjection = Math.Round(projectionTimer.ElapsedTicks / totalTicks * 100, 2);
                double percentageOptimization = Math.Round(optimizationTimer.ElapsedTicks / totalTicks * 100, 2);
                tr.Info("Percentage spent on projection to surface: " + percentageProjection.ToString(NumberFormatInfo.InvariantInfo) + "%");
                tr.Info("Percentage spent on optimization: " + percentageOptimization.ToString(NumberFormatInfo.InvariantInfo) + "%");

                return result;
            }
        }

        public RefElement RefElement {
            get {
                return tracker.GridDat.Grid.RefElements[0];
            }
        }

        #endregion

        private void SwitchOrder(int order) {
            currentOrder = order;

            int noOfEdges = tracker.GridDat.Grid.RefElements[0].NoOfFaces;

            phiBasis = new DivergenceFreeFaceBasis[noOfEdges];
            for (int e = 0; e < noOfEdges; e++) {
                phiBasis[e] = new DivergenceFreeFaceBasis(tracker.GridDat, this.RefElement, order, e);
            }

            int noOfPhis = GetNumberOfPhis();

            int minOrder = 1;
            double safetyFactor = 1.6;
            while (RefElement.FaceRefElement.GetQuadratureRule(minOrder).NoOfNodes < safetyFactor * noOfPhis) {
                minOrder += 1;
            }

            baseRule = new CellBoundaryFromEdgeRuleFactory<CellBoundaryQuadRule>(
                tracker.GridDat,
                RefElement,
                new FixedRuleFactory<QuadRule>(RefElement.FaceRefElement.GetQuadratureRule(minOrder))).
                GetQuadRuleSet(new CellMask(tracker.GridDat, Chunk.GetSingleElementChunk(0)), -1).
                First().Rule;
        }

        protected CellBoundaryQuadRule GetOptimizedRule(int cell, int order) {
            using (var tr = new FuncTrace()) {
                int numberOfPhis = GetNumberOfPhis();
                int noOfEdges = tracker.GridDat.Grid.RefElements[0].NoOfFaces;

                PhiQuadrature quadrature = new PhiQuadrature(tracker, this, order, cell);
                quadrature.Execute();
                double[,] quadResults = quadrature.Results;

                CellBoundaryQuadRule optimizedRule = new CellBoundaryQuadRule() {
                    NumbersOfNodesPerFace = baseRule.NumbersOfNodesPerFace.CloneAs()
                };

                for (int e = 0; e < noOfEdges; e++) {
                    bool edgeIsCut = false;
                    for (int i = 0; i < numberOfPhis; i++) {
                        edgeIsCut |= quadResults[e, i].Abs() > 1e-9;
                    }

                    if (!edgeIsCut) {
                        optimizedRule.NumbersOfNodesPerFace[e] = 0;
                    }
                }

                optimizedRule.Weights = MultidimensionalArray.Create(
                    optimizedRule.NumbersOfNodesPerFace.Sum());
                optimizedRule.Nodes = new NodeSet(
                    this.RefElement,
                    optimizedRule.Weights.Length,
                    tracker.GridDat.SpatialDimension);
                //optimizedRule.Nodes.LockForever();

                if (optimizedRule.NoOfNodes == 0) {
                    var emptyRule = CellBoundaryQuadRule.CreateEmpty(
                        RefElement, 1, tracker.GridDat.SpatialDimension, noOfEdges);
                    emptyRule.Nodes.LockForever();
                    return emptyRule;
                }

                for (int e = 0; e < noOfEdges; e++) {
                    int noOfNodesOnEdge = optimizedRule.NumbersOfNodesPerFace[e];
                    if (noOfNodesOnEdge == 0) {
                        continue;
                    }

                    CopyNodes(baseRule, optimizedRule, e);
                }
                optimizedRule.Nodes.LockForever();

                if (tracker.GridDat.Cells.Cells2Edges[cell].Length != noOfEdges) {
                    throw new NotImplementedException("Hanging nodes not yet supported");
                }

                int noOfProcessedNodes = 0;
                for (int e = 0; e < noOfEdges; e++) {
                    int noOfNodesOnEdge = optimizedRule.NumbersOfNodesPerFace[e];
                    if (noOfNodesOnEdge == 0) {
                        continue;
                    }

                    MultidimensionalArray normals = EvaluateRefNormalsOnEdge(this.LevelSetData, cell, optimizedRule, e);
                    //MultidimensionalArray metrics = GetMetricTermsOnEdge(tracker, levelSetIndex, optimizedRule, cell, e);

                    //lh = tracker.GridDat.NSC.LockNodeSetFamily(tracker.GridDat.NSC.CreateContainer(
                    //    optimizedRule.Nodes.ExtractSubArrayShallow(
                    //        new int[] { noOfProcessedNodes, 0 },
                    //        new int[] { noOfProcessedNodes + noOfNodesOnEdge - 1, optimizedRule.SpatialDim - 1 }),
                    //    0,
                    //    0.0));
                    NodeSet irgendwelcheNodes = new NodeSet(this.RefElement,
                        optimizedRule.Nodes.ExtractSubArrayShallow(
                           new int[] { noOfProcessedNodes, 0 },
                           new int[] { noOfProcessedNodes + noOfNodesOnEdge - 1, optimizedRule.SpatialDim - 1 }));
                    MultidimensionalArray phiValues = EvaluatePhis(irgendwelcheNodes, cell, e);

                    double[] matrix = new double[numberOfPhis * noOfNodesOnEdge];
                    // Additional space required by Fortran routine
                    double[] rhs = new double[Math.Max(noOfNodesOnEdge, numberOfPhis)];

                    for (int j = 0; j < numberOfPhis; j++) {
                        rhs[j] = quadResults[e, j];

                        for (int i = 0; i < noOfNodesOnEdge; i++) {
                            // Beware of Fortran order!
                            int index = i * numberOfPhis + j;

                            for (int d = 0; d < tracker.GridDat.SpatialDimension; d++) {
                                matrix[index] += phiValues[i, j, d] * normals[i, d];
                            }
                        }
                    }

                    LAPACK.F77_LAPACK.DGELSY(numberOfPhis, noOfNodesOnEdge, matrix, rhs, 1, 1e-12);

                    int edge = Math.Abs(tracker.GridDat.Cells.Cells2Edges[cell][e]) - 1;
                    double maxWeight = 0.0;
                    for (int i = 0; i < noOfNodesOnEdge; i++) {
                        //optimizedRule.Weights[noOfProcessedNodes + i] = rhs[i] / metrics[i];
                        optimizedRule.Weights[noOfProcessedNodes + i] = rhs[i];
                        maxWeight = Math.Max(optimizedRule.Weights[noOfProcessedNodes + i].Abs(), maxWeight);
                    }
                    noOfProcessedNodes += noOfNodesOnEdge;

                    if (maxWeight > 4.0 * tracker.GridDat.Edges.h_max_Edge[edge]) {
                        tr.Info(String.Format(
                            "Warning: Abnormally large integration weight detected"
                            + " for level set edge surface integral on edge {0} of cell {1} "
                            + " (|w| = {2}). This may indicate a loss of"
                            + " integration accuracy.",
                            e,
                            cell,
                            maxWeight));
                    }
                }

                return optimizedRule;
            }
        }

        public static MultidimensionalArray EvaluateRefNormalsOnEdge(LevelSetTracker.LevelSetData lsData, int cell, CellBoundaryQuadRule rule, int localEdge) {
            if (rule.NumbersOfNodesPerFace[localEdge] == 0) {
                return null;
            }

            MultidimensionalArray gradients = EvaluateLevelSetReferenceGradientsOnEdge(lsData, cell, rule, localEdge);
            int noOfNodes = gradients.GetLength(0);
            int D = lsData.Tracker.GridDat.SpatialDimension;

            MultidimensionalArray normals = MultidimensionalArray.Create(noOfNodes, D);
            for (int i = 0; i < noOfNodes; i++) {
                double norm = 0.0;
                for (int d = 0; d < D; d++) {
                    norm += gradients[i, d] * gradients[i, d];
                    normals[i, d] = gradients[i, d];
                }

                if (norm == 0.0) {
                    continue;
                }

                norm = Math.Sqrt(norm);
                for (int d = 0; d < D; d++) {
                    normals[i, d] /= norm;
                }
            }

            return normals;
        }

        private static MultidimensionalArray EvaluateLevelSetGradientsOnEdge(LevelSetTracker.LevelSetData lsData, int cell, CellBoundaryQuadRule rule, int localEdge) {
            var gdat = lsData.Tracker.GridDat;
            int edge = Math.Abs(gdat.Cells.Cells2Edges[cell][localEdge]) - 1;
            int D = gdat.SpatialDimension;
            MultidimensionalArray volumeGradients = lsData.GetLevelSetGradients(rule.Nodes, cell, 1);

            int noOfNodes = rule.NumbersOfNodesPerFace[localEdge];
            if (noOfNodes == 0) {
                return null;
            }

            int offset = rule.NumbersOfNodesPerFace.Take(localEdge).Sum();
            MultidimensionalArray normals = MultidimensionalArray.Create(noOfNodes, D);
            var NormalsForAffine = gdat.Edges.NormalsForAffine;
            for (int i = 0; i < noOfNodes; i++) {
                double normalComponent = 0.0;
                for (int d = 0; d < D; d++) {
                    normalComponent += volumeGradients[0, offset + i, d] * NormalsForAffine[edge, d];
                }

                // Subtract normal component
                for (int d = 0; d < D; d++) {
                    normals[i, d] = volumeGradients[0, offset + i, d]
                        - normalComponent * NormalsForAffine[edge, d];
                }
            }

            return normals;
        }

        private static MultidimensionalArray EvaluateLevelSetReferenceGradientsOnEdge(LevelSetTracker.LevelSetData LevelSetData, int cell, CellBoundaryQuadRule rule, int localEdge) {
            int D = LevelSetData.Tracker.GridDat.SpatialDimension;
            MultidimensionalArray volumeGradients = LevelSetData.GetLevelSetReferenceGradients(rule.Nodes, cell, 1);
            MultidimensionalArray edgeNormals = LevelSetData.Tracker.GridDat.Grid.RefElements[0].FaceNormals;

            int noOfNodes = rule.NumbersOfNodesPerFace[localEdge];
            if (noOfNodes == 0) {
                return null;
            }

            int offset = rule.NumbersOfNodesPerFace.Take(localEdge).Sum();
            MultidimensionalArray normals = MultidimensionalArray.Create(noOfNodes, D);
            for (int i = 0; i < noOfNodes; i++) {
                double normalComponent = 0.0;
                for (int d = 0; d < D; d++) {
                    normalComponent += volumeGradients[
                        0, offset + i, d] * edgeNormals[localEdge, d];
                }

                // Subtract normal component
                for (int d = 0; d < D; d++) {
                    normals[i, d] = volumeGradients[0, offset + i, d] - normalComponent * edgeNormals[localEdge, d];
                }
            }

            return normals;
        }

        public static MultidimensionalArray GetMetricTermsOnEdge(LevelSetTracker.LevelSetData LevelSetData, int levSetIndex, CellBoundaryQuadRule rule, int cell, int localEdge) {
            if (rule.NumbersOfNodesPerFace[localEdge] == 0) {
                return null;
            }

            MultidimensionalArray physGradients = EvaluateLevelSetGradientsOnEdge(LevelSetData, cell, rule, localEdge);
            MultidimensionalArray refGradients = EvaluateLevelSetReferenceGradientsOnEdge(LevelSetData, cell, rule, localEdge);
            int noOfNodes = physGradients.GetLength(0);
            int D = LevelSetData.Tracker.GridDat.SpatialDimension;
            MultidimensionalArray result = MultidimensionalArray.Create(noOfNodes);

            int edge = Math.Abs(LevelSetData.Tracker.GridDat.Cells.Cells2Edges[cell][localEdge]) - 1;

            var JacobiDet = LevelSetData.Tracker.GridDat.Cells.JacobiDet;
            var SqrtGramian = LevelSetData.Tracker.GridDat.Edges.SqrtGramian;
            for (int j = 0; j < noOfNodes; j++) {
                double normPhys = 0.0;
                double normRef = 0.0;

                for (int d = 0; d < D; d++) {
                    normPhys += physGradients[j, d] * physGradients[j, d];
                    normRef += refGradients[j, d] * refGradients[j, d];
                }

                result[j] = Math.Sqrt(normRef / normPhys) / SqrtGramian[edge] *
                    Math.Sqrt(JacobiDet[cell]);
                //    result[j] = Math.Sqrt(normRef / normPhys) / tracker.GridDat.Edges.SqrtGramian[edge] / tracker.Ctx.GridDat.OneOverSqrt_AbsDetTransformation[cell];
            }

            return result;
        }

        private static void CopyNodes(CellBoundaryQuadRule source, CellBoundaryQuadRule target, int localEdge) {
            Debug.Assert(
                source.NumbersOfNodesPerFace[localEdge] == target.NumbersOfNodesPerFace[localEdge],
                "Source and target are incompatible");

            int offsetSource = 0;
            int offsetTarget = 0;
            for (int e = 0; e < localEdge; e++) {
                offsetSource += source.NumbersOfNodesPerFace[e];
                offsetTarget += target.NumbersOfNodesPerFace[e];
            }

            for (int i = 0; i < source.NumbersOfNodesPerFace[localEdge]; i++) {
                for (int d = 0; d < source.SpatialDim; d++) {
                    target.Nodes[offsetTarget + i, d] = source.Nodes[offsetSource + i, d];
                }
            }
        }

        private MultidimensionalArray EvaluatePhis(NodeSet NS, int cell, int localEdge) {
            int noOfNodes = NS.NoOfNodes;
            int D = tracker.GridDat.SpatialDimension;
            int noOfPhis = GetNumberOfPhis();

            MultidimensionalArray phiValues = phiBasis[localEdge].Evaluate(NS);
            MultidimensionalArray transformedPhiValues = MultidimensionalArray.Create(
                noOfNodes * noOfPhis, D);
            tracker.GridDat.Grid.RefElements[0].TransformFaceVectors(
                localEdge,
                phiValues.ResizeShallow(noOfNodes * noOfPhis, RefElement.FaceRefElement.SpatialDimension),
                transformedPhiValues);

            return transformedPhiValues.ResizeShallow(noOfNodes, noOfPhis, D);
        }

        private int GetNumberOfPhis() {
            return phiBasis[0].Count / RefElement.FaceRefElement.SpatialDimension;
        }

        class PhiQuadrature : CellBoundaryQuadrature<CellEdgeBoundaryQuadRule> {

            private LevelSetTracker tracker;

            private LevelSetEdgeSurfaceQuadRuleFactory owner;

            public double[,] Results;

            public PhiQuadrature(LevelSetTracker tracker, LevelSetEdgeSurfaceQuadRuleFactory owner, int maxPhiDegree, int element)
                : base(
                    new int[] { owner.GetNumberOfPhis() },
                    tracker.GridDat,
                    new CellEdgeBoundaryQuadratureScheme(
                        false,
                        owner.edgeRuleFactory,
                        new CellMask(tracker.GridDat, Chunk.GetSingleElementChunk(element)))
                .Compile(tracker.GridDat, maxPhiDegree),
                CoordinateSystem.Reference) {
                this.tracker = tracker;
                this.owner = owner;
            }

            public override void Execute() {
                Results = new double[gridData.iGeomCells.RefElements[0].NoOfFaces, IntegralCompDim[0]];
                base.Execute();
            }

            /*
            protected override NodeSetController.NodeSetContainer[] CreateNodeSetFamily(MultidimensionalArray NodesUntransformed, int iKref) {
                int noOfRelevantEdges = CurrentRule.NumbersOfNodesPerFace.Count(n => n > 0);
                NodeSetController.NodeSetContainer[] family =
                    new NodeSetController.NodeSetContainer[1 + noOfRelevantEdges];

                family[0] = gridData.NSC.CreateContainer(NodesUntransformed, iKref, 0.0);

                int noOfProcessedNodes = 0;
                int nodeSetIndex = 1;
                for (int e = 0; e < tracker.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                    int noOfNodesOnEdge = CurrentRule.NumbersOfNodesPerFace[e];
                    if (noOfNodesOnEdge == 0) {
                        continue;
                    }

                    family[nodeSetIndex] = gridData.NSC.CreateContainer(
                        CurrentRule.Nodes.ExtractSubArrayShallow(
                            new int[] { noOfProcessedNodes, 0 },
                            new int[] { noOfProcessedNodes + noOfNodesOnEdge - 1, CurrentRule.SpatialDim - 1 }),
                        iKref,
                        0.0);
                    noOfProcessedNodes += noOfNodesOnEdge;
                    nodeSetIndex++;
                }

                return family;
            }

             */

            NodeSet[] family;

            protected override void QuadNodesChanged(NodeSet newNodes) {
                Debug.Assert(object.ReferenceEquals(newNodes, this.CurrentRule.Nodes));

                int noOfRelevantEdges = CurrentRule.NumbersOfNodesPerFace.Count(n => n > 0);
                family = new NodeSet[1 + noOfRelevantEdges];

                family[0] = newNodes;

                int noOfProcessedNodes = 0;
                int nodeSetIndex = 1;
                for (int e = 0; e < tracker.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                    int noOfNodesOnEdge = CurrentRule.NumbersOfNodesPerFace[e];
                    if (noOfNodesOnEdge == 0) {
                        continue;
                    }

                    family[nodeSetIndex] = new NodeSet(
                        newNodes.RefElement,
                        CurrentRule.Nodes.ExtractSubArrayShallow(
                            new int[] { noOfProcessedNodes, 0 },
                            new int[] { noOfProcessedNodes + noOfNodesOnEdge - 1, CurrentRule.SpatialDim - 1 }));
                    noOfProcessedNodes += noOfNodesOnEdge;
                    nodeSetIndex++;
                }
            }



            protected override void Evaluate(int i0, int Length, CellEdgeBoundaryQuadRule qr, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = qr.Nodes;
                RefElement Kref = this.GridDat.iGeomCells.GetRefElement(i0);
                Debug.Assert(object.ReferenceEquals(Kref, owner.RefElement));

                int noOfFaces = Kref.NoOfFaces;
                int noOfFacesOfFace = Kref.FaceRefElement.NoOfFaces;
                int D = gridData.SpatialDimension;
                MultidimensionalArray edgeOfEdgeNormals = MultidimensionalArray.Create(
                    1, owner.RefElement.FaceRefElement.SpatialDimension);
                MultidimensionalArray edgeNormals = MultidimensionalArray.Create(1, D);

                if (Length > 1 || D != 3) {
                    throw new NotSupportedException();
                }

                Debug.Assert(object.ReferenceEquals(base.CurrentRule.Nodes, QuadNodes));
                Debug.Assert(object.ReferenceEquals(this.family[0], QuadNodes));

                Debug.Assert(
                    owner.RefElement.FaceRefElement.FaceRefElement is Line,
                    "Edge of edge _must_ be a line");

                for (int i = 0; i < Length; i++) {
                    int nodeIndex = -1;

                    int nodeSetIndex = 1;
                    for (int e = 0; e < noOfFaces; e++) {
                        int noOfNodesOnEdge = CurrentRule.NumbersOfNodesPerFace[e];
                        if (noOfNodesOnEdge == 0) {
                            continue;
                        }

                        //int edge = Math.Abs(gridData.Cells.Cells2Edges[i0 + i][e]) - 1;

                        MultidimensionalArray phiValues = owner.EvaluatePhis(this.family[nodeSetIndex], i0, e);

                        int noOfProcessedNodesOnEdge = 0;
                        for (int ee = 0; ee < noOfFacesOfFace; ee++) {
                            for (int d = 0; d < owner.RefElement.FaceRefElement.SpatialDimension; d++) {
                                edgeOfEdgeNormals[0, d] = owner.RefElement.FaceRefElement.FaceNormals[ee, d];
                            }

                            owner.tracker.GridDat.Grid.RefElements[0].TransformFaceVectors(
                                e, edgeOfEdgeNormals, edgeNormals);

                            int noOfNodesOnEdgeOfEdge = CurrentRule.NumbersOfNodesPerFaceOfFace[e, ee];
                            for (int j = 0; j < noOfNodesOnEdgeOfEdge; j++) {
                                nodeIndex++;

                                for (int k = 0; k < IntegralCompDim[0]; k++) {
                                    for (int d = 0; d < D; d++) {
                                        EvalResult[i, nodeIndex, k] +=
                                            phiValues[noOfProcessedNodesOnEdge + j, k, d] * edgeNormals[0, d];
                                    }
                                }
                            }

                            noOfProcessedNodesOnEdge += noOfNodesOnEdgeOfEdge;
                        }

                        nodeSetIndex++;
                    }
                }
            }

            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int e = 0; e < tracker.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                    for (int k = 0; k < IntegralCompDim[0]; k++) {
                        Results[e, k] = ResultsOfIntegration[0, e, k];
                    }
                }
            }

            protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {
                throw new NotSupportedException("special treatment for this class");
            }
        }
    }
}
