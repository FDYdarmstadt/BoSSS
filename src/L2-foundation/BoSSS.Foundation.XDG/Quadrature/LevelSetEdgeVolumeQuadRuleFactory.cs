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

    /// <summary>
    /// Uses the hierarchical moment-fitting (HMF) strategy to compute
    /// two-dimensional volume integrals over the boundaries of cells that
    /// are intersected by the level set. In other words: If the (planar) edge
    /// of a three-dimensional element is intersected by the level set, this
    /// factory applies the HMF strategy in order to create accurate quadrature
    /// rules for the integration over the negative/positive (depending on the
    /// jump type, see
    /// <see cref="LevelSetEdgeVolumeQuadRuleFactory.LevelSetEdgeVolumeQuadRuleFactory"/>)
    /// level set region.
    /// </summary>
    /// <remarks>
    /// For details about the algorithm, see
    /// Mueller, Kummer &amp; Oberlack: Highly accurate surface and volume
    /// integration on implicit domains by means of moment-fitting (2013)
    /// </remarks>
    public class LevelSetEdgeVolumeQuadRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {

        /// <summary>
        /// See
        /// <see cref="LevelSetEdgeVolumeQuadRuleFactory.LevelSetEdgeVolumeQuadRuleFactory"/>
        /// </summary>
        private JumpTypes jumpType;

        ///// <summary>
        ///// Tracks the level set
        ///// </summary>
        //private LevelSetTracker tracker;

        /// <summary>
        /// Factory for one-dimensional integration rules for the edges of the
        /// edges of the three-dimensional volume element
        /// </summary>
        private IQuadRuleFactory<CellEdgeBoundaryQuadRule> CoFaceQuadRuleFactory;

        /// <summary>
        /// Factory for the one-dimensional integration over the intersections
        /// of the zero iso-contour of the level set function with the edges of
        /// the three-dimensional volume element
        /// </summary>
        private LevelSetEdgeSurfaceQuadRuleFactory edgeSurfaceRuleFactory;

        /// <summary>
        /// Base rule for the optimization procedure within the HMF procedure.
        /// In essence, provides the (fixed) integration nodes
        /// </summary>
        private CellBoundaryQuadRule baseRule;

        /// <summary>
        /// Rule used to determine the sign of the level set function for
        /// <b>uncut</b> edges.
        /// </summary>
        private CellBoundaryQuadRule signTestRule;

        /// <summary>
        /// Values of <see cref="lambdaBasis"/> for a single edge.
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Polynomial index</item>
        /// </list>
        /// </summary>
        private MultidimensionalArray basisValuesEdge;

        /// <summary>
        /// Vector-valued moment-fitting basis constructed from the
        /// 'anti-derivatives' of a standard basis. Here, the 'anti-derivative'
        /// of a scalar polynomial \f$ p\f$  is defined as
        /// a vector-valued polynomial
        /// \f$ \vec{\Lambda}\f$  such that
        /// \f$ 
        /// p = \nabla \cdot \vec{\Lambda}
        /// \f$ 
        /// </summary>
        private PolynomialList lambdaBasis;

        /// <summary>
        /// The sub-grid for the current execution mask
        /// </summary>
        private SubGrid subGrid;

        /// <summary>
        /// The requested polynomial order in the last call of
        /// <see cref="GetQuadRuleSet"/>; used for caching
        /// </summary>
        private int lastOrder = -1;

        /// <summary>
        /// Cache for quadrature rules
        /// <list type="bullet">
        ///     <item>Key: Local cell index</item>
        ///     <item>Value: Quadrature rule</item>
        /// </list>
        /// </summary>
        private Dictionary<int, CellBoundaryQuadRule> cache =
            new Dictionary<int, CellBoundaryQuadRule>();

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return cache.Keys.ToArray();
        }

        /// <summary>
        /// Index of the considered level set for <see cref="tracker"/>
        /// </summary>
        private int levelSetIndex;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="lsData"></param>
        /// <param name="rootFindingAlgorithm">
        /// One-dimensional root-finding algorithm for determining roots of the
        /// level set function on edges of the edges of the volume simplex.
        /// Default is <see cref="LineSegment.DefaultRootFindingAlgorithm"/>
        /// </param>
        /// <param name="jumpType">
        /// Determines the level set region to be integrated over (negative
        /// level set values, positive level set values, or both)
        /// </param>
        public LevelSetEdgeVolumeQuadRuleFactory(
            LevelSetTracker.LevelSetData lsData, LineSegment.IRootFindingAlgorithm rootFindingAlgorithm = null, JumpTypes jumpType = JumpTypes.Heaviside) {
            if (lsData.GridDat.Cells.RefElements.Length > 1) {
                throw new NotImplementedException(
                    "Multiple reference elements currently not supported");
            }
            this.LevelSetData = lsData;
            this.jumpType = jumpType;
            this.levelSetIndex = lsData.LevelSetIndex;
            CoFaceQuadRuleFactory = new CutLineOnEdgeQuadRuleFactory(lsData, rootFindingAlgorithm, jumpType);
            edgeSurfaceRuleFactory = new LevelSetEdgeSurfaceQuadRuleFactory(lsData, CoFaceQuadRuleFactory, jumpType);

            // Use vertices; Since it is only used on edges that are considered
            // uncut, they _should_ all have the same sign
            RefElement simplex = LevelSetData.GridDat.Grid.RefElements[0];
            RefElement edgeSimplex = simplex.FaceRefElement;
            QuadRule signEdgeRule = new QuadRule() {
                Nodes = edgeSimplex.Vertices,
                Weights = MultidimensionalArray.Create(edgeSimplex.NoOfVertices)
            };

            signTestRule = new CellBoundaryFromEdgeRuleFactory<CellBoundaryQuadRule>(
                LevelSetData.GridDat, simplex, new FixedRuleFactory<QuadRule>(signEdgeRule)).
                GetQuadRuleSet(new CellMask(LevelSetData.GridDat, Chunk.GetSingleElementChunk(0)), -1).
                First().Rule;
        }

        internal LevelSetTracker.LevelSetData LevelSetData {
            get;
            private set;
        }

        private int[] localCellIndex2SubgridIndex;

        #region IQuadRuleFactory<CellBoundaryQuadRule> Members

        /// <summary>
        /// Constructs the quadrature rules all edges of all cells in
        /// <paramref name="mask"/>. For edges that are not intersected by the
        /// zero iso-contour, standard Gaussian quadrature rules of
        /// sufficiently high order will be used.
        /// </summary>
        /// <param name="mask">
        /// Cells for which quadrature rules shall be created
        /// </param>
        /// <param name="order">
        /// Desired order of the moment-fitting system. Assuming that
        /// <see cref="edgeSurfaceRuleFactory"/> integrates the basis
        /// polynomials exactly over the zero iso-contour (which it usually
        /// doesn't!), the resulting quadrature rules will be exact up to this
        /// order.
        /// </param>
        /// <returns>A set of quadrature rules</returns>
        /// <remarks>
        /// Since the selected level set is generally discontinuous across cell
        /// boundaries, this method does not make use of the fact that
        /// neighboring cells share edges. That is, the optimization will be
        /// performed twice for each inner edge in <paramref name="mask"/>.
        /// </remarks>
        public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            using (var tr = new FuncTrace()) {
                if (!(mask is CellMask)) {
                    throw new ArgumentException("CellMask required", "mask");
                }
#if DEBUG 
                CellMask differingCells = ((CellMask)mask).Except(this.LevelSetData.Region.GetCutCellMask4LevSet(this.levelSetIndex));
                if (differingCells.NoOfItemsLocally > 0) {
                    throw new ArgumentException("The provided mask has to be a sub-set of the cut cells. Cells " + differingCells.GetSummary() +  " are not in the CutCellMaks of this tracker.");
                }
#endif


                int noOfEdges = LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces;
                subGrid = new SubGrid((CellMask)mask);
                localCellIndex2SubgridIndex = subGrid.LocalCellIndex2SubgridIndex;

                if (order != lastOrder) {
                    cache.Clear();
                    SwitchOrder(order);
                }

                var result = new List<ChunkRulePair<CellBoundaryQuadRule>>(mask.NoOfItemsLocally);
                CellBoundaryQuadRule[] optimizedRules = GetOptimizedRules((CellMask)mask, order);
                int n = 0;
                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        if (cache.ContainsKey(cell)) {
                            result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                Chunk.GetSingleElementChunk(cell), cache[cell]));
                        } else {
                            cache.Add(cell, optimizedRules[n]);
                            result.Add(new ChunkRulePair<CellBoundaryQuadRule>(
                                Chunk.GetSingleElementChunk(cell), optimizedRules[n]));
                        }

                        n++;
                    }
                }

                return result;
            }
        }

        #endregion

        /// <summary>
        /// The first reference element of the grid
        /// </summary>
        public RefElement RefElement {
            get {
                return LevelSetData.GridDat.Grid.RefElements[0];
            }
        }

        /// <summary>
        /// Uses a moment-fitting basis of order <paramref name="order"/> for
        /// determining optimal weights of a quadrature rule for a sub-region
        /// of each edge of each cell in the given <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// Cells for which rules shall be created
        /// </param>
        /// <param name="order">
        /// Desired order of the moment-fitting system. Assuming that
        /// <see cref="edgeSurfaceRuleFactory"/> integrates the basis
        /// polynomials exactly over the zero iso-contour (which it usually
        /// doesn't), the resulting quadrature rules will be exact up to this
        /// order.
        /// </param>
        /// <returns>
        /// A set of quadrature rules
        /// </returns>
        private CellBoundaryQuadRule[] GetOptimizedRules(CellMask mask, int order) {
            using (var tr = new FuncTrace()) {
                int maxLambdaDegree = order + 1;
                int noOfLambdas = GetNumberOfLambdas();
                int noOfFaces = LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces;
                int D = LevelSetData.GridDat.SpatialDimension;
                int noOfNodesPerEdge = baseRule.NumbersOfNodesPerFace[0];
                CellBoundaryQuadRule[] optimizedRules = new CellBoundaryQuadRule[mask.NoOfItemsLocally];

                Debug.Assert(
                    baseRule.NumbersOfNodesPerFace.Distinct().Count() == 1,
                    "Assumption violated: Number of nodes varies from edge to edge.");
                Debug.Assert(noOfLambdas < noOfNodesPerEdge, "Not enough integration points");

                LambdaEdgeBoundaryQuadrature cellBoundaryQuadrature =
                    new LambdaEdgeBoundaryQuadrature(this, CoFaceQuadRuleFactory, maxLambdaDegree, mask);
                cellBoundaryQuadrature.Execute();
                double[, ,] boundaryResults = cellBoundaryQuadrature.Results;

                LambdaLevelSetSurfaceQuadrature surfaceQuadrature =
                    new LambdaLevelSetSurfaceQuadrature(this, edgeSurfaceRuleFactory, maxLambdaDegree, mask);
                surfaceQuadrature.Execute();
                double[, ,] surfaceResults = surfaceQuadrature.Results;

                int noOfRhs = 0;
                int[,] rhsIndexMap = new int[mask.NoOfItemsLocally, noOfFaces];
                foreach (Chunk chunk in mask) {
                    MultidimensionalArray levelSetValues = LevelSetData.GetLevSetValues(signTestRule.Nodes, chunk.i0, chunk.Len);
                    
                    for (int i = 0; i < chunk.Len; i++) {
                        int cell = i + chunk.i0;
                        int iSubGrid = localCellIndex2SubgridIndex[cell];

                        optimizedRules[iSubGrid] = new CellBoundaryQuadRule() {
                            Nodes = baseRule.Nodes,
                            Weights = MultidimensionalArray.Create(baseRule.NoOfNodes),
                            NumbersOfNodesPerFace = baseRule.NumbersOfNodesPerFace.CloneAs()
                        };

                        int noOfProcessedNodes = 0;
                        for (int e = 0; e < noOfFaces; e++) {
                            int noOfNodesOnEdge = baseRule.NumbersOfNodesPerFace[e];

                            bool faceIsCut = false;
                            for (int k = 0; k < noOfLambdas; k++) {
                                faceIsCut |= surfaceResults[iSubGrid, e, k].Abs() > 1e-9;
                            }

                            //edgeIsCut = tracker.EdgeIsCut[cell, e];

                            if (!faceIsCut) {
                                // Sign is checked in multiple points to avoid
                                // some weird edge cases. Sign with most 
                                // occurrences on the test nodes wins
                                int numNeg = 0;
                                int numPos = 0;
                                int offset = signTestRule.NumbersOfNodesPerFace.Take(e).Sum();
                                for (int j = 0; j < signTestRule.NumbersOfNodesPerFace[e]; j++) {
                                    double val = levelSetValues[i, offset + j];
                                    if (val < 0.0) {
                                        numNeg++;
                                    } else if (val > 0.0) {
                                        numPos++;
                                    }
                                }

                                int sign = numPos - numNeg;

                                if (sign == 0) {
                                    throw new Exception(String.Format(
                                        "Could not determine sign of face {0} of cell {1}", e, cell));
                                }

                                switch (jumpType) {
                                    case JumpTypes.Heaviside:
                                        if (sign > 0) {
                                            CopyWeights(baseRule, optimizedRules[iSubGrid], e, 1.0);
                                        }
                                        break;

                                    case JumpTypes.OneMinusHeaviside:
                                        if (sign < 0) {
                                            CopyWeights(baseRule, optimizedRules[iSubGrid], e, -1.0);
                                        }
                                        break;

                                    case JumpTypes.Sign:
                                        CopyWeights(baseRule, optimizedRules[iSubGrid], e, levelSetValues[i, e].Sign());
                                        break;

                                    default:
                                        throw new NotImplementedException();
                                }

                                rhsIndexMap[iSubGrid, e] = -1;
                                noOfProcessedNodes += noOfNodesOnEdge;
                                continue;
                            }

                            rhsIndexMap[iSubGrid, e] = noOfRhs;
                            noOfRhs++;
                            noOfProcessedNodes += noOfNodesOnEdge;
                        }
                    }
                }

                // Leading dimension of B (rhs); required by DGELSY
                int LDB = Math.Max(noOfLambdas, noOfNodesPerEdge);
                double[] rhs = new double[LDB * noOfRhs];
                foreach (Chunk chunk in mask) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int cell = i + chunk.i0;
                        int iSubGrid = localCellIndex2SubgridIndex[cell];

                        for (int e = 0; e < noOfFaces; e++) {
                            int rhsIndex = rhsIndexMap[iSubGrid, e];
                            if (rhsIndex < 0) {
                                continue;
                            }

                            switch (jumpType) {
                                case JumpTypes.Heaviside:
                                    for (int k = 0; k < noOfLambdas; k++) {
                                        rhs[rhsIndex * LDB + k] = boundaryResults[iSubGrid, e, k] - surfaceResults[iSubGrid, e, k];
                                    }
                                    break;

                                case JumpTypes.OneMinusHeaviside:
                                    for (int k = 0; k < noOfLambdas; k++) {
                                        rhs[rhsIndex * LDB + k] = boundaryResults[iSubGrid, e, k] + surfaceResults[iSubGrid, e, k];
                                    }
                                    break;

                                case JumpTypes.Sign:
                                    for (int k = 0; k < noOfLambdas; k++) {
                                        rhs[rhsIndex * LDB + k] = boundaryResults[iSubGrid, e, k] - 2.0 * surfaceResults[iSubGrid, e, k];
                                    }
                                    break;

                                default:
                                    throw new NotImplementedException();
                            }
                        }
                    }
                }

                if (rhs.Length > 0) {
                    LAPACK.F77_LAPACK.DGELSY(noOfLambdas, noOfNodesPerEdge, basisValuesEdge.Storage, rhs, noOfRhs, 1e-12);
                }

                foreach (Chunk chunk in mask) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int cell = i + chunk.i0;
                        int iSubGrid = localCellIndex2SubgridIndex[cell];

                        int noOfProcessedNodes = 0;
                        for (int e = 0; e < noOfFaces; e++) {
                            int noOfNodesOnEdge = baseRule.NumbersOfNodesPerFace[e];
                            int rhsIndex = rhsIndexMap[iSubGrid, e];
                            if (rhsIndex < 0) {
                                noOfProcessedNodes += noOfNodesOnEdge;
                                continue;
                            }

                            for (int j = 0; j < noOfNodesOnEdge; j++) {
                                optimizedRules[iSubGrid].Weights[noOfProcessedNodes + j] = rhs[rhsIndex * LDB + j];
                            }

                            noOfProcessedNodes += noOfNodesOnEdge;
                        }


                        double max = optimizedRules[iSubGrid].Weights.Max(d => d.Abs());
                        if (max > 2.0 * RefElement.Volume) {
                            tr.Info(String.Format(
                                "Warning: Abnormally large integration weight detected"
                                + " for level set edge volume integral in cell {0}"
                                + " (|w| = {1}). This may indicate a loss of"
                                + " integration accuracy.",
                                cell,
                                max));
                        }
                    }
                }

                return optimizedRules;
            }
        }

        /// <summary>
        /// Updates <see cref="lambdaBasis"/>, <see cref="baseRule"/> and
        /// <see cref="basisValuesEdge"/> every time a different order is
        /// requested in <see cref="GetQuadRuleSet"/>
        /// </summary>
        /// <param name="order">
        /// The new order of the moment-fitting basis
        /// </param>
        private void SwitchOrder(int order) {
            int iKref = this.LevelSetData.GridDat.Cells.RefElements.IndexOf(this.RefElement, (A, B) => object.ReferenceEquals(A, B));
            int noOfFaces = LevelSetData.GridDat.Grid.RefElements[iKref].NoOfFaces;
            int D = RefElement.FaceRefElement.SpatialDimension;
            
            Polynomial[] basePolynomials = RefElement.FaceRefElement.GetOrthonormalPolynomials(order).ToArray();
            Polynomial[] antiderivativePolynomials =
                new Polynomial[basePolynomials.Length * D];
            for (int i = 0; i < basePolynomials.Length; i++) {
                Polynomial p = basePolynomials[i];
                for (int d = 0; d < D; d++) {
                    Polynomial pNew = p.CloneAs();
                    for (int j = 0; j < p.Coeff.Length; j++) {
                        pNew.Exponents[j, d]++;
                        pNew.Coeff[j] /= pNew.Exponents[j, d];
                        // Make sure divergence is Phi again
                        pNew.Coeff[j] /= D;
                    }
                    antiderivativePolynomials[i * D + d] = pNew;
                }
            }
            lambdaBasis = new PolynomialList(antiderivativePolynomials);

            int minNoOfPoints = GetNumberOfLambdas();

            int minOrder = 1;
            double safetyFactor = 1.6;
            while (RefElement.FaceRefElement.GetQuadratureRule(minOrder).NoOfNodes < safetyFactor * minNoOfPoints) {
                minOrder += 1;
            }

            QuadRule singleEdgeRule = RefElement.FaceRefElement.GetQuadratureRule(minOrder);
            baseRule = new CellBoundaryFromEdgeRuleFactory<CellBoundaryQuadRule>(
                LevelSetData.GridDat,
                LevelSetData.GridDat.Grid.RefElements[0],
                new FixedRuleFactory<QuadRule>(singleEdgeRule)).
                GetQuadRuleSet(new CellMask(LevelSetData.GridDat, Chunk.GetSingleElementChunk(0)), -1).
                First().Rule;

            //Basis singleEdgeBasis = new Basis(tracker.GridDat, order, RefElement.FaceRefElement);
            PolynomialList singleEdgeBasis = new PolynomialList(RefElement.FaceRefElement.GetOrthonormalPolynomials(order));
            //MultidimensionalArray edgeNodes = singleEdgeRule.Nodes.CloneAs();
            //basisValuesEdge = MultidimensionalArray.Create(
            //    singleEdgeRule.NoOfNodes, singleEdgeBasis.Count);
            //MultidimensionalArray monomials = Polynomial.GetMonomials(
            //    edgeNodes, RefElement.FaceRefElement.SpatialDimension, singleEdgeBasis.Degree);
            //for (int j = 0; j < singleEdgeBasis.MinimalLength; j++) {
            //    singleEdgeBasis.Polynomials[iKref, j].Evaluate(
            //        basisValuesEdge.ExtractSubArrayShallow(-1, j),
            //        edgeNodes,
            //        monomials);
            //}
            this.basisValuesEdge = singleEdgeBasis.Values.GetValues(singleEdgeRule.Nodes);
        }

        /// <summary>
        /// Copies the (scaled) weights for the local edge referenced by
        /// <paramref name="localEdge"/> from <paramref name="source"/> to
        /// <paramref name="target"/>
        /// </summary>
        /// <param name="source">Source quad rule</param>
        /// <param name="target">Target quad rule</param>
        /// <param name="localEdge">
        /// The considered edge of the volume element (cf.
        /// <see cref="RefElement"/>)
        /// </param>
        /// <param name="scaling">
        /// Scaling applied to the weights before saving them in
        /// <paramref name="target"/>.
        /// </param>
        private void CopyWeights(CellBoundaryQuadRule source, CellBoundaryQuadRule target, int localEdge, double scaling) {
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
                target.Weights[offsetTarget + i] = scaling * source.Weights[offsetSource + i];
            }
        }

        /// <summary>
        /// Computes the number of vector-valued basis function (cf.
        /// <see cref="lambdaBasis"/>)
        /// </summary>
        /// <returns></returns>
        private int GetNumberOfLambdas() {
            return lambdaBasis.Count / RefElement.FaceRefElement.SpatialDimension;
        }

        /// <summary>
        /// Evaluates the vector-valued anti-derivatives defining the
        /// moment-fitting basis (cf. <see cref="lambdaBasis"/>) in each node
        /// of <paramref name="rule"/> in cell <paramref name="cell"/>
        /// </summary>
        /// <param name="cell">
        /// Local cell index
        /// </param>
        /// <param name="rule">
        /// Evaluation nodes
        /// </param>
        /// <returns>
        /// The values of <see cref="lambdaBasis"/> in each node of
        /// <paramref name="rule"/>
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Basis function index</item>
        ///     <item>3rd index: Spatial dimension</item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// This method does not evaluate a Lambda which is constant since the
        /// derivative of such function is zero. As a result, the integral over
        /// this function is zero, too, which makes it useless for the
        /// construction of a quadrature rule
        /// </remarks>
        private MultidimensionalArray EvaluateLambdas(int cell, CellBoundaryQuadRule rule) {
            int noOfLambdas = GetNumberOfLambdas();
            int noOfNodes = rule.NoOfNodes;
            int D = LevelSetData.GridDat.SpatialDimension;
            MultidimensionalArray allNodes = rule.Nodes.CloneAs();

            MultidimensionalArray result = MultidimensionalArray.Create(
                noOfNodes, noOfLambdas, D);
            int noOfProcessedNodes = 0;
            for (int e = 0; e < LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                int noOfNodesOnEdge = rule.NumbersOfNodesPerFace[e];

                if (noOfNodesOnEdge == 0) {
                    continue;
                }

                MultidimensionalArray nodes = allNodes.ExtractSubArrayShallow(
                    new int[] { noOfProcessedNodes, 0 },
                    new int[] { noOfProcessedNodes + noOfNodesOnEdge - 1, D - 1 });
                NodeSet edgeNodes = new NodeSet(RefElement.FaceRefElement, RefElement.GetInverseFaceTrafo(e).Transform(nodes));

                //MultidimensionalArray monomials = Polynomial.GetMonomials(
                //    edgeNodes, RefElement.FaceRefElement.SpatialDimension, lambdaBasis.Degree);
                //MultidimensionalArray monomials = Caching.MonomialCache.Instance.GetMonomials(edgeNodes, lambdaBasis.MaxAbsoluteDegree);
                MultidimensionalArray lambdaValuesEdge = MultidimensionalArray.Create(
                    noOfNodesOnEdge, noOfLambdas, RefElement.FaceRefElement.SpatialDimension);
                Debug.Assert(noOfLambdas == lambdaBasis.Count / RefElement.FaceRefElement.SpatialDimension);

                //for (int j = 0; j < noOfLambdas; j++) {
                //    for (int d = 0; d < RefElement.FaceRefElement.SpatialDimension; d++) {
                //        lambdaBasis[j * RefElement.FaceRefElement.SpatialDimension + d].Evaluate(
                //            lambdaValuesEdge.ExtractSubArrayShallow(-1, j, d),
                //            edgeNodes,
                //            monomials);
                //    }
                //}
                lambdaValuesEdge = lambdaBasis.Values.GetValues(edgeNodes);
                
                MultidimensionalArray resultOnEdge = MultidimensionalArray.Create(
                    noOfNodesOnEdge, noOfLambdas, D);
                this.RefElement.TransformFaceCoordinates(
                    e,
                    lambdaValuesEdge.ResizeShallow(
                        noOfLambdas * noOfNodesOnEdge, RefElement.FaceRefElement.SpatialDimension),
                    resultOnEdge.ResizeShallow(noOfLambdas * noOfNodesOnEdge, D));

                for (int i = 0; i < noOfNodesOnEdge; i++) {
                    for (int j = 0; j < noOfLambdas; j++) {
                        for (int d = 0; d < D; d++) {
                            result[noOfProcessedNodes + i, j, d] = resultOnEdge[i, j, d];
                        }
                    }
                }

                noOfProcessedNodes += noOfNodesOnEdge;
            }

            return result;
        }

        /// <summary>
        /// Used for the computation of the integral of the basis functions
        /// (see <see cref="lambdaBasis"/>) over the boundary of an edge. More
        /// specifically, the boundary of the edge is given by line elements
        /// that are potentially cut by the zero iso-contour of the level set,
        /// which is why this requires a special treatment; see
        /// <see cref="CoFaceQuadRuleFactory"/>.
        /// </summary>
        private class LambdaEdgeBoundaryQuadrature : CellBoundaryQuadrature<CellEdgeBoundaryQuadRule> {

            /// <summary>
            /// Owner of this object
            /// </summary>
            private LevelSetEdgeVolumeQuadRuleFactory owner;

            /// <summary>
            /// Number of items in the current execution mask
            /// </summary>
            private int noOfItemsLocally;

            /// <summary>
            /// Integration results for the last call to <see cref="Execute"/>
            /// <list type="bullet">
            ///     <item>
            ///     1st index: Cell index within
            ///     <see cref="LevelSetEdgeVolumeQuadRuleFactory.subGrid"/>
            ///     </item>
            ///     <item>
            ///     2nd index: Local edge index
            ///     </item>
            ///     <item>
            ///     3rd index: Basis function (see
            ///     <see cref="LevelSetEdgeVolumeQuadRuleFactory.EvaluateLambdas"/>
            ///     </item>
            /// </list>
            /// </summary>
            public double[, ,] Results;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="edgeRuleFactory"></param>
            /// <param name="maxLambdaDegree"></param>
            /// <param name="mask"></param>
            public LambdaEdgeBoundaryQuadrature(
                LevelSetEdgeVolumeQuadRuleFactory owner,
                IQuadRuleFactory<CellEdgeBoundaryQuadRule> edgeRuleFactory,
                int maxLambdaDegree,
                CellMask mask)
                : base(
                    new int[] { owner.GetNumberOfLambdas() },
                    owner.LevelSetData.GridDat,
                    (new CellEdgeBoundaryQuadratureScheme(false, edgeRuleFactory, mask)).Compile(owner.LevelSetData.GridDat, maxLambdaDegree),
                    CoordinateSystem.Reference) {
                this.owner = owner;
                noOfItemsLocally = mask.NoOfItemsLocally;
            }

            /// <summary>
            /// Computes the integrals while overwriting <see cref="Results"/>
            /// </summary>
            public override void Execute() {
                Results = new double[
                    noOfItemsLocally,
                    gridData.iGeomCells.RefElements[0].NoOfFaces,
                    IntegralCompDim[0]];
                base.Execute();
            }

            
            /// <summary>
            /// Not implemented.
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="L"></param>
            /// <returns></returns>
            protected override MultidimensionalArray GetScalingsForLinearElements(int i0, int L) {
                throw new NotSupportedException("special treatment for this class");
            }

            /// <summary>
            /// For each edge \f$ E\f$  of each cell in
            /// the given range and for each
            /// \f$ \vec{\Lambda}\f$  in
            /// <see cref="LevelSetEdgeVolumeQuadRuleFactory.lambdaBasis"/>:
            /// Computes
            /// \f$ 
            /// \int \limits_{\partial E} \vec{\Lambda} \cdot \vec{n} H(\varphi) \;ds,
            /// \f$ 
            /// where \f$ \vec{n}\f$  denotes the outer
            /// unit normal vector on \f$ \partial E\f$ 
            /// (<b>not</b> \f$ E\f$ !). Moreover,
            /// \f$ H\f$  a weight function that depends
            /// on the level set function \f$ \varphi\f$ .
            /// Typically, \f$ H\f$  is given by the
            /// Heaviside function. For more details, see
            /// <see cref="CutLineQuadRuleFactory"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="NoOfNodes"></param>
            /// <param name="EvalResult"></param>
            protected override void Evaluate(int i0, int Length, CellEdgeBoundaryQuadRule CQR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = CQR.Nodes;
                int noOfEdges = owner.RefElement.NoOfFaces;
                int noOfEdgesOfEdge = owner.RefElement.FaceRefElement.NoOfFaces;
                int D = gridData.SpatialDimension;
                MultidimensionalArray edgeOfEdgeNormals = MultidimensionalArray.Create(
                    1, owner.RefElement.FaceRefElement.SpatialDimension);
                MultidimensionalArray edgeNormals = MultidimensionalArray.Create(1, D);

                Debug.Assert(
                    owner.RefElement.FaceRefElement.FaceRefElement is Line,
                    "Edge of edge _must_ be a line");

                for (int i = 0; i < Length; i++) {
                    MultidimensionalArray lambdaValues = owner.EvaluateLambdas(i0 + i, CurrentRule);

                    int nodeIndex = -1;

                    for (int e = 0; e < noOfEdges; e++) {
                        for (int ee = 0; ee < noOfEdgesOfEdge; ee++) {

                            for (int d = 0; d < owner.RefElement.FaceRefElement.SpatialDimension; d++) {
                                edgeOfEdgeNormals[0, d] = owner.RefElement.FaceRefElement.FaceNormals[ee, d];
                            }

                            owner.LevelSetData.GridDat.Grid.RefElements[0].TransformFaceVectors(
                                e, edgeOfEdgeNormals, edgeNormals);

                            for (int j = 0; j < CurrentRule.NumbersOfNodesPerFaceOfFace[e, ee]; j++) {
                                nodeIndex++;

                                for (int k = 0; k < IntegralCompDim[0]; k++) {
                                    for (int d = 0; d < D; d++) {
                                        EvalResult[i, nodeIndex, k] +=
                                            lambdaValues[nodeIndex, k, d] * edgeNormals[0, d];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// Saves the integration results to <see cref="Results"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    int cell = i0 + i;
                    int iSubGrid = owner.localCellIndex2SubgridIndex[cell];

                    for (int e = 0; e < owner.LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                        for (int k = 0; k < IntegralCompDim[0]; k++) {
                            Results[iSubGrid, e, k] = ResultsOfIntegration[i, e, k];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Used for the computation of the integral of the basis functions
        /// (see <see cref="lambdaBasis"/>) over the intersection of the zero
        /// iso-contour of the level set with an edge. This requires a special
        /// treatment; see <see cref="edgeSurfaceRuleFactory"/>.
        /// </summary>
        private class LambdaLevelSetSurfaceQuadrature : CellBoundaryQuadrature<CellBoundaryQuadRule> {

            /// <summary>
            /// Owner of this object
            /// </summary>
            private LevelSetEdgeVolumeQuadRuleFactory owner;

            /// <summary>
            /// Number of items in the current execution mask
            /// </summary>
            private int noOfItemsLocally;

            /// <summary>
            /// Integration results for the last call to <see cref="Execute"/>
            /// <list type="bullet">
            ///     <item>
            ///     1st index: Cell index within
            ///     <see cref="LevelSetEdgeVolumeQuadRuleFactory.subGrid"/>
            ///     </item>
            ///     <item>
            ///     2nd index: Local edge index
            ///     </item>
            ///     <item>
            ///     3rd index: Basis function (see
            ///     <see cref="LevelSetEdgeVolumeQuadRuleFactory.EvaluateLambdas"/>
            ///     </item>
            /// </list>
            /// </summary>
            public double[, ,] Results;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="surfaceRuleFactory"></param>
            /// <param name="maxLambdaDegree"></param>
            /// <param name="mask"></param>
            public LambdaLevelSetSurfaceQuadrature(LevelSetEdgeVolumeQuadRuleFactory owner, IQuadRuleFactory<CellBoundaryQuadRule> surfaceRuleFactory, int maxLambdaDegree, CellMask mask)
                : base(new int[] { owner.GetNumberOfLambdas() },
                       owner.LevelSetData.GridDat,
                       (new CellBoundaryQuadratureScheme(surfaceRuleFactory, mask)).Compile(owner.LevelSetData.GridDat, maxLambdaDegree),
                       CoordinateSystem.Reference) //
            {
                this.owner = owner;
                noOfItemsLocally = mask.NoOfItemsLocally;
            }



            /// <summary>
            /// Computes the integrals while overwriting <see cref="Results"/>
            /// </summary>
            public override void Execute() {
                Results = new double[
                    noOfItemsLocally,
                    gridData.iGeomCells.RefElements[0].NoOfFaces,
                    IntegralCompDim[0]];
                base.Execute();
            }

            

            /// <summary>
            /// For each edge \f$ E\f$  of each cell in
            /// the given range and for each
            /// \f$ \vec{\Lambda}\f$  in
            /// <see cref="LevelSetEdgeVolumeQuadRuleFactory.lambdaBasis"/>:
            /// Computes
            /// \f$ 
            /// \int_{ \{ \vec{x}; \varphi(\vec{x}) = 0 \}  \cap E} \vec{\Lambda} \cdot \vec{n}_I \;ds,
            /// \f$ 
            /// where \f$ \varphi\f$  is the level set
            /// function and \f$ \vec{n}_I\f$  denotes the
            /// unit normal vector on
            /// \f$ \varphi \cap E\f$ 
            /// (<b>not</b> \f$ E\f$ !)
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="NoOfNodes"></param>
            /// <param name="EvalResult"></param>
            protected override void Evaluate(int i0, int Length, CellBoundaryQuadRule CBQR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = CBQR.Nodes;
                int D = gridData.SpatialDimension;

                for (int i = 0; i < Length; i++) {
                    MultidimensionalArray lambdaValues = owner.EvaluateLambdas(i0 + i, CurrentRule);
                    int nodeIndex = -1;

                    for (int e = 0; e < gridData.iGeomCells.RefElements[0].NoOfFaces; e++) {
                        MultidimensionalArray levelSetNormals =
                            LevelSetEdgeSurfaceQuadRuleFactory.EvaluateRefNormalsOnEdge(this.owner.LevelSetData, i0 + i, CurrentRule, e);
                        //MultidimensionalArray metrics = LevelSetEdgeSurfaceQuadRuleFactory.GetMetricTermsOnEdge(
                        //    this.owner.LevelSetData, this.owner.levelSetIndex, CurrentRule, i0 + i, e);

                        for (int j = 0; j < CurrentRule.NumbersOfNodesPerFace[e]; j++) {
                            nodeIndex++;

                            for (int k = 0; k < IntegralCompDim[0]; k++) {
                                for (int d = 0; d < D; d++) {
                                    EvalResult[i, nodeIndex, k] += lambdaValues[nodeIndex, k, d] * levelSetNormals[j, d];
                                }

                                //EvalResult[i, nodeIndex, k] *= metrics[j];
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// Saves the integration results to <see cref="Results"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    int cell = i0 + i;
                    int iSubGrid = owner.localCellIndex2SubgridIndex[cell];

                    for (int e = 0; e < owner.LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                        for (int k = 0; k < IntegralCompDim[0]; k++) {
                            Results[iSubGrid, e, k] = ResultsOfIntegration[i, e, k];
                        }
                    }
                }
            }
        }
    }
}
