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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform.Utils.Geom;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// This factory produces quadrature rules which are, for each cell
    /// \f$ K\f$  in a volume mask, capable of computing
    /// (an approximation of)
    /// \f[ 
    ///    \int\limits_{\{ \vec{x}; \varphi(\vec{x}) {\leq \atop \geq} 0 \} \cap K}  f \ d \vec{x},
    /// \f]
    /// where \f$ \varphi\f$  denotes the level set function.
    /// </summary>
    public class LevelSetVolumeQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        /// <summary>
        /// Determines the level set region to be integrated over (negative
        /// level set values, positive level set values, or both)
        /// </summary>
        private JumpTypes jumpType;

        /// <summary>
        /// Tracker containing the level set to be integrated over
        /// </summary>
        private LevelSetTracker tracker;

        /// <summary>
        /// Index of the considered level set for <see cref="tracker"/>
        /// </summary>
        private int levelSetIndex;

        /// <summary>
        /// Factory for one- or two-dimensional integration rules for the edges
        /// of the of the volume element.
        /// </summary>
        private IQuadRuleFactory<CellBoundaryQuadRule> edgeRuleFactory;

        /// <summary>
        /// Factory for the one- or two-dimensional integration over the
        /// interface that is defined as the intersection of the zero 
        /// iso-contour of the level set function with the three-dimensional
        /// volume element
        /// </summary>
        private IQuadRuleFactory<QuadRule> surfaceRuleFactory;

        private int[] localCellIndex2SubgridIndex;

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
        /// Safety factor for the number of the nodes with respect to the
        /// number of modes. A factor of 1 (default) indicates that at least
        /// as many nodes as modes are used; a higher value may increase
        /// robustness and runtime.
        /// </summary>
        public static double NodeCountSafetyFactor = 1.0;

        public static bool RestrictNodes;

        public static bool UseGaussNodes = true;

        private const double RCOND = 1e-11;

        private Dictionary<int, ChunkRulePair<QuadRule>[]> cachedRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        private AffineTrafo[] trafos;

        /// <summary>
        /// constructor.
        /// </summary>
        /// <param name="lsData">
        /// </param>
        /// <param name="edgeRuleFactory">
        /// Some factory that provides quadrature rules for the integration 
        /// over 
        /// \f[ 
        ///  \partial K \cap \{ \vec{x}; \varphi(\vec{x}) {\leq \atop \geq} 0 \}.
        /// \f]
        /// Here, \f$ \partial K\f$  the boundary of
        /// some cell \f$ K\f$  and 
        /// \f$ \varphi\f$  denotes the level set function
        /// </param>
        /// <param name="surfaceRuleFactory">
        /// Some factory that provides quadrature rules for the integration 
        /// over the zero level set, i.e.
        /// \f[ 
        ///   \{ \vec{x}; \varphi(\vec{x}) = 0 \} \cap K.
        /// \f]
        /// (ere, \f$ \partial K\f$  the boundary of some
        /// cell \f$ K\f$  and 
        /// \f$ \varphi\f$  denotes the level set function
        /// </param>
        /// <param name="levSetIndex">
        /// Index of the considered level set in <paramref name="tracker"/>
        /// </param>

        public LevelSetVolumeQuadRuleFactory(
            LevelSetTracker.LevelSetData lsData,
            IQuadRuleFactory<CellBoundaryQuadRule> edgeRuleFactory,
            IQuadRuleFactory<QuadRule> surfaceRuleFactory,
            JumpTypes jumpType = JumpTypes.Heaviside) {

            if (jumpType == JumpTypes.Implicit) {
                throw new NotSupportedException();
            }

            this.tracker = lsData.Tracker;
            this.levelSetIndex = lsData.LevelSetIndex;
            this.edgeRuleFactory = edgeRuleFactory;
            this.surfaceRuleFactory = surfaceRuleFactory;
            this.jumpType = jumpType;
            this.LevelSetData = lsData;
        }

        LevelSetTracker.LevelSetData LevelSetData;

        #region IQuadRuleFactory<QuadRule> Members

        /// <summary>
        /// The reference element of the grid
        /// </summary>
        public RefElement RefElement {
            get {
                return tracker.GridDat.Grid.RefElements[0];
            }
        }

        /// <summary>
        /// Constructs suitable quadrature rules cells in
        /// <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// Cells for which quadrature rules shall be created
        /// </param>
        /// <param name="order">
        /// Desired order of the moment-fitting system. Assuming that
        /// <see cref="surfaceRuleFactory"/> integrates the basis polynomials
        /// exactly over the zero iso-contour (which it usually
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
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            using (var tr = new FuncTrace()) {
                CellMask cellMask = mask as CellMask;
                if (cellMask == null) {
                    throw new ArgumentException("Mask must be a volume mask", "mask");
                }

                // Note: This is a parallel call, so do this early to avoid parallel confusion
                localCellIndex2SubgridIndex = new SubGrid(cellMask).LocalCellIndex2SubgridIndex;

                int maxLambdaDegree = order + 1;
                int noOfLambdas = GetNumberOfLambdas(maxLambdaDegree);
                int noOfEdges = tracker.GridDat.Grid.RefElements[0].NoOfFaces;
                int D = RefElement.SpatialDimension;

                // Get the basis polynomials and integrate them analytically
                Polynomial[] basePolynomials = RefElement.GetOrthonormalPolynomials(order).ToArray();
                Polynomial[] polynomials = new Polynomial[basePolynomials.Length * D];
                for (int i = 0; i < basePolynomials.Length; i++) {
                    Polynomial p = basePolynomials[i];

                    for (int d = 0; d < D; d++) {
                        Polynomial pNew = p.CloneAs();
                        for (int j = 0; j < p.Coeff.Length; j++) {
                            pNew.Exponents[j, d]++;
                            pNew.Coeff[j] /= pNew.Exponents[j, d];
                            pNew.Coeff[j] /= D; // Make sure divergence is Phi again
                        }
                        polynomials[i * D + d] = pNew;
                    }
                }

                // basePolynomials[i] == div(polynomials[i*D], ... , polynomials[i*D + D - 1])
                lambdaBasis = new PolynomialList(polynomials);


                if (RestrictNodes) {
                    trafos = new AffineTrafo[mask.NoOfItemsLocally];

                    foreach (Chunk chunk in mask) {
                        foreach (var cell in chunk.Elements.AsSmartEnumerable()) {
                            CellMask singleElementMask = new CellMask(
                                tracker.GridDat, Chunk.GetSingleElementChunk(cell.Value));

                            LineAndPointQuadratureFactory.LineQRF lineFactory = this.edgeRuleFactory as LineAndPointQuadratureFactory.LineQRF;
                            if (lineFactory == null) {
                                throw new Exception();
                            }
                            var lineRule = lineFactory.GetQuadRuleSet(singleElementMask, order).Single().Rule;
                            var pointRule = lineFactory.m_Owner.GetPointFactory().GetQuadRuleSet(singleElementMask, order).Single().Rule;

                            // Also add point rule points since line rule points
                            // are constructed from Gauss rules that do not include
                            // the end points
                            BoundingBox box = new BoundingBox(lineRule.Nodes);
                            box.AddPoints(pointRule.Nodes);

                            int noOfRoots = pointRule.Nodes.GetLength(0);
                            if (noOfRoots <= 1) {
                                // Cell is considered cut because the level set
                                // is very close, but actually isn't. Note that
                                // we can NOT omit the cell (as in the surface
                                // case) as it will be missing in the list of
                                // uncut cells, i.e. this cell would be ignored
                                // completely
                                trafos[localCellIndex2SubgridIndex[cell.Value]] =
                                    AffineTrafo.Identity(RefElement.SpatialDimension);
                                continue;
                            } else if (noOfRoots == 2) {
                                // Go a bit into the direction of the normal
                                // from the center between the nodes in order
                                // not to miss regions with strong curvature
                                double[] center = box.Min.CloneAs();
                                center.AccV(1.0, box.Max);
                                center.ScaleV(0.5);
                                NodeSet centerNode = new NodeSet(RefElement, center);
                                centerNode.LockForever();

                                MultidimensionalArray normal = LevelSetData.GetLevelSetReferenceNormals(centerNode, cell.Value, 1);
                                MultidimensionalArray dist = LevelSetData.GetLevSetValues(centerNode, cell.Value, 1);

                                double scaling = Math.Sqrt(tracker.GridDat.Cells.JacobiDet[cell.Value]);

                                double[] newPoint = new double[D];
                                for (int d = 0; d < D; d++) {
                                    newPoint[d] = center[d] - normal[0, 0, d] * dist[0, 0] / scaling;
                                }

                                box.AddPoint(newPoint);

                                // Make sure points stay in box
                                for (int d = 0; d < D; d++) {
                                    box.Min[d] = Math.Max(box.Min[d], -1);
                                    box.Max[d] = Math.Min(box.Max[d], 1);
                                }
                            }

                            MultidimensionalArray preImage = RefElement.Vertices.ExtractSubArrayShallow(
                                new int[] { 0, 0 }, new int[] { D, D - 1 });

                            MultidimensionalArray image = MultidimensionalArray.Create(D + 1, D);
                            image[0, 0] = box.Min[0]; // Top left
                            image[0, 1] = box.Max[1];
                            image[1, 0] = box.Max[0]; // Top right
                            image[1, 1] = box.Max[1];
                            image[2, 0] = box.Min[0]; // Bottom left;
                            image[2, 1] = box.Min[1];

                            AffineTrafo trafo = AffineTrafo.FromPoints(preImage, image);
                            trafos[localCellIndex2SubgridIndex[cell.Value]] = trafo;
                        }
                    }
                }

                LambdaCellBoundaryQuadrature cellBoundaryQuadrature =
                    new LambdaCellBoundaryQuadrature(this, edgeRuleFactory, cellMask);
                cellBoundaryQuadrature.Execute();

                LambdaLevelSetSurfaceQuadrature surfaceQuadrature =
                    new LambdaLevelSetSurfaceQuadrature(this, surfaceRuleFactory, cellMask);
                surfaceQuadrature.Execute();

                // Must happen _after_ all parallel calls (e.g., definition of
                // the sub-grid or quadrature) in order to avoid problems in
                // parallel runs
                if (mask.NoOfItemsLocally == 0) {
                    var empty = new ChunkRulePair<QuadRule>[0];
                    return empty;
                }

                if (cachedRules.ContainsKey(order)) {
                    order = cachedRules.Keys.Where(cachedOrder => cachedOrder >= order).Min();
                    CellMask cachedMask = new CellMask(mask.GridData, cachedRules[order].Select(p => p.Chunk).ToArray());

                    if (cachedMask.Equals(mask)) {
                        return cachedRules[order];
                    } else {
                        throw new NotImplementedException(
                            "Case not yet covered yet in combination with caching; deactivate caching to get rid of this message");
                    }
                }

                double[,] quadResults = cellBoundaryQuadrature.Results;
                foreach (Chunk chunk in mask) {
                    for (int i = 0; i < chunk.Len; i++) {
                        int iSubGrid = localCellIndex2SubgridIndex[chunk.i0 + i];

                        switch (jumpType) {
                            case JumpTypes.Heaviside:
                                for (int k = 0; k < noOfLambdas; k++) {
                                    quadResults[iSubGrid, k] -= surfaceQuadrature.Results[iSubGrid, k];
                                }
                                break;

                            case JumpTypes.OneMinusHeaviside:
                                for (int k = 0; k < noOfLambdas; k++) {
                                    quadResults[iSubGrid, k] += surfaceQuadrature.Results[iSubGrid, k];
                                }
                                break;

                            case JumpTypes.Sign:
                                for (int k = 0; k < noOfLambdas; k++) {
                                    quadResults[iSubGrid, k] -= 2.0 * surfaceQuadrature.Results[iSubGrid, k];
                                }
                                break;

                            default:
                                throw new NotImplementedException();
                        }
                    }
                }

                BitArray voidCellsArray = new BitArray(tracker.GridDat.Cells.NoOfLocalUpdatedCells);
                BitArray fullCellsArray = new BitArray(tracker.GridDat.Cells.NoOfLocalUpdatedCells);
                foreach (Chunk chunk in cellMask) {
                    foreach (var cell in chunk.Elements) {
                        double rhsL2Norm = 0.0;
                        for (int k = 0; k < noOfLambdas; k++) {
                            double entry = quadResults[localCellIndex2SubgridIndex[cell], k];
                            rhsL2Norm += entry * entry;
                        }

                        if (rhsL2Norm < 1e-14) {
                            // All integrals are zero => cell not really cut
                            // (level set is tangent) and fully in void region
                            voidCellsArray[cell] = true;
                            continue;
                        }

                        double l2NormFirstIntegral = quadResults[localCellIndex2SubgridIndex[cell], 0];
                        l2NormFirstIntegral *= l2NormFirstIntegral;
                        double rhsL2NormWithoutFirst = rhsL2Norm - l2NormFirstIntegral;
                        
                        // Beware: This check is only sensible if basis is orthonormal on RefElement!
                        if (rhsL2NormWithoutFirst < 1e-14 && 
                            Math.Abs(l2NormFirstIntegral - RefElement.Volume) < 1e-14) {
                            // All integrals are zero except integral over first integrand
                            // If basis is orthonormal, this implies that cell is uncut and
                            // fully in non-void region since then
                            // \int_K \Phi_i dV = \int_A \Phi_i dV = \delta_{0,i}
                            // However, we have to compare RefElement.Volume since
                            // integration is performed in reference coordinates!
                            fullCellsArray[cell] = true;
                        }
                    }
                }

                var result = new List<ChunkRulePair<QuadRule>>(cellMask.NoOfItemsLocally);

                CellMask emptyCells = new CellMask(tracker.GridDat, voidCellsArray);
                foreach (Chunk chunk in emptyCells) {
                    foreach (int cell in chunk.Elements) {
                        QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, RefElement.SpatialDimension);
                        emptyRule.Nodes.LockForever();
                        result.Add(new ChunkRulePair<QuadRule>(
                            Chunk.GetSingleElementChunk(cell), emptyRule));
                    }
                }

                CellMask fullCells = new CellMask(tracker.GridDat, fullCellsArray);
                foreach (Chunk chunk in fullCells) {
                    foreach (int cell in chunk.Elements) {
                        QuadRule fullRule = RefElement.GetQuadratureRule(order);
                        result.Add(new ChunkRulePair<QuadRule>(
                            Chunk.GetSingleElementChunk(cell), fullRule));
                    }
                }

                CellMask realCutCells = cellMask.Except(emptyCells).Except(fullCells);
                if (RestrictNodes) {
                    foreach (Chunk chunk in realCutCells) {
                        foreach (int cell in chunk.Elements) {
                            CellMask singleElementMask = new CellMask(
                                tracker.GridDat, Chunk.GetSingleElementChunk(cell));

                            AffineTrafo trafo = trafos[localCellIndex2SubgridIndex[cell]];
                            Debug.Assert(Math.Abs(trafo.Matrix.Determinant()) > 1e-10);

                            NodeSet nodes = GetNodes(noOfLambdas).CloneAs();
                            NodeSet mappedNodes = new NodeSet(RefElement, trafo.Transform(nodes));
                            mappedNodes.LockForever();

                            // Remove nodes in negative part
                            MultidimensionalArray levelSetValues = LevelSetData.GetLevSetValues(mappedNodes, cell, 1);
                            List<int> nodesToBeCopied = new List<int>(mappedNodes.GetLength(0));
                            for (int n = 0; n < nodes.GetLength(0); n++) {
                                if (levelSetValues[0, n] >= 0.0) {
                                    nodesToBeCopied.Add(n);
                                }
                            }

                            NodeSet reducedNodes = new NodeSet(
                                this.RefElement, nodesToBeCopied.Count, D);
                            for (int n = 0; n < nodesToBeCopied.Count; n++) {
                                for (int d = 0; d < D; d++) {
                                    reducedNodes[n, d] = mappedNodes[nodesToBeCopied[n], d];
                                }
                            }
                            reducedNodes.LockForever();

                            QuadRule optimizedRule = GetOptimizedRule(
                                cell,
                                trafo,
                                reducedNodes,
                                quadResults,
                                order);

                            result.Add(new ChunkRulePair<QuadRule>(
                                singleElementMask.Single(), optimizedRule));
                        }
                    }
                } else {
                    // Use same nodes in all cells
                        QuadRule[] optimizedRules = GetOptimizedRules(
                            realCutCells, GetNodes(noOfLambdas), quadResults, order);                 
                    int ruleIndex = 0;
                    foreach (Chunk chunk in realCutCells) {
                        foreach (var cell in chunk.Elements) {
                            result.Add(new ChunkRulePair<QuadRule>(
                                Chunk.GetSingleElementChunk(cell), optimizedRules[ruleIndex]));
                            ruleIndex++;
                        }
                    }
                }

                cachedRules[order] = result.OrderBy(p => p.Chunk.i0).ToArray();
                return cachedRules[order];
            }

        }

        private NodeSet GetNodes(int noOfLambdas) {
            bool gaussianRuleFound = false;

            NodeSet nodes = null;
            if (UseGaussNodes) {
                int minOrder = 1;
                gaussianRuleFound = true;
                while (true) {
                    if (minOrder > RefElement.HighestKnownOrder) {
                        gaussianRuleFound = false;
                        break;
                    }

                    if (RefElement.GetQuadratureRule(minOrder).NoOfNodes < NodeCountSafetyFactor * noOfLambdas) {
                        minOrder++;
                    } else {
                        break;
                    }
                }

                if (gaussianRuleFound) {
                    nodes = RefElement.GetQuadratureRule(minOrder).Nodes;
                    if (nodes.GetLength(0) < noOfLambdas) {
                        throw new Exception();
                    }
                }
            }

            if (!gaussianRuleFound) {
                int D = RefElement.SpatialDimension;
                double targetNumber = NodeCountSafetyFactor * noOfLambdas;
                int noOfNodesPerDirection = (int)Math.Ceiling(Math.Pow(targetNumber, 1.0 / D));
                double[] linearNodes = GenericBlas.Linspace(-1.0, 1.0, noOfNodesPerDirection);

                nodes = new NodeSet(RefElement, noOfNodesPerDirection * noOfNodesPerDirection, D);
                int node = 0;
                for (int i = 0; i < noOfNodesPerDirection; i++) {
                    for (int j = 0; j < noOfNodesPerDirection; j++) {
                        nodes[node, 0] = linearNodes[i];
                        nodes[node, 1] = linearNodes[j];
                        node++;
                    }
                }
                nodes.LockForever();
            }

            return nodes;
        }

        /// <summary>
        /// If there are any cached rules, this method returns their order.
        /// </summary>
        public int[] GetCachedRuleOrders() {
            return cachedRules.Keys.ToArray();
        }

        /// <summary>
        /// Uses a moment-fitting basis of order <paramref name="order"/> for
        /// determining optimal weights of a quadrature rule for a sub-region
        /// of each cell in the given <paramref name="mask"/>.
        /// </summary>
        /// <param name="mask">
        /// Cells for which rules shall be created
        /// </param>
        /// <param name="nodes">Integration nodes</param>
        /// <param name="quadResults">
        /// Integration results for each function in the
        /// <see cref="lambdaBasis"/>
        /// </param>
        /// <param name="order">
        /// Desired order of the moment-fitting system. Assuming that
        /// <see cref="surfaceRuleFactory"/> integrates the basis
        /// polynomials exactly over the zero iso-contour (which it usually
        /// doesn't), the resulting quadrature rules will be exact up to this
        /// order.
        /// </param>
        /// <returns>
        /// A set of quadrature rules
        /// </returns>
        private QuadRule[] GetOptimizedRules(CellMask mask, NodeSet nodes, double[,] quadResults, int order) {
            using (var tr = new FuncTrace()) {
                int maxLambdaDegree = order + 1;
                int noOfLambdas = GetNumberOfLambdas(maxLambdaDegree);
                int noOfNodes = nodes.GetLength(0);

                QuadRule[] optimizedRules = new QuadRule[mask.NoOfItemsLocally];

                if (mask.NoOfItemsLocally ==0)
                        return optimizedRules;
                // Leading dimension of B (rhs); required by DGELSY
                int LDB = Math.Max(noOfLambdas, noOfNodes);
                double[] rhs = new double[LDB * mask.NoOfItemsLocally];

                Basis basis = new Basis(tracker.GridDat, order);
                MultidimensionalArray basisValues = basis.Evaluate(nodes);

                int n = 0;
                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        int iSubGrid = localCellIndex2SubgridIndex[cell];

                        for (int k = 0; k < noOfLambdas; k++) {
                            rhs[n * LDB + k] += quadResults[iSubGrid, k];
                        }

                        n++;
                    }
                }
                
                double[] matrix;
                if (basisValues.IsContinious) {
                    matrix = basisValues.Storage;
                } else {
                    matrix = new double[basisValues.Length];
                    basisValues.CopyTo(matrix, true, 0);
                }

                
                LAPACK.F77_LAPACK.DGELSY(noOfLambdas, noOfNodes, matrix, rhs, mask.NoOfItemsLocally, RCOND);

                n = 0;
                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        optimizedRules[n] = new QuadRule() {
                            Nodes = nodes,
                            Weights = MultidimensionalArray.Create(noOfNodes),
                            OrderOfPrecision = order
                        };

                        for (int j = 0; j < noOfNodes; j++) {
                            optimizedRules[n].Weights[j] = rhs[n * LDB + j];
                        }

                        double max = optimizedRules[n].Weights.Max(d => d.Abs());
                        if (max > 2.0 * RefElement.Volume) {
                            tr.Info(String.Format(
                                "Warning: Abnormally large integration weight detected"
                                + " for level set volume integral in cell {0}"
                                + " (|w| = {1}). This may indicate a loss of"
                                + " integration accuracy.",
                                cell,
                                max));
                        }

                        n++;
                    }
                }

                return optimizedRules;
            }
        }

        /// <summary>
        /// Non-vectorized reference implementation of
        /// <see cref="GetOptimizedRule(int, AffineTrafo, NodeSet, double[,], int)"/>
        /// </summary>
        /// <param name="cell"></param>
        /// <param name="trafo"></param>
        /// <param name="nodes"></param>
        /// <param name="quadResults"></param>
        /// <param name="order"></param>
        /// <returns></returns>
        private QuadRule GetOptimizedRule(int cell, AffineTrafo trafo, NodeSet nodes, double[,] quadResults, int order) {
            int maxLambdaDegree = order + 1;
            int noOfLambdas = GetNumberOfLambdas(maxLambdaDegree);
            int noOfNodes = nodes.GetLength(0);

            // Leading dimension of B (rhs); required by DGELSY
            int LDB = Math.Max(noOfLambdas, noOfNodes);
            double[] rhs = new double[LDB];

            AffineTrafo inverseTrafo = trafo.Invert();
            NodeSet trafoNodes = new NodeSet(RefElement, inverseTrafo.Transform(nodes));
            trafoNodes.LockForever();

            Basis basis = new Basis(tracker.GridDat, order);
            MultidimensionalArray basisValues = basis.Evaluate(trafoNodes);

            int iSubGrid = localCellIndex2SubgridIndex[cell];
            for (int k = 0; k < noOfLambdas; k++) {
                rhs[k] += quadResults[iSubGrid, k];
            }

            LAPACK.F77_LAPACK.DGELSY(noOfLambdas, noOfNodes, basisValues.Storage, rhs, 1, RCOND);

            QuadRule optimizedRule = new QuadRule() {
                Nodes = nodes,
                Weights = MultidimensionalArray.Create(noOfNodes),
                OrderOfPrecision = order
            };

            for (int j = 0; j < noOfNodes; j++) {
                optimizedRule.Weights[j] = rhs[j];
            }

            return optimizedRule;
        }

        /// <summary>
        /// Computes the number of vector-valued basis function (cf.
        /// <see cref="lambdaBasis"/>)
        /// </summary>
        /// <returns></returns>
        private int GetNumberOfLambdas(int maxLambdaDegree) {
            switch (tracker.GridDat.SpatialDimension) {
                case 1:
                    return maxLambdaDegree;

                case 2:
                    return maxLambdaDegree * (maxLambdaDegree + 1) / 2;

                case 3:
                    return maxLambdaDegree * (maxLambdaDegree + 1) * (maxLambdaDegree + 2) / 6;

                default:
                    throw new ApplicationException("Impossible");
            }
        }

        /// <summary>
        /// Evaluates the vector-valued anti-derivatives defining the
        /// moment-fitting basis (cf. <see cref="lambdaBasis"/>) in each node
        /// of all cells in the given range
        /// </summary>
        /// <param name="i0">
        /// First cell in range
        /// </param>
        /// <param name="length">
        /// Number of cells
        /// </param>
        /// <returns>
        /// The values of <see cref="lambdaBasis"/> in each node
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
        /// <param name="NS">
        /// Nodes at which to evaluate.
        /// </param>
        private MultidimensionalArray EvaluateLambdas(int cell, NodeSet NS) {
            int D = tracker.GridDat.SpatialDimension;
            Debug.Assert(
                lambdaBasis.Count % D == 0,
                "Number of polynomials in basis should be divisible by D = " + D);

            int noOfLambdas = lambdaBasis.Count / D;
            int noOfNodes = NS.NoOfNodes;

            if (RestrictNodes) {
                AffineTrafo trafo = trafos[localCellIndex2SubgridIndex[cell]];

                AffineTrafo inverse = trafo.Invert();
                NS = new NodeSet(RefElement, inverse.Transform(NS));
                NS.LockForever();

                MultidimensionalArray lambdaValues = lambdaBasis.Values.GetValues(NS);
                lambdaValues = lambdaValues.ResizeShallow(noOfNodes, noOfLambdas, D);

                for (int i = 0; i < noOfNodes; i++) {
                    for (int j = 0; j < noOfLambdas; j++) {
                        for (int d = 0; d < D; d++) {

                            // Bounding box transformation is assumed to just a
                            // stretching, i.e. off-diagonals are zero
                            lambdaValues[i, j, d] *= trafo.Matrix[d, d];
                        }
                    }
                }

                return lambdaValues;
            } else {
                MultidimensionalArray lambdaValues = lambdaBasis.Values.GetValues(NS);
                return lambdaValues.ResizeShallow(noOfNodes, noOfLambdas, D);
            }
        }

        public void InvalidateCache() {
            cachedRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();
        }

        #endregion

        /// <summary>
        /// Used for the computation of the integral of the basis functions
        /// (see <see cref="lambdaBasis"/>) over the boundary of a cell. More
        /// specifically, the boundary of the cell is given by
        /// (D-1)-dimensional elements that are potentially cut by the zero
        /// iso-contour of the level set, which is why this requires a special
        /// treatment; see <see cref="surfaceRuleFactory"/>.
        /// </summary>
        private class LambdaCellBoundaryQuadrature : CellBoundaryQuadrature<CellBoundaryQuadRule> {

            /// <summary>
            /// Owner of this object
            /// </summary>
            private LevelSetVolumeQuadRuleFactory owner;

            /// <summary>
            /// Number of items in the current execution mask
            /// </summary>
            private int NoOfItemsLocally;

            /// <summary>
            /// <see cref="SubGrid.LocalCellIndex2SubgridIndex"/>
            /// </summary>
            private int[] localCellIndex2SubgridIndex;

            /// <summary>
            /// Integration results for the last call to <see cref="Execute"/>
            /// <list type="bullet">
            ///     <item>
            ///     1st index: Cell index within
            ///     <see cref="localCellIndex2SubgridIndex"/>
            ///     </item>
            ///     <item>
            ///     3rd index: Basis function (see
            ///     <see cref="LevelSetVolumeQuadRuleFactory.EvaluateLambdas"/>
            ///     </item>
            /// </list>
            /// </summary>
            public double[,] Results {
                get;
                private set;
            }

            /// <summary>
            /// 
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="edgeRuleFactory"></param>
            /// <param name="mask"></param>
            public LambdaCellBoundaryQuadrature(
                LevelSetVolumeQuadRuleFactory owner, IQuadRuleFactory<CellBoundaryQuadRule> edgeRuleFactory, CellMask mask)
                : base(new int[] { owner.GetNumberOfLambdas(owner.lambdaBasis.MaxAbsoluteDegree) },
                       owner.tracker.GridDat,
                       (new CellBoundaryQuadratureScheme(edgeRuleFactory, mask)).Compile(owner.tracker.GridDat, owner.lambdaBasis.MaxAbsoluteDegree + 1),
                       CoordinateSystem.Reference) //
            {
                this.owner = owner;
                NoOfItemsLocally = mask.NoOfItemsLocally;
            }

            /// <summary>
            /// Computes the integrals while overwriting <see cref="Results"/>
            /// </summary>
            public override void Execute() {
                Results = new double[NoOfItemsLocally, m_TotalNoOfIntegralsPerItem];
                localCellIndex2SubgridIndex = owner.localCellIndex2SubgridIndex;
                base.Execute();
            }

            /// <summary>
            /// For each cell \f$ K\f$  in the given
            /// range and for each \f$ \vec{\Lambda}\f$  
            /// in <see cref="LevelSetEdgeVolumeQuadRuleFactory.lambdaBasis"/>:
            /// Computes
            /// \f$ 
            /// \int \limits_{\partial K} \vec{\Lambda} \cdot \vec{n} H(\varphi) \;ds,
            /// \f$ 
            /// where \f$ \vec{n}\f$  denotes the outer
            /// unit normal vector on \f$ \partial K\f$ .
            /// Moreover, \f$ H\f$  a weight function that
            /// depends on the level set function
            /// \f$ \varphi\f$ . Typically,
            /// \f$ H\f$  is given by the Heaviside
            /// function. For more details, see
            /// <see cref="CutLineQuadRuleFactory"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="QR"></param>
            /// <param name="EvalResult"></param>
            protected override void Evaluate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                int noOfEdges = owner.tracker.GridDat.Grid.RefElements[CurrentRuleRefElementIndex].NoOfFaces;
                int D = gridData.SpatialDimension;

                for (int i = 0; i < Length; i++) {
                    int nodeIndex = -1;
                    MultidimensionalArray lambdaValues = owner.EvaluateLambdas(i0 + i, QuadNodes);

                    for (int e = 0; e < noOfEdges; e++) {
                        for (int j = 0; j < CurrentRule.NumbersOfNodesPerFace[e]; j++) {
                            nodeIndex++;

                            for (int k = 0; k < m_TotalNoOfIntegralsPerItem; k++) {
                                for (int d = 0; d < D; d++) {
                                    EvalResult[i, nodeIndex, k] += lambdaValues[nodeIndex, k, d] *
                                        owner.RefElement.FaceNormals[e, d];
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
                int noOfFaces = owner.tracker.GridDat.Grid.RefElements[CurrentRuleRefElementIndex].NoOfFaces;

                for (int i = 0; i < Length; i++) {
                    int cell = i0 + i;
                    int iSubGrid = localCellIndex2SubgridIndex[cell];

                    for (int e = 0; e < noOfFaces; e++) {
                        for (int k = 0; k < m_TotalNoOfIntegralsPerItem; k++) {
                            Results[iSubGrid, k] += ResultsOfIntegration[i, e, k];
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Used for the computation of the integral of the basis functions
        /// (see <see cref="lambdaBasis"/>) over the intersection of the zero
        /// iso-contour of the level set with an edge. This requires a special
        /// treatment; see <see cref="surfaceRuleFactory"/>.
        /// </summary>
        private class LambdaLevelSetSurfaceQuadrature : CellQuadrature {

            /// <summary>
            /// Owner of this object
            /// </summary>
            private LevelSetVolumeQuadRuleFactory owner;

            /// <summary>
            /// Number of items in the current sub-grid
            /// </summary>
            private int NoOfItemsLocally;

            /// <summary>
            /// <see cref="SubGrid.LocalCellIndex2SubgridIndex"/>
            /// </summary>
            private int[] localCellIndex2SubgridIndex;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="surfaceRuleFactory"></param>
            /// <param name="mask"></param>
            public LambdaLevelSetSurfaceQuadrature(
                LevelSetVolumeQuadRuleFactory owner, IQuadRuleFactory<QuadRule> surfaceRuleFactory, CellMask mask)
                : base(
                    new int[] { owner.GetNumberOfLambdas(owner.lambdaBasis.MaxAbsoluteDegree) },
                    owner.tracker.GridDat,
                    new CellQuadratureScheme(surfaceRuleFactory, mask).Compile(owner.tracker.GridDat, owner.lambdaBasis.MaxAbsoluteDegree),
                    CoordinateSystem.Reference) //
            {
                this.owner = owner;
                NoOfItemsLocally = mask.NoOfItemsLocally;
            }

            /// <summary>
            /// Results of the last integration
            /// <list type="bullet">
            ///     <item>
            ///     1st index: Cell index within the current sub-grid
            ///     </item>
            ///     <item>
            ///     2nd index: Basis function index
            ///     </item>
            /// </list>
            /// </summary>
            public double[,] Results {
                get;
                private set;
            }

            /// <summary>
            /// Performs the integration while overwriting old results in
            /// <see cref="Results"/>
            /// </summary>
            public override void Execute() {
                Results = new double[NoOfItemsLocally, m_TotalNoOfIntegralsPerItem];
                localCellIndex2SubgridIndex = owner.localCellIndex2SubgridIndex;
                base.Execute();
            }

            /// <summary>
            /// For each cell \f$ K\f$   in the given
            /// range and for each \f$ \vec{\Lambda}\f$ 
            /// in <see cref="LevelSetVolumeQuadRuleFactory.lambdaBasis"/>:
            /// Computes
            /// \f$ 
            /// \int \limits_{\{\vec{x}; \varphi(\vec{x}) = 0 \}  \cap K} \vec{\Lambda} \cdot \vec{n}_I \;ds,
            /// \f$ 
            /// where \f$ \varphi\f$  is the level set
            /// function and \f$ \vec{n}_I\f$  denotes the
            /// unit normal vector on \f$ \varphi\f$ 
            /// </summary>
            protected override void Evaluate(int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                int D = gridData.SpatialDimension;
                int NoOfNodes = QuadNodes.NoOfNodes;
                
                MultidimensionalArray levelSetNormals = owner.tracker.GetLevelSetReferenceNormals(owner.levelSetIndex, QuadNodes, i0, Length);
                MultidimensionalArray metrics = owner.tracker.GetLevelSetNormalReferenceToPhysicalMetrics(
                    owner.levelSetIndex, QuadNodes, i0, Length);

                for (int i = 0; i < Length; i++) {
                    MultidimensionalArray lambdaValues = owner.EvaluateLambdas(i0 + i, QuadNodes);

                    for (int j = 0; j < NoOfNodes; j++) {
                        for (int k = 0; k < m_TotalNoOfIntegralsPerItem; k++) {
                            for (int d = 0; d < D; d++) {
                                EvalResult[i, j, k] += lambdaValues[j, k, d] * levelSetNormals[i, j, d];
                            }

                            EvalResult[i, j, k] *= metrics[i, j];
                        }
                    }
                }
            }

            /// <summary>
            /// Saves the results of the integration to <see cref="Results"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    int cell = i0 + i;
                    int iSubGrid = localCellIndex2SubgridIndex[cell];

                    for (int k = 0; k < m_TotalNoOfIntegralsPerItem; k++) {
                        Results[iSubGrid, k] = ResultsOfIntegration[i, k];
                    }
                }
            }
        }
    }
}
