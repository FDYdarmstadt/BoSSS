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
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.XDG.Quadrature.HMF {

    /// <summary>
    /// This factory uses the hierarchical moment-fitting (HMF) strategy in
    /// order to produce quadrature rules which are, for each cell
    /// \f$ K\f$  in a volume mask, capable of computing
    /// (an approximation of)
    /// \f[ 
    ///    \oint\limits_{\partial K \cap \{ \vec{x}; \varphi(\vec{x}) = 0 \} }  f \;dS,
    /// \f]
    /// where \f$ \varphi\f$  denotes the level set
    /// function.
    /// </summary>
    /// <remarks>
    /// For details about the algorithm, see
    /// Mueller, Kummer &amp; Oberlack: Highly accurate surface and volume
    /// integration on implicit domains by means of moment-fitting (2013)
    /// </remarks>
    public class LevelSetSurfaceQuadRuleFactory : IQuadRuleFactory<QuadRule> {

        ///// <summary>
        ///// Tracker containing the level set to be integrated over
        ///// </summary>
        //private LevelSetTracker tracker;

        /// <summary>
        /// Quadrature rule factory for the integration over the edges of the
        /// volume simplex, cf. <see cref="PhiQuadrature"/>
        /// </summary>
        private IQuadRuleFactory<CellBoundaryQuadRule> edgeRuleFactory;

        /// <summary>
        /// Mapping between local cell index and indices into current cell mask
        /// (cf. <see cref="GetQuadRuleSet(ExecutionMask, int)"/>)
        /// </summary>
        private int[] localCellIndex2SubgridIndex;

        /// <summary>
        /// A vector-valued basis consisting of divergence-free vectors of
        /// polynomials. These are required for the construction of the
        /// moment-fitting system, cf. <see cref="GetOptimizedRules"/>
        /// </summary>
        private DivergenceFreeBasis phiBasis;

        /// <summary>
        /// Index of the considered level set within <see cref="tracker"/>.
        /// </summary>
        private int levelSetIndex;

        /// <summary>
        /// Transformations from the volume coordinate system of a cell to
        /// the bounding box defined by the vertices with positive level set
        /// values and the roots of the level set
        /// </summary>
        private AffineTrafo[] trafosToBoundingBox;

        /// <summary>
        /// Indicates whether the integration nodes should be projected onto
        /// the zero level set, or not. The latter approach is much more
        /// efficient (since the projection is expensive), but also less
        /// accurate.
        /// </summary>
        public static bool UseNodesOnLevset;

        /// <summary>
        /// If <see cref="UseNodesOnLevset"/> is true, this switch
        /// indicates whether projected nodes that land outside of the
        /// corresponding cell should be kept (false) or discarded (true).
        /// If <see cref="UseNodesOnLevset"/> is false, this switch indicates
        /// whether nodes located outside of the sub-cell of interest should
        /// be kept (false) or discarded (true)
        /// </summary>
        public static bool RestrictNodes = false;

        /// <summary>
        /// 
        /// </summary>
        public static bool UseGaussNodes = true;

        private const double RCOND = 1e-11;

        /// <summary>
        /// Cache for quad rules
        /// </summary>
        private Dictionary<int, ChunkRulePair<QuadRule>[]> cachedRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="levelSetData">
        /// </param>
        /// <param name="edgeRuleFactory">
        /// Quadrature rule factory for the integration over the edges of the
        /// volume simplex. Note that the related quadrature rules must be able
        /// to cope with discontinuous integrands. For example, see
        /// <see cref="CutLineQuadRuleFactory"/> and
        /// <see cref="LevelSetEdgeVolumeQuadRuleFactory"/>
        /// </param>
        public LevelSetSurfaceQuadRuleFactory(
            LevelSetTracker.LevelSetData levelSetData,
            IQuadRuleFactory<CellBoundaryQuadRule> edgeRuleFactory) {

            if(levelSetData.GridDat.Cells.RefElements.Length > 1)
                throw new NotSupportedException();

            this.LevelSetData = levelSetData;
            if (!object.ReferenceEquals(RefElement, edgeRuleFactory.RefElement)) {
                throw new ArgumentException();
            }
            this.edgeRuleFactory = edgeRuleFactory;
                        
            this.levelSetIndex = levelSetData.LevelSetIndex;
        }

        LevelSetTracker.LevelSetData LevelSetData;

        #region IQuadRuleFactory<QuadRule>

        /// <summary>
        /// The first reference element of the current grid.
        /// </summary>
        public RefElement RefElement {
            get {
                return this.LevelSetData.GridDat.Grid.RefElements[0];
            }
        }

        /// <summary>
        /// For each cell in the given <paramref name="mask"/>: Constructs a
        /// quadrature rule that can be used to evaluate
        /// \f[ 
        ///    \oint\limits_{\partial K \cap \{ \vec{x}; \varphi(\vec{x}) = 0 \} }  f \;dS,
        /// \f]
        /// with relatively high accuracy (depending on the selected integration
        /// <paramref name="order"/>)
        /// </summary>
        /// <param name="mask">
        /// Cells to be integrated over
        /// </param>
        /// <param name="order">
        /// The desired order of the moment-fitting system. Note that the
        /// integration result will not be exact, even for polynomial functions
        /// f with degree &lt; <paramref name="order"/>. However, increasing
        /// the will typically strongly increase the integration accuracy
        /// </param>
        /// <returns>
        /// A set of quadrature rules. Note that the quadrature nodes will not
        /// necessarily coincide with the zero level set, depending on the
        /// arguments passed to the constructor.
        /// </returns>
        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
            using (var tr = new FuncTrace()) {
                // Ansatz order for moment-fitting system is one order higher
                // than requested order because the moment-fitting is based on
                // 'anti-derivatives' of the real integrands
                order += 1;

               

                CellMask cellMask = mask as CellMask;
                if (mask == null) {
                    throw new ArgumentException("Cell mask required", "mask");
                }
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");


                localCellIndex2SubgridIndex = cellMask.GetIndex2MaskItemMap();

                phiBasis = new DivergenceFreeBasis(LevelSetData.GridDat, this.RefElement, order);
                int noOfPhis = this.NumberOfMoments;
                int D = RefElement.SpatialDimension;

                if (RestrictNodes) {
                    trafosToBoundingBox = ConstructTransformationsToBoundingBox(cellMask, order);
                }

                PhiQuadrature edgeQuadrature = new PhiQuadrature(this, cellMask);
                edgeQuadrature.Execute();

                // Must happen _after_ all parallel calls (e.g., definition of
                // the sub-grid or quadrature) in order to avoid problems in
                // parallel runs
                if (mask.NoOfItemsLocally == 0) {
                    ChunkRulePair<QuadRule>[] empty = new ChunkRulePair<QuadRule>[0];
                    return empty;
                }

                // Uncommented for FSI after speaking to Florian
                //______________________________________
                /*
                if (cachedRules.Count > 0 && cachedRules.Keys.Max() >= order) {
                    order = cachedRules.Keys.Where(cachedOrder => cachedOrder >= order).Min();
                    CellMask cachedMask = new CellMask(mask.GridData, cachedRules[order].Select(p => p.Chunk).ToArray());

                    if (cachedMask.Equals(mask)) {
                        Console.WriteLine("Cached rule: " + order);
                        return cachedRules[order];
                    } else {
                        throw new NotImplementedException(
                            "Case not yet covered yet in combination with caching; deactivate caching to get rid of this message");
                    }
                }
                */

                var result = new List<ChunkRulePair<QuadRule>>(mask.NoOfItemsLocally);
                if (UseNodesOnLevset) {
                    // use only quadrature Nodes on the zero-level-set (by closest-point finding)
                    // => different nodes for each cell
                    result.AddRange(GetOptimizedQuadRulesWithNodesOnLevelSet(
                        cellMask, order, edgeQuadrature.IntegrationResults));
                } else {
                    // => quadrature nodes are NOT on the zero-level-set
                    if (RestrictNodes) {
                        // => Different nodes in all cells
                        result.AddRange(GetOptimizedQuadRulesWithNodesInPositiveRegion(
                            cellMask, order, edgeQuadrature.IntegrationResults));
                    } else {
                        // => Same nodes in all cells
                        NodeSet baseNodes = GetSeedNodes(noOfPhis);
                        foreach (Chunk chunk in mask) {
                            QuadRule[] optimizedRules = GetOptimizedRules(
                                baseNodes, order, chunk.i0, chunk.Len, edgeQuadrature.IntegrationResults, levelSetIndex);

                            for (int i = 0; i < chunk.Len; i++) {
                                result.Add(new ChunkRulePair<QuadRule>(
                                    Chunk.GetSingleElementChunk(i + chunk.i0), optimizedRules[i]));
                            }
                        }
                    }
                }

                cachedRules[order] = result.ToArray();


                return result;
            }
        }

        /// <summary>
        /// Minimum desired ratio between number of nodes and number of moment
        /// equations. If not set, a value that is determined from numerical
        /// experiments will be used. The value for the 2D case has been
        /// verified independently by Martin Gehrke. The 3D value is still not
        /// fully settled.
        /// </summary>
        private double nodeCountSafetyFactor {
            get {
                if (RestrictNodes) {
                    switch (RefElement.SpatialDimension) {
                        case 2:
                            return 5.0;

                        case 3:
                            throw new NotImplementedException();

                        default:
                            throw new NotImplementedException();
                    }
                } else {
                    switch (RefElement.SpatialDimension) {
                        case 2:
                            return 1.6;

                        case 3:
                            return 3.0;

                        default:
                            throw new NotImplementedException();
                    }
                }
            }
        }

        private NodeSet GetSeedNodes(int noOfPhis) {
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

                    if (RefElement.GetQuadratureRule(minOrder).NoOfNodes < nodeCountSafetyFactor * noOfPhis) {
                        minOrder++;
                    } else {
                        break;
                    }
                }

                if (gaussianRuleFound) {
                    nodes = RefElement.GetQuadratureRule(minOrder).Nodes;
                    if (nodes.GetLength(0) < noOfPhis) {
                        throw new Exception();
                    }
                }
            }

            if (!gaussianRuleFound) {
                int D = RefElement.SpatialDimension;
                double targetNumber = nodeCountSafetyFactor * noOfPhis;
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
        /// Projects the given <paramref name="Nodes"/> onto the zero level-set
        /// of cell <paramref name="jCell"/>
        /// </summary>
        /// <param name="jCell">Local cell index</param>
        /// <param name="Nodes">
        /// The nodes to be projected
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Spatial dimension</item>
        /// </list>
        /// </param>
        /// <returns>
        /// The projection of all nodes in <paramref name="Nodes"/> onto the
        /// zero level-set
        /// </returns>
        private NodeSet ProjectNodesOntoLevelSet(int jCell, NodeSet Nodes) {
            int D = Nodes.GetLength(1);
            int NoOfNodes = Nodes.GetLength(0);
            var m_Context = this.LevelSetData.GridDat;
            LevelSet LevSet = (LevelSet)( this.LevelSetData.LevelSet);
            RefElement Kref = m_Context.Cells.GetRefElement(jCell);

            MultidimensionalArray LevSetValues = MultidimensionalArray.Create(1, NoOfNodes);
            MultidimensionalArray LevSetGrad = MultidimensionalArray.Create(1, NoOfNodes, D);

            MultidimensionalArray x0_i_Local = MultidimensionalArray.Create(1, NoOfNodes, D);
            MultidimensionalArray x0_i_Global = MultidimensionalArray.Create(1, NoOfNodes, D); // quadrature nodes in global coordinates

            MultidimensionalArray x0_ip1_Local = MultidimensionalArray.Create(1, NoOfNodes, D);
            MultidimensionalArray x0_ip1_Global = MultidimensionalArray.Create(NoOfNodes, D);

            // set initial value;
            x0_i_Local.SetSubArray(Nodes, 0, -1, -1);

            int NN = NoOfNodes;
            for (int i = 0; i < 10; i++) {
                double radiusError = 0;
                int j = jCell;

                LevSet.Evaluate(j, 1, Nodes, LevSetValues, 0, 0.0);
                LevSet.EvaluateGradient(j, 1, Nodes, LevSetGrad);

                m_Context.TransformLocal2Global(new NodeSet(this.RefElement, x0_i_Local.ExtractSubArrayShallow(0, -1, -1)), j, 1, x0_i_Global, 0);

                for (int nn = 0; nn < NN; nn++) {

                    double sc = 0;
                    for (int d = 0; d < D; d++) {
                        sc += LevSetGrad[0, nn, d].Pow2();
                    }


                    for (int d = 0; d < D; d++) {
                        double xd = x0_i_Global[0, nn, d]
                            - LevSetGrad[0, nn, d] * LevSetValues[0, nn] / sc;
                        x0_ip1_Global[nn, d] = xd;
                    }

                    radiusError += Math.Abs(LevSetValues[0, nn]);

                }

                m_Context.TransformGlobal2Local(x0_ip1_Global, x0_ip1_Local, j, 1, 0);

                //Console.WriteLine("iter #" + i + ", radiusError = " + radiusError.ToStringDot());
                //tr.Info("iter #" + i + ", radiusError = " + radiusError.ToStringDot());

                // next iter: x0_i <- x0_{i+1}
                x0_i_Local.Set(x0_ip1_Local);
                Nodes = new NodeSet(Kref, x0_i_Local.ExtractSubArrayShallow(0, -1, -1));
            }


            if (RestrictNodes) {
                var allPoints = x0_i_Local.ExtractSubArrayShallow(0, -1, -1);
                bool[] InOrOut = new bool[allPoints.GetLength(0)];
                int Incnt = 0;

                for (int k = 0; k < allPoints.GetLength(0); k++) {

                    double[] pt = new double[D];
                    for (int d = 0; d < D; d++) {
                        pt[d] = allPoints[k, d];
                    }

                    InOrOut[k] = this.RefElement.IsWithin(pt);
                    if (InOrOut[k])
                        Incnt++;
                }

                var ret = MultidimensionalArray.Create(Incnt, D);
                int i = 0;
                for (int k = 0; k < allPoints.GetLength(0); k++) {
                    if (InOrOut[k]) {
                        for (int d = 0; d < D; d++) {
                            ret[i, d] = allPoints[k, d];
                        }
                        i++;
                    }
                }

                return new NodeSet(Kref, ret);
            } else {
                return new NodeSet(Kref, x0_i_Local.ExtractSubArrayShallow(0, -1, -1));
            }
        }

        #endregion

        /// <summary>
        /// Determines the (in some sense) optimal weights for the given
        /// quadrature <paramref name="nodes"/> by solving the moment-fitting
        /// system of degree <paramref name="phiDegree"/> in each cell in the
        /// given range of cells, where right-hand side follows from the given
        /// <paramref name="quadResults"/>.
        /// </summary>
        /// <param name="nodes">Integration nodes</param>
        /// <param name="phiDegree">Order of the moment-fitting system</param>
        /// <param name="i0">First cell index</param>
        /// <param name="length">Number of cells</param>
        /// <param name="quadResults">
        /// Integration results for each function in the divergence-free basis
        /// <see cref="phiBasis"/>
        /// <list type="bullet">
        ///     <item>
        ///     1st index: Cell index within the current <see cref="subGrid"/>
        ///     </item>
        ///     <item>
        ///     2nd index: Local edge index of <see cref="RefElement"/>
        ///     </item>
        ///     <item>
        ///     3rd index: Basis function index
        ///     </item>
        /// </list>
        /// </param>
        /// <param name="levSetIndex"></param>
        /// <returns></returns>
        protected QuadRule[] GetOptimizedRules(
            NodeSet nodes,
            int phiDegree,
            int i0,
            int length,
            MultidimensionalArray quadResults,
            int levSetIndex) {

            using (var tr = new FuncTrace()) {
                int noOfNodes = nodes.GetLength(0);
                int noOfPhis = NumberOfMoments;

                MultidimensionalArray phis = EvaluatePhis(i0, length, nodes);
                MultidimensionalArray normals =
                    LevelSetData.GetLevelSetReferenceNormals(nodes, i0, length);
                MultidimensionalArray metrics =
                    LevelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(nodes, i0, length);

                if (!phis.IsContinious || !quadResults.IsContinious || !normals.IsContinious || !metrics.IsContinious) {
                    throw new NotImplementedException(
                        String.Format(
                            "This method assumes that all input arrays have a continuous memory layout, but we have"
                                + " phis.IsContinious={0}, quadResults.IsContinious={1}, normals.IsContinious={2}, metrics.IsContinious={3}",
                            phis.IsContinious,
                            quadResults.IsContinious,
                            normals.IsContinious,
                            metrics.IsContinious));
                }

                // Additional space required by Fortran routine
                double[] rhs = new double[Math.Max(noOfNodes, noOfPhis)];
                var _matrix = MultidimensionalArray.Create(noOfPhis, noOfNodes);
                double[] matrix = _matrix.Storage;
                Debug.Assert(_matrix.IsContinious && _matrix.Index(0, 0) == 0);
                QuadRule[] optimizedRules = new QuadRule[length];
                var nodesArray = nodes;
                var NoOfFaces = this.RefElement.NoOfFaces;
            
                unsafe
                {
                    fixed (double* pRhs = &rhs[0],
                        pMatrix = &matrix[0],
                        pQuad = &quadResults.Storage[0],
                        pPhis = &phis.Storage[0],
                        pNormals = &normals.Storage[0],
                        pMetrics = &metrics.Storage[0]) {

                        for (int i = 0; i < length; i++) {
                            int cell = i0 + i;
                            int iSubGrid = localCellIndex2SubgridIndex[cell];

                            Array.Clear(matrix, 0, matrix.Length);
                            Array.Clear(rhs, 0, rhs.Length);

                            double* pRhsCur = pRhs;
                            double rhsL2Norm = 0.0;
                            for (int k = 0; k < noOfPhis; k++) {
                                for (int e = 0; e < NoOfFaces; e++) {
                                    *pRhsCur += *(pQuad + quadResults.Index(iSubGrid, e, k));
                                }
                                rhsL2Norm += *pRhsCur * *pRhsCur;
                                pRhsCur++;
                            }

                            if (rhsL2Norm < 1e-14) {
                                // All integrals are zero => cell not really cut
                                // (happens e.g. if level set is tangent)
                                QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, RefElement.SpatialDimension);
                                emptyRule.Nodes.LockForever();
                                optimizedRules[i] = emptyRule;
                                continue;
                            }

                            double* pPhisCur = pPhis;
                            for (int j = 0; j < noOfNodes; j++) {
                                for (int k = 0; k < noOfPhis; k++) {
                                    double* pNormalsCur = pNormals + normals.Index(i, j, 0);
                                    // Beware of Fortran order!
                                    double* pMatrixCur = pMatrix + j * noOfPhis + k;

                                    for (int d = 0; d < RefElement.SpatialDimension; d++) {
                                        *pMatrixCur += *(pPhisCur++) * *(pNormalsCur++);
                                    }
                                }
                            }

                            LAPACK.F77_LAPACK.DGELSY(noOfPhis, noOfNodes, matrix, rhs, 1, RCOND);

                            QuadRule optimizedRule = new QuadRule() {
                                Nodes = nodesArray,
                                Weights = MultidimensionalArray.Create(noOfNodes),
                                OrderOfPrecision = phiDegree
                            };

                            double maxWeight = 0.0;
                            pRhsCur = pRhs;
                            double* pMetricsCur = pMetrics + metrics.Index(i, 0);
                            for (int j = 0; j < noOfNodes; j++) {
                                optimizedRule.Weights[j] = *(pRhsCur++) / *(pMetricsCur++);
                                maxWeight = Math.Max(optimizedRule.Weights[j].Abs(), maxWeight);
                            }

                            if (maxWeight > 5.0 * RefElement.Volume) {
                                tr.Info(String.Format(
                                    "Warning: Abnormally large, ugly integration weight detected"
                                    + " for level set surface integral in cell {0}"
                                    + " (|w| = {1}). This may indicate a loss of"
                                    + " integration accuracy.",
                                    cell,
                                    maxWeight));
                            }

                            optimizedRules[i] = optimizedRule;
                        }
                    }
                }

                return optimizedRules;
            }
        }

        /// <summary>
        /// Non-vectorized reference implementation of
        /// <see cref="GetOptimizedRules(NodeSet, int, int, int, MultidimensionalArray, int)"/>
        /// which is used for optimizing the weights in single cells
        /// </summary>
        /// <param name="nodes"></param>
        /// <param name="phiDegree"></param>
        /// <param name="jCell"></param>
        /// <param name="quadResults"></param>
        /// <param name="levSetIndex"></param>
        /// <returns></returns>
        protected QuadRule GetOptimizedRule(
            NodeSet nodes, int phiDegree, int jCell, MultidimensionalArray quadResults, int levSetIndex) {

            using (var tr = new FuncTrace()) {
                int noOfNodes = nodes.GetLength(0);
                int noOfPhis = NumberOfMoments;

                MultidimensionalArray phis = EvaluatePhis(jCell, nodes);
                MultidimensionalArray normals =
                    LevelSetData.GetLevelSetReferenceNormals(nodes, jCell, 1);
                MultidimensionalArray metrics =
                    LevelSetData.GetLevelSetNormalReferenceToPhysicalMetrics(nodes, jCell, 1);

                // Additional space required by Fortran routine
                double[] rhs = new double[Math.Max(noOfNodes, noOfPhis)];
                var _matrix = MultidimensionalArray.Create(noOfPhis, noOfNodes);
                double[] matrix = _matrix.Storage;
                Debug.Assert(_matrix.IsContinious && _matrix.Index(0, 0) == 0);
                var nodesArray = nodes;
                var NoOfFaces = this.RefElement.NoOfFaces;

                // Reference implementation
                int iSubGrid = localCellIndex2SubgridIndex[jCell];

                Array.Clear(matrix, 0, matrix.Length);
                Array.Clear(rhs, 0, rhs.Length);

                for (int e = 0; e < LevelSetData.GridDat.Grid.RefElements[0].NoOfFaces; e++) {
                    for (int k = 0; k < noOfPhis; k++) {
                        rhs[k] += quadResults[iSubGrid, e, k];
                    }
                }

                if (rhs.L2NormPow2() < 1e-14) {
                    // All integrals are zero => cell not really cut
                    // (happens e.g. if level set is tangent)
                    QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, RefElement.SpatialDimension);
                    emptyRule.Nodes.LockForever();
                    return emptyRule;
                }


                for (int j = 0; j < noOfNodes; j++) {
                    for (int k = 0; k < noOfPhis; k++) {
                        // Beware of Fortran order!
                        int index = j * noOfPhis + k;

                        for (int d = 0; d < RefElement.SpatialDimension; d++) {
                            matrix[index] += phis[j, k, d] * normals[0, j, d];
                        }
                    }
                }

                LAPACK.F77_LAPACK.DGELSY(noOfPhis, noOfNodes, matrix, rhs, 1, RCOND);

                QuadRule optimizedRule = new QuadRule() {
                    Nodes = nodesArray,
                    Weights = MultidimensionalArray.Create(noOfNodes)
                };

                for (int j = 0; j < noOfNodes; j++) {
                    optimizedRule.Weights[j] = rhs[j] / metrics[0, j];
                }

                double max = optimizedRule.Weights.Max(d => d.Abs());
                if (max > 5.0 * RefElement.Volume) {
                    tr.Info(String.Format(
                        //Console.WriteLine(String.Format(
                        "Warning: Abnormally large integration weight detected"
                        + " for level set surface integral in cell {0}"
                        + " (|w| = {1}). This may indicate a loss of"
                        + " integration accuracy.",
                        jCell,
                        max));
                }

                return optimizedRule;
            }
        }

        /// <summary>
        /// Evaluates <see cref="phiBasis"/> in all cells within the given
        /// range
        /// </summary>
        /// <param name="i0">First cell index</param>
        /// <param name="length">Number of cells</param>
        /// <returns>
        /// The vector-valued basis values
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Basis function index</item>
        ///     <item>3rd index: Spatial dimension</item>
        /// </list>
        /// </returns>
        /// <param name="NS">
        /// Node Set.
        /// </param>
        private MultidimensionalArray EvaluatePhis(int i0, int length, NodeSet NS) {
            int noOfNodes = NS.NoOfNodes;
            int D = LevelSetData.GridDat.Grid.SpatialDimension;
            int noOfPhis = phiBasis.Count / D;

            MultidimensionalArray phiValues = phiBasis.Values.GetValues(NS);
            return phiValues.ResizeShallow(noOfNodes, noOfPhis, D);
        }

        private MultidimensionalArray EvaluatePhis(int cell, NodeSet NS) {
            int noOfNodes = NS.NoOfNodes;
            int D = LevelSetData.GridDat.Grid.SpatialDimension;
            int noOfPhis = phiBasis.Count / D;

            if (RestrictNodes) {
                AffineTrafo trafo = trafosToBoundingBox[localCellIndex2SubgridIndex[cell]];

                if (trafo == null) {
                    return MultidimensionalArray.Create(noOfNodes, noOfPhis, D);
                }

                AffineTrafo inverse = trafo.Invert();
                NS = new NodeSet(RefElement, inverse.Transform(NS));
                NS.LockForever();

                MultidimensionalArray phiValues = phiBasis.Values.GetValues(NS);
                phiValues = phiValues.ResizeShallow(noOfNodes, noOfPhis, D);

                for (int i = 0; i < noOfNodes; i++) {
                    for (int j = 0; j < noOfPhis; j++) {
                        for (int d = 0; d < D; d++) {

                            // Bounding box transformation is assumed to just a
                            // stretching, i.e. off-diagonals are zero
                            phiValues[i, j, d] *= trafo.Matrix[d, d];
                        }
                    }
                }

                return phiValues;
            } else {
                MultidimensionalArray phiValues = phiBasis.Values.GetValues(NS);
                return phiValues.ResizeShallow(noOfNodes, noOfPhis, D);
            }
        }

        /// <summary>
        /// Determines the number of basis functions represented by
        /// <see cref="phiBasis"/>
        /// </summary>
        /// <returns></returns>
        private int NumberOfMoments {
            get {
                return phiBasis.Count / LevelSetData.GridDat.Grid.SpatialDimension;
            }
        }

        public void InvalidateCache() {
            cachedRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();
        }
        
        private IEnumerable<ChunkRulePair<QuadRule>> GetOptimizedQuadRulesWithNodesInPositiveRegion(
            CellMask mask, int order, MultidimensionalArray integrationResults) {

            var result = new List<ChunkRulePair<QuadRule>>(mask.NoOfItemsLocally);
            int noOfPhis = NumberOfMoments;
            int D = RefElement.SpatialDimension;

            // Use same nodes in all cells where possible, but use
            // nodes in positive part, only
            foreach (Chunk chunk in mask) {
                for (int jCell = chunk.i0; jCell < chunk.JE; jCell++) {
                    CellMask singleElementMask = new CellMask(
                        LevelSetData.GridDat, Chunk.GetSingleElementChunk(jCell));

                    NodeSet nodes = GetSeedNodes(noOfPhis).CloneAs();
                    AffineTrafo trafo = trafosToBoundingBox[localCellIndex2SubgridIndex[jCell]];

                    if (trafo == null) {
                        QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, 2);
                        emptyRule.Nodes.LockForever();
                        result.Add(new ChunkRulePair<QuadRule>(
                            singleElementMask.Single(), emptyRule));
                        continue;
                    }

                    NodeSet mappedNodes = new NodeSet(RefElement, trafo.Transform(nodes));
                    mappedNodes.LockForever();

                    // Remove nodes in negative part
                    MultidimensionalArray levelSetValues = LevelSetData.GetLevSetValues(mappedNodes, jCell, 1);
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

                    if (nodesToBeCopied.Count == 0) {
                        QuadRule emptyRule = QuadRule.CreateEmpty(RefElement, 1, 2);
                        emptyRule.Nodes.LockForever();
                        result.Add(new ChunkRulePair<QuadRule>(
                            singleElementMask.Single(), emptyRule));
                        continue;
                    }

                    QuadRule optimizedRule = GetOptimizedRule(
                        reducedNodes,
                        order,
                        jCell,
                        integrationResults,
                        levelSetIndex);

                    result.Add(new ChunkRulePair<QuadRule>(singleElementMask.Single(), optimizedRule));
                }
            }

            return result;
        }

        private List<ChunkRulePair<QuadRule>> GetOptimizedQuadRulesWithNodesOnLevelSet(CellMask mask, int order, MultidimensionalArray integrationResults) {
            var result = new List<ChunkRulePair<QuadRule>>(mask.NoOfItemsLocally);
            int noOfPhis = NumberOfMoments;

            foreach (Chunk chunk in mask) {
                for (int jCell = chunk.i0; jCell < chunk.JE; jCell++) {

                    Chunk singleElement = new Chunk() {
                        i0 = jCell,
                        Len = 1
                    };

                    int iteration = 1;
                    int maxNodeIncrementIterations = 8;
                    NodeSet projectedNodes;

                    do {
                        NodeSet baseNodes = GetSeedNodes(iteration * noOfPhis);
                        projectedNodes = ProjectNodesOntoLevelSet(jCell, baseNodes);
                        iteration++;
                        if (iteration > maxNodeIncrementIterations) {
                            Console.WriteLine("WARNING (surface HMF): Max iterations in SurfaceRuleFactory for cell {0} reached, only {1} of recommend {2} nodes",
                                jCell,
                                projectedNodes.GetLength(0),
                                nodeCountSafetyFactor * noOfPhis);
                            break;
                        }
                    } while (projectedNodes.GetLength(0) < nodeCountSafetyFactor * noOfPhis);

                    QuadRule optimizedRule = GetOptimizedRule(
                        projectedNodes,
                        order,
                        jCell,
                        integrationResults,
                        levelSetIndex);
                    
                    result.Add(new ChunkRulePair<QuadRule>(singleElement, optimizedRule));
                }
            }

            return result;
        }

        private AffineTrafo[] ConstructTransformationsToBoundingBox(CellMask mask, int order) {
            AffineTrafo[] trafos = new AffineTrafo[mask.NoOfItemsLocally];
            int D = RefElement.SpatialDimension;

            foreach (Chunk chunk in mask) {
                foreach (var cell in chunk.Elements.AsSmartEnumerable()) {
                    CellMask singleElementMask = new CellMask(
                        LevelSetData.GridDat, new[] { Chunk.GetSingleElementChunk(cell.Value) }, MaskType.Geometrical);

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
                    if (noOfRoots == 1) {
                        trafos[localCellIndex2SubgridIndex[cell.Value]] = null;
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

                        double scaling = Math.Sqrt(LevelSetData.GridDat.Cells.JacobiDet[cell.Value]);

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

            return trafos;
        }

        /// <summary>
        /// Used for the integration over the divergence-free basis functions
        /// (cf. <see cref="subGrid"/>) over the boundary of the current element.
        /// </summary>
        private class PhiQuadrature : CellBoundaryQuadrature<CellBoundaryQuadRule> {

            /// <summary>
            /// Owner of this object.
            /// </summary>
            private LevelSetSurfaceQuadRuleFactory owner;

            /// <summary>
            /// Number of cells in the current sub-grid
            /// </summary>
            private int NoOfItemsLocally;

            /// <summary>
            /// Constructor
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="mask"></param>
            public PhiQuadrature(LevelSetSurfaceQuadRuleFactory owner, CellMask mask)
                : base(
                    new int[] { owner.NumberOfMoments },
                     owner.LevelSetData.GridDat,
                     (new CellBoundaryQuadratureScheme(owner.edgeRuleFactory, mask)).Compile(owner.LevelSetData.GridDat, GetQuadratureDegree(owner)),
                     CoordinateSystem.Reference) {
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
            ///     2nd index: Local edge of the reference element
            ///     </item>
            ///     <item>
            ///     3rd index: Basis function index
            ///     </item>
            /// </list>
            /// </summary>
            public MultidimensionalArray IntegrationResults {
                get;
                private set;
            }

            /// <summary>
            /// Performs the integration while overwriting old results in
            /// <see cref="IntegrationResults"/>
            /// </summary>
            public override void Execute() {
                IntegrationResults = MultidimensionalArray.Create(
                    NoOfItemsLocally, owner.RefElement.NoOfFaces, IntegralCompDim[0]);
                base.Execute();
            }



            /// <summary>
            /// For each cell \f$ K\f$  and for each
            /// divergence-free basis polynomial
            /// \f$ \vec{\Phi}\f$ : Computes the
            /// integral
            /// \f$ 
            ///     \int \limits_{\partial K} \vec{\Phi} H(\varphi) \;dS,
            /// \f$ 
            /// where \f$ \varphi\f$  is the level set
            /// function and \f$ H\f$  is an indicator
            /// function that restricts the integration domain to positive
            /// level set values (Heaviside) or negative level set values
            /// (One minus Heaviside). This choice is implicitly given by
            /// the edge rule factory supplied to
            /// <see cref="LevelSetSurfaceQuadRuleFactory.LevelSetSurfaceQuadRuleFactory"/>
            /// </summary>
            protected override void Evaluate(int i0, int Length, CellBoundaryQuadRule QR, MultidimensionalArray EvalResult) {
                NodeSet QuadNodes = QR.Nodes;
                int noOfEdges = owner.RefElement.NoOfFaces;
                int D = gridData.SpatialDimension;

                for (int i = 0; i < Length; i++) {
                    MultidimensionalArray phiValues = owner.EvaluatePhis(i0 + i, QuadNodes);
                    int nodeIndex = 0;

                    // loop over edges of the reference element
                    for (int e = 0; e < noOfEdges; e++) {
                        // loop over nodes per edge
                        for (int j = 0; j < CurrentRule.NumbersOfNodesPerFace[e]; j++) {
                            // loop over (vector) test functions
                            for (int k = 0; k < m_TotalNoOfIntegralsPerItem; k++) {
                                double acc = 0;
                                // sum over vector components
                                for (int d = 0; d < D; d++) {
                                    acc += phiValues[nodeIndex, k, d] * owner.RefElement.FaceNormals[e, d];
                                }
                                EvalResult[i, nodeIndex, k] += acc;
                            }

                            nodeIndex++;
                        }
                    }
                }
            }

            /// <summary>
            /// Saves the results of the integration to
            /// <see cref="IntegrationResults"/>
            /// </summary>
            /// <param name="i0"></param>
            /// <param name="Length"></param>
            /// <param name="ResultsOfIntegration"></param>
            protected override void SaveIntegrationResults(int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                for (int i = 0; i < Length; i++) {
                    int iSubGrid = owner.localCellIndex2SubgridIndex[i + i0];
                    for (int e = 0; e < owner.RefElement.NoOfFaces; e++) {
                        for (int j = 0; j < IntegralCompDim[0]; j++) {
                            IntegrationResults[iSubGrid, e, j] = ResultsOfIntegration[i, e, j];
                        }
                    }
                }
            }

            /// <summary>
            /// Determines the required degree of the boundary integration
            /// </summary>
            /// <param name="owner">
            /// The owner of this object.
            /// </param>
            /// <returns>
            /// <paramref name="owner"/>.phiBasis.Degree in 2D, but
            /// <paramref name="owner"/>.phiBasis.Degree + 1 in 3D.
            /// </returns>
            /// <remarks>
            /// In 2D, it is given by the degree p of the basis we want to
            /// integrate, and, as the boundary integral can be evaluated
            /// <b>exactly</b> guarantees a convergence order h^(p+1). In 3D,
            /// however, the situation is a bit different: If we use order p on
            /// the boundary, the related error converges with h^(p+1), but the
            /// surface case is based on the divergence theorem, i.e. the
            /// surface integral only converges with h^p. Since we want h^(p+1)
            /// in both cases, we thus use p+1 for the boundary integral.
            /// </remarks>
            private static int GetQuadratureDegree(LevelSetSurfaceQuadRuleFactory owner) {
                switch (owner.RefElement.SpatialDimension) {
                    case 2:
                        return owner.phiBasis.MaxAbsoluteDegree;

                    case 3:
                        return owner.phiBasis.MaxAbsoluteDegree + 1;

                    default:
                        throw new ArgumentException("Invalid spatial dimension", "owner");
                }
            }
        }
    }
}
