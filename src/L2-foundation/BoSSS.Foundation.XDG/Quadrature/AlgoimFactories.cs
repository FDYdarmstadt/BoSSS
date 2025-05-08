// Ignore Spelling: Algoim

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.Linq;


namespace BoSSS.Foundation.XDG.Quadrature.Algoim {

	/// <summary>
	/// Provides a factory configuration for processing based on specified parameters.
	/// </summary>
	public class AlgoimFactories {

		/// <summary>
		/// This would return a factory object with the configuration 
		/// </summary>
		/// <param name="ls">level set data</param>
		/// <param name="e">reference element</param>
		/// <param name="negativeLevelSet">are the negative level sets considered (default positive, i.e., ls > 0) </param>
		/// <param name="callSurfaceAndVolumeAtOnce">if enabled, performs calculation of volume and surface rules at once to reduce computational cost. 
		/// (it can cause redundancy if both are not required at the same time)</param>
		public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e, bool negativeLevelSet = false, bool callSurfaceAndVolumeAtOnce = true) {
			if (!(e is Grid.RefElements.Cube ||
	              e is Grid.RefElements.Square ||
	              e is Grid.RefElements.Line ||
	              e is Grid.RefElements.Point))
				throw new NotSupportedException("Algoim is only supported for Cube, Square, Line and Point reference elements");

			refElement = e;
            lsData = ls;
            VolumeSign = negativeLevelSet;
            SurfaceAndVolumeAtOnce = callSurfaceAndVolumeAtOnce;
            MaxGrid = lsData.GridDat.Cells.GetCells4Refelement(refElement).Intersect(lsData.Region.GetCutCellMask4LevSet(lsData.LevelSetIndex).ToGeometicalMask());
        }

        LevelSetTracker.LevelSetData lsData;

        RefElement refElement;

        CellMask MaxGrid;

        /// <summary>
        /// Enables combined surface and volume quadrature rule calls with caching to avoid repeated polynomial interpolations. This speeds up computations when both rules are needed but may cause unnecessary calculations if only one rule is required.
        /// </summary>
        bool SurfaceAndVolumeAtOnce;

        /// <summary>
        /// A boolean to change the sign of level set. true: negative level set, false: positive level set (ls > 0).
        /// </summary>
        bool VolumeSign;

        /// <summary>
        /// - key: quadrature order 
        /// - value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_SurfaceRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        /// <summary>
        /// - key: quadrature order 
        /// - value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_CellBoundaryRules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();

        /// <summary>
        /// - key: quadrature order 
        /// - value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_EdgeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        /// <summary>
        /// - key: quadrature order 
        /// - value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_VolumeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            var factory = new AlgoimFactory() {
                m_Owner = this,
                m_Rules = this.m_SurfaceRules
            };
            factory.useMetrics = true;
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetSurfaceQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            var factory = new AlgoimFactory() {
                m_Owner = this,
                m_Rules = this.m_VolumeRules
            };
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetVolumeQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundarySurfaceFactory() {
            var factory = new AlgoimCellBoundaryFactory() {
                m_Owner = this,
                m_Rules = this.m_CellBoundaryRules
            };
            //factory.useMetrics = true;
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetSurfaceQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundaryVolumeFactory() {
            var factory = new AlgoimCellBoundaryFactory() {
                m_Owner = this,
                m_Rules = this.m_CellBoundaryRules
            };
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetVolumeQuadratureRules;
            return factory;
        }

        public IQuadRuleFactory<QuadRule> GetEdgeVolumeFactoryOnEdge() {
            var factory = new AlgoimEdgeRuleFactoryOnEdge() {
                m_Owner = this,
                m_Rules = this.m_EdgeRules
            };
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetVolumeQuadratureRules;
            return factory;
        }

        class AlgoimFactory : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<QuadRule>[]> m_Rules;

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            bool VolumeSign => m_Owner.VolumeSign;

            CellMask m_MaxGridMask => m_Owner.MaxGrid;

            bool m_SurfaceAndVolumeAtOnce => spaceDim > 1 ? m_Owner.SurfaceAndVolumeAtOnce : false;

            public bool useMetrics = false;

            // only valid for single quadrature rule calculations (either surface or volume)
            internal GetQuadratureRule m_CalculateQuadRule;
            internal delegate UnsafeAlgoim.QuadScheme GetQuadratureRule(int dim, int p, int q, int[] lengths, double[] x, double[] y);

            public LevelSetTracker.LevelSetData lsData => m_Owner.lsData;

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            /// <summary>
            /// Returns the quadrature rule of given order for each cell in mask.
            /// </summary>
            /// <param name="mask">Cell mask for which quadrature nodes are requested</param>
            /// <param name="RequestedOrder">Order of quadrature</param>
            /// <returns>An array of quadrature rules for cells in mask</returns>
            /// <exception cref="ArgumentException"></exception>
            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");

                if (!m_Rules.ContainsKey(RequestedOrder)) {
                    if (m_SurfaceAndVolumeAtOnce)
                        CalculateQuadRuleSetCombo(RequestedOrder);
                    else
                        CalculateQuadRuleSetSingle(RequestedOrder);
                }

                // check if all the mask or a submask is requested
                if (mask.NoOfItemsLocally == m_MaxGridMask.NoOfItemsLocally) {
                    return m_Rules[RequestedOrder];
                } else {
                    var Rule = m_Rules[RequestedOrder];
                    int localLength = mask.NoOfItemsLocally, totalLength = Rule.Length;
                    var Ret = new ChunkRulePair<QuadRule>[localLength];
                    int t = 0;
                    int jsub = 0;
                    foreach (int jCell in mask.ItemEnum) {
                        Debug.Assert(Rule[t].Chunk.Len == 1);
                        while (jCell > Rule[t].Chunk.i0) {
                            t++;
                        }

                        Debug.Assert(jCell == Rule[t].Chunk.i0);
                        Ret[jsub] = Rule[t];
#if DEBUG
                        Ret[jsub].Rule.Weights.CheckForNanOrInf();
                        Ret[jsub].Rule.Nodes.CheckForNanOrInf();
#endif
                        jsub++;
                    }
                    Debug.Assert(jsub == localLength);

                    return Ret;
                }
            }

            private void CalculateQuadRuleSetSingle(int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (int cell in m_MaxGridMask.ItemEnum) {
                    var quadRule = GetNodesAndWeights(cell, RequestedOrder);
                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));
                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
            }

            /// <summary>
            /// Combo rule calculation for surface and volume rules at the same time. Returns the values of either but stores both in cache once calculated
            /// </summary>
            /// <param name="RequestedOrder"></param>
            /// <returns></returns>
            private void CalculateQuadRuleSetCombo(int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> retSurf = new List<ChunkRulePair<QuadRule>>();
                List<ChunkRulePair<QuadRule>> retVol = new List<ChunkRulePair<QuadRule>>();

                foreach (int cell in m_MaxGridMask.ItemEnum) {
                    var quadRule = GetNodesAndWeightsCombo(cell, RequestedOrder);
                    retSurf.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule[0]));
                    retVol.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule[1]));
                }

                m_Owner.m_SurfaceRules.Add(RequestedOrder, retSurf.ToArray());
                m_Owner.m_VolumeRules.Add(RequestedOrder, retVol.ToArray());
            }

            /// <summary>
            /// Get the quadrature nodes and weights for a cell. (either surface or volume, declared at the instantiation of the factory)
            /// </summary>
            /// <param name="jCell">local index of the cell</param>
            /// <param name="RequestedOrder">requested order for quadrature</param>
            /// <returns></returns>
            private QuadRule GetNodesAndWeights(int jCell, int RequestedOrder) {
                //number of nodes in 1d for level set interpolation = degree + 1
                int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

                //create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
                double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

                //Cartesian pair product for the points
                int numberOfCombinations = (int)Math.Pow(points.Length, spaceDim); 
                MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(points, spaceDim);

                //Get level set values on combinations
                NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, jCell, 1);

                //Convert to 1D Array for the Wrapper
                double[] y = new double[ret.Length];

                if (VolumeSign) {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = -ret[0, i];
                } else {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = ret[0, i];
                }

                // create the double array for the coordinates that level set are queried (1D version, for details see AlgoimWrapper)
                double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

                // create the double array for the number of points in each coordinate
                int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

                // get the quadrature rule from the wrapper
                UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(spaceDim, n, RequestedOrder, sizes, x, y);
 
                // If quadrature rule is empty, return.
                if (qs.length < 1) 
                    return CreateEmptyQuadRule();
                
                // Create quadrature rule and copy from the scheme
                QuadRule quadRule = QuadRule.CreateBlank(RefElement, qs.length, qs.dimension);

                for (int row = 0; row < qs.length; row++) {
                    quadRule.Weights[row] = qs.weights[row];
                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension +d;
                        quadRule.Nodes[row, d] = qs.nodes[ind];
                    }
                }
                quadRule.Nodes.LockForever();
                quadRule.OrderOfPrecision = RequestedOrder;

                //// In order to calculate surface and volume needs in the same loop, there is a "hack", in which surface nodes multiplied their transformation coefficients but divided by determinants for volume.
                //if (useMetrics)
                //    ApplyMetrics(quadRule, qs, jCell);

                return quadRule;
            }

            /// <summary>
            /// GetNodesAndWeightsCombo (both for surface and volume at the same time)
            /// </summary>
            /// <param name="jCell">local index of the cell</param>
            /// <param name="RequestedOrder">requested order for quadrature</param>
            /// <returns></returns>
            private QuadRule[] GetNodesAndWeightsCombo(int jCell, int RequestedOrder) {
                //number of nodes in 1d for level set interpolation = degree + 1
                int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

                //create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
                double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

                //Cartesian pair product for the points
                int numberOfCombinations = (int)Math.Pow(points.Length, spaceDim);
                MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(points, spaceDim);

                //Get level set values on combinations
                NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, jCell, 1);

                //Convert to 1D Array for the Wrapper
                double[] y = new double[ret.Length];

                if (VolumeSign) {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = -ret[0, i];
                } else {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = ret[0, i];
                }

                // create the double array for the coordinates that level set are queried (1D version, for details see AlgoimWrapper)
                double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

                // create the double array for the number of points in each coordinate
                int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

                // get the quadrature rule from the wrapper
                UnsafeAlgoim.QuadSchemeCombo qs = ilPSP.Utils.Algoim.GetComboQuadratureRules(spaceDim, n, RequestedOrder, sizes, x, y);
                
                QuadRule[] quadRules = new QuadRule[2];

                // If the surface quadrature rule is empty.
                if (qs.lengthVol + qs.lengthSurf < 1) {
                    quadRules[0] = CreateEmptyQuadRule(); quadRules[0].OrderOfPrecision = RequestedOrder;
                    quadRules[1] = CreateEmptyQuadRule(); quadRules[1].OrderOfPrecision = RequestedOrder;
                    return quadRules;
                }

                // Create surface quadrature rule and copy from the scheme
                QuadRule quadRuleSurf = QuadRule.CreateBlank(RefElement, qs.lengthSurf, qs.dimension);
                for(int row = 0; row < qs.lengthSurf; row++) {
                    quadRuleSurf.Weights[row] = qs.weights[row];
                    for(int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension + d;
                        quadRuleSurf.Nodes[row, d] = qs.nodes[ind];
                    }
                }
                quadRuleSurf.Nodes.LockForever();
                quadRules[0] = quadRuleSurf;
                quadRules[0].OrderOfPrecision = RequestedOrder;

                // Create volume quadrature rule and copy from the scheme
                QuadRule quadRuleVol = QuadRule.CreateBlank(RefElement, qs.lengthVol, qs.dimension);
                for (int row = qs.lengthSurf, rowVol=0; rowVol < qs.lengthVol; rowVol++,row++) { //lengthSurf+lengthVol = total length
                    quadRuleVol.Weights[rowVol] = qs.weights[row];
                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension + d;
                        quadRuleVol.Nodes[rowVol, d] = qs.nodes[ind];
                    }
                }
                quadRuleVol.Nodes.LockForever();
                quadRules[1] = quadRuleVol;
                quadRules[1].OrderOfPrecision = RequestedOrder;

                //// Apply metrics if required (applied only to surface rule)
                //ApplyMetrics(quadRuleSurf, qs, jCell);

                return quadRules;
            }

            
            //// In order to calculate surface and volume needs in the same loop, there is a "hack", in which surface nodes multiplied their transformation coefficients but divided by determinants for volume.
            //private void ApplyMetrics(QuadRule quadRule, UnsafeAlgoim.QuadScheme qs, int jCell) {
            //    var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);
            //    for (int k = 0; k < qs.length; k++) {
            //        quadRule.Weights[k] /= metrics[0, k];
            //    }
            //}

            //private void ApplyMetrics(QuadRule quadRule, UnsafeAlgoim.QuadSchemeCombo qs, int jCell) {
            //    var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);
            //    for (int k = 0; k < qs.lengthSurf; k++) {
            //        quadRule.Weights[k] /= metrics[0, k];
            //    }
            //}

            private QuadRule CreateEmptyQuadRule() {
                QuadRule quadRuleEmpty = QuadRule.CreateBlank(RefElement, 1, spaceDim);
                quadRuleEmpty.Nodes.LockForever();
                return quadRuleEmpty;
            }

        }
        
        class AlgoimCellBoundaryFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_Rules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            // Edge rule is D-1 for the original space D
            int ruleDim => spaceDim - 1;

            bool m_SurfaceAndVolumeAtOnce => false;

            bool VolumeSign => m_Owner.VolumeSign;

            //public bool useMetrics = false;

            internal GetQuadratureRule m_CalculateQuadRule;

            internal delegate UnsafeAlgoim.QuadScheme GetQuadratureRule(int dim, int p, int q, int[] lengths, double[] x, double[] y);

            public LevelSetTracker.LevelSetData lsData => m_Owner.lsData;

            public GridData GridDat => lsData.GridDat;

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask cellMask))
                    throw new ArgumentException("Expecting a cell mask.");

				if (m_Rules.ContainsKey(RequestedOrder))
                    return m_Rules[RequestedOrder];

                if (m_SurfaceAndVolumeAtOnce)
                    return CalculateQuadRuleSetCombo(cellMask, RequestedOrder);
                else
                    return CalculateQuadRuleSetSingle(cellMask, RequestedOrder);
            }

            /// <summary>
            /// Combines different edge rules into a single CellBoundaryQuadRule rule.
            /// </summary>
            /// <param name="rules">An array of edge rules to be combined.</param>
            /// <returns>A combined <c>CellBoundaryQuadRule</c> with arranged <c>NumbersOfNodesPerFace</c>.</returns>
            CellBoundaryQuadRule CombineEdgeRules(QuadRule[] rules) {
                int numberOfNodes = 0;
                foreach (QuadRule rule in rules) {
                    numberOfNodes += rule.NoOfNodes;
                }
                CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(RefElement, numberOfNodes, spaceDim, RefElement.NoOfFaces);
                int subarrayPointer = 0;
                for (int i = 0; i < rules.Length; ++i) {
                    int subNumberOfNodes = rules[i].NoOfNodes;
                    combinedRule.NumbersOfNodesPerFace[i] = subNumberOfNodes;
                    for (int j = 0; j < subNumberOfNodes; ++j) {
                        for (int d = 0; d < spaceDim; ++d) {
                            combinedRule.Nodes[subarrayPointer + j, d] = rules[i].Nodes[j, d];
                        }
                        combinedRule.Weights[subarrayPointer + j] = rules[i].Weights[j];
                    }
                    subarrayPointer += subNumberOfNodes;
                }
                combinedRule.Nodes.LockForever();
                combinedRule.OrderOfPrecision = rules.Select(r => r.OrderOfPrecision).Min();
                return combinedRule;
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetSingle(CellMask mask, int RequestedOrder) {
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

                List<ChunkRulePair<CellBoundaryQuadRule>> ret = new List<ChunkRulePair<CellBoundaryQuadRule>>();

                int noOfFaces = RefElement.NoOfFaces;

                foreach (int cell in mask.ItemEnum) {
                    QuadRule[] edgeRules = new QuadRule[noOfFaces];
                    for (int e = 0; e < noOfFaces; e++) { 
                        edgeRules[e] = GetNodesAndWeights(cell, RequestedOrder, e);
                    }

                    var combinedEdgeRules = CombineEdgeRules(edgeRules);
                    ret.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(cell), combinedEdgeRules));

                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetCombo(CellMask mask, int requestedOrder) {
                throw new NotImplementedException();
            }

            //Map from bosss-cube indices to edges
            //RowNumber: index of edge in bosss-cube 
            //Each Row: first entry: index of axis, second entry: coordinate in the axis
            static int[,] edge2CubeMap = new int[6, 2]{
            { 0, -1 },
            { 0,  1 },
            { 1,  1 },
            { 1, -1},
            { 2, 1 },
            { 2, -1 }};

            private CellBoundaryQuadRule GetNodesAndWeights(int jCell, int RequestedOrder, int faceIndex) {

                //var scalings =  this.GridDat.Edges.SqrtGramian;
                //int iEdge = this.GridDat.GetEdgesForFace(jCell, faceIndex, out _, out var furtherEdges);
                //if(furtherEdges != null && furtherEdges.Length > 0)
                //    throw new NotImplementedException(); 


                //int[] Cells2Edge = .Cells.Cells2Edges[jCell];
                //int iEdge = this.GridDat.GetEdgesForFace(int jCell, int iFace, out int InOrOut, out int[] MoreEdges)

                //number of nodes in 1d for level set interpolation = degree + 1
                int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

                //create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
                double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

                int numberOfCombinations = points.Length.Pow(ruleDim); //Cartesian pair product for the points
                MultidimensionalArray combinationsOnEdge = MultidimensionalArray.CreateCartesianPairProduct(points, ruleDim);

                // Transform edge based coordinates to cell based coordinates (to query level sets)
                MultidimensionalArray combinations  = MultidimensionalArray.Create(numberOfCombinations, spaceDim);  // to transform edge coordinates to cell, create an array in the original space dim
				RefElement.TransformFaceCoordinates(faceIndex, combinationsOnEdge, combinations);

				NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, jCell, 1);

                // create data for Algoim wrapper and query the quad rule
                double[] y = new double[ret.Length];
                if (VolumeSign) {
					for (int i = 0; i < ret.Length; i++)
                        y[i] = -ret[0, i];
                } else {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = ret[0, i];
                }

				// notice that on Algoim side, this is still an edge rule (e.g., a plane on ruleDim=spaceDim-1 for a cube)
				double[] x = Enumerable.Repeat(points, ruleDim).SelectMany(i => i).ToArray(); 
                int[] sizes = Enumerable.Repeat(points.Length, ruleDim).ToArray();

                // calculate the quadrature rule
                UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(ruleDim, n, RequestedOrder, sizes, x, y);

                if (qs.length < 1) {
                    CellBoundaryQuadRule quadRuleEmpty = CellBoundaryQuadRule.CreateEmpty(RefElement, 1, spaceDim, RefElement.NoOfFaces);
                    quadRuleEmpty.Nodes.LockForever();
                    quadRuleEmpty.OrderOfPrecision = RequestedOrder;
                    return quadRuleEmpty;
                }

				QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement.FaceRefElement, qs.length, qs.dimension);

                for (int row = 0; row < qs.length; row++) {
                    quadRuleOnEdge.Weights[row] = qs.weights[row]; // * (1.0 / scalings[])

                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension + d;
                        quadRuleOnEdge.Nodes[row, d] = qs.nodes[ind];
                    }
                }
                
				// Algoim returns an edge based rule, it must be converted to a CellBoundaryQuadRule (cell based coordinates)
				CellBoundaryQuadRule quadRule = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, spaceDim, RefElement.NoOfFaces);
                quadRule.Weights = quadRuleOnEdge.Weights;
                RefElement.TransformFaceCoordinates(faceIndex, quadRuleOnEdge.Nodes, quadRule.Nodes); // to transform edge rule back to cell coordinates (since we are creating CellBoundaryQuadRule, cell based rule)
				quadRule.Nodes.LockForever();
                quadRule.OrderOfPrecision = RequestedOrder;

                //if(useMetrics) {
                //    var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);
                //    for(int k = 0; k < qs.length; k++)
                //        quadRule.Weights[k] /= metrics[0, k];
                //}

                return quadRule;
            }

        }

		#region Edge rules on edge (Obsolete)
		// Obsolete, one needs to arrange how to couple quadrature nodes for an edge, since an edge is shared by two cells and can lead to different nodes and values
		class AlgoimEdgeRuleFactoryOnEdge : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<QuadRule>[]> m_Rules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            // Edge rule is D-1 for the original space D
            int ruleDim => spaceDim - 1;

            bool m_SurfaceAndVolumeAtOnce => false;

            bool VolumeSign => m_Owner.VolumeSign;

            public bool useMetrics = false;

            internal GetQuadratureRule m_CalculateQuadRule;

            internal delegate UnsafeAlgoim.QuadScheme GetQuadratureRule(int dim, int p, int q, int[] lengths, double[] x, double[] y);

            public LevelSetTracker.LevelSetData lsData => m_Owner.lsData;

            public GridData GridDat => lsData.GridDat;

            private int GetEdges(RefElement refElement) {
                // helper vars
                int D = this.GridDat.SpatialDimension;
                var CellToEdge = this.GridDat.iLogicalCells.Cells2Edges;
                var EdgeToCell = this.GridDat.iLogicalEdges.CellIndices;
                var FaceIndices = this.GridDat.iGeomEdges.FaceIndices;
                var TrafoIdx = this.GridDat.iGeomEdges.Edge2CellTrafoIndex;
                var EdgeToCellTrafos = this.GridDat.iGeomEdges.Edge2CellTrafos;
                int NoOfFaces = RefElement.NoOfFaces;
                var Scalings = this.GridDat.iGeomEdges.Edge2CellTrafos_SqrtGramian;
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                var faceRefElement = refElement.FaceRefElement;
                return 0;
            }

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is EdgeMask edgeMask))
                    throw new ArgumentException("Expecting an edge mask.");

                if (m_Rules.ContainsKey(RequestedOrder))
                    return m_Rules[RequestedOrder];

                if (m_SurfaceAndVolumeAtOnce)
                    return CalculateQuadRuleSetCombo(edgeMask, RequestedOrder);
                else
                    return CalculateQuadRuleSetSingle(edgeMask, RequestedOrder);
            }

            public IEnumerable<IChunkRulePair<QuadRule>> CalculateQuadRuleSetSingle(EdgeMask mask, int RequestedOrder) {
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (int edge in mask.ItemEnum) {

                    QuadRule edgeRule = GetNodesAndWeightsOnEdge(edge, RequestedOrder);

                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(edge), edgeRule));

                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            // to do
            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetCombo(ExecutionMask mask, int requestedOrder) {
                throw new NotImplementedException();
            }

            private QuadRule GetNodesAndWeightsOnEdge(int edge, int requestedOrder) {
                if (edge >= 0)
                    throw new NotImplementedException();

				LevelSet levelSet = (LevelSet)lsData.LevelSet;

				//number of nodes in 1d for level set interpolation = degree + 1
				int n = levelSet.Basis.Degree + 1;

				//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
				double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

				int numberOfCombinations = (int)Math.Pow(points.Length, ruleDim);        //Cartesian pair product for the points
				MultidimensionalArray combinationsOnEdge = MultidimensionalArray.CreateCartesianPairProduct(points, ruleDim);

				NodeSet NS = new NodeSet(RefElement.FaceRefElement, combinationsOnEdge, false);

				//create a new array but the first index is reserved for node index, which would be 0 always in this scope
				int[] lengthsModified = new int[NS.Lengths.Length+1];
                lengthsModified[0] = 1; //we have only one node
				Array.Copy(NS.Lengths, 0, lengthsModified, 1, NS.Lengths.Length);   // copy for the rest excluding the first index   

				MultidimensionalArray LSinValues = MultidimensionalArray.Create(lengthsModified);
				MultidimensionalArray LSoutValues = MultidimensionalArray.Create(lengthsModified);

				levelSet.EvaluateEdge(edge, 1, NS, LSinValues, LSoutValues); //calculate the level set values on edge, notice that it is designed for multiple nodes

                var ret = LSinValues;
				double[] y = new double[ret.Length];

				if (VolumeSign) {
					for (int i = 0; i < ret.Length; i++)
						y[i] = -ret[0, i];
				} else {
					for (int i = 0; i < ret.Length; i++)
						y[i] = ret[0, i];
				}

				double[] x = Enumerable.Repeat(points, ruleDim).SelectMany(i => i).ToArray();

				int[] sizes = Enumerable.Repeat(points.Length, ruleDim).ToArray();

				UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(ruleDim, n, requestedOrder, sizes, x, y);



				QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement, qs.length, spaceDim);


				quadRuleOnEdge.Nodes.LockForever();

                //if (useMetrics) {
                //	var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRuleOnEdge.Nodes, jCell, 1);
                //	//useMetrics.SaveToTextFile("useMetrics.txt");

                //	////quadRule.Weights.To1DArray().ToList().SaveToTextFileDebugUnsteady("beforeWeights");

                //	for (int k = 0; k < qs.length; k++)
                //		quadRuleOnEdge.Weights[k] /= metrics[0, k];
                //}


                //quadRuleGlobal.OutputQuadratureRuleAsVtpXML($"bosssj{jCell}e{edgeIndex}.vtp");
                return quadRuleOnEdge;




                //var Edg2Cel = this.GridDat.iGeomEdges.CellIndices;
                //var Edg2Fac = this.GridDat.iGeomEdges.FaceIndices;
                //int J = this.GridDat.Cells.NoOfLocalUpdatedCells;
                //QuadRule DefaultRule = this.RefElement.GetQuadratureRule(requestedOrder); ;

                //int myIKrfeEdge = this.GridDat.Edges.EdgeRefElements.IndexOf(this.RefElement, (a, b) => object.ReferenceEquals(a, b));
                //if (myIKrfeEdge < 0)
                //    throw new ApplicationException("fatal error");


                //int[] EdgeIndices = mask.ItemEnum.ToArray();
                //int NoEdg = EdgeIndices.Length;

                //// return value
                //ChunkRulePair<QuadRule>[] Ret = new ChunkRulePair<QuadRule>[NoEdg];

                //                // find cells
                //                // ==========

                //                BitArray CellBitMask = new BitArray(J);
                //                int[] Cells = new int[NoEdg]; // mapping: Edge Index --> Cell Index (both geometrical)
                //                int[] Faces = new int[NoEdg]; // mapping: Edge Index --> Face 
                //                BitArray MaxDomainMask = m_maxDomain.GetBitMask();
                //                {
                //                    for (int i = 0; i < NoEdg; i++) {

                //                        int iEdge = EdgeIndices[i];
                //                        if (this.grd.Edges.GetRefElementIndex(iEdge) != myIKrfeEdge)
                //                            throw new ArgumentException("illegal edge mask");

                //                        if (!(grd.Edges.IsEdgeConformalWithCell1(iEdge) || grd.Edges.IsEdgeConformalWithCell2(iEdge))) {
                //                            throw new NotSupportedException("For an edge that is not conformal with at least one cell, no edge rule can be created from a cell boundary rule.");
                //                        }

                //                        int jCell0 = Edg2Cel[iEdge, 0];
                //                        int jCell1 = Edg2Cel[iEdge, 1];
                //                        bool conf0 = grd.Edges.IsEdgeConformalWithCell1(iEdge);
                //                        bool conf1 = grd.Edges.IsEdgeConformalWithCell2(iEdge);

                //                        // this gives no errors for surface elements in 3D
                //                        bool Allow0 = MaxDomainMask[jCell0];
                //                        bool Allow1 = (jCell1 >= 0 && jCell1 < J) ? MaxDomainMask[jCell1] : false;

                //                        // //this is required for MPI parallel calculations
                //                        //bool Allow0 = true;// AllowedCells[jCell0];
                //                        //bool Allow1 = (jCell1 >= 0 && jCell1 < J);// ? AllowedCells[jCell1] : false;

                //                        if (!Allow0 && !Allow1) {
                //                            // fallback onto default rule, if allowed

                //                            //if (this.m_DefaultRuleFallbackAllowed) {
                //                            //    Cells[i] = -1; // by a negative index, we mark that we take the default rule
                //                            //    Faces[i] = -1;
                //                            //    Ret[i] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(EdgeIndices[i]), DefaultRule);


                //                            //} else {
                //                            throw new ArgumentException("unable to find a cell from which the edge rule can be taken.");
                //                            //}
                //                        } else {
                //                            Debug.Assert(Allow0 || Allow1);

                //                            if (conf0 && Allow0) {
                //                                // cell 0 is allowed and conformal:
                //                                // take this, it won't get better

                //                                CellBitMask[jCell0] = true;
                //                                Faces[i] = Edg2Fac[iEdge, 0] + 333; // arbitrary shift, as we need the sign to distinguish a case later on, otherwise Face = 0 is undefined and may lead to wrong results
                //                                Cells[i] = jCell0;
                //                            } else if (conf1 && Allow1) {
                //                                // cell 1 is allowed and conformal:
                //                                // take this, it won't get better

                //                                CellBitMask[jCell1] = true;
                //                                Faces[i] = Edg2Fac[iEdge, 1] + 333;
                //                                Cells[i] = jCell1;

                //                            } else if (Allow0) {
                //                                // cell 0 is allowed, but NOT conformal

                //                                CellBitMask[jCell0] = true;
                //                                Faces[i] = -(Edg2Fac[iEdge, 0] + 333); // by a negative index, we mark a non-conformal edge
                //                                Cells[i] = jCell0;
                //                            } else if (Allow1) {
                //                                // cell 1 is allowed, but NOT conformal

                //                                CellBitMask[jCell1] = true;
                //                                Faces[i] = -(Edg2Fac[iEdge, 1] + 333); // by a negative index, we mark a non-conformal edge
                //                                Cells[i] = jCell1;
                //                            }

                //                        }
                //                    }
                //                }


                //                // get cell boundary rule
                //                // ======================
                //                var CellMask = new CellMask(this.grd, CellBitMask, MaskType.Geometrical);


                //                IChunkRulePair<CellBoundaryQuadRule>[] cellBndRule = this.m_cellBndQF.GetQuadRuleSet(CellMask, order).ToArray();
                //                int[] jCell2PairIdx = new int[J];
                //                for (int i = 0; i < cellBndRule.Length; i++) {
                //                    var chk = cellBndRule[i].Chunk; // cell chunk
                //                    for (int jCell = chk.JE - 1; jCell >= chk.i0; jCell--) {
                //                        jCell2PairIdx[jCell] = i + 555;
                //                    }
                //                }

                //                int[] iChunk = new int[NoEdg]; // which chunk for which edge?
                //                for (int i = 0; i < NoEdg; i++) {
                //                    if (Cells[i] >= 0) {
                //                        iChunk[i] = jCell2PairIdx[Cells[i]] - 555;
                //                    } else {
                //                        iChunk[i] = int.MinValue;
                //                    }
                //                }
                //                // build rule
                //                // ==========
                //                {
                //                    for (int i = 0; i < NoEdg; i++) { // loop over edges
                //                                                      //if (MaxDomainMask[Cells[i]] == false)
                //                                                      //    Debugger.Break();

                //                        //if (Cells[i] >= 0) {
                //                        var CellBndR = cellBndRule[iChunk[i]].Rule;
                //                        QuadRule qrEdge = null;

                //                        if (Faces[i] >= 0)
                //                            qrEdge = this.CombineQr(null, CellBndR, Faces[i] - 333);
                //                        else {
                //                            qrEdge = this.CombineQrNonConformal(null, CellBndR, -Faces[i] - 333, EdgeIndices[i], Cells[i]);
                //                        }

                //                        qrEdge.Nodes.LockForever();
                //                        qrEdge.Weights.LockForever();
                //                        Ret[i] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(EdgeIndices[i]), qrEdge);
                //                        //} else {
                //                        //    Debug.Assert(Ret[i] != null);
                //                        //}
                //                    }
                //                }

                //                // return
                //                // ======

                //#if DEBUG
                //                for (int i = 0; i < Ret.Length; i++) {
                //                    Chunk c = Ret[i].Chunk;
                //                    for (int e = c.i0; e < c.JE; e++) {
                //                        Debug.Assert(maskBitMask[e] == true);
                //                    }
                //                }
                //#endif

                //                return Ret;


            }

        }
		#endregion
	}





}


