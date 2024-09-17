// Ignore Spelling: Algoim

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
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
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineAndPointQuadratureFactory;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;
using static BoSSS.Foundation.XDG.Quadrature.SayeSquare;

namespace BoSSS.Foundation.XDG.Quadrature {


    public class AlgoimFactories {


        public IQuadRuleFactory<QuadRule> GetSurfaceFactory() {
            var factory = new Factory() {
                m_Owner = this,
                m_Rules = this.m_SurfaceRules
            };
            factory.useMetrics = true;
            factory.m_CalculateQuadRule = Algoim.GetSurfaceQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            var factory = new Factory() {
                m_Owner = this,
                m_Rules = this.m_SurfaceRules
            };
            factory.m_CalculateQuadRule = Algoim.GetVolumeQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundarySurfaceFactory() {
            var factory = new CellBoundaryFactory() {
                m_Owner = this,
                m_Rules = this.m_CellBoundaryRules
            };
            factory.useMetrics = true;
            factory.m_CalculateQuadRule = Algoim.GetSurfaceQuadratureRules; 
            return factory;
        }

        public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundaryVolumeFactory() {
            var factory = new CellBoundaryFactory() {
                m_Owner = this,
                m_Rules = this.m_CellBoundaryRules
            };
            factory.m_CalculateQuadRule = Algoim.GetVolumeQuadratureRules;
            return factory;
        }

        public IQuadRuleFactory<QuadRule> GetEdgeVolumeFactoryOnEdge() {
            var factory = new EdgeRuleFactoryOnEdge() {
                m_Owner = this,
                m_Rules = this.m_EdgeRules
            };
            factory.m_CalculateQuadRule = Algoim.GetVolumeQuadratureRules;
            return factory;
        }

        //This would return a factory object with the configuration 
        //input: level set, tolerance etc.
        public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e, bool negativeLevelSet = false, bool callSurfaceAndVolumeAtOnce = false) {
            refElement = e;
            lsData = ls;
            VolumeSign = negativeLevelSet;
            SurfaceAndVolumeAtOnce = callSurfaceAndVolumeAtOnce;
        }

        LevelSetTracker.LevelSetData lsData;

        RefElement refElement;

        /// <summary>
        /// Enables combined surface and volume quadrature rule calls with caching to avoid repeated polynomial interpolations. This speeds up computations when both rules are needed but may cause unnecessary calculations if only one rule is required.
        /// </summary>
        bool SurfaceAndVolumeAtOnce;

        /// <summary>
        /// A boolean to change the sign of level set. true: negative level set, false: positive level set (ls > 0).
        /// </summary>
        bool VolumeSign;

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_SurfaceRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_CellBoundaryRules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_EdgeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[]> m_VolumeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

        #region Edge rules

        class Factory : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<QuadRule>[]> m_Rules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            bool VolumeSign => m_Owner.VolumeSign;

            public bool useMetrics = false;

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

                if (m_Rules.ContainsKey(RequestedOrder))
                    return m_Rules[RequestedOrder];

                if (m_Owner.SurfaceAndVolumeAtOnce)
                    return CalculateQuadRuleSetCombo(mask, RequestedOrder);
                else
                    return CalculateQuadRuleSetSingle(mask, RequestedOrder);
            }

            public IEnumerable<IChunkRulePair<QuadRule>> CalculateQuadRuleSetSingle(ExecutionMask mask, int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (int cell in mask.ItemEnum) {
                    var quadRule = GetNodesAndWeights(cell, RequestedOrder);
                        ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));
                    
                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            // to do
            public IEnumerable<IChunkRulePair<QuadRule>> CalculateQuadRuleSetCombo(ExecutionMask mask, int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (int cell in mask.ItemEnum) {
                    var quadRule = GetNodesAndWeights(cell, RequestedOrder);
                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));
                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

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

                UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(spaceDim, n, RequestedOrder, sizes, x, y);
 
                // If quadrature rule is empty, return.
                if (qs.length < 1) {
                    QuadRule quadRuleEmpty = QuadRule.CreateEmpty(RefElement, 1, spaceDim);
                    quadRuleEmpty.Nodes.LockForever();
                    return quadRuleEmpty;
                }

                // Create quadrature rule and copy from the scheme
                QuadRule quadRule = QuadRule.CreateEmpty(RefElement, qs.length, qs.dimension);

                for (int row = 0; row < qs.length; row++) {
                    quadRule.Weights[row] = qs.weights[row];

                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension +d;
                        quadRule.Nodes[row, d] = qs.nodes[ind];
                    }
                }

                quadRule.Nodes.LockForever();

                // In order to calculate surface and volume needs in the same loop, there is a "hack", in which surface nodes multiplied their transformation coefficients but divided by determinants for volume.
                if (useMetrics) { 
                    var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);

                    for (int k = 0; k < qs.length; k++)
                        quadRule.Weights[k] /= metrics[0, k];
                }

                return quadRule;
            }

        }

        class CellBoundaryFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]> m_Rules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[]>();

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            // Edge rule is D-1 for the original space D
            int ruleDim => spaceDim - 1;

            bool VolumeSign => m_Owner.VolumeSign;

            public bool useMetrics = false;

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

                if (m_Owner.SurfaceAndVolumeAtOnce)
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
                return combinedRule;
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetSingle(CellMask mask, int RequestedOrder) {
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

                List<ChunkRulePair<CellBoundaryQuadRule>> ret = new List<ChunkRulePair<CellBoundaryQuadRule>>();

                int noOfEdges = RefElement.NoOfFaces;

                foreach (int cell in mask.ItemEnum) {
                    QuadRule[] edgeRules = new QuadRule[noOfEdges];
                    for (int e = 0; e < noOfEdges; e++) {
                        edgeRules[e] = GetNodesAndWeights(cell, RequestedOrder, e);
                    }

                    var combinedEdgeRules = CombineEdgeRules(edgeRules);
                    ret.Add(new ChunkRulePair<CellBoundaryQuadRule>(Chunk.GetSingleElementChunk(cell), combinedEdgeRules));

                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            // to do
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

            private CellBoundaryQuadRule GetNodesAndWeights(int jCell, int RequestedOrder, int edgeIndex) {
                //number of nodes in 1d for level set interpolation = degree + 1
                int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

                //create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
                double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

                int numberOfCombinations = (int)Math.Pow(points.Length, ruleDim); //Cartesian pair product for the points

                MultidimensionalArray combinationsOnEdge = MultidimensionalArray.CreateCartesianPairProduct(points, ruleDim);
                MultidimensionalArray combinations  = MultidimensionalArray.Create(numberOfCombinations, spaceDim);
    
                // This operation basically insert a column at 'edgeMap[0]' into matrix combinationsOnEdge with all values equal to 'edgeMap[1]'
                int dimOffset = 0;
                for (int d = 0; d < spaceDim; d++) {
                    if (d == edge2CubeMap[edgeIndex, 0]) {
                        for (int i = 0; i < numberOfCombinations; i++) {
                            //Assign the coordinate from mapping (e.g. if left edge: -1)
                            combinations[i, d] = edge2CubeMap[edgeIndex, 1];
                        }
                        dimOffset = 1;
                    } else {
                        for (int i = 0; i < numberOfCombinations; i++) {
                            //Assign from combinationsOnEdge
                            combinations[i, d] = combinationsOnEdge[i,d-dimOffset];
                        }
                    }
                }

                //combinations.SaveToTextFile($"combsj{jCell}e{edgeIndex}.txt");

                NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, jCell, 1);

                //if boundary take in value
                //in out cells for discontinous (take average if not same )
                //(lsData.LevelSet as LevelSet).EvaluateEdge()

                //Console.WriteLine(ret.ToString());
                //ret.SaveToTextFile("ret.txt");

                //Console.WriteLine("Ret dim: " + ret.Dimension);
                //foreach (var l in ret.Lengths)
                //    Console.WriteLine("Length in:  " + l);

                //Console.WriteLine("Total number of elements: " + ret.Length);

                double[] y = new double[ret.Length];

                if (VolumeSign) {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = -ret[0, i];
                } else {
                    for (int i = 0; i < ret.Length; i++)
                        y[i] = ret[0, i];
                }

                double[] x = Enumerable.Repeat(points, ruleDim).SelectMany(i => i).ToArray();

                //new double[points.Length + points.Length];
                //Array.Copy(points, 0, x, 0, points.Length);
                //Array.Copy(points, 0, x, points.Length, points.Length);
                //int[] sizes = new int[] { 3, 3 };
                int[] sizes = Enumerable.Repeat(points.Length, ruleDim).ToArray();

                UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(ruleDim, n, RequestedOrder, sizes, x, y);
                //qs.OutputQuadratureRuleAsVtpXML("bosssj" + jCell + ".vtp");
                //Console.WriteLine("qs.length = " + qs.length);

                if (qs.length < 1) {
                    CellBoundaryQuadRule quadRuleEmpty = CellBoundaryQuadRule.CreateEmpty(RefElement, 1, spaceDim, RefElement.NoOfFaces);
                    quadRuleEmpty.Nodes.LockForever();
                    return quadRuleEmpty;
                }


                CellBoundaryQuadRule quadRuleOnEdge = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, spaceDim, RefElement.NoOfFaces);

                //int ind = 0;
                for (int row = 0; row < qs.length; row++) {
                    quadRuleOnEdge.Weights[row] = qs.weights[row];

                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension + d;
                        quadRuleOnEdge.Nodes[row, d] = qs.nodes[ind];
                        //Console.WriteLine($"quadRule.Nodes[{row}, {d}] = qs.nodes[{ind}] = {qs.nodes[ind]} ");
                    }
                }

                CellBoundaryQuadRule quadRule = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, spaceDim, RefElement.NoOfFaces);
                quadRule.Weights = quadRuleOnEdge.Weights;
                dimOffset = 0;
                for (int d = 0; d < spaceDim; d++) {
                    if (d == edge2CubeMap[edgeIndex, 0]) {
                        for (int i = 0; i < qs.length; i++) {
                            //Assign the coordinate from mapping (e.g. if left edge: -1)
                            quadRule.Nodes[i, d] = edge2CubeMap[edgeIndex, 1];
                        }
                        dimOffset = 1;
                    } else {
                        for (int i = 0; i < qs.length; i++) {
                            //Assign from combinationsOnEdge
                            quadRule.Nodes[i, d] = quadRuleOnEdge.Nodes[i, d - dimOffset];
                        }
                    }
                }


                //quadRule.OutputQuadratureRuleAsVtpXML("quadRule" + jCell + ".vtp");
                quadRule.Nodes.LockForever();
                //var quadRuleGlobal = quadRule.CloneAs();
                //quadRuleGlobal.TransformLocal2Global(lsData.GridDat.Grid, jCell);

                //var LevelSetNormals = lsData.GetLevelSetReferenceNormals(quadRule.Nodes, jCell, 1).ExtractSubArrayShallow(0, -1, -1);
                //LevelSetNormals.SaveToTextFile("lsnormals.txt");
                if (useMetrics) {
                    var metrics = lsData.GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);
                    //useMetrics.SaveToTextFile("useMetrics.txt");

                    ////quadRule.Weights.To1DArray().ToList().SaveToTextFileDebugUnsteady("beforeWeights");

                    for (int k = 0; k < qs.length; k++)
                        quadRule.Weights[k] /= metrics[0, k];
                }
                //quadRule.Weights.MatMatMul(useMetrics);

                //quadRule.Weights.To1DArray().ToList().SaveToTextFileDebugUnsteady("afterWeights");

                //quadRuleGlobal.OutputQuadratureRuleAsVtpXML($"bosssj{jCell}e{edgeIndex}.vtp");
                return quadRule;
            }

        }

        class EdgeRuleFactoryOnEdge : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            internal Dictionary<int, ChunkRulePair<QuadRule>[]> m_Rules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;

            // Edge rule is D-1 for the original space D
            int ruleDim => spaceDim - 1;

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

                if (m_Owner.SurfaceAndVolumeAtOnce)
                    return CalculateQuadRuleSetCombo(edgeMask, RequestedOrder);
                else
                    return CalculateQuadRuleSetSingle(edgeMask, RequestedOrder);
            }

            public IEnumerable<IChunkRulePair<QuadRule>> CalculateQuadRuleSetSingle(EdgeMask mask, int RequestedOrder) {
                if (mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("Expecting a geometrical mask.");

                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                int noOfEdges = RefElement.NoOfFaces;

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

                QuadRule ret = new QuadRule();
                ret.OrderOfPrecision = int.MaxValue - 1;
                ret.Nodes = new NodeSet(this.RefElement, 1, Math.Max(1, ruleDim), false);
                ret.Weights = MultidimensionalArray.Create(1);  // this is an empty rule, since the weight is zero!
                                                                // (rules with zero nodes may cause problems at various places.)
                return ret;

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
    }






    #endregion




    //QuadRule RuleToRuleThemAll = QuadRule.CreateEmpty(RefElement, count, spatialDim);
    //RuleToRuleThemAll.Nodes = new NodeSet(RefElement, nodes, true);
    //RuleToRuleThemAll.Weights = weights;
    //        return RuleToRuleThemAll;





}


