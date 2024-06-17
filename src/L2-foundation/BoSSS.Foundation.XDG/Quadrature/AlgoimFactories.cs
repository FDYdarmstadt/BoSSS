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

        public IQuadRuleFactory<CellBoundaryQuadRule> GetEdgeVolumeFactory() {
            var factory = new EdgeRuleFactory() {
                m_Owner = this,
                m_Rules = this.m_CellBoundaryRules
            };
            factory.m_CalculateQuadRule = Algoim.GetVolumeQuadratureRules;
            return factory;
        }

        //This would return a factory object with the configuration 
        //input: level set, tolerance etc.
        public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e, bool negativeVolume = false, bool callSurfaceAndVolumeAtOnce = false) {
            refElement = e;
            lsData = ls;
            VolumeSign = negativeVolume;
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
                    QuadRule quadRuleEmpty = QuadRule.CreateEmpty(RefElement, 1, 1);
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

        class EdgeRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
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

            private int GetEdges(RefElement refElement) {
                // helper vars
                int D = this.GridDat.SpatialDimension;
                var CellToEdge = this.GridDat.iLogicalCells.Cells2Edges;
                var EdgeToCell = this.GridDat.iLogicalEdges.CellIndices;
                var FaceIndices = this.GridDat.iGeomEdges.FaceIndices;
                var TrafoIdx = this.GridDat.iGeomEdges.Edge2CellTrafoIndex;
                var EdgeToCellTrafos = this.GridDat.iGeomEdges.Edge2CellTrafos;
                int NoOfFaces = RefElement.NoOfFaces;
                var Scalings = GridDat.iGeomEdges.Edge2CellTrafos_SqrtGramian;
                int J = this.GridDat.iLogicalCells.NoOfLocalUpdatedCells;

                var faceRefElement = refElement.FaceRefElement;
                return 0;
            }

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");

                if (m_Rules.ContainsKey(RequestedOrder))
                    return m_Rules[RequestedOrder];

                if (m_Owner.SurfaceAndVolumeAtOnce)
                    return CalculateQuadRuleSetCombo(mask, RequestedOrder);
                else
                    return CalculateQuadRuleSetSingle(mask, RequestedOrder);
            }

            CellBoundaryQuadRule CombineEdgeRules(QuadRule[] rules) {
                int numberOfNodes = 0;
                foreach (QuadRule rule in rules) {
                    numberOfNodes += rule.NoOfNodes;
                }
                CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(RefElement, numberOfNodes, ruleDim, RefElement.NoOfFaces);
                int subarrayPointer = 0;
                for (int i = 0; i < rules.Length; ++i) {
                    int subNumberOfNodes = rules[i].NoOfNodes;
                    combinedRule.NumbersOfNodesPerFace[i] = subNumberOfNodes;
                    for (int j = 0; j < subNumberOfNodes; ++j) {
                        for (int d = 0; d < ruleDim; ++d) {
                            combinedRule.Nodes[subarrayPointer + j, d] = rules[i].Nodes[j, d];
                        }
                        combinedRule.Weights[subarrayPointer + j] = rules[i].Weights[j];
                    }
                    subarrayPointer += subNumberOfNodes;
                }
                combinedRule.Nodes.LockForever();
                return combinedRule;
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetSingle(ExecutionMask mask, int RequestedOrder) {
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
            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> CalculateQuadRuleSetCombo(ExecutionMask mask, int RequestedOrder) {
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

                combinations.SaveToTextFile("combs.txt");


                NodeSet NS = new NodeSet(RefElement.FaceRefElement, combinations, false);
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
                    CellBoundaryQuadRule quadRuleEmpty = CellBoundaryQuadRule.CreateEmpty(RefElement, 1, 1, RefElement.NoOfFaces);
                    quadRuleEmpty.Nodes.LockForever();
                    return quadRuleEmpty;
                }


                CellBoundaryQuadRule quadRule = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, qs.dimension, RefElement.NoOfFaces);

                //int ind = 0;
                for (int row = 0; row < qs.length; row++) {
                    quadRule.Weights[row] = qs.weights[row];

                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension + d;
                        quadRule.Nodes[row, d] = qs.nodes[ind];
                        //Console.WriteLine($"quadRule.Nodes[{row}, {d}] = qs.nodes[{ind}] = {qs.nodes[ind]} ");
                    }
                }
                //quadRule.OutputQuadratureRuleAsVtpXML("quadRule" + jCell + ".vtp");
                quadRule.Nodes.LockForever();
                var quadRuleGlobal = quadRule.CloneAs();
                quadRuleGlobal.TransformLocal2Global(lsData.GridDat.Grid, jCell);

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

                quadRuleGlobal.OutputQuadratureRuleAsVtpXML("bosssj" + jCell  + ".vtp");
                return quadRule;
            }

        }
    }






    #endregion




    //QuadRule RuleToRuleThemAll = QuadRule.CreateEmpty(RefElement, count, spatialDim);
    //RuleToRuleThemAll.Nodes = new NodeSet(RefElement, nodes, true);
    //RuleToRuleThemAll.Weights = weights;
    //        return RuleToRuleThemAll;





}


