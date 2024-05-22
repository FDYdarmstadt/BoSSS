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

        //This would return a factory object with the configuration 
        //input: level set, tolerance etc.
        public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e, bool negativeVolume = false, bool callSurfaceAndVolumeAtOnce = true) {
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
                
                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        var quadRule = CalculateLevSetOnCell(cell, RequestedOrder);
                        ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));
                    }
                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            // to do
            public IEnumerable<IChunkRulePair<QuadRule>> CalculateQuadRuleSetCombo(ExecutionMask mask, int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        var quadRule = CalculateLevSetOnCell(cell, RequestedOrder);
                        //Algoim.GetVolumeQuadratureRules()
                        ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));

                    }
                }

                m_Rules.Add(RequestedOrder, ret.ToArray());
                return ret.ToArray();
            }

            private IEnumerable<IChunkRulePair<QuadRule>> GetChunkRulePair(Chunk chunk, int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                //QuadRule quadRule = new QuadRule();


                foreach (int cell in chunk.Elements) {
                    var quadRule = CalculateLevSetOnCell(cell, RequestedOrder);
                    //Algoim.GetVolumeQuadratureRules()
                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(chunk.i0), quadRule));

                }

                return ret;  ;

            }

            private QuadRule CalculateLevSetOnCell(int jCell, int RequestedOrder) {
                int pOrder = 3; //Math.Max(3, RequestedOrder); //ensure the polynomial interpolation is at least degree of 3.
                int n = (pOrder + 1) ;
                double[] points = GenericBlas.ChebyshevNodes(-1.0, 1.0, Math.Max(n, 3));
                //double[] points = GenericBlas.Linspace(-1.0, 1.0, Math.Max( n,3)); //testing

                //double[] points = ChebyshevPoints.Generate(n);


                int numberOfCombinations = (int)Math.Pow(points.Length, spaceDim); //Cartesian pair product for the points

                MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(points, spaceDim);

                //combinations.SaveToTextFile("combs.txt");


                NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, jCell, 1);
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


                double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

                    //new double[points.Length + points.Length];
                //Array.Copy(points, 0, x, 0, points.Length);
                //Array.Copy(points, 0, x, points.Length, points.Length);
                //int[] sizes = new int[] { 3, 3 };
                int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

                UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRule(spaceDim, pOrder, RequestedOrder, sizes, x, y);
                //qs.OutputQuadratureRuleAsVtpXML("bosssj" + jCell + ".vtp");
                //Console.WriteLine("qs.length = " + qs.length);
                QuadRule quadRule = QuadRule.CreateEmpty(RefElement, qs.length, qs.dimension);

                //for (int i = 0; i < q.length; i++) {
                //    writer.WriteString($"{q.nodes[i * dim]} {q.nodes[i * dim + 1]} {(dim == 3 ? q.nodes[i * dim + 2] : 0.0)}\n");
                //}

                //int ind = 0;
                for (int row = 0; row < qs.length; row++) {
                    quadRule.Weights[row] = qs.weights[row];

                    for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                        int ind = row * qs.dimension +d;
                        quadRule.Nodes[row, d] = qs.nodes[ind];
                        //Console.WriteLine($"quadRule.Nodes[{row}, {d}] = qs.nodes[{ind}] = {qs.nodes[ind]} ");
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

                //quadRuleGlobal.OutputQuadratureRuleAsVtpXML("bosssj" + jCell + ".vtp");
                return quadRule;
            }

        }

        class AlgoimEdgeRuleFactory : IQuadRuleFactory<CellBoundaryQuadRule> {
            public RefElement RefElement => throw new NotImplementedException();


            public AlgoimEdgeRuleFactory() { }

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<CellBoundaryQuadRule>> GetQuadRuleSet(ExecutionMask mask, int order) {
                throw new NotImplementedException();
            }



        }
    }






    #endregion




    //QuadRule RuleToRuleThemAll = QuadRule.CreateEmpty(RefElement, count, spatialDim);
    //RuleToRuleThemAll.Nodes = new NodeSet(RefElement, nodes, true);
    //RuleToRuleThemAll.Weights = weights;
    //        return RuleToRuleThemAll;





}


