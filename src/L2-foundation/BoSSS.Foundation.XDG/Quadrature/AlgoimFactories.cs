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
            return new Factory() {
                m_Owner = this
            };
        }

        public IQuadRuleFactory<QuadRule> GetVolumeFactory() {
            return new Factory() {
                m_Owner = this
            };
        }

        //This would return an factory object with the configuration 
        //input: level set, tolerance etc.
        public AlgoimFactories(LevelSetTracker.LevelSetData ls, RefElement e) {
            refElement = e;
            lsData = ls;
        }

        LevelSetTracker.LevelSetData lsData;
        RefElement refElement;

        #region Edge rules


        class Factory : IQuadRuleFactory<QuadRule> {
            internal AlgoimFactories m_Owner;

            public RefElement RefElement => m_Owner.refElement;

            int spaceDim => RefElement.SpatialDimension;


            public LevelSetTracker.LevelSetData lsData => m_Owner.lsData;

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");
                //if (mask.MaskType != MaskType.Geometrical)
                //    throw new ArgumentException("Expecting a geometrical mask.");

                //CellMask _mask = mask as CellMask;
                //SubGrid sgrd = new SubGrid(_mask);

                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                foreach (Chunk chunk in mask) {
                    foreach (int cell in chunk.Elements) {
                        var quadRule = calculateLSonCell(cell, RequestedOrder);
                        //Algoim.GetVolumeQuadratureRules()
                        ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(cell), quadRule));

                    }
                }

                //    quadRule.Nodes =;
                //quadRule.Weights;


                //VolumeRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), VolRule);

                return ret;

            }

            private IEnumerable<IChunkRulePair<QuadRule>> GetChunkRulePair(Chunk chunk, int RequestedOrder) {
                List<ChunkRulePair<QuadRule>> ret = new List<ChunkRulePair<QuadRule>>();

                //QuadRule quadRule = new QuadRule();


                foreach (int cell in chunk.Elements) {
                    var quadRule = calculateLSonCell(cell, RequestedOrder);
                    //Algoim.GetVolumeQuadratureRules()
                    ret.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(chunk.i0), quadRule));

                }

                return ret;  ;

            }

            public class ChebyshevPoints {
                public static double[] Generate(int n) {
                    if (n < 3) n = 3; // Ensure there are at least three points
                    double[] points = new double[n];
                    for (int i = 0; i < n; i++) {
                        points[i] = Math.Cos(Math.PI * (2 * i + 1) / (2 * n));
                    }
                    return points;
                }
            }

            private QuadRule calculateLSonCell(int j0, int RequestedOrder) {
                int n = (RequestedOrder + 3) ;
                double[] points = GenericBlas.Linspace(-1, 1, Math.Max( n,3)); //testing

                //double[] points = ChebyshevPoints.Generate(n);


                int numberOfCombinations = (int)Math.Pow(points.Length, spaceDim); //Cartesian pair product for the points

                MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(points, spaceDim);

                //combinations.SaveToTextFile("combs.txt");


                NodeSet NS = new NodeSet(RefElement, combinations, false);
                var ret = lsData.GetLevSetValues(NS, j0, 1);
                //Console.WriteLine(ret.ToString());
                //ret.SaveToTextFile("ret.txt");

                //Console.WriteLine("Ret dim: " + ret.Dimension);
                //foreach (var l in ret.Lengths)
                //    Console.WriteLine("Length in:  " + l);

                //Console.WriteLine("Total number of elements: " + ret.Length);

                double[] y = new double[ret.Length];

                for (int i = 0; i < ret.Length; i++) 
                    y[i] = ret[0,i];

                double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

                    //new double[points.Length + points.Length];
                //Array.Copy(points, 0, x, 0, points.Length);
                //Array.Copy(points, 0, x, points.Length, points.Length);
                //int[] sizes = new int[] { 3, 3 };
                int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

                var qs = Algoim.GetVolumeQuadratureRules(spaceDim, RequestedOrder, sizes, x, y);
                //qs.OutputQuadratureRuleAsVtpXML("bosssj" + j0 + ".vtp");
                //Console.WriteLine("qs.length = " + qs.length);
                QuadRule quadRule = QuadRule.CreateEmpty(RefElement,qs.length, qs.dimension);

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
                //quadRule.OutputQuadratureRuleAsVtpXML("quadRule" + j0 + ".vtp");
                quadRule.Nodes.LockForever();
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


