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
                    ret.Add(GetChunkRulePair(chunk));
                }

                //    quadRule.Nodes =;
                //quadRule.Weights;


                //VolumeRule[jSub] = new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), VolRule);

                return ret;

            }

            private ChunkRulePair<QuadRule> GetChunkRulePair(Chunk chunk) {
                QuadRule quadRule = new QuadRule();


                foreach (int cell in chunk.Elements) {
                    var q = calculateLSonCell(cell);
                    //Algoim.GetVolumeQuadratureRules()

                    quadRule = q;

                }

                return new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(chunk.i0), quadRule);

            }

            private QuadRule calculateLSonCell(int j0) {
                double[] points = GenericBlas.Linspace(-1, 1, 3); //testing

                MultidimensionalArray combinations = MultidimensionalArray.Create(points.Length * points.Length, 2);
                combinations.SetAll(-1);

                // Fill the array with all combinations of the points array
                int index = 0;
                for (int j = 0; j < points.Length; j++) {
                    for (int i = 0; i < points.Length; i++) {
                        combinations[index, 0] = points[i];
                        combinations[index, 1] = points[j];
                        //Console.WriteLine("combinations[" + index + ", 0]: " + combinations[index, 0]);
                        //Console.WriteLine("combinations[" + index + ", 1]: " + combinations[index, 1]);

                        index++;
                    }
                }
                combinations.SaveToTextFile("combds.txt");

                double[] x = new double[points.Length+ points.Length];


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


                Array.Copy(points, 0, x, 0, points.Length);
                Array.Copy(points, 0, x, points.Length, points.Length);
                int[] sizes = new int[] { 3, 3 };

                var qsExample = Algoim.GetVolumeQuadratureRulesTest();

                var qs = Algoim.GetVolumeQuadratureRules(2, 5, sizes, x, y);
                //Console.WriteLine("qs.length = " + qs.length);
                QuadRule quadRule = QuadRule.CreateEmpty(RefElement,qs.length, qs.dimension);

                int ind = 0;
                for (int row = 0; row < qs.length; row++) {
                    quadRule.Weights[row] = qs.weights[row];

                    for (int d = 0; d < qs.dimension; d++, ind++) // map 1d array back to 2d
                    quadRule.Nodes[row, d] = qs.nodes[ind];
                }
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


