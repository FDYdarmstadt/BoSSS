using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using log4net.Core;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Remoting.Messaging;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.LevelSetTracker;
using static BoSSS.Foundation.XDG.XQuadFactoryHelper;

namespace BoSSS.Foundation.XDG.Quadrature
{
    class CombinedLevelSetQuadrature
    {
        CombinedLevelSetData combinedLsDat;

        XQuadFactoryHelper combinedLevSetQuadFactories;

        public CombinedLevelSetQuadrature(LevelSetData[] lsDatas)
        {
            combinedLsDat = new CombinedLevelSetData(lsDatas, 3);
            combinedLevSetQuadFactories = new XQuadFactoryHelper(combinedLsDat.GetCombinedLevelSetDatas(), MomentFittingVariants.Saye);
        }

        public IQuadRuleFactory<QuadRule> GetEdgeRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol)
        {
            int combinedLevSetIndex = combinedLsDat.GetCombinedLevelSetIndice(levSetIndex0, jmp0, levSetIndex1, jmp1);
            return combinedLevSetQuadFactories.GetEdgeRuleFactory(combinedLevSetIndex, JumpTypes.Heaviside, KrefVol);
        }

        public IQuadRuleFactory<QuadRule> GetSurfaceFactory(int levSetIndex0, 
            int levSetIndex1, 
            JumpTypes jmp1, 
            RefElement KrefVol, 
            IQuadRuleFactory<QuadRule> levSet1SurfaceRule)
        {
            int combinedLevSetIndex0 = combinedLsDat.GetCombinedLevelSetIndice(levSetIndex0, JumpTypes.Heaviside, levSetIndex1, jmp1);
            int combinedLevSetIndex1 = combinedLsDat.GetCombinedLevelSetIndice(levSetIndex0, JumpTypes.OneMinusHeaviside, levSetIndex1, jmp1);
            //Get Rule for all 3 and then  (0 + 1 - original) / 2

            var rule0 = combinedLevSetQuadFactories.GetSurfaceFactory(combinedLevSetIndex0, KrefVol);
            var rule1 = combinedLevSetQuadFactories.GetSurfaceFactory(combinedLevSetIndex1, KrefVol);

            var rule0plus1 = QuadRuleFactoryArithmetic.Add(rule0, rule1);
            var rule0plus1minusLevSet1Rule = QuadRuleFactoryArithmetic.Add(QuadRuleFactoryArithmetic.Scale(-1, levSet1SurfaceRule), rule0plus1);
            var finalRule = QuadRuleFactoryArithmetic.Scale(0.5, rule0plus1minusLevSet1Rule);

            return finalRule;
        }

        public IQuadRuleFactory<QuadRule> GetVolRuleFactory(int levSetIndex0, JumpTypes jmp0, int levSetIndex1, JumpTypes jmp1, RefElement KrefVol)
        {
            int combinedLevSetIndex = combinedLsDat.GetCombinedLevelSetIndice(levSetIndex0, jmp0, levSetIndex1, jmp1);
            return combinedLevSetQuadFactories.GetVolRuleFactory(combinedLevSetIndex, JumpTypes.Heaviside, KrefVol);
        }
    }

    static class QuadRuleFactoryArithmetic
    {
        public static IQuadRuleFactory<QuadRule> Add(IQuadRuleFactory<QuadRule> a, IQuadRuleFactory<QuadRule> b) 
        {
            return new QuadRuleFactorySum(a, b);
        }

        public static IQuadRuleFactory<QuadRule> Scale(double scale, IQuadRuleFactory<QuadRule> b)
        {
            return new QuadRuleFactoryScale(b, scale);
        }

    }

    class QuadRuleFactorySum : IQuadRuleFactory<QuadRule>
    {
        public RefElement RefElement => a.RefElement;

        IQuadRuleFactory<QuadRule> a;

        IQuadRuleFactory<QuadRule> b;

        public QuadRuleFactorySum(IQuadRuleFactory<QuadRule> a, IQuadRuleFactory<QuadRule> b)
        {
            if(a.RefElement != b.RefElement)
            {
                throw new ArgumentException("Ref Elements must be equal");
            }
            this.a = a;
            this.b = b;
        }

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            IEnumerator<IChunkRulePair<QuadRule>> ruleA = a.GetQuadRuleSet(mask, order).GetEnumerator();
            IEnumerator<IChunkRulePair<QuadRule>> ruleB = b.GetQuadRuleSet(mask, order).GetEnumerator();
            List<IChunkRulePair<QuadRule>> ruleSum = new List<IChunkRulePair<QuadRule>>();
            while(ruleA.MoveNext() && ruleB.MoveNext())
            {
                IChunkRulePair<QuadRule> currentA = ruleA.Current;
                IChunkRulePair<QuadRule> currentB = ruleB.Current;
                ChunkRulePair<QuadRule> sum = Add(currentA, currentB);
                ruleSum.Add(sum);
            }

            return ruleSum;
        }

        static ChunkRulePair<QuadRule> Add(IChunkRulePair<QuadRule> currentA, IChunkRulePair<QuadRule> currentB)
        {
            if (!currentA.Chunk.Equals(currentB.Chunk))
            {
                throw new Exception("Chunks are not equal");
            };
            QuadRule sum = Add(currentA.Rule, currentB.Rule);
            ChunkRulePair<QuadRule> pairSum = new ChunkRulePair<QuadRule>(currentA.Chunk, sum);
            return pairSum;
        }

        static QuadRule Add(QuadRule A, QuadRule B)
        {
            int noOfNodes = A.NoOfNodes + B.NoOfNodes;
            int D = A.SpatialDim;
            QuadRule sum = QuadRule.CreateEmpty(A.RefElement, noOfNodes, D);

            sum.Nodes.SetSubArray(A.Nodes, new int[] { 0,0}, new int[] {A.NoOfNodes - 1, D - 1});
            sum.Weights.SetSubArray(A.Weights, new int[] { 0 }, new int[] { A.NoOfNodes - 1});

            sum.Nodes.SetSubArray(B.Nodes, new int[] { A.NoOfNodes, 0 }, new int[] { noOfNodes - 1, D - 1});
            sum.Weights.SetSubArray(B.Weights, new int[] { A.NoOfNodes }, new int[] { noOfNodes - 1 });

            sum.Nodes.LockForever();
            return sum;
        }
    }

    class QuadRuleFactoryScale : IQuadRuleFactory<QuadRule>
    {
        public RefElement RefElement => a.RefElement;

        IQuadRuleFactory<QuadRule> a;

        double scale;

        public QuadRuleFactoryScale(IQuadRuleFactory<QuadRule> a, double scale)
        {
            this.a = a;
            this.scale = scale;
        }

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            IEnumerator<IChunkRulePair<QuadRule>> ruleA = a.GetQuadRuleSet(mask, order).GetEnumerator();
            List<IChunkRulePair<QuadRule>> ruleSum = new List<IChunkRulePair<QuadRule>>();
            while (ruleA.MoveNext())
            {
                IChunkRulePair<QuadRule> currentA = ruleA.Current;
                ChunkRulePair<QuadRule> sum = Scale(currentA, scale);
                ruleSum.Add(sum);
            }
            return ruleSum;
        }

        static ChunkRulePair<QuadRule> Scale(IChunkRulePair<QuadRule> currentA, double scale)
        {
            QuadRule scaled = Scale(currentA.Rule, scale);
            ChunkRulePair<QuadRule> scaledPair = new ChunkRulePair<QuadRule>(currentA.Chunk, scaled);
            return scaledPair;
        }

        static QuadRule Scale(QuadRule A, double scale)
        {
            QuadRule scaled = QuadRule.CreateEmpty(A.RefElement, A.NoOfNodes, A.SpatialDim);
            scaled.Nodes.Set(A.Nodes);
            scaled.Nodes.LockForever();
            scaled.Weights.Set(A.Weights);
            scaled.Weights.Scale(scale);
            return scaled;
        }
    }
}
