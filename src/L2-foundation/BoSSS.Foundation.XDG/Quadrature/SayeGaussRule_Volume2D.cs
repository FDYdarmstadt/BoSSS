using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;

namespace BoSSS.Foundation.XDG.Quadrature
{

    /// <summary>
    /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
    /// </summary>
    public class SayeGaussRule_Volume2D : IQuadRuleFactory<QuadRule>
    {
        LevelSetTracker tracker;
        IRootFindingAlgorithm rootFinder;

        SayeGaussRule_Volume2D(LevelSetTracker Tracker, IRootFindingAlgorithm RootFinder)
        {
            tracker = Tracker;
            rootFinder = RootFinder;
        }

        #region IQaudRuleFactory<QuadRule>

        public RefElement RefElement {
            get {
                return Square.Instance;
            }
        }

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(order);
            var result = new List<ChunkRulePair<QuadRule>>();

            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    Integrand integrand = CreateSayeIntegrand( gaussRule_1D);
                    QuadRule sayeRule = integrand.Evaluate();
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(chunk, sayeRule);
                    result.Add(sayePair);
                }
            }
            return result;
        }

#endregion

        Integrand CreateSayeIntegrand(QuadRule Rule)
        {
            throw new NotImplementedException();
        }

        class Integrand
        {
            QuadRule rule;
            TreeNode<Argument> recursionInfo;

            Integrand()
            {

            }

            public QuadRule Evaluate()
            {
                throw new NotImplementedException();
            }
        }
    }
}
