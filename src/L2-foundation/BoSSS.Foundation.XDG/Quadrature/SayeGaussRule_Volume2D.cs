using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using static BoSSS.Foundation.XDG.Quadrature.HMF.LineSegment;
using System.Collections;

namespace BoSSS.Foundation.XDG.Quadrature
{

    /// <summary>
    /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
    /// </summary>
    public class SayeGaussRule_Volume2D :
        SayeIntegrand<Square>,
        IQuadRuleFactory<QuadRule>
    {
        LevelSetTracker.LevelSetData lsData;
        IRootFindingAlgorithm rootFinder;

        public SayeGaussRule_Volume2D(
            LevelSetTracker.LevelSetData _lsData, 
            IRootFindingAlgorithm RootFinder)
        {
            this.lsData = _lsData;
            this.rootFinder = RootFinder;

            //var ns = new NodeSet(this.RefElement, 3, 1);
            //ns[2, 1] = 0.3;
            //ns.LockForever();

            //lsData.GetLevelSetReferenceGradients(RefElement.Center, 0, 1);

        }

        #region IQaudRuleFactory<QuadRule>

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            QuadRule gaussRule_1D = Line.Instance.GetQuadratureRule(order);
            var result = new List<ChunkRulePair<QuadRule>>();

            //Find quadrature nodes and weights in each cell/chunk
            foreach (Chunk chunk in mask)
            {
                foreach (int cell in chunk.Elements)
                {
                    QuadRule sayeRule = this.Evaluate(cell, gaussRule_1D);
                    ChunkRulePair<QuadRule> sayePair = new ChunkRulePair<QuadRule>(chunk, sayeRule);
                    result.Add(sayePair);
                }
            }
            return result;
        }

        public RefElement RefElement {
            get {
                return Square.Instance;
            }
        }

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        #endregion

        #region SayeIntegrand

        protected override void SetCellCenter(SayeArgument<Square> arg)
        {
            throw new NotImplementedException();
        }

        protected override double EvaluateBounds(Square psi)
        {
            throw new NotImplementedException();
        }

        protected override Dim FindPromisingHeightDirection()
        {
            throw new NotImplementedException();
        }

        protected override bool HeightDirectionIsSuitable()
        {
            throw new NotImplementedException();
        }

        protected override TreeNode<SayeArgument<Square>>[] 
            SubdivideNode(TreeNode<SayeArgument<Square>> node)
        {
            throw new NotImplementedException();
        }

        protected override bool SubdivideSuitable(SayeArgument<Square> arg)
        {
            throw new NotImplementedException();
        }

        #endregion

    }

    public class Square : 
        IPsi
    {
        Square()
        {

        }

        public double EvaluateAt(double[] point)
        {
            throw new NotImplementedException();
        }
    }

}

