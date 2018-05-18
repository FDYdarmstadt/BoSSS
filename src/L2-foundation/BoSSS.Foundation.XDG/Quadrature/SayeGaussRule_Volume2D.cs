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
        SayeGaussRule_Volume2D(LevelSetTracker tracker, IRootFindingAlgorithm rootFinder)
        {
            

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
            Line.Instance.GetQuadratureRule(order);
            //Find quadrature nodes and weights in each cell/chunk

                //Find Points and weights with saye

                
            throw new NotImplementedException();
        }

#endregion

        private object SayeRecursion()
        {
            return new object(); 
        }
    }
}
