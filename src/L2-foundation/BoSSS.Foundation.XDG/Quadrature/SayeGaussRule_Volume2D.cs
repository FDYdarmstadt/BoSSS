using BoSSS.Foundation.Quadrature;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.XDG.Quadrature
{

    /// <summary>
    /// Gauss rules for \f$ \int_{\frakA \cap K_j } \ldots \dV \f$ in the 2D case
    /// </summary>
    public class SayeGaussRule_Volume2D : IQuadRuleFactory<QuadRule>
    {
        public RefElement RefElement => throw new NotImplementedException();

        public int[] GetCachedRuleOrders()
        {
            throw new NotImplementedException();
        }

        public IEnumerable<IChunkRulePair<QuadRule>> GetQuadRuleSet(ExecutionMask mask, int order)
        {
            throw new NotImplementedException();
        }
    }
}
