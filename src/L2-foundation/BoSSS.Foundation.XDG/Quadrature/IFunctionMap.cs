using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    public interface IFunctionMap {

        RefElement Domain { get; }

        HyperRectangle Codomain(int j);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="rule">
        /// quadrature rule within the co-domain returned by <see cref="Codomain(int)"/>
        /// </param>
        /// <param name="j">
        /// cell or edge index
        /// </param>
        /// <returns></returns>
        QuadRule MapFromCodomainToDomain(QuadratureRule rule, int j);

        IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int j);
    }

}
