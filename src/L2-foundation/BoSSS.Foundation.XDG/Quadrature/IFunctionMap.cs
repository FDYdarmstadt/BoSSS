using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    interface IFunctionMap {

        RefElement Domain { get; }

        HyperRectangle Codomain(int j);

        QuadRule MapFromCodomainToDomain(QuadratureRule rule, int j);

        IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int j);
    }

}
