using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using System;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature {

    class GlobalCellMapping : IFunctionMap {

        int dim;
        GridData grid;

        public RefElement Domain { get; private set; }

        public GlobalCellMapping(GridData grid) {
            this.grid = grid;
            dim = grid.SpatialDimension;

            switch(dim) {
                case 1:
                Domain = Line.Instance;
                break;
                case 2:
                Domain = Square.Instance;
                break;
                case 3:
                Domain = Cube.Instance;
                break;
                default:
                throw new NotImplementedException();
            }
        }

        public HyperRectangle Codomain(int jCell) {
            HyperRectangle codomain = new HyperRectangle(grid.SpatialDimension);
            codomain.Dimension = grid.SpatialDimension;
            for(int i = 0; i < codomain.Dimension; ++i) {
                codomain.Center[i] = 0.0;
                codomain.Diameters[i] = 2.0;
            }
            return codomain;
        }

        public QuadRule MapFromCodomainToDomain(QuadratureRule rule, int jCell) {
            QuadRule q = QuadRule.CreateBlank(Domain, rule.Count, dim);
            for(int i = 0; i < rule.Count; ++i) {
                for(int d = 0; d < dim; ++d) {
                    q.Nodes[i, d] = rule[i].Point[d];
                }
            }
            q.Nodes.LockForever();

            for(int i = 0; i < rule.Count; ++i) {
                q.Weights[i] = rule[i].Weight;
            }
            return q;
        }

        public IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int jCell) {
            return new GlobalCellFunction(levelSet, jCell);
        }
    }
}
