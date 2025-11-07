using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform.LinAlg;
using ilPSP;
using IntersectingQuadrature;
using IntersectingQuadrature.Tensor;
using System;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    class GlobalEdgeMapping : IFunctionMap {

        class Selection {

            int except;

            public Selection(int except) {
                this.except = except;
            }

            public int FromSubIndexToIndex(int i) {
                if(i < except) {
                    return i;
                } else {
                    return ++i;
                }
            }
        }

        GridData grid;

        public RefElement Domain { get; private set; }

        RefElement CellDomain;

        public GlobalEdgeMapping(GridData grid) {
            this.grid = grid;
            int dim = grid.SpatialDimension;

            switch(dim - 1) {
                case 1:
                Domain = Line.Instance;
                CellDomain = Square.Instance;
                break;
                case 2:
                Domain = Square.Instance;
                CellDomain = Cube.Instance;
                break;
                default:
                throw new NotImplementedException();
            }
            //center = MultidimensionalArray.Create(1, dim);
        }

        //MultidimensionalArray center;

        int EmptyDim(MultidimensionalArray center) {
            double max = double.MinValue;
            int emptyDim = -1;
            for(int i = 0; i < center.GetLength(1); ++i) {
                if(Math.Abs(center[0, i]) > max) {
                    emptyDim = i;
                    max = Math.Abs(center[0, i]);
                }
            }
            return emptyDim;
        }

        public HyperRectangle Codomain(int jEdge) {
            HyperRectangle codomain = new HyperRectangle(Domain.SpatialDimension);
            codomain.Dimension = Domain.SpatialDimension;

            //int jCell = grid.Edges.CellIndices[jEdge, 0];
            //byte iFace = grid.Edges.FaceIndices[jEdge, 0];
            //NodeSet faceCenter = CellDomain.GetFaceCenter(iFace);
            //int emptyDim = EmptyDim(faceCenter);

            //grid.TransformLocal2Global(faceCenter, center, jCell);
            //Selection edgeToCell = new Selection(emptyDim);

            for(int i = 0; i < Domain.SpatialDimension; ++i) {
                //int j = edgeToCell.FromSubIndexToIndex(i);
                codomain.Center[i] = 0.0; // center[0, j];
                codomain.Diameters[i] = 2.0; // grid.Cells.Transformation[jCell, j, j] * 2;
            }
            return codomain;
        }

        public QuadRule MapFromCodomainToDomain(QuadratureRule rule, int jEdge) {

            MultidimensionalArray globalNodes = MultidimensionalArray.Create(rule.Count, Domain.SpatialDimension);
            for(int i = 0; i < rule.Count; ++i) {
                for(int d = 0; d < Domain.SpatialDimension; ++d) {
                    globalNodes[i, d] = rule[i].Point[d];
                }
            }

            AffineTrafo t = FromCodomainToDomain(jEdge);
            QuadRule q = QuadRule.CreateBlank(Domain, rule.Count, Domain.SpatialDimension);
            t.Transform(globalNodes, q.Nodes);

            for(int i = 0; i < rule.Count; ++i) {
                q.Weights[i] = rule[i].Weight;
            }
            q.Nodes.LockForever();
            return q;
        }

        AffineTrafo FromCodomainToDomain(int jEdge) {

            byte iFace = grid.Edges.FaceIndices[jEdge, 0];
            NodeSet faceCenter = CellDomain.GetFaceCenter(iFace);
            int emptyDim = EmptyDim(faceCenter);
            Selection edgeToCell = new Selection(emptyDim);

            AffineTrafo T = grid.Edges.Edge2CellTrafos[grid.Edges.Edge2CellTrafoIndex[jEdge, 0]];
            AffineTrafo t = new AffineTrafo(Domain.SpatialDimension);
            HyperRectangle codomain = Codomain(jEdge);

            for(int i = 0; i < Domain.SpatialDimension; ++i) {
                int k = edgeToCell.FromSubIndexToIndex(i);
                for(int j = 0; j < Domain.SpatialDimension; ++j) {
                    t.Affine[i] += T.Matrix[k, j] * -codomain.Center[j] * 2.0 / codomain.Diameters[j];
                    t.Matrix[i, j] += T.Matrix[k, j] * 2.0 / codomain.Diameters[j];
                }
            }
            return t;
        }

        public IScalarFunction MapFromDomainToCodomain(LevelSetData levelSet, int jEdge) {
            int jCell = grid.Edges.CellIndices[jEdge, 0];
            IScalarFunction globalLevelSet = new GlobalCellFunction(levelSet, jCell);
            IVectorFunction T = FromCodomainToEmbeddedCodomain(jEdge);
            IScalarFunction alpha = new ScalarComposition(globalLevelSet, T);
            return alpha;
        }

        IVectorFunction FromCodomainToEmbeddedCodomain(int jEdge) {
            int jCell = grid.Edges.CellIndices[jEdge, 0];
            byte iFace = grid.Edges.FaceIndices[jEdge, 0];
            NodeSet faceCenter = CellDomain.GetFaceCenter(iFace);
            int emptyDim = EmptyDim(faceCenter);

            //grid.TransformLocal2Global(faceCenter, center, jCell);
            Selection edgeToCell = new Selection(emptyDim);

            Tensor2 m = Tensor2.Zeros(CellDomain.SpatialDimension, Domain.SpatialDimension);
            for(int i = 0; i < Domain.SpatialDimension; ++i) {
                int j = edgeToCell.FromSubIndexToIndex(i);
                m[j, i] = 1.0;
            }

            Tensor1 affine = Tensor1.Zeros(CellDomain.SpatialDimension);
            //affine[emptyDim] = center[0, emptyDim];
            affine[emptyDim] = faceCenter[0, emptyDim];

            LinearVectorPolynomial t = new LinearVectorPolynomial(affine, m);
            return t;
        }
    }

}
