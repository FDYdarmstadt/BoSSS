using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using IntersectingQuadrature.Tensor;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {
    class GlobalCellFunction : IScalarFunction {

        readonly LevelSetData levelSet;
        readonly RefElement Kref;
        readonly int jCell;

        public GlobalCellFunction(LevelSetData levelSet, int jCell) {
            this.jCell = jCell;
            this.M = levelSet.GridDat.SpatialDimension;
            this.Kref = levelSet.GridDat.iGeomCells.GetRefElement(jCell);
            this.levelSet = levelSet;
        }

        //Spatial Dimension
        public int M { get; private set; }

        public double Evaluate(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            return _evaluation[0, 0];
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            var _gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            return (_evaluation[0, 0], ToTensor1(_gradient));
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            X.LockForever();
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            var _gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            var _hessian = levelSet.GetLevelSetReferenceHessian(X, jCell, 1);
            return (_evaluation[0, 0], ToTensor1(_gradient), ToTensor2(_hessian));
        }

        NodeSet FromGlobalToReferenceNodeSet(Tensor1 x) {
            NodeSet n = new NodeSet(Kref, 1, x.M, false);
            for(int i = 0; i < x.M; ++i) {
                n[0, i] = x[i];
            }
            n.LockForever();
            return n;
        }

        static Tensor1 ToTensor1(MultidimensionalArray x) {
            Tensor1 X = Tensor1.Zeros(x.Length);
            for(int i = 0; i < X.M; ++i) {
                X[i] = x[0, 0, i];
            }
            return X;
        }

        static Tensor2 ToTensor2(MultidimensionalArray x) {
            Tensor2 X = Tensor2.Zeros(x.GetLength(x.Dimension - 1));
            for(int i = 0; i < X.M; ++i) {
                for(int j = 0; j < X.M; ++j) {
                    X[i, j] = x[0, 0, i, j];
                }
            }
            return X;
        }
    }
}
