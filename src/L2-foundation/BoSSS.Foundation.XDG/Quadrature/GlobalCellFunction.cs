using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using IntersectingQuadrature.Tensor;
using System.Diagnostics;
using static BoSSS.Foundation.XDG.LevelSetTracker;

namespace BoSSS.Foundation.XDG.Quadrature.Intersecting {



    class GlobalCellFunction : IScalarFunction {

        readonly LevelSetData levelSet;
        readonly RefElement Kref;
        readonly int jCell;
        /*
        static public int EvalCounterV = 0;
        static public int EvalCounterG = 0;
        static public int EvalCounterH = 0;
        static public Stopwatch stwV = new Stopwatch();
        static public Stopwatch stwG = new Stopwatch();
        static public Stopwatch stwH = new Stopwatch();
        */

        public GlobalCellFunction(LevelSetData levelSet, int jCell) {
            this.jCell = jCell;
            this.M = levelSet.GridDat.SpatialDimension;
            this.Kref = levelSet.GridDat.iGeomCells.GetRefElement(jCell);
            this.levelSet = levelSet;
        }

        //Spatial Dimension
        public int M { get; private set; }

        public double Evaluate(Tensor1 x) {
            //stwV.Start();
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            //EvalCounterV++;
            //stwV.Stop();
            return _evaluation[0, 0];
        }

        public (double evaluation, Tensor1 gradient) EvaluateAndGradient(Tensor1 x) {
            //stwG.Start();
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            var _gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            //EvalCounterG++;
            //stwG.Stop();
            return (_evaluation[0, 0], ToTensor1(_gradient));
        }

        public (double evaluation, Tensor1 gradient, Tensor2 hessian) EvaluateAndGradientAndHessian(Tensor1 x) {
            //stwH.Start(); 
            NodeSet X = FromGlobalToReferenceNodeSet(x);
            X.LockForever();
            var _evaluation = levelSet.GetLevSetValues(X, jCell, 1);
            var _gradient = levelSet.GetLevelSetReferenceGradients(X, jCell, 1);
            var _hessian = levelSet.GetLevelSetReferenceHessian(X, jCell, 1);
            //EvalCounterH++;
            //stwH.Stop();
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
