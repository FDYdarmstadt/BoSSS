using System;
using BoSSS.Foundation.IO;
using BoSSS.Solution;

namespace TwoPhaseSimple {

    class MyClass : Application {

        static void Main(string[] args) {
            Application._Main(args, true, null, delegate() {
                return new MyClass();
            });
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            throw new NotImplementedException();
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L)
        {
            throw new NotImplementedException();
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0)
        {
            throw new NotImplementedException();
        }

    }
}
