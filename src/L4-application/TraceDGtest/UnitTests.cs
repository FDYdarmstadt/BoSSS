using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.TraceDGtest {
    public static class UnitTests {


        [Test]
        public static void CartesianWithLevelSetAMR([Values(2, 3, 4)] int degree,
                                                    [Values(false, true)] bool useAMR) {
            var c = ControlExamples.SteadystateWithLevelSetAMR(degree, useAMR);
            c.ImmediatePlotPeriod = 1;


            using(var s = new TraceDGtestMain()) {
                s.Init(c);
                s.RunSolverMode();

                foreach(var f in s.CurrentResidual.Fields) {
                    var f_L2 = f.L2Norm();
                    Console.WriteLine($"  L2-Norm of {f}: {f_L2:0.####e-00}");
                }
                foreach(var f in s.CurrentResidual.Fields) {
                    var f_L2 = f.L2Norm();
                    Assert.LessOrEqual(f_L2, 1.0e-8,  $"  L2-Norm of {f} is to large.");
                }
            }

        }



    }
}
