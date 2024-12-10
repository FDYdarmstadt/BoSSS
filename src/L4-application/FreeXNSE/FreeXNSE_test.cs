using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using ilPSP.LinSolvers;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Utils;
using NUnit.Framework;

namespace FreeXNSE.Test {

    [TestFixture]
    public static class FreeXNSE_test {


        [Test]
        public static void FreeXNSE_ChannelTest() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.ChannelFlow();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();

            var Solution = ((XDGField)solver.CurrentState.Fields.Where(f => f.Identification == "VelocityX").Single()).GetSpeciesShadowField("A");
            var ExactSolution = new SinglePhaseField(Solution.Basis);
            ExactSolution.ProjectField(NonVectorizedScalarFunction.Vectorize(X => 1.0 - X[1] * X[1]));
            double err = Solution.L2Error(ExactSolution);

            Assert.IsTrue(err < 1e-10);
        }

        [Test]
        public static void FreeXNSE_CircleTest() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.Ellipse();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_EllipseTest() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.Ellipse(ecc: 0.8);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_EllipseTestParameterized() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.EllipseParameterized(ecc: 0.8);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_TestParameterizedAdvection() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.ParameterizedAdvection();
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }

        [Test]
        public static void FreeXNSE_EllipsoidTest() {
            var solver = new FreeXNSE();
            var C = FreeXNSE_controlfiles.Ellipsoid(ecc: 0.8);
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;
            solver.Init(C);
            solver.RunSolverMode();
        }
    }
}
