using BoSSS.Foundation.XDG;

namespace IntersectingLevelSetTest {
    class Program {
        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;
            //AllUpTest.AllUp(2, XQuadFactoryHelper.MomentFittingVariants.Saye);
            //AllUpTest.LocalTestWithPlotting(2, XQuadFactoryHelper.MomentFittingVariants.Saye);

            //AllUpTest.ParabolaTest(3);
            //AllUpTest.TwoStraightTest(3);
            //AllUpTest.TransformTest(3);
            //AllUpTest.Convergence2DTest(1);


            //AllUpTest.Rotation3DTest(3);
            //AllUpTest.Transform3DTest(2);
            //AllUpTest.MovingSphere3DTest(2);
            //AllUpTest.Convergence3DTest(3);

            //AllUpTest.RandomTest(3);
            //AllUpTest.RealTransformTest(3);


            //BoSSS.Solution.Application<PlotControl>._Main(
            //    args,
            //    true,
            //    () => new ZwoLsCoupledSolver<PlotControl>() {MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye });
        }
    }
}
