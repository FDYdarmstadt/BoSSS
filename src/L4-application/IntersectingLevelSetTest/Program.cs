using BoSSS.Foundation.XDG;

namespace IntersectingLevelSetTest {
    class Program {
        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;
            //AllUpTest.AllUp(2, XQuadFactoryHelper.MomentFittingVariants.Saye);
            //AllUpTest.LocalTestWithPlotting(2, XQuadFactoryHelper.MomentFittingVariants.Saye);

            AllUpTest.ParabolaTest(2);
            AllUpTest.TwoStraightTest(2);
            AllUpTest.TransformTest(2);


            AllUpTest.Rotation3DTest(2);
            AllUpTest.Transform3DTest(2);
            AllUpTest.MovingSphere3DTest(2);

            //AllUpTest.RandomTest(2);
            //AllUpTest.RealTransformTest(2);


            //BoSSS.Solution.Application<PlotControl>._Main(
            //    args,
            //    true,
            //    () => new ZwoLsCoupledSolver<PlotControl>() {MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye });
        }
    }
}
