using BoSSS.Foundation.XDG;

namespace IntersectingLevelSetTest
{
    class Program
    {
        static void Main(string[] args)
        {            
            XQuadFactoryHelper.CheckQuadRules = true;
            //AllUpTest.AllUp(0, XQuadFactoryHelper.MomentFittingVariants.Saye);
            AllUpTest.LocalTestWithPlotting(2, XQuadFactoryHelper.MomentFittingVariants.Saye);

            //BoSSS.Solution.Application<PlotControl>._Main(
            //    args,
            //    true,
            //    () => new ZwoLsCoupledSolver<PlotControl>() {MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye });
        }
    }
}
