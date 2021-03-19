using BoSSS.Foundation.XDG;

namespace IntersectingLevelSetTest
{
    class Program
    {
        static void Main(string[] args)
        {
            XQuadFactoryHelper.CheckQuadRules = true;
            
            BoSSS.Solution.Application<PlotControl>._Main(
                args,
                true,
                () => new ZwoLsSinglePhaseSolver<PlotControl>() {MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye });
        }
    }
}
