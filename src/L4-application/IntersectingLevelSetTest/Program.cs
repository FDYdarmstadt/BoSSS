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
                () => new ZwoLsSolver<PlotControl>() { DEGREE = 3, MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye });
        }
    }
}
