using BoSSS.Foundation.XDG;

namespace IntersectingLevelSetTest {
    class Program {
        static void Main(string[] args) {
            XQuadFactoryHelper.CheckQuadRules = true;

            IntegrationTests.SquareTest(3);
            AllUpTest.ParabolaTestSaye(3);


            //BoSSS.Solution.Application<PlotControl>._Main(
            //    args,
            //    false,
            //    () => new ZwoLsCoupledSolver<PlotControl>() { MomentFittingVariant = CutCellQuadratureMethod.Saye });
        }
    }
}
