/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace IntersectingLevelSetTest {

    [TestFixture]
    static public class AllUpTest {

        [OneTimeSetUp]
        static public void SetUp()
        {
            BoSSS.Solution.Application.InitMPI();
        }

        [Test]
        static public void AllUp(
            [Values(1, 2, 3)] int DGdegree,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants quadVariant)            
         {
            var C = new PlotControl();
            ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl> p = null;
            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                delegate() {
                    p = new ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl>();
                    p.DEGREE = DGdegree;
                    p.MomentFittingVariant = quadVariant;
                    return p;
                });
        }

        static public void LocalTestWithPlotting(
            [Values(1, 2, 3)] int DGdegree,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants quadVariant)
        {
            BoSSS.Solution.Application.InitMPI();
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            var C = new PlotControl();
            var p = new ZwoLsSolver<PlotControl>();
            p.DEGREE= DGdegree;    
            p.MomentFittingVariant = quadVariant;
            C.SuperSampling = 5;
            p.Init(C);
            p.RunSolverMode();
        }
        //static public void OneCellTwoIntersectingLevelSets()
        //{
        //    var C = new PlotControl();
        //    C.SetDGdegree(0);

        //    C.SetGrid(Grid2D.Cartesian2DGrid(new double[] {0,1}, new double[] { 0, 2}));

        //    ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl> p = null;
        //    BoSSS.Solution.Application._Main(
        //        new string[0],
        //        true,
        //        delegate () {
        //            p = new ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl>();
        //            p.DEGREE = 0;
        //            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //            return p;
        //        });
        //}
    }
}
