using BoSSS.Foundation.XDG;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.CutCellQuadratureScaling {

    [TestFixture]
    static public class AllTests {


        [Test]
        public static void OneLevelSet_2D(
#if DEBUG
            [Values(3, 6, 7)] int quadOrder,
#else
            [Values(3, 4, 7, 8, 9, 10)] int quadOrder,
#endif
            [Values(
            XQuadFactoryHelper.MomentFittingVariants.Classic, XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants cutCellQuadType
            ) {
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_2D


            using(var Ref = new TestSetupSingleLevset2D(1.0, quadOrder, cutCellQuadType)) {
                Ref.Init();
                Ref.RunSolverMode();

                Ref.CompareTotalSurface();
                Ref.CompareTotalVolume();


                using (var Test = new TestSetupSingleLevset2D(0.5, quadOrder, cutCellQuadType)) {
                    Test.Init();
                    Test.RunSolverMode();

                    Test.CompareTotalSurface();
                    Test.CompareTotalVolume();
                    Test.CompareSurfaceTo(Ref);
                    Test.CompareVolumeTo(Ref);
                    Test.CompareEdgeAreaTo(Ref);
                    Test.CompareCutLineTo(Ref);
                }

            }

        }


        [Test]
        public static void OneLevelSet_3D(
#if DEBUG
            [Values(3, 6, 7)] int quadOrder,
#else
            [Values(3, 4, 7, 8, 9, 10)] int quadOrder,
#endif
            [Values(XQuadFactoryHelper.MomentFittingVariants.Classic, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants cutCellQuadType
            ) {
            using(var Ref = new TestSetupSingleLevset3D(1.0, quadOrder, cutCellQuadType)) {
                Ref.Init();
                Ref.RunSolverMode();

                Ref.CompareTotalSurface();
                Ref.CompareTotalVolume();

                using (var Test = new TestSetupSingleLevset3D(0.5, quadOrder, cutCellQuadType)) {
                    Test.Init();
                    Test.RunSolverMode();

                    Test.CompareTotalSurface();
                    Test.CompareTotalVolume();
                    Test.CompareSurfaceTo(Ref);
                    Test.CompareVolumeTo(Ref);
                    Test.CompareEdgeAreaTo(Ref);
                    Test.CompareCutLineTo(Ref);
                }

            }

        }

        /// <summary>
        /// The 3D test setting is a prismatic extension of the 2D case, i.e., all properties are constant in z-direction.
        /// Therefore, e.g., the 2D level-set-area in the 3D test could be obtained by multiplying the 1D level-set-length from the 2D test with the mesh with in z-direction.
        /// </summary>
        [Test]
        public static void OneLevelSet_2Dvs3D(
#if DEBUG
            [Values(3, 6, 7)] int quadOrder,
#else
            [Values(3, 4, 7, 8, 9, 10)] int quadOrder,
#endif
            [Values(XQuadFactoryHelper.MomentFittingVariants.Classic, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants cutCellQuadType
            ) {
            using(var Ref = new TestSetupSingleLevset2D(1.0, quadOrder, cutCellQuadType)) {
                Ref.Init();
                Ref.RunSolverMode();

                var surf = Ref.latestCCM.InterfaceArea[Ref.LsTrk.GetSpeciesId("A")].To1DArray().Sum();
                


                using(var Test = new TestSetupSingleLevset3D(1, quadOrder, cutCellQuadType)) {
                    Test.Init();
                    Test.RunSolverMode();

                    var surf3D = Test.latestCCM.InterfaceArea[Ref.LsTrk.GetSpeciesId("A")].To1DArray().Sum();



                    Test.CompareSurfaceTo2D(Ref);
                }

            }

        }



    }
}
