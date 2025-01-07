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
            [Values(3, 4, 7, 8, 9, 10)] int quadOrder,
            [Values(CutCellQuadratureMethod.Classic, CutCellQuadratureMethod.OneStepGauss, CutCellQuadratureMethod.OneStepGaussAndStokes, CutCellQuadratureMethod.Saye, CutCellQuadratureMethod.Algoim)] CutCellQuadratureMethod cutCellQuadType
            ) {
            //BoSSS.Application.CutCellQuadratureScaling.AllTests.OneLevelSet_2D


            using(var Ref = new TestSetupSingleLevset2D(1.0, quadOrder)) {
                Ref.Init();
                Ref.RunSolverMode();

                using(var Test = new TestSetupSingleLevset2D(0.5, quadOrder)) {
                    Test.Init();
                    Test.RunSolverMode();

                    Test.CompareSurfaceTo(Ref);
                    Test.CompareVolumeTo(Ref);
                    Test.CompareEdgeAreaTo(Ref);
                }

            }

        }


        [Test]
        public static void OneLevelSet_3D(
            [Values(CutCellQuadratureMethod.Classic, CutCellQuadratureMethod.Saye, CutCellQuadratureMethod.Algoim)] CutCellQuadratureMethod cutCellQuadType
            ) {


        }

        /// <summary>
        /// The 3D test setting is a prismatic extension of the 2D case, i.e., all properties are constant in z-direction.
        /// Therefore, e.g., the 2D level-set-area in the 3D test could be obtained by multiplying the 1D level-set-length from the 2D test with the mesh with in z-direction.
        /// </summary>
        [Test]
        public static void OneLevelSet_2Dvs3D() {


        }



    }
}
