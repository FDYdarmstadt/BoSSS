using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// collection of tests
    /// </summary>
    [TestFixture]
    public static class VariousTests {


        /// <summary>
        /// Tests that multiple calls to <see cref="LevelSetTracker.UpdateTracker"/> do no harm.
        /// </summary>
        [Test]
        public static void MultipleTrackerUpdateCalls([Values(0, 1, 2)] int NearRegionWidht) {

            var Nodes = GenericBlas.Linspace(-2, 2, 20);
            GridCommons grd = Grid2D.Cartesian2DGrid(Nodes, Nodes);
            var gdat = grd.GridData;
            var phi = new LevelSet(new Basis(gdat, 2), "Phi");
            var trk = new LevelSetTracker(gdat, XQuadFactoryHelper.MomentFittingVariants.Saye, NearRegionWidht, new[] { "A", "B" }, phi);

            Func<double[], double> UA = X => X[0] * X[1];
            Func<double[], double> UB = X => X[0].Pow2();

            XDGField u = new XDGField(new XDGBasis(trk, 2), "u");
            // check that no one changes the default value accidentally or intentionally.
            // If you really want to change this, talk to Florian.
            Assert.AreEqual(u.UpdateBehaviour, BehaveUnder_LevSetMoovement.PreserveMemory, "Default update behavior has changed");
            u.GetSpeciesShadowField("A").ProjectField(UA);
            u.GetSpeciesShadowField("B").ProjectField(UB);

            void CheckError() {
                var Amask = trk.Regions.GetSpeciesMask("A");
                var Bmask = trk.Regions.GetSpeciesMask("B");
                double uAerr = u.GetSpeciesShadowField("A").L2Error(UA.Vectorize(), new CellQuadratureScheme(true, Amask));
                double uBerr = u.GetSpeciesShadowField("B").L2Error(UA.Vectorize(), new CellQuadratureScheme(true, Bmask));

                Assert.LessOrEqual(uAerr, 1.0e-10, "error for A to high");
                Assert.LessOrEqual(uAerr, 1.0e-10, "error for B to high");
            }

            for(int i = 0; i < 3; i++) {

                // check before update
                CheckError();

                // after one update
                trk.UpdateTracker(0.0);
                CheckError();

                // after more updates
                trk.UpdateTracker(0.0);
                trk.UpdateTracker(0.0);
                CheckError();

                // Push the stacks
                trk.PushStacks();
                trk.IncreaseHistoryLength(i+1);
            }
        }


    }
}
