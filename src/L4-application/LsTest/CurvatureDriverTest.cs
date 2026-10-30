using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Linq;

namespace BoSSS.Application.LsTest {

    /// <summary>
    /// Tests of <see cref="CurvatureAlgorithms.CurvatureDriver"/>.
    /// </summary>
    [TestFixture]
    public static class CurvatureDriverTest {

        /// <summary>
        /// Curvature of a circle of radius R (exact: 1/R) from <see cref="CurvatureAlgorithms.CurvatureDriver"/> in a Laplace-Beltrami mode,
        /// i.e. in the branch where the driver computes the Hessian itself, with the recommended filter and with / without filtering the Hessian.
        /// There, both out parameters of <c>ComputeHessian</c> used to be the same variable, so that the curvature was zero with a filtered Hessian.
        /// </summary>
        [Test]
        public static void CurvatureOfCircle_LaplaceBeltramiMode([Values(false, true)] bool useFiltLevSetHess) {
            const double R = 0.6;
            const int deg = 3;
            var Nodes = GenericBlas.Linspace(-1, 1, 17);
            var gdat = Grid2D.Cartesian2DGrid(Nodes, Nodes).GridData;
            var phi = new LevelSet(new Basis(gdat, deg), "Phi");
            phi.ProjectField(X => X[0] * X[0] + X[1] * X[1] - R * R);
            var trk = new LevelSetTracker(gdat, CutCellQuadratureMethod.Saye, 1, new[] { "A", "B" }, phi);
            trk.UpdateTracker(0.0);

            var config = CurvatureAlgorithms.FilterConfiguration.Default;
            config.useFiltLevSetHess = useFiltLevSetHess;
            var curvature = new SinglePhaseField(new Basis(gdat, deg), "Curvature");
            CurvatureAlgorithms.CurvatureDriver(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, config, curvature, out _, trk, 2 * deg + 1, phi);

            // curvature on the circle, at points away from the grid lines
            const int N = 24;
            double[] kappa = Enumerable.Range(0, N).Select(i => {
                double t = 2 * Math.PI * (i + 0.37) / N;
                return curvature.ProbeAt(R * Math.Cos(t), R * Math.Sin(t));
            }).ToArray();
            double meanErr = Math.Abs(kappa.Average() * R - 1.0), maxErr = kappa.Max(k => Math.Abs(k * R - 1.0));
            Console.WriteLine($"useFiltLevSetHess = {useFiltLevSetHess}: curvature on the circle, mean {kappa.Average()}, relative error of the mean {meanErr}, maximum relative error {maxErr} (exact {1.0 / R})");
            // the quadratic level set is exactly represented, the error is about 1e-4 (with the bug: 1, i.e. zero curvature)
            Assert.Less(meanErr, 0.01, "mean curvature on the circle");
            Assert.Less(maxErr, 0.05, "pointwise curvature on the circle");
        }
    }
}
