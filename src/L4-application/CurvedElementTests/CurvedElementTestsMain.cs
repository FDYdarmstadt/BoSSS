using System;
using BoSSS.Solution;

namespace BoSSS.Application.CurvedElementTests {

    /// <summary>
    /// Entry point for running the tests in this assembly without the test runner, e.g. for debugging.
    /// For CI, the NUnit tests in <see cref="TaylorCouetteTests"/> are discovered by the <c>PublicTestRunner</c>.
    /// </summary>
    public static class CurvedElementTestsMain {

        static void Main(string[] args) {
            BoSSS.Solution.Application.InitMPI(args);

            // call the test(s) you want to debug here:
            foreach(int resolution in new[] { 1, 2 }) {
                for(int degree = 1; degree <= 3; degree++) {
                    TaylorCouetteTests.SteadyStokesOnPeriodicWedge(degree, resolution);
                    TaylorCouetteTests.SteadyNavierStokesOnPeriodicWedge(degree, resolution);
                }
            }

            BoSSS.Solution.Application.FinalizeMPI();
        }
    }
}
