using ilPSP.Utils;
using PublicTestRunner;
using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using ilPSP;
using BoSSS.Application.LsTest;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Application.BoSSSpad;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.IO;
using System.Threading;
using BoSSS.Solution.Gnuplot;
using BoSSS.Foundation.Quadrature;

namespace ValidationTestRunner {
    /// <summary>
    /// Partial class, managing long(er) running tests for the various level set algorithms.
    /// Maintainer: Matthias Rieckmann
    /// </summary>
    /// <remarks>
    /// All these tests here are intended to be run at the local MS windows HPC cluster (aka. FDYcluster) at Chair of Fluid Dynamics (FDY)
    /// </remarks>
    [TestFixture]
    [NUnitNumThreads(1)]
    static public partial class WorksheetTests_Local_LevelSet {

        // each setup as a single test, so when deploying via 'runjobmanager' all cases are executed simultaneuos
        #region Swirling Flow   
        
        [Test]
        [NUnitFileToCopyHack("levelset/SwirlingFlow_TemporalConvergence_Run.ipynb", "levelset/SwirlingFlow_TemporalConvergence_Evaluate.ipynb")]
        public static void SwirlingFlow_TemporalConvergence() {
            string really = System.Environment.GetEnvironmentVariable("RUN_LEVELSET");
            if(really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping level-set validation ");
                return;
            } else {
                Console.WriteLine("RUN_LEVELSET = " + really);
            }
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
               "SwirlingFlow_TemporalConvergence",
               "SwirlingFlow_TemporalConvergence*",
               "delete_SwirlingFlow_TemporalConvergence", new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("levelset/SwirlingFlow_TemporalConvergence_Run.ipynb");
            ValidationTestRunnerMain.RunWorksheet("levelset/SwirlingFlow_TemporalConvergence_Evaluate.ipynb");

            Console.WriteLine("SwirlingFlow - Temporal Convergence @ FDYcluster");

        }

        [Test]
        [NUnitFileToCopyHack("levelset/SwirlingFlow_SpatialConvergence_Run.ipynb", "levelset/SwirlingFlow_SpatialConvergence_Evaluate.ipynb")]
        public static void SwirlingFlow_SpatialConvergence() {
            string really = System.Environment.GetEnvironmentVariable("RUN_LEVELSET");
            if(really.IsEmptyOrWhite()) {
                Console.WriteLine("skipping level-set validation ");
                return;
            } else {
                Console.WriteLine("RUN_LEVELSET = " + really);
            }
            ValidationTestRunnerMain.DeleteDatabaseAndDeploymentsWhenOld(
               "SwirlingFlow_SpatialConvergence",
               "SwirlingFlow_SpatialConvergence*",
               "delete_SwirlingFlow_SpatialConvergence", new TimeSpan(days: 25, hours: 1, minutes: 0, seconds: 1));

            ValidationTestRunnerMain.RunWorksheet("levelset/SwirlingFlow_SpatialConvergence_Run.ipynb");
            ValidationTestRunnerMain.RunWorksheet("levelset/SwirlingFlow_SpatialConvergence_Evaluate.ipynb");

            Console.WriteLine("SwirlingFlow - Spatial Convergence @ FDYcluster");

        }       

        #endregion
    }
}
