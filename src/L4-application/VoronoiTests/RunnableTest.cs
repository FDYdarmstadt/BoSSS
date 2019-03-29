using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using BoSSS.Solution;
using NUnit.Framework;
using BoSSS.Solution.Control;

namespace VoronoiTests
{
    /// <summary>
    ///  Base class for tests that run an BoSSS application.
    /// </summary>
    [TestFixture]
    class RunnableTest
    {
        /// <summary>
        /// Performs bootstrapping.
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp()
        {
            ilPSP.Environment.Bootstrap(
                new string[0],
                Application.GetBoSSSInstallDir(),
                out bool dummy);
            Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CurrentCulture;
        }

        public static void RunApplication(IApplication app, AppControl ctrl)
        {
            app.Init(ctrl);
            app.RunSolverMode();
        }
    }
}
