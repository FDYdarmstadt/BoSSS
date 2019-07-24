using BoSSS.Solution;
using BoSSS.Solution.Control;
using System.Threading;

namespace VoronoiTests
{
    /// <summary>
    ///  Base class for tests that run an BoSSS application.
    /// </summary>
    public class BoSSSTestBench : TestBench
    {
        /// <summary>
        /// Performs bootstrapping.
        /// </summary>
        public override void SetUp()
        {
            ilPSP.Environment.Bootstrap(
                new string[0],
                Application.GetBoSSSInstallDir(),
                out bool dummy);
            Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CurrentCulture;
        }

        protected static void RunApplication(IApplication app, AppControl ctrl)
        {
            app.Init(ctrl);
            app.RunSolverMode();
            app.Dispose();
        }

        public override void Run(){}

        public override void TearDown(){}
    }
}
