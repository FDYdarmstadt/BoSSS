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
            BoSSS.Solution.Application.InitMPI();
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
