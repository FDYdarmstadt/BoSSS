using NUnit.Framework;
using MPI.Wrappers;

namespace AdvancedSolverTests {
    [TestFixture]
    public abstract class TestBench {
        [TestFixtureSetUp]
        public virtual void SetUp() {
            bool MpiInit;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }

        public abstract void Run();

        [TestFixtureTearDown]
        public virtual void TearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}
