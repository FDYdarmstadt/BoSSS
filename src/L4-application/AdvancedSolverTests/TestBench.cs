using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;
using System;

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

        /// <summary>
        /// 
        /// </summary>
        [TestFixtureTearDown]
        public virtual void TearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}
