using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;
using System;


namespace AdvancedSolverTests {
    [TestFixture]
    public abstract class TestBench {
        [OneTimeSetUp]
        public virtual void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }

        public abstract void Run();

        /// <summary>
        /// 
        /// </summary>
        [OneTimeTearDown]
        public virtual void TearDown() {
        }


        

    }
}
