using System;
using System.Globalization;
using System.Threading;
using NUnit.Framework;
using MPI.Wrappers;
using System.Diagnostics;

namespace BoSSS.Application.DatabaseTests
{
    [TestFixture]
    class MPITest
    {
        [OneTimeSetUp]
        public static void InitOnce()
        {
            BoSSS.Solution.Application.InitMPI();
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            //  removed MPI shutdown, this causes the test to crash without result (for some reason, this method is called multiple times)
            //  tests seem to work without MPI shutdown anyway
            csMPI.Raw.mpiFinalize();
        }
    }
}
