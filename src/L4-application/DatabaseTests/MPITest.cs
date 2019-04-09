using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;
using System.Threading;
using NUnit.Framework;
using MPI.Wrappers;

namespace BoSSS.Application.DatabaseTests
{
    [TestFixture]
    class MPITest
    {
        [TestFixtureSetUp]
        public static void InitOnce()
        {
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out bool dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        [TestFixtureTearDown]
        public static void TearDown()
        {
            csMPI.Raw.mpiFinalize();
        }
    }
}
