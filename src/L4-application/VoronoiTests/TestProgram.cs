using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using BoSSS.Solution;
using NUnit.Framework;

namespace VoronoiTests
{
    /// <summary>
    /// Abstract base class for tests.
    /// </summary>
    [TestFixture]
    public abstract class TestProgram
    {
        /// <summary>
        /// Performs bootstrapping.
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp()
        {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.CurrentCulture;
        }
    }
}
