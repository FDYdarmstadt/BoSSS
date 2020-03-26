using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using BoSSS.Solution;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;

namespace AdvancedSolverTests {
    
    [TestFixture]
    class AdvancedSolverMain {

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        public static void Main() {
            SetUp();
            Test();
            TestFixtureTearDown();
        }

        public static void Test() {
            //SubBlockTests.LocalIndexTest(XDGusage.all,2);
            //SubBlockTests.ExternalIndexTest(XDGusage.all, 2);
            //SubBlockTests.MapConsistencyTest(XDGusage.all, 2);
            //SubBlockTests.WriteOutTestMatrices();
            //SubBlockTests.SubMatrixExtractionWithCoupling(XDGusage.all, 2,MatrixShape.diagonal_var_spec);
            SubBlockTests.SplitVectorOperation(XDGusage.all, 2, MatrixShape.diagonal);
        }
    }
}
