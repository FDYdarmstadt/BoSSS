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

            Console.WriteLine("TOTAL Time spend in matrix operations: " + " sec.");


        }

        [Test]
        public static void Test() {
            //SubBlockTests.ExtractDiagonalBlocks(XDGusage.all, 1,MatrixShape.diagonal);
            //SubBlockTests.ExtractDiagonalBlocks(XDGusage.all, 1, MatrixShape.diagonal_var);
            //SubBlockTests.ExtractDiagonalBlocks(XDGusage.none, 1, MatrixShape.diagnoal_var_spec);
            SubBlockTests.ExtractDiagonalBlocks(XDGusage.all, 1, MatrixShape.diagnoal_spec);
        }

    }
}
