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
using AdvancedSolverTests.SubBlocking;

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
            System.Threading.Thread.Sleep(10000);
            //SubBlockTests.LocalIndexTest(XDGusage.all,2);
            //SubBlockTests.ExternalIndexTest(XDGusage.all, 2);
            //SubBlockTests.MapConsistencyTest(XDGusage.all, 2);
            //SubBlockTests.WriteOutTestMatrices();
            //SubBlockTests.SubMatrixExtractionWithCoupling(XDGusage.all, 2,MatrixShape.diagonal_var_spec);
            //SubBlockTests.SplitVectorOperations(XDGusage.none, 2, MatrixShape.diagonal_var);
            //SubBlockTests.SubSelection(SelectionType.species);
            //LocalTests.CellwiseSubSelection(SelectionType.species);
            //ExternalTests.ExternalIndexTest(XDGusage.all,2);
            //ExternalTests.ExternalIndexTest(XDGusage.all, 2,4);
            //ExternalTests.GetExternalRowsTest(XDGusage.all, 2,4);
            //ExternalTests.FastSubMatrixExtraction(XDGusage.all, 2,MatrixShape.laplace,4);
            //ExternalTests.SubMatrixIgnoreCoupling(XDGusage.all, 2, MatrixShape.diagonal_var_spec,4);
            //ExternalTests.SubMatrixExtraction(XDGusage.all, 2, MatrixShape.full_var_spec,4);
            //ExternalTests.SubBlockExtraction(XDGusage.all, 2, MatrixShape.laplace, 4);
            ExternalTests.VectorCellwiseOperation(XDGusage.all, 2, MatrixShape.diagonal_var_spec, 4);
        }
    }
}
