using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;
using System;

namespace BoSSS.Application.LoadBalancingTest {

    /// <summary>
    /// Complete Test for the load balancing.
    /// </summary>
    [TestFixture]
    static public class AllUpTest {

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool MpiInit;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }

       
        /*
        /// <summary>
        /// Da Test!
        /// </summary>
        [Test]
        static public void RuntimeCostDynamicBalanceTest(
            [Values(1, 2)] int DGdegree) {
            LoadBalancingTestMain p = null;
            
            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                null,
                delegate () {
                    p = new LoadBalancingTestMain();
                    p.DynamicBalance = true;
                    p.DEGREE = DGdegree;
                    p.cellCostEstimatorFactory = CellCostEstimatorLibrary.MeasureCostOfWhateverFlorianThinksIsRelevantButIDontUnderstand;
                    return p;
                });
        }
        */

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}
