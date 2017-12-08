using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.RefineAndLoadBal {

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
            Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }

        /// <summary>
        /// Da Test!
        /// </summary>
        [Test]
        static public void LoadBalancing([Values(1, 2)] int DGdegree) {
            RefineAndLoadBalTestMain p = null;
           

            BoSSS.Solution.Application<RefineAndLoadBalControl>._Main(
                new string[0],
                true,
                null,
                delegate() {
                    p = new RefineAndLoadBalTestMain();
                    p.DEGREE = DGdegree;
                    return p;
                });
        }

        /// <summary>
        /// MPI shutdown.
        /// </summary>
        [TestFixtureTearDown]
        public static void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }
    }
}
