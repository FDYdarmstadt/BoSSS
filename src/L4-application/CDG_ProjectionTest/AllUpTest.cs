using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.CDG_ProjectionTest {

    /// <summary>
    /// Unit test for the continuous L2-projection
    /// </summary>
    [TestFixture]
    public class AllUpTest {

        /// <summary>
        /// MPI finalization
        /// </summary>
        [TestFixtureTearDown]
        public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public void TestFixtureSetUp() {
            bool MpiInit;
            Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }


        [Test]
        public void AllUp() {

            CDGprojectionMain p = null;

            System.Threading.Thread.Sleep(10000);

            BoSSS.Solution.Application._Main(new string[0], true, null, delegate () {
                p = new CDGprojectionMain();
                return p;
            });

           Assert.IsTrue(p.passed);

        }
    }



}
