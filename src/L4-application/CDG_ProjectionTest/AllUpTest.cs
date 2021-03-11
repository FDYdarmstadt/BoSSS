using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.CDG_ProjectionTest {

    /// <summary>
    /// Unit test for the continuous L2-projection
    /// </summary>
    [TestFixture]
    public class AllUpTest {


        [Test]
        public void AllUp(
            [Values(2)] int dimension,
            [Values(2,3,4)] int degree,
            [Values(2,4,8)] int gridResolution
            ) {

            CDGprojectionMain p = null;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new CDGprojectionMain();
                p.dimension = dimension;
                p.degree = degree;
                p.gridResolution = gridResolution;
                return p;
            });

           Assert.IsTrue(p.passed);

        }

    }



}
