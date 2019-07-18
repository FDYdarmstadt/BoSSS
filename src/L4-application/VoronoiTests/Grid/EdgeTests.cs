using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using BoSSS.Foundation.Grid.Voronoi;

namespace VoronoiTests.Grid
{
    class EdgeTests : BoSSSTestBench
    {
        public override void Run()
        {
            TestRotationVelocity();
        }

        [Test]
        static void TestRotationVelocity()
        {
            double rotation = VoronoiEdge.RotationVelocity(
                new double[] { 1, 0 },
                new double[] { 0, 1 },
                new double[] { -1, 0 },
                new double[] { 0, -1 },
                new double[] { 0, -1 });
            double expectedRotation = 1;
            Assert.IsTrue(Math.Abs(rotation - expectedRotation) < 1e-12,
                "Speed of rotation is not correct.");
        }
    }
}
