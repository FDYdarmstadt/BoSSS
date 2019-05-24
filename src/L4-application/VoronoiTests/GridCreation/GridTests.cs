using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Platform.LinAlg;
using NUnit.Framework;

namespace VoronoiTests.GridCreation
{
    class GridTests : BoSSSTestBench
    {
        public override void Run()
        {
            LShapedGrid();
        }

        [Test]
        public void LShapedGrid()
        {
            var LShapePolygon = LShape();
            int[] NodeSeedsNumbers = { 200 };
            for(int i = 0; i < NodeSeedsNumbers.Length; ++i)
            {
                int ammountOfNodeSeeds = NodeSeedsNumbers[i];
                IGrid grid = VoronoiGrid2D.Polygonal(LShapePolygon, 5, ammountOfNodeSeeds);
            }
        }

        Vector[] LShape()
        {
            double a = 1;
            Vector[] LShapedPolygon = new[] {
                    new Vector(-a,a),
                    new Vector(a,a),
                    new Vector(a,-a),
                    new Vector(0,-a),
                    new Vector(0,0),
                    new Vector(-a,0)
                };
            return LShapedPolygon;
        }
    }
}
