using NUnit.Framework;
using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Platform.LinAlg;

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
            int[] NodeSeedsNumbers = { 100, 1000, 10000, 1000000 };
            for(int i = 0; i < NodeSeedsNumbers.Length; ++i)
            {
                int ammountOfNodeSeeds = NodeSeedsNumbers[i];
                VoronoiGrid2D.FromPolygonalDomain(LShapePolygon, 5, ammountOfNodeSeeds);
            }
        }

        Vector[] LShape()
        {
            Vector[] LShapedPolygon = new[] {
                    new Vector(-1,1),
                    new Vector(1,1),
                    new Vector(1,-1),
                    new Vector(0,-1),
                    new Vector(0,0),
                    new Vector(-1,0)
                };
            return LShapedPolygon;
        }
    }
}
