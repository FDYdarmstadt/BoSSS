using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Platform.LinAlg;
using ilPSP;
using NUnit.Framework;

namespace VoronoiTests.Grid
{
    class GridCreationTests : BoSSSTestBench
    {
        public override void Run()
        {
            FiveNodesInRectangle();
        }

        [Test]
        public void FiveNodesInRectangle()
        {
            var rectangle = Rectangle(2, 2);
            MultidimensionalArray nodes = MultidimensionalArray.Create(5, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5});
            nodes.SetRow(1, new double[] { -0.8, -0.4});
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 0.8, 0.4 });
            nodes.SetRow(4, new double[] { 0.5, -0.5 });

            IGrid grid = VoronoiGrid2D.Polygonal( nodes,rectangle, 10, 0);
        }

        [Test]
        public void LShapedGrid()
        {
            var LShapePolygon = LShape();
            int[] NodeSeedsNumbers = { 20 };
            for(int i = 0; i < NodeSeedsNumbers.Length; ++i)
            {
                int ammountOfNodeSeeds = NodeSeedsNumbers[i];
                IGrid grid = VoronoiGrid2D.Polygonal(LShapePolygon, 5, ammountOfNodeSeeds);
            }
        }

        static Vector[] LShape()
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

        static Vector[] Rectangle(double width, double height)
        {
            Vector[] polygonBoundary = new Vector[]
            {
                new Vector(-width / 2, height / 2),
                new Vector(width / 2, height / 2),
                new Vector(width / 2, -height / 2),
                new Vector(-width / 2, -height / 2)
            };
            return polygonBoundary;
        }

        [Test]
        public void ToDo()
        {
            MultidimensionalArray nodes1 = MultidimensionalArray.Create(2, 2);
            nodes1.SetRowPt(0, new Vector(-1, 1));
            nodes1.SetRowPt(1, new Vector(1, 1));
            //klappt nicht bei den 2 Zellen, the fuck?!
            VoronoiGrid grid1 = VoronoiGrid2D.Polygonal(nodes1, Rectangle(2,2), 0, 0);
        }
    }
}
