using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using BoSSS.Platform.LinAlg;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace VoronoiTests.Grid
{
    class GridCreationTests : BoSSSTestBench
    {
        public override void Run()
        {
            PeriodicBoundaryPair();
        }

        [Test]
        public void PeriodicBoundaryPair()
        {
            byte[] tags = { 1, 155, 1, 155 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(1);
            tagNames.Add(1, "A");
            tagNames.Add(155, "B");

            VoronoiBoundary rectangle = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            MultidimensionalArray nodes = MultidimensionalArray.Create(5, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5 });
            nodes.SetRow(1, new double[] { -0.8, -0.4 });
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 0.8, 0.4 });
            nodes.SetRow(4, new double[] { 0.5, -0.5 });

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, rectangle, 10, 0);
        }

        [Test]
        public void FiveNodesInRectangle()
        {
            var rectangle = GridShapes.Rectangle(2, 2);
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
            var LShapePolygon = GridShapes.LShape();
            int[] NodeSeedsNumbers = { 20 };
            for(int i = 0; i < NodeSeedsNumbers.Length; ++i)
            {
                int ammountOfNodeSeeds = NodeSeedsNumbers[i];
                IGrid grid = VoronoiGrid2D.Polygonal(LShapePolygon, 5, ammountOfNodeSeeds);
            }
        }

        [Test]
        public void ToDo()
        {
            MultidimensionalArray nodes1 = MultidimensionalArray.Create(2, 2);
            nodes1.SetRowPt(0, new Vector(-1, 1));
            nodes1.SetRowPt(1, new Vector(1, 1));
            //klappt nicht bei den 2 Zellen, the fuck?!
            VoronoiGrid grid1 = VoronoiGrid2D.Polygonal(nodes1, GridShapes.Rectangle(2,2), 0, 0);
        }
    }
}
