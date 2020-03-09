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
            OnlyEdgesOnBoundary();
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
        public void CornerSpeziale()
        {
            MultidimensionalArray nodes1 = MultidimensionalArray.Create(3, 2);
            nodes1.SetRowPt(0, new Vector(-1, 1));
            nodes1.SetRowPt(1, new Vector(1, 1));
            nodes1.SetRowPt(2, new Vector(-1, -1));
            VoronoiGrid grid1 = VoronoiGrid2D.Polygonal(nodes1, GridShapes.Rectangle(2,2), 0, 0);
        }

        [Test]
        public void SomeEdgesOnBoundary()
        {
            MultidimensionalArray nodes1 = MultidimensionalArray.Create(4, 2);
            nodes1.SetRowPt(0, new Vector(0, 2));
            nodes1.SetRowPt(1, new Vector(0, 0));
            nodes1.SetRowPt(2, new Vector(-1, 0));
            nodes1.SetRowPt(3, new Vector(-1, -1));
            VoronoiGrid grid1 = VoronoiGrid2D.Polygonal(nodes1, GridShapes.Rectangle(2, 2), 0, 0);
        }

        [Test]
        public void OnlyEdgesOnBoundary()
        {
            MultidimensionalArray nodes1 = MultidimensionalArray.Create(8, 2);
            nodes1.SetRowPt(0, new Vector(0, 1.5));
            nodes1.SetRowPt(1, new Vector(-0.5, 0));
            nodes1.SetRowPt(2, new Vector(0, -0.5));
            nodes1.SetRowPt(3, new Vector(0.5, 0));
            nodes1.SetRowPt(4, new Vector(0, 0.5));
            nodes1.SetRowPt(5, new Vector(1.5, 0));
            nodes1.SetRowPt(6, new Vector(-1.5, 0));
            nodes1.SetRowPt(7, new Vector(0, -1.5));
            VoronoiGrid grid1 = VoronoiGrid2D.Polygonal(nodes1, GridShapes.Rectangle(2, 2), 0, 1);
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
    }
}
