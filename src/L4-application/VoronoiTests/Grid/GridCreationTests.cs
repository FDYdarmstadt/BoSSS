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
            Large();
        }

        [Test]
        public void Large()
        {
            var rectangle = GridShapes.Rectangle(2, 2);

            IGrid grid = VoronoiGrid2D.Polygonal(rectangle, 0, 10000);
            Plotter.Plot(grid);
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
        public void EquidistandGrid()
        {
            var rectangle = GridShapes.Rectangle(2, 2);
            MultidimensionalArray nodes = MultidimensionalArray.Create(4, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5 });
            nodes.SetRow(1, new double[] { -0.5, -0.5 });
            nodes.SetRow(2, new double[] { 0.5, 0.5 });
            nodes.SetRow(3, new double[] { 0.5, -0.5 });

            IGrid grid = VoronoiGrid2D.Polygonal(nodes, rectangle, 0, 0);
        }

        [Test]
        public void EquidistandGridLarge()
        {
            int nodesPerDimension = 100;
            var rectangle = GridShapes.Rectangle(2, 2);
            
            MultidimensionalArray nodes = MultidimensionalArray.Create(nodesPerDimension * nodesPerDimension, 2);
            for(int i = 0; i < nodesPerDimension; ++i)
            {
                for(int j = 0; j < nodesPerDimension; ++j)
                {
                    nodes[i * nodesPerDimension +j, 0] = 2.0/ (nodesPerDimension - 1) * i -1.0 ;
                    nodes[i * nodesPerDimension +j, 1] = 2.0 / (nodesPerDimension - 1) * j -1.0;
                }
            }
            IGrid grid = VoronoiGrid2D.Polygonal(nodes, rectangle, 0, 0);
            //Plotter.Plot(grid);
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
    }
}
