using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using ilPSP;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace VoronoiTests.Grid
{
    class BoundaryTests : BoSSSTestBench
    {
        public override void Run()
        {
            AllPeriodicBoundaries();
        }

        [Test]
        public void SetBoundaryTags()
        {
            byte[] tags = { 1, 1, 1, 1 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(1);
            tagNames.Add(1, "A");

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
            bool tagsAreCorrect = CheckAllBoundaryTagsAre(1, grid.iGridData);
            Assert.IsTrue(tagsAreCorrect);
        }

        static bool CheckAllBoundaryTagsAre(byte tag, IGridData gridData)
        {
            int D = 2;
            double[] x = new double[D];
            MultidimensionalArray GlobalVerticesOut = MultidimensionalArray.CreateWrapper(x, 1, D);
            for (int iEdge = 0; iEdge < gridData.iGeomEdges.Count; ++iEdge)
            {
                if (gridData.iGeomEdges.IsEdgeBoundaryEdge(iEdge))
                {
                    int jCell = gridData.iGeomEdges.CellIndices[iEdge, 0];
                    int iFace = gridData.iGeomEdges.FaceIndices[iEdge, 0];
                    var KRef = gridData.iGeomCells.GetRefElement(jCell);

                    gridData.TransformLocal2Global(KRef.GetFaceCenter(iFace), GlobalVerticesOut, jCell);

                    if (gridData.iGeomEdges.EdgeTags[iEdge] != tag)
                    {
                        Console.WriteLine("Did not find boundary with midpoint: " + x[0] + " " + x[1]);
                        return false;
                    }
                    else
                    {
                        Console.WriteLine("Found boundary with midpoint: " + x[0] + " " + x[1]);
                    }
                }
            }
            return true;
        }

        [Test]
        public void PeriodicBoundaryPair()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            MultidimensionalArray nodes = MultidimensionalArray.Create(6, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5 });
            nodes.SetRow(1, new double[] { -0.8, -0.4 });
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 0.8, 0.4 });
            nodes.SetRow(4, new double[] { 0.9, 0.0 });
            nodes.SetRow(5, new double[] { 0.5, -0.5 });

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
        }

        [Test]
        public void PeriodicBoundaryPairSpeziale()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            MultidimensionalArray nodes = MultidimensionalArray.Create(5, 2);
            nodes.SetRow(0, new double[] { 0.8, 0.5 });
            nodes.SetRow(1, new double[] { 0.8, -0.5 });
            nodes.SetRow(2, new double[] { 0.6, 0 });
            nodes.SetRow(3, new double[] { -0.1, -0.6 });
            nodes.SetRow(4, new double[] { -0.1, 0.6 });


            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
        }

        [Test]
        public void PeriodicBoundaryPairLarge()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(gridBoundary, 0, 500);
        }

        [Test]
        public void AllPeriodicBoundaries()
        {
            byte[] tags = { 182, 181, 182, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 182, "Periodic-Y" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            MultidimensionalArray nodes = MultidimensionalArray.Create(6, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5 });
            nodes.SetRow(1, new double[] { -0.8, -0.4 });
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 0.8, 0.8 });
            nodes.SetRow(4, new double[] { 0.9, 0.2 });
            nodes.SetRow(5, new double[] { 0.5, -0.6 });

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
        }

        [Test]
        public void AllPeriodicBoundariesLarge()
        {
            byte[] tags = { 182, 181, 182, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 182, "Periodic-Y" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(gridBoundary, 0, 500);
        }

        [Test]
        public void LShapePeriodicBoundaries()
        {

            byte[] tags = { 1, 1, 181, 1, 1, 181};
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "Dirichlet" },
                { 181, "Periodic" }
            };

            MultidimensionalArray nodes = MultidimensionalArray.Create(6, 2);
            nodes.SetRow(0, new double[] { -0.5, 0.5 });
            nodes.SetRow(1, new double[] { -0.8, -0.4 });
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 0.8, 0.8 });
            nodes.SetRow(4, new double[] { 0.9, 0.2 });
            nodes.SetRow(5, new double[] { 0.5, -0.6 });

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.LShape(),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 1);
        }

        [Test]
        public void LShapePeriodicBoundariesLarge()
        {

            byte[] tags = { 1, 1, 181, 1, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "Dirichlet" },
                { 181, "Periodic" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.LShape(),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(gridBoundary, 10, 500);
        }

        public void FShape()
        {
            byte[] tags = { 1, 181, 1, 1, 1, 182, 1, 1, 181, 1, 1, 182, 1, 1};
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "Dirichlet" },
                { 181, "Periodic" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.FShape(),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(gridBoundary, 40, 500);
        }
    }
}
