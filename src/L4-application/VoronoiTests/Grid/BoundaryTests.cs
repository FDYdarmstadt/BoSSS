using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Voronoi;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;

namespace VoronoiTests.Grid
{
    class BoundaryTests : BoSSSTestBench
    {
        public override void Run()
        {
            PeriodicPairSkewCheckerBoard();
        }

        [Test]
        public void Shift()
        {
            double offset = 1e-7;
            MultidimensionalArray nodes = MultidimensionalArray.Create(6, 2);
            nodes.SetRowPt(0, new Vector(0.5 + offset, 0.5));
            nodes.SetRowPt(1, new Vector(-0.5 + offset, 0.5));
            nodes.SetRowPt(2, new Vector(0.5 - offset, -0.5));
            nodes.SetRowPt(3, new Vector(-0.5 - offset, -0.5));
            nodes.SetRowPt(4, new Vector(0, -0.5));
            nodes.SetRowPt(5, new Vector(0, 0.5));

            // ### Grid
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "AdiabaticSlipWall"},
                { 181, "Periodic-X"}
            };
            var gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
        }

        [Test]
        public void PeriodicPairCheckerBoard()
        {
            // ### Grid
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 1, "AdiabaticSlipWall"},
                { 181, "Periodic-X"}
            };
            var gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(1, 0.5 ),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            double[] xTics = GenericBlas.Linspace(-0.5, 0.5, 3);
            double[] yTics = GenericBlas.Linspace(-0.25, 0.25, 5);
            for (int i = 0; i < yTics.Length; ++i)
            {
                yTics[i] *= -1;
            }
            MultidimensionalArray nodes = Checkerize(xTics, yTics);

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
            Plotter.Plot(grid);

        }

        MultidimensionalArray Checkerize(double[] a, double[] b)
        {
            MultidimensionalArray abT = MultidimensionalArray.Create((a.Length - 1) * (b.Length - 1), 2);
            for (int i = 0; i < a.Length - 1; ++i)
            {
                for (int j = 0; j < b.Length - 1; ++j)
                {
                    abT[i * (b.Length - 1) + j, 0] = (a[i] + a[i + 1]) / 2.0;
                    abT[i * (b.Length - 1) + j, 1] = (b[j] + b[j + 1]) / 2.0;
                }
            }
            return abT;
        }

        [Test] 
        public void PeriodicPairSkewCheckerBoard()
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
            nodes.SetRow(0, new double[] { -0.71, 0.7 });
            nodes.SetRow(1, new double[] { 0.8, 0.7 });
            nodes.SetRow(2, new double[] { 0, 0.5 });
            nodes.SetRow(3, new double[] { -0.7, -0.7 });
            nodes.SetRow(4, new double[] { 0.7, -0.7 });
            nodes.SetRow(5, new double[] { 0, -0.5 });

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
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
        public void PeriodicBoundaryPairBoundaryOnEdge()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "Dirichlet" }
            };
            MultidimensionalArray nodes = MultidimensionalArray.Create(6, 2);
            nodes.SetRowPt(0, new Vector(-0.8, 0.6));
            nodes.SetRowPt(1, new Vector(-0.8, -0.6));
            nodes.SetRowPt(2, new Vector(-0.2, 0.0));
            nodes.SetRowPt(3, new Vector(0.2, 0.0));
            nodes.SetRowPt(4, new Vector(0.8, 0.6));
            nodes.SetRowPt(5, new Vector(0.8, -0.6));
            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2, 2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };
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
                Polygon = GridShapes.Rectangle(8, 8),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };

            Random random = new Random(10);
            MultidimensionalArray nodes = default(MultidimensionalArray);
            for (int i = 0; i < 10; i += 1)
            {
                Console.WriteLine($"Roll number{i}");
                nodes = RandomNodesInSquare(4.09, 4.0, 300, random);
                nodes[0, 0] = -1 + 1e-5;
                nodes[0, 1] = 1 - 1e-5;
                try
                {
                    VoronoiGrid garid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                }
            }
            //VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
            //Plotter.Plot(grid);
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
        public void AllPeriodicBoundariesWithOustideStartNodes()
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
            nodes.SetRow(4, new double[] { 1.2, 0.2 });
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
            Random random = new Random(1);
            MultidimensionalArray nodes = default(MultidimensionalArray);
            for (int i = 0; i < 48; ++i)
            {
                Console.WriteLine($"Roll number {i}");
                nodes = RandomNodesInSquare(1.0, 1.0, 300, random);
                nodes[0, 0] = -1 + 1e-6;
                nodes[0, 1] = 1 - 1e-6;
                try
                {
                    //VoronoiGrid garid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
                    //Plotter.Plot(grid);
                }
                catch (Exception e)
                {
                    Console.WriteLine(e);
                }
            }
            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
        }

        MultidimensionalArray RandomNodesInSquare(double height, double width, int number, Random random = null)
        {
            if(random == null)
            {
                random = new Random();
            }
            MultidimensionalArray nodes = MultidimensionalArray.Create(number, 2);
            
            for (int i = 0; i < number; ++i)
            {
                nodes[i, 0] = width * 2 * (random.NextDouble() -0.5);  
                nodes[i, 1] = height * 2 * (random.NextDouble() - 0.5);
            }
            return nodes;
        }

        [Test]
        public void AllPeriodicBoundariesOutside()
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
            nodes.SetRow(0, new double[] { -1.1, 1.1 });
            nodes.SetRow(1, new double[] { 0.6, 1.4 });
            nodes.SetRow(2, new double[] { 0, 0 });
            nodes.SetRow(3, new double[] { 1, -0.9 });
            nodes.SetRow(4, new double[] { -0.9, -1.2 });
            nodes.SetRow(5, new double[] { -1.1, 0.2 });

            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 10, 0);
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

        [Test]
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

        [Test]
        public void Remap()
        {
            byte[] tags = { 1, 181, 1, 181 };
            SortedList<byte, string> tagNames = new SortedList<byte, string>(2)
            {
                { 181, "Periodic-X" },
                { 1, "bondary" }
            };

            VoronoiBoundary gridBoundary = new VoronoiBoundary
            {
                Polygon = GridShapes.Rectangle(2,2),
                EdgeTags = tags,
                EdgeTagNames = tagNames
            };
            double[] positions = new double[] 
            {
                -0.64874644688322713,
                0.83818004313993111,
                -0.39475947428553138,
                0.23663302374998896,
                0.58922918492853482,
                0.83854511946848365,
                -0.811382461267156,
                -0.4159610860057516,
                -0.19666215667077264,
                -0.24376388607043981,
                0.3385324063754323,
                0.086134041417832763,
                0.80498108279434089,
                0.22350558445791927,
                -0.68131747598283521,
                -0.87257764806623139,
                0.48863086193005234,
                -0.51183362983054159,
                0.99309783411349173,
                -0.10141430352239808,
            };
            MultidimensionalArray nodes = MultidimensionalArray.CreateWrapper(positions, 10, 2);
            for(int i = 0; i < 10; ++i) 
            {
                nodes[i, 0] += 0.01;
            }
            
            VoronoiGrid grid = VoronoiGrid2D.Polygonal(nodes, gridBoundary, 0, 0);
            Plotter.Plot(grid);
        }
    }
}
