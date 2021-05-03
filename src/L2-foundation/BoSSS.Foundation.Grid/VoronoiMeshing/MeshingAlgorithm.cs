using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class MeshingAlgorithm
    {
        public class State
        {
            public Vector[] BoundingBox;

            public Vector[] Boundary;

            public PeriodicMap PeriodicMap = null;

            public int NumberOfLloydIterations = 10;

            public int FirstCellNodeIndice = 0;
        }

        static void AssertCorrectness<T>(State settings, IList<T> nodes)
        {
            Debug.Assert(settings.BoundingBox.Length == 4);
            Debug.Assert(settings.Boundary.Length > 2);
            Debug.Assert(settings.FirstCellNodeIndice > -1 && settings.FirstCellNodeIndice < nodes.Count);
        }

        public static Domain<T> ComputeMesh<T>(IList<T> nodes, State settings)
            where T : ICloneable<T>, new()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            AssertCorrectness(settings, nodes);
            MeshGenerator<T> voronoiMesher = new MeshGenerator<T>(settings);
            Domain<T> domain = voronoiMesher.Generate(nodes, settings.FirstCellNodeIndice);
            for (int iLloyd = 0; iLloyd < settings.NumberOfLloydIterations; ++iLloyd)
            {
                if (iLloyd != settings.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(domain.Cells);
                }
                domain = voronoiMesher.Generate(domain.Nodes, domain.Boundary.FirstCorner);
            }
            stopwatch.Stop();
            Console.WriteLine(stopwatch.ElapsedMilliseconds);

            //MatlabPlotter.Plot(domain.Mesh);
            Debug.Assert(InnerEdgesAlign(domain.Mesh));
            
            return domain;
        }

        static void MoveNodesTowardsCellCenter<T>(IReadOnlyList<MeshCell<T>> cells)
            where T : ILocatable
        {
            for (int i = 0; i < cells.Count; ++i)
            {
                MeshCell<T> cell = cells[i];
                Vector centerOfGravity = CenterOf(cell);
                double relaxValue = 0.3;
                centerOfGravity = centerOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);
                cell.Node.Position = centerOfGravity;
            }
        }

        static Vector VertexMean<T>(MeshCell<T> cell)
        {
            Vector centerOfGravity = new Vector(0, 0);
            foreach (Vertex vertex in cell.Vertices)
            {
                centerOfGravity += vertex.Position;
            }
            centerOfGravity.ScaleInPlace(1.0 / cell.Vertices.Length);
            return centerOfGravity;
        }

        static Vector CenterOf<T>(MeshCell<T> cell)
        {
            Vector centerOfGravity = new Vector(0, 0);
            Vector root = cell.Vertices[0].Position;
            double cellArea = 0;

            for (int i = 1; i < cell.Vertices.Length - 1; ++i)
            {
                Vector a = cell.Vertices[i].Position;
                Vector b = cell.Vertices[i + 1].Position;
                double area = AreaOfTriangle(root, a, b);
                centerOfGravity += (root + a + b) / 3 * area;
                cellArea += area;
            }
            centerOfGravity.ScaleInPlace(1.0 / cellArea);
            return centerOfGravity;
        }

        static double AreaOfTriangle(Vector a, Vector b, Vector c)
        {
            double area = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
            return area / 2;
        }

        static bool InnerEdgesAlign<T>(Mesh<T> mesh)
        {
            foreach(MeshCell<T> cell in mesh.Cells)
            {
                for(int i = 0; i < cell.Edges.Length; ++i)
                {
                    Edge<T> edge = cell.Edges[i];
                    if (!edge.IsBoundary)
                    {
                        if((edge.Start.ID != edge.Twin.End.ID) || (edge.End.ID != edge.Twin.Start.ID))
                        {
                            return false;
                        }
                        if(edge.Start.ID == edge.End.ID)
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
    }
}