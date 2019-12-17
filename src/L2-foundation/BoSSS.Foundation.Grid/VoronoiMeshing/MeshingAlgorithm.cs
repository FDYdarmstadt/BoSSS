using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class MeshingAlgorithm
    {
        public class Settings
        {
            public Vector[] BoundingBox;

            public Vector[] Boundary;

            public PeriodicMap PeriodicMap = null;

            public int NumberOfLloydIterations = 10;

            public int FirstCellNodeIndice = 0;
        }

        static void AssertCorrectness<T>(Settings settings, IList<T> nodes)
        {
            Debug.Assert(settings.BoundingBox.Length == 4);
            Debug.Assert(settings.Boundary.Length > 2);
            Debug.Assert(settings.FirstCellNodeIndice > -1 && settings.FirstCellNodeIndice < nodes.Count);
        }

        public static IMesh<T> ComputeMesh<T>(IList<T> nodes, Settings settings)
            where T : ICloneable<T>, new()
        {
            AssertCorrectness(settings, nodes);
            MeshGenerator<T> voronoiMesher = new MeshGenerator<T>(settings);
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
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
            return domain;
        }

        static void MoveNodesTowardsCellCenter<T>(IReadOnlyList<MeshCell<T>> cells)
            where T : ILocatable
        {
            for (int i = 0; i < cells.Count; ++i)
            {
                MeshCell<T> cell = cells[i];
                Vector centerOfGravity = CenterOf(cell);
                double relaxValue = 0.1;
                centerOfGravity = centerOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);
                cell.Node.Position = centerOfGravity;
            }
        }

        static Vector CenterOf<T>(MeshCell<T> cell)
        {
            Vector centerOfGravity = new Vector(0, 0);
            foreach (Vertex vertex in cell.Vertices)
            {
                centerOfGravity += vertex.Position;
            }
            centerOfGravity.Scale(1.0 / cell.Vertices.Length);
            return centerOfGravity;
        }
    }
}