using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public interface ILocatable
    {
        Vector Position { get; set; }
    }

    static class MeshingAlgorithm
    {
        public class Settings
        {
            public Vector[] BoundingBox;

            public Vector[] Boundary;

            public IDictionary<int, int> PeriodicBoundaryMap = null;

            public int NumberOfLloydIterations = 10;

            public int FirstCellNode_indice = 0;
        }

        static void AssertCorrectness<T>(Settings settings, IList<T> nodes)
        {
            Debug.Assert(settings.BoundingBox.Length == 4);
            Debug.Assert(settings.Boundary.Length > 2);
            Debug.Assert(settings.FirstCellNode_indice > -1 && settings.FirstCellNode_indice < nodes.Count);
        }

        public static Mesh<T> ComputeMesh<T>(IList<T> nodes, Settings settings)
            where T : ILocatable, new()
        {
            AssertCorrectness(settings, nodes);
            Mesh<T> mesh = null;
            MeshGenerator<T> voronoiMesher = new MeshGenerator<T>(settings);

            for (int iLloyd = 0; iLloyd <= settings.NumberOfLloydIterations; ++iLloyd)
            {
                mesh = voronoiMesher.Generate(nodes);

                // Lloyds algorithm (Voronoi relaxation)
                if (iLloyd != settings.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(mesh.Cells, ref settings.FirstCellNode_indice);
                }
                nodes = mesh.Nodes;
            }
            return mesh;
        }

        static void MoveNodesTowardsCellCenter<T>(IReadOnlyList<MeshCell<T>> Cells, ref int FirstCellNode_indice)
            where T : ILocatable, new()
        {
            for (int i = 0; i < Cells.Count; ++i)
            {
                MeshCell<T> cell = Cells[i];
                Vector centerOfGravity = CenterOf(cell);
                double relaxValue = 0.1;
                centerOfGravity = centerOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);
                cell.Node.Position = centerOfGravity;

                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = i;
                }
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