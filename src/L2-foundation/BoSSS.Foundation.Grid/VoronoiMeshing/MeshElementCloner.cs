using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class MeshElementCloner
    {
        public static MeshCell<T>[] Clone<T>(IList<MeshCell<T>> cells)
            where T : ILocatable, new()
        {
            MeshCell<T>[] clones = new MeshCell<T>[cells.Count];
            for (int i = 0; i < cells.Count; ++i)
            {
                clones[i] = Clone(cells[i]);
            }
            return clones;
        }

        public static List<MeshCell<T>> Clone<T>(List<MeshCell<T>> cells)
            where T : ILocatable, new()
        {
            List<MeshCell<T>> clones = new List<MeshCell<T>>(cells.Count);
            for (int i = 0; i < cells.Count; ++i)
            {
                clones.Add(Clone(cells[i]));
            }
            return clones;
        }

        public static MeshCell<T> Clone<T>(MeshCell<T> cell)
            where T : ILocatable, new()
        {
            MeshCell<T> clone = new MeshCell<T>()
            {
                Node = Clone(cell.Node),
                Type = cell.Type,
                ID = cell.ID,
                Vertices = Clone(cell.Vertices),
                Edges = Clone(cell.Edges)
            };
            FuseEdgeVertices(clone.Edges, clone.Vertices);
            return clone;
        }

        static void FuseEdgeVertices<T>(Edge<T>[] edges, Vertex[] vertices)
        {
            for(int i = 0; i < edges.Length; ++i)
            {
                edges[i].Start = vertices[i];
                edges[i].End = vertices[( i + 1 ) % edges.Length];
            }
        }

        public static T Clone<T>(T node)
            where T : ILocatable, new()
        {
            T clone = new T
            {
                Position = node.Position
            };
            return clone;
        }

        public static Edge<T>[] Clone<T>(IList<Edge<T>> edges)
        {
            Edge<T>[] clones = new Edge<T>[edges.Count];
            for (int i = 0; i < edges.Count; ++i)
            {
                clones[i] = Clone(edges[i]);
            }
            return clones;
        }

        public static Edge<T> Clone<T>(Edge<T> edge)
        {
            Edge<T> clone = new Edge<T>()
            {
                IsBoundary = edge.IsBoundary,
                Cell = edge.Cell,
                BoundaryEdgeNumber = edge.BoundaryEdgeNumber,
                Start = Clone(edge.Start),
                End = Clone(edge.End),
                Twin = edge.Twin
            };
            return clone;
        }

        static Vertex[] Clone(IList<Vertex> vertices)
        {
            Vertex[] clones = new Vertex[vertices.Count];
            for (int i = 0; i < vertices.Count; ++i)
            {
                clones[i] = Clone(vertices[i]);
            }
            return clones;
        }

        static Vertex Clone(Vertex vertex)
        {
            Vertex clone = new Vertex
            {
                Position = vertex.Position,
                ID = vertex.ID
            };
            return clone;
        }
    }
}
