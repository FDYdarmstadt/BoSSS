using ilPSP;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class Vertex
    {
        public Vector Position { get; set; }

        public int ID { get; set; }

        public static explicit operator Vector(Vertex vtx)
        {
            return vtx.Position;
        }
    }

    class Line
    {
        public Vertex Start { get; set; }

        public Vertex End { get; set; }
    }

    class Edge<T> : Line
    {
        public bool IsBoundary = false;

        public int BoundaryEdgeNumber = -1;

        public Edge<T> Twin { get; set; }

        public MeshCell<T> Cell { get; set; }
    }

    class EdgeComparer<T> : EqualityComparer<Edge<T>>
    {
        public override bool Equals(Edge<T> x, Edge<T> y)
        {
            return ((x.Start.ID == y.Start.ID && x.End.ID == y.End.ID) || (x.Start.ID == y.End.ID && x.End.ID == y.Start.ID));
        }

        public override int GetHashCode(Edge<T> obj)
        {
            //http://eternallyconfuzzled.com/tuts/algorithms/jsw_tut_hashing.aspx 
            //If x == y hash must be hash(x) = hash(y)
            int start = obj.End.ID > obj.Start.ID ? obj.Start.ID : obj.End.ID;
            int end = obj.End.ID < obj.Start.ID ? obj.Start.ID : obj.End.ID;

            int hash = 17;
            hash = hash * 23 + start.GetHashCode();
            hash = hash * 23 + end.GetHashCode();

            return hash;
        }
    }

    enum MeshCellType { Inside, NotDetermined, Outside }

    class MeshCell<T>
    {
        public MeshCellType Type = MeshCellType.NotDetermined;

        public T Node;
        public int ID { get; set; }
        public Vertex[] Vertices { get; set; }
        public Edge<T>[] Edges { get; set; }
        public int IntersectionVertex { get; set; }
    }
}
