using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class Boundary<T>
    {
        public MeshCell<T> FirstCorner;

        public BoundaryLine[] BoundaryLines;
    }

    class Domain<T> : IMesh<T>
    {
        public IDMesh<T> Mesh { get; set; }

        public Boundary<T> Boundary { get; set; }

        public T CornerNode => Boundary.FirstCorner.Node;

        public List<MeshCell<T>> Cells => Mesh.Cells;

        public List<T> Nodes => Mesh.Nodes;
    }
}
