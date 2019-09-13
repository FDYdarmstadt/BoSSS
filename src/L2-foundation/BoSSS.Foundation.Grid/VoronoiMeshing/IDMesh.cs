using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class IDMesh<T> : Mesh<T>
    {
        public List<Vertex> Vertices;

        public int AddCell(MeshCell<T> cell)
        {
            cell.ID = Cells.Count;
            Cells.Add(cell);
            return cell.ID;
        }

        public int AddVertex(Vertex vert)
        {
            vert.ID = Vertices.Count;
            Vertices.Add(vert);
            return vert.ID;
        }
    }

    class Mesh<T>
    {
        public List<MeshCell<T>> Cells;

        public List<T> Nodes;
    }
}
