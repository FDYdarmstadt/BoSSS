using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class IDMesh<T>
    {
        public List<MeshCell<T>> Cells;

        public List<Vertex> Vertices;

        public List<T> Nodes;

        public IDMesh()
        {
            Cells = new List<MeshCell<T>>();
            Vertices = new List<Vertex>();
            Nodes = new List<T>();
        }

        public IDMesh(List<MeshCell<T>> Cells, List<Vertex> Vertices)
        {
            this.Cells = Cells;
            this.Vertices = Vertices;
            InitializeNodes();
        }

        void InitializeNodes()
        {
            Nodes = new List<T>(Cells.Count);
            foreach(MeshCell<T> cell in Cells)
            {
                Nodes.Add(cell.Node);
            }
        }

        public int AddCell(MeshCell<T> cell)
        {
            cell.ID = Cells.Count;
            Cells.Add(cell);
            Nodes.Add(cell.Node);
            return cell.ID;
        }

        public MeshCell<T> GetCell(int ID)
        {
            return Cells[ID];
        }

        public int AddVertex(Vertex vert)
        {
            vert.ID = Vertices.Count;
            Vertices.Add(vert);
            return vert.ID;
        }

        public Vertex GetVertex(int ID)
        {
            return Vertices[ID];
        }
    }
}
