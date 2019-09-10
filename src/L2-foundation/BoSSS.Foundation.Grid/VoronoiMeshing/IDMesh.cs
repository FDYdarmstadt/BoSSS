using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class IDMesh<T>
    {
        /// <summary>
        /// Cell Id is Position in list Cells 
        /// </summary>
        public List<MeshCell<T>> Cells;

        /// <summary>
        /// Vertex Id is Position inidce in list Vertices
        /// </summary>
        public List<Vertex> Vertices;

        public List<T> Nodes;

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
}
