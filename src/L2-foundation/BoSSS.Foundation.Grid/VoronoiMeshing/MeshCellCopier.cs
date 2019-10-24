using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshCellCopier<T>
    {
        readonly IDMesh<T> mesh;

        public MeshCellCopier(IDMesh<T> mesh)
        {
            this.mesh = mesh;
        }

        public MeshCell<T>[] Copy(IList<MeshCell<T>> cells)
        {
            MeshCell<T>[] copies = MeshCellCloner.Clone(cells);
            foreach(MeshCell<T> cell in copies)
            {
                RegisterToMesh(cell);
            }
            return copies;
        }

        void RegisterToMesh(MeshCell<T> cell)
        {
            mesh.AddCell(cell);
            foreach(Vertex node in cell.Vertices)
            {
                mesh.AddVertex(node);
            }
        }

        public MeshCell<T> Copy(MeshCell<T> cell)
        {
            MeshCell<T> copy = MeshCellCloner.Clone(cell);
            RegisterToMesh(copy);
            return copy;
        }


    }
}
