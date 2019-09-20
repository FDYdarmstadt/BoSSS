using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryCellRemover<T>
        where T : IMesherNode
    {
        readonly Mesh<T> mesh;

        readonly InsideCellEnumerator<T> insideCells;

        public BoundaryCellRemover(Mesh<T> mesh, int firstCellNodeIndice)
        {
            this.mesh = mesh;
            insideCells = new InsideCellEnumerator<T>(mesh, firstCellNodeIndice);
        }

        public void Remove(IList<MeshCell<T>> cells)
        {
            foreach(MeshCell<T> cell in cells)
            {
                SetAllEdgeTwinsToBoundaryEdge(cell.Edges);
            }
            RemoveCellsFromMesh();
        }

        void RemoveCellsFromMesh()
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach(MeshCell<T> cell in insideCells.Cells())
            {
                cells.Add(cell);
            }
            mesh.Cells = cells;
            mesh.Nodes.Clear();
            foreach (MeshCell<T> cell in mesh.Cells)
            {
                mesh.Nodes.Add(cell.Node);
            }
        }

        void SetAllEdgeTwinsToBoundaryEdge(Edge<T>[] edges)
        {
            foreach (Edge<T> edge in edges)
            {
                if (edge.Twin != null)
                {
                    edge.Twin.IsBoundary = true;
                }
            }
        }
    }
}
