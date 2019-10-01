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

        public void SetAsBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
        {
            foreach (MeshCell<T> cell in cells)
            {
                SetAllEdgesToBoundary(cell.Edges, boundaryEdgeNumber);
            }
        }

        void SetAllEdgesToBoundary(Edge<T>[] edges, int boundaryEdgeNumber)
        {
            foreach (Edge<T> edge in edges)
            {
                edge.IsBoundary = true;
                if (edge.Twin != null)
                {
                    edge.Twin.BoundaryEdgeNumber = boundaryEdgeNumber;
                    edge.Twin.IsBoundary = true;
                }
            }
        }

        public void RemoveOuterCellsFromMesh()
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach (MeshCell<T> cell in insideCells.Cells())
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
    }
}
