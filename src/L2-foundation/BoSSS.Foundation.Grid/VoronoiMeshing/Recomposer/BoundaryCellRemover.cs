using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryCellRemover<T>
        where T : ILocatable
    {
        readonly Mesh<T> mesh;

        readonly InsideCellEnumerator<T> insideCells;

        readonly Queue<(IList<MeshCell<T>>, int)> meshCellsToRemove;

        public BoundaryCellRemover(Mesh<T> mesh, int firstCellNodeIndice)
        {
            this.mesh = mesh;
            insideCells = new InsideCellEnumerator<T>(mesh, firstCellNodeIndice);
            meshCellsToRemove = new Queue<(IList<MeshCell<T>>, int)>();
        }

        public void EnqueueForRemoval(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
        {
            meshCellsToRemove.Enqueue((cells, boundaryEdgeNumber));
        }

        public void SetCellsAsBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
        {
            SetAsBoundary(cells, boundaryEdgeNumber);
        }

        void SetAsBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
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
                edge.BoundaryEdgeNumber = boundaryEdgeNumber;
                if (edge.Twin != null)
                {
                    if(edge.Twin.Start.ID == edge.End.ID || edge.Twin.End.ID == edge.Start.ID)
                    {
                        edge.Twin.BoundaryEdgeNumber = boundaryEdgeNumber;
                        edge.Twin.IsBoundary = true;
                    }
                }
            }
        }

        public void RemoveQueuedCells()
        {
            while(meshCellsToRemove.Count > 0)
            {
                (IList<MeshCell<T>> cells, int boundaryNumber) = meshCellsToRemove.Dequeue();
                SetAsBoundary(cells, boundaryNumber);
            }
            RemoveOuterCellsFromMesh();
        }

        void RemoveOuterCellsFromMesh()
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach (MeshCell<T> cell in insideCells.EnumerateCellsInConcentricCircles())
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
