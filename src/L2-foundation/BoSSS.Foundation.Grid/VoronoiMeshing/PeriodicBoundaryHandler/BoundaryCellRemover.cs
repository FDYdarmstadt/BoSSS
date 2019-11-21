﻿using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryCellRemover<T>
        where T : ILocatable
    {
        readonly Mesh<T> mesh;

        readonly InsideCellEnumerator<T> insideCells;

        readonly Queue<(IList<MeshCell<T>>, int, int)> meshCellsToRemove;

        public BoundaryCellRemover(Mesh<T> mesh, int firstCellNodeIndice)
        {
            this.mesh = mesh;
            insideCells = new InsideCellEnumerator<T>(mesh, firstCellNodeIndice);
            meshCellsToRemove = new Queue<(IList<MeshCell<T>>, int, int)>();
        }

        public void EnqueueForRemoval(IList<MeshCell<T>> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            meshCellsToRemove.Enqueue((cells, boundaryEdgeNumber, pairedBoundaryEdgeNumber));
        }

        void SetAsBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber, int pairedBoundaryNumber)
        {
            foreach (MeshCell<T> cell in cells)
            {
                SetAllEdgesToBoundary(cell.Edges, boundaryEdgeNumber, pairedBoundaryNumber);
            }
        }

        void SetAllEdgesToBoundary(Edge<T>[] edges, int boundaryEdgeNumber, int pairedBoundaryNumber)
        {
            foreach (Edge<T> edge in edges)
            {
                edge.IsBoundary = true;
                edge.BoundaryEdgeNumber = pairedBoundaryNumber;
                if (edge.Twin != null)
                {
                    edge.Twin.BoundaryEdgeNumber = boundaryEdgeNumber;
                    edge.Twin.IsBoundary = true;
                }
            }
        }

        public void RemoveQueuedCells()
        {
            while(meshCellsToRemove.Count > 0)
            {
                (IList<MeshCell<T>> cells, int boundaryNumber, int pairedBoundaryNumber) = meshCellsToRemove.Dequeue();
                SetAsBoundary(cells, boundaryNumber, pairedBoundaryNumber);
            }
            RemoveOuterCellsFromMesh();
        }

        public void SetQueuedCellsAsOuter()
        {
            while (meshCellsToRemove.Count > 0)
            {
                (IList<MeshCell<T>> cells, int boundaryNumber, int pairedBoundaryNumber) = meshCellsToRemove.Dequeue();
                SetAsBoundary(cells, boundaryNumber, pairedBoundaryNumber);
            }
        }

        public void RemoveOuterCellsFromMesh()
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
