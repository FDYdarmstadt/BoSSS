using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryMover<T>
    {
        readonly Mesh<T> mesh;

        readonly InsideCellEnumerator<T> insideCells;

        readonly MeshEdgeDivider<T> edgeDivider;

        public BoundaryMover(IDMesh<T> mesh, int firstCellNodeIndice)
        {
            this.mesh = mesh;
            insideCells = new InsideCellEnumerator<T>(mesh, firstCellNodeIndice);
            edgeDivider = new MeshEdgeDivider<T>(mesh);
        }

        public void MoveBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            foreach (MeshCell<T> cell in cells)
            {
                cell.type = MeshCellType.Outside;
                SetAllEdgesToBoundary(cell.Edges, boundaryEdgeNumber, pairedBoundaryEdgeNumber);
            }
        }

        void SetAllEdgesToBoundary(Edge<T>[] edges, int boundaryEdgeNumber, int pairedBoundaryNumber)
        {
            foreach (Edge<T> edge in edges)
            {
                SwitchBoundary(edge, pairedBoundaryNumber);
                if (edge.Twin != null)
                {
                    SwitchBoundary(edge.Twin, boundaryEdgeNumber);
                }
            }
        }

        void SwitchBoundary(Edge<T> edge, int boundaryEdgeNumber)
        {
            if (edge.IsBoundary)
            {
                edge.IsBoundary = false;
                edge.BoundaryEdgeNumber = -1;
            }
            else
            {
                edge.IsBoundary = true;
                edge.BoundaryEdgeNumber = boundaryEdgeNumber;
            }
        }

        public void DivideBoundary(IList<MeshCell<T>> cells)
        {
            foreach(MeshCell<T> cell in cells)
            {
                for(int i = 0; i < cell.Edges.Length; ++i)
                {
                    Edge<T> edge = cell.Edges[i];
                    if (edge.IsBoundary)
                    {
                        edgeDivider.DivideEdge(edge, i);
                    }
                }
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
