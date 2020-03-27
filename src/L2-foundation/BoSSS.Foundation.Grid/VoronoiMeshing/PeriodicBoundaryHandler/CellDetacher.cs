using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class CellDetacher<T>
    {
        readonly MeshEdgeDivider<T> edgeDivider;

        readonly Boundary<T> boundary;

        readonly BoundaryElementEnumerator<T> boundaryEnumerator;

        public CellDetacher(Domain<T> mesh)
        {
            edgeDivider = new MeshEdgeDivider<T>(mesh.Mesh);
            boundary = mesh.Boundary;
            boundaryEnumerator = new BoundaryElementEnumerator<T>(mesh);
        }

        public void DetachCells(IList<MeshCell<T>> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            MoveBoundary(cells, boundaryEdgeNumber, pairedBoundaryEdgeNumber);
            edgeDivider.DivideBoundary(cells);
        }

        public void MoveBoundary(IList<MeshCell<T>> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            foreach (MeshCell<T> cell in cells)
            {
                cell.Type = MeshCellType.Outside;
            }
            if (boundary.FirstCorner.Type == MeshCellType.Outside)
            {
                MoveBoundaryCorner();
            }
            foreach (MeshCell<T> cell in cells)
            {
                SetAllEdgesToBoundary(cell.Edges, boundaryEdgeNumber, pairedBoundaryEdgeNumber);
            }
        }

        void MoveBoundaryCorner()
        {
            foreach (MeshCell<T> cell in boundaryEnumerator.CycleCells())
            {
                if (cell.Type != MeshCellType.Outside)
                {
                    boundary.FirstCorner = cell;
                    break;
                }
            }
        }

        void SetAllEdgesToBoundary(Edge<T>[] edges, int boundaryEdgeNumber, int pairedBoundaryNumber)
        {
            foreach (Edge<T> edge in edges)
            {
                SwitchBoundary(edge, pairedBoundaryNumber, boundaryEdgeNumber);
                if (edge.Twin != null)
                {
                    SwitchBoundary(edge.Twin, boundaryEdgeNumber, pairedBoundaryNumber);
                }
            }
        }

        bool IsNotPeriodic(Edge<T> edge)
        {
            if (edge.Twin != null)
            {
                return ((edge.Twin.Start.ID == edge.End.ID) 
                    || (edge.Twin.End.ID == edge.Start.ID));
            }
            else
            {
                return true;
            }
        }

        void SwitchBoundary(Edge<T> edge, int innerEdgeNumber, int pairedEdgeNumber)
        {
            if (edge.IsBoundary)
            {
                if(edge.BoundaryEdgeNumber == innerEdgeNumber || edge.BoundaryEdgeNumber == pairedEdgeNumber)
                {
                    edge.IsBoundary = false;
                    edge.BoundaryEdgeNumber = -1;
                }
            }
            else
            {
                edge.IsBoundary = true;
                edge.BoundaryEdgeNumber = innerEdgeNumber;
            }
        }
    }
}
