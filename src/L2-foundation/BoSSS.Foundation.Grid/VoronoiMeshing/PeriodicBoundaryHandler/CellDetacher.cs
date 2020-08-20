using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class TwinEdgeComparer<T> : EqualityComparer<Edge<T>>
    {
        public override bool Equals(Edge<T> x, Edge<T> y)
        {
            return ((x.Start.ID == y.Start.ID && x.End.ID == y.End.ID) || (x.Start.ID == y.Twin.Start.ID && x.End.ID == y.Twin.End.ID));
        }

        public override int GetHashCode(Edge<T> edge)
        {
            //http://eternallyconfuzzled.com/tuts/algorithms/jsw_tut_hashing.aspx 
            //If x == y hash must be hash(x) = hash(y)
            int start = Math.Max(edge.Start.ID, edge.Twin.Start.ID);
            int end = Math.Max(edge.End.ID, edge.Twin.End.ID);

            int hash = 17;
            hash = hash * 23 + start.GetHashCode();
            hash = hash * 23 + end.GetHashCode();

            return hash;
        }
    }

    class CellDetacher<T>
    {
        readonly MeshEdgeDivider<T> edgeDivider;

        readonly Boundary<T> boundary;

        readonly BoundaryElementEnumerator<T> boundaryEnumerator;

        readonly PeriodicMap map;

        readonly Dictionary<Edge<T>, (int, int, int)> edges;

        public CellDetacher(Domain<T> mesh, PeriodicMap map)
        {
            edgeDivider = new MeshEdgeDivider<T>(mesh.Mesh);
            boundary = mesh.Boundary;
            boundaryEnumerator = new BoundaryElementEnumerator<T>(mesh);
            this.map = map;
            edges = new Dictionary<Edge<T>, (int, int, int)>(new TwinEdgeComparer<T>());
        }

        public void DetachCells(IList<(MeshCell<T>, bool)> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            MoveBoundary(cells, boundaryEdgeNumber, pairedBoundaryEdgeNumber);
            edgeDivider.DivideBoundary(cells);
        }

        public void MoveBoundary(IList<(MeshCell<T>, bool)> cells, int boundaryEdgeNumber, int pairedBoundaryEdgeNumber)
        {
            foreach ((MeshCell<T> cell, bool isSplit) cell in cells)
            {
                if (cell.isSplit)
                {
                    cell.cell.Type = MeshCellType.Outside;
                }
            }
            if (boundary.FirstCorner.Type == MeshCellType.Outside)
            {
                MoveBoundaryCorner();
            }
            foreach ((MeshCell<T> cell, bool isSplit) cell in cells)
            {
                if (cell.isSplit)
                {
                    MoveBoundary(cell.cell.Edges, boundaryEdgeNumber, pairedBoundaryEdgeNumber);
                }
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

        void MoveBoundary(Edge<T>[] edges, int boundaryEdgeNumber, int pairedBoundaryNumber)
        {
            foreach (Edge<T> edge in edges)
            {
                SwitchBoundary(edge, boundaryEdgeNumber, pairedBoundaryNumber);
            }
        }

        void SwitchBoundary(Edge<T> edge, int boundaryNumber, int pairedBoundaryNumber)
        {
            if (edge.IsBoundary)
            {
                if (edge.BoundaryEdgeNumber == boundaryNumber || edge.BoundaryEdgeNumber == pairedBoundaryNumber)
                {
                    edge.IsBoundary = false;
                    edge.BoundaryEdgeNumber = -1;
                    if (edge.Twin != null)
                    {
                        edge.Twin.IsBoundary = false;
                        edge.Twin.BoundaryEdgeNumber = -1;
                    }
                }
                else
                {
                    if (edge.BoundaryEdgeNumber != edge.Twin.BoundaryEdgeNumber)
                    {
                        SwitchCornerEdges(edge, boundaryNumber);
                    }
                }
            }
            else
            {
                edge.IsBoundary = true;
                edge.BoundaryEdgeNumber = pairedBoundaryNumber;
                if (edge.Twin != null)
                {
                    edge.Twin.IsBoundary = true;
                    edge.Twin.BoundaryEdgeNumber = boundaryNumber;
                }
            }
        }

        void SwitchCornerEdges(Edge<T> edge, int boundaryNumber)
        {
            if(edges.TryGetValue(edge, out(int, int, int) a))
            {
                if(edge.Start.ID == a.Item1)
                {
                    edge.BoundaryEdgeNumber = a.Item2;
                    edge.Twin.BoundaryEdgeNumber = a.Item3;
                }
                else
                {
                    edge.BoundaryEdgeNumber = a.Item3;
                    edge.Twin.BoundaryEdgeNumber = a.Item2;
                }
            }
            else
            {
                edges.Add(edge, (edge.Start.ID, edge.BoundaryEdgeNumber, edge.Twin.BoundaryEdgeNumber));
                Corner newCorner = new Corner
                {
                    FirstEdge = boundaryNumber,
                    SecondEdge = edge.Twin.BoundaryEdgeNumber
                };
                int cornerBoundary = map.PeriodicCornerCorrelation[newCorner];
                edge.BoundaryEdgeNumber = map.PeriodicBoundaryCorrelation[cornerBoundary];
                edge.Twin.BoundaryEdgeNumber = cornerBoundary;
            }
        }
    }
}
