using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerBoundaryIdentifier<T>
    {
        struct EdgeID : IEquatable<EdgeID>
        {
            int startID;

            int endID;

            int twinStartID;

            int twinEndID;

            public EdgeID(Edge<T> edge)
            {
                startID = edge.Start.ID;
                endID = edge.End.ID;
                twinStartID = edge.Twin.Start.ID;
                twinEndID = edge.Twin.End.ID;
            }

            public bool Equals(EdgeID other)
            {
                bool so = (startID == other.startID && endID == other.endID)
                    || (twinStartID == other.twinStartID && twinEndID == other.twinEndID);

                bool oderSo = (startID == other.twinStartID && endID == other.twinEndID)
                    || (twinStartID == other.startID && twinEndID == other.endID);

                return so || oderSo;
            }

            public override bool Equals(object other)
            {
                if (other is EdgeID otherID)
                {
                    return Equals(otherID);
                }
                else
                {
                    return false;
                }
            }

            public override int GetHashCode()
            {
                return (startID + endID + twinStartID + twinEndID).GetHashCode();
            }
        }

        int firstEdgeID;

        bool firstCircle;

        readonly ICollection<EdgeID> visitedEdges;

        readonly PeriodicMap map;

        readonly Queue<Edge<T>> periodicEdges;

        Corner periodicCorner;

        PeriodicCornerBoundaryAssigner<T> boundaryAssigner;

        public PeriodicCornerBoundaryIdentifier(PeriodicMap map)
        {
            this.map = map;
            visitedEdges = new LinkedList<EdgeID>();
            periodicEdges = new Queue<Edge<T>>();
            boundaryAssigner = new PeriodicCornerBoundaryAssigner<T>(map);
            periodicCorner = new Corner
            {
                FirstEdge = -1,
                SecondEdge = -1
            };
        }

        public void SetEdges(MeshCell<T> cornerCell)
        {
            FindWronglyAssignedBoundaries(cornerCell);
            boundaryAssigner.AssignBoundariesOfPeriodicCorners(periodicEdges, periodicCorner);
            ClearData();
        }

        void ClearData()
        {
            periodicCorner.FirstEdge = -1;
            periodicCorner.SecondEdge = -1;
            visitedEdges.Clear();
            periodicEdges.Clear();
        }

        void FindWronglyAssignedBoundaries(MeshCell<T> cornerCell)
        {
            Edge<T> firstEdge = GetFirstEdgeOfBoundary(cornerCell);
            IEnumerable<Edge<T>> firstEdges = PositiveRotationOfBoundaryEdgesBeginningWith(firstEdge);
            Edge<T> edge = FindFirstBridge(firstEdges);

            firstCircle = true;
            while (firstCircle)
            {
                Edge<T> twin = edge.Twin;
                IEnumerable<Edge<T>> edges = NegativeRotationOfBoundaryEdgesBeginningWith(twin);
                edge = FindBridge(edges);
            }
        }

        Edge<T> GetFirstEdgeOfBoundary(MeshCell<T> cornerCell)
        {
            foreach (var edge in new Convolution<Edge<T>>(cornerCell.Edges))
            {
                if (edge.Current.IsBoundary && !edge.Previous.IsBoundary)
                {
                    return edge.Current;
                }
            }
            throw new Exception("Cell does not have a boundary.");
        }

        IEnumerable<Edge<T>> PositiveRotationOfBoundaryEdgesBeginningWith(Edge<T> first)
        {
            IEnumerable<Edge<T>> edges = PositiveEdgeRotatationOfCellAfter(first, 0);
            while (true)
            {
                foreach (Edge<T> edge in edges)
                {
                    if (edge.IsBoundary)
                    {
                        yield return edge;
                    }
                    else
                    {
                        edges = PositiveEdgeRotatationOfCellAfter(edge.Twin, 1);
                        break;
                    }
                }
            }
        }

        CyclicArray<Edge<T>> PositiveEdgeRotatationOfCellAfter(Edge<T> first, int offset)
        {
            MeshCell<T> cell = first.Cell;
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(cell.Edges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(cell.Edges, firstEdgeIndice + offset);
            return edges;
        }

        Edge<T> FindFirstBridge(IEnumerable<Edge<T>> firstEdges)
        {
            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(firstEdges);
            foreach (var edgePair in edges)
            {
                if (IsBridgeInPositveRotation(edgePair, out Edge<T> current))
                {
                    firstEdgeID = current.Start.ID;
                    return current;
                }
            }
            throw new Exception("Did not find exit edge.");
        }

        bool IsBridgeInPositveRotation(Pair<Edge<T>> edgePair, out Edge<T> edge)
        {
            edge = null;
            if ((edgePair.Current.BoundaryEdgeNumber != edgePair.Previous.BoundaryEdgeNumber)
                || (edgePair.Current.Twin.BoundaryEdgeNumber != edgePair.Previous.Twin.BoundaryEdgeNumber))
            {
                map.PeriodicBoundaryCorrelation.TryGetValue(edgePair.Previous.BoundaryEdgeNumber, out int pairedBoundary);
                if (edgePair.Previous.Twin.BoundaryEdgeNumber == pairedBoundary)
                {
                    edge = edgePair.Previous;
                    return true;
                }
            }
            return false;
        }

        static int FindIndiceOfEdgeInItsCell(Edge<T>[] edges, Edge<T> first)
        {
            for (int i = 0; i < edges.Length; ++i)
            {
                if (first.Start.ID == edges[i].Start.ID)
                {
                    return i;
                }
            }
            throw new Exception("Edge not found.");
        }

        IEnumerable<Edge<T>> NegativeRotationOfBoundaryEdgesBeginningWith(Edge<T> first)
        {
            IList<Edge<T>> edges = NegativeEdgeRotationOfCellAfter(first, 0);
            while (true)
            {
                foreach (Edge<T> edge in edges)
                {
                    if (edge.IsBoundary)
                    {
                        yield return edge;
                    }
                    else
                    {
                        edges = NegativeEdgeRotationOfCellAfter(edge.Twin, 1);
                        break;
                    }
                }
            }
        }

        CyclicArray<Edge<T>> NegativeEdgeRotationOfCellAfter(Edge<T> first, int offset)
        {
            Edge<T>[] cell = first.Cell.Edges;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(cell);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice + offset);
            return edges;
        }

        Edge<T> FindBridge(IEnumerable<Edge<T>> boundaryEdges)
        {
            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(boundaryEdges);

            Edge<T> current = null;
            foreach (var edgePair in edges)
            {
                current = edgePair.Current;
                Edge<T> previous = edgePair.Previous;
                if (current.Start.ID == firstEdgeID)
                {
                    firstCircle = false;
                    return current;
                }
                else
                {
                    if (AlreadyVisited(current))
                    {
                        if (periodicCorner.FirstEdge == periodicCorner.SecondEdge)
                        {
                            periodicCorner = FromCell(current.Cell);
                        }
                        periodicEdges.Enqueue(current);
                    }
                    else
                    {
                        if (IsBridgeInNegativeRotation(edgePair, out Edge<T> exitEdge))
                        {
                            return exitEdge;
                        }
                    };
                }
            }
            throw new Exception("Cell does not have a boundary.");
        }

        bool IsBridgeInNegativeRotation(Pair<Edge<T>> edgePair, out Edge<T> edge)
        {
            edge = null;
            if ((edgePair.Current.BoundaryEdgeNumber != edgePair.Previous.BoundaryEdgeNumber)
                || (edgePair.Current.Twin.BoundaryEdgeNumber != edgePair.Previous.Twin.BoundaryEdgeNumber))
            {
                map.PeriodicBoundaryCorrelation.TryGetValue(edgePair.Current.BoundaryEdgeNumber, out int pairedBoundary);
                if (edgePair.Current.Twin.BoundaryEdgeNumber == pairedBoundary)
                {
                    edge = edgePair.Current;
                    return true;
                }
            }
            return false;
        }

        bool AlreadyVisited(Edge<T> current)
        {
            EdgeID id = new EdgeID(current);
            if (visitedEdges.Contains(id))
            {
                return true;
            }
            else
            {
                visitedEdges.Add(id);
                return false;
            }
        }

        Corner FromCell(MeshCell<T> cell)
        {
            int firstEdge = int.MaxValue;
            int secondEdge = int.MinValue;

            foreach (Edge<T> edge in cell.Edges)
            {
                int edgeNumber = edge.BoundaryEdgeNumber;
                if (edgeNumber > -1)
                {
                    firstEdge = Math.Min(edgeNumber, firstEdge);
                    secondEdge = Math.Max(edgeNumber, secondEdge);
                }
            }

            Corner cellCorner = new Corner
            {
                FirstEdge = firstEdge,
                SecondEdge = secondEdge
            };
            return cellCorner;
        }
    }
}
