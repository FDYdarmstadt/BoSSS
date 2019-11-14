using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerBoundaryAssigner<T>
    {
        int firstEdgeID;

        bool firstCircle;

        readonly ICollection<EdgeID> visitedEdges;

        readonly PeriodicMap map;

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
                bool so = startID == other.startID
                    && endID == other.endID
                    && twinStartID == other.twinStartID
                    && twinEndID == other.twinEndID;

                bool oderSo = startID == other.twinStartID
                    && endID == other.twinEndID
                    && twinStartID == other.startID
                    && twinEndID == other.endID;

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

        readonly Queue<(Corner, Edge<T>)> periodicCorners;

        BoundaryAssigner<T> boundaryAssigner;

        public PeriodicCornerBoundaryAssigner(PeriodicMap map)
        {
            this.map = map;
            visitedEdges = new LinkedList<EdgeID>();
            periodicCorners = new Queue<(Corner, Edge<T>)>();
            boundaryAssigner = new BoundaryAssigner<T>(map);
        }

        public void SetEdges(MeshCell<T> cornerCell)
        {
            visitedEdges.Clear();
            FindWronglyAssignedBoundaries(cornerCell);
            boundaryAssigner.AssignBoundariesOfPeriodicCorners(periodicCorners);
        }

        Edge<T> FindFirstExitEdge(MeshCell<T> cornerCell)
        {
            Edge<T> firstEdge = GetArbitraryInnerEdge(cornerCell);
            IEnumerable<Edge<T>> firstEdges = PositiveRotationBoundaryEdgesBeginningWith(firstEdge);

            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(firstEdges);
            foreach (var edgePair in edges)
            {
                if (FindExitEdgeInPositveRotation(edgePair, out Edge<T> current))
                {
                    if (!current.IsBoundary)
                    {
                        current = FindFirstExitEdge(current.Twin.Cell);
                    }
                    return current;
                }
            }
            throw new Exception("Did not find exit edge.");
        }

        Edge<T> GetArbitraryInnerEdge(MeshCell<T> cornerCell)
        {
            foreach (var edge in cornerCell.Edges)
            {
                if (!edge.IsBoundary)
                {
                    return edge;
                }
            }
            throw new Exception("Cell does not have a boundary.");
        }

        bool FindExitEdgeInPositveRotation(Pair<Edge<T>> edgePair, out Edge<T> edge)
        {
            edge = null;
            if (edgePair.Current.BoundaryEdgeNumber != edgePair.Previous.BoundaryEdgeNumber)
            {
                map.PeriodicBoundaryCorrelation.TryGetValue(edgePair.Current.BoundaryEdgeNumber, out int pairedBoundary);
                if (edgePair.Current.Twin.BoundaryEdgeNumber == pairedBoundary)
                {
                    edge = edgePair.Current;
                    return true;
                }
                else
                {
                    if (!edgePair.Previous.IsBoundary)
                    {
                        edge = edgePair.Previous;
                        return true;
                    }
                }
            }
            return false;
        }

        bool FindExitEdgeInNegativeRotation(Pair<Edge<T>> edgePair, out Edge<T> edge)
        {
            edge = null;
            if (edgePair.Current.BoundaryEdgeNumber != edgePair.Previous.BoundaryEdgeNumber)
            {
                if (edgePair.Current.IsBoundary)
                {
                    map.PeriodicBoundaryCorrelation.TryGetValue(edgePair.Current.BoundaryEdgeNumber, out int pairedBoundary);
                    if (edgePair.Current.Twin.BoundaryEdgeNumber == pairedBoundary)
                    {
                        edge = edgePair.Current;
                        return true;
                    }
                }
                else
                {
                    edge = edgePair.Current;
                    return true;
                }

            }
            return false;
        }

        IEnumerable<Edge<T>> PositiveRotationBoundaryEdgesBeginningWith(Edge<T> first)
        {
            MeshCell<T> cell = first.Cell;
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(cell.Edges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(cell.Edges, firstEdgeIndice);
            foreach (Edge<T> edge in edges)
            {
                yield return edge;
            }
        }

        IEnumerable<Edge<T>> NegativeRotationBoundaryEdgesBeginningWith(Edge<T> first)
        {
            MeshCell<T> cell = first.Cell;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(cell.Edges);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice);
            foreach (Edge<T> edge in edges)
            {
                yield return edge;
            }
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

        void FindWronglyAssignedBoundaries(MeshCell<T> cornerCell)
        {
            Edge<T> edge = FindFirstExitEdge(cornerCell);
            firstEdgeID = edge.Start.ID;
            firstCircle = true;
            while (firstCircle)
            {
                Edge<T> twin = edge.Twin;
                IEnumerable<Edge<T>> edges = NegativeRotationBoundaryEdgesBeginningWith(twin);
                edge = CheckForCorners(edges);
            }
        }

        Edge<T> CheckForCorners(IEnumerable<Edge<T>> boundaryEdges)
        {
            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(boundaryEdges);

            Edge<T> current = null;
            foreach (var edgePair in edges)
            {
                current = edgePair.Current;
                Edge<T> previous = edgePair.Previous;
                if (current.Start.ID != firstEdgeID)
                {
                    if (FindExitEdgeInNegativeRotation(edgePair, out Edge<T> exitEdge))
                    {
                        return exitEdge;
                    }
                    if (AlreadyVisited(current))
                    {
                        Corner periodicCorner = FromCell(current.Cell);
                        periodicCorners.Enqueue((periodicCorner, current));
                    };
                }
                else
                {
                    firstCircle = false;
                    return current;
                }
            }
            throw new Exception("Cell does not have a boundary.");
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
                firstEdge = Math.Min(Math.Max(edgeNumber, 0), firstEdge);
                secondEdge = Math.Max(edgeNumber, secondEdge);
            }

            Corner cellCorner = new Corner
            {
                FirstEdge = firstEdge,
                SecondEdge = secondEdge
            };
            return cellCorner;
        }
    }

    class BoundaryAssigner<T>
    {
        readonly PeriodicMap map;

        public BoundaryAssigner(PeriodicMap map)
        {
            this.map = map;
        }

        public void AssignBoundariesOfPeriodicCorners(Queue<(Corner, Edge<T>)> periodicCorners)
        {
            while (periodicCorners.Count > 0)
            {
                (Corner periodicCorner, Edge<T> current) = periodicCorners.Dequeue();
                AssignEdge(current, periodicCorner);
            }
        }

        void AssignEdge(Edge<T> edge, Corner corner)
        {
            Corner twin = CreateTwinOf(corner);

            map.PeriodicCornerCorrelation.TryGetValue(corner, out int boundary);
            edge.BoundaryEdgeNumber = boundary;
            map.PeriodicCornerCorrelation.TryGetValue(twin, out int twinBoundary);
            edge.Twin.BoundaryEdgeNumber = twinBoundary;

            //Debug.Assert(IsEdgeTransformation(map.PeriodicBoundaryTransformations[boundary], edge));
        }

        Corner CreateTwinOf(Corner corner)
        {
            map.PeriodicBoundaryCorrelation.TryGetValue(corner.FirstEdge, out int cornerFirstTwin);
            map.PeriodicBoundaryCorrelation.TryGetValue(corner.SecondEdge, out int cornerSecondTwin);
            Corner currentTwinCorner = new Corner
            {
                FirstEdge = cornerFirstTwin,
                SecondEdge = cornerSecondTwin
            };

            Debug.Assert(AreRegisteredInMap(corner, currentTwinCorner));

            return currentTwinCorner;
        }

        bool AreRegisteredInMap(Corner From, Corner To)
        {
            bool areCorners = map.PeriodicCornerCorrelation.ContainsKey(From)
                && map.PeriodicCornerCorrelation.ContainsKey(To);
            return areCorners;
        }

        static bool IsEdgeTransformation(Transformation transformation, Edge<T> edge)
        {
            Vector transformedStart = transformation.Transform(edge.Start.Position);
            double distance = Vector.Dist(transformedStart, edge.Twin.End.Position);
            return distance <= 1e-13;
        }
    }
}
