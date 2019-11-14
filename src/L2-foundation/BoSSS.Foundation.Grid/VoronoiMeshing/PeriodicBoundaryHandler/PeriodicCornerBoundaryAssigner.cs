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

        public PeriodicCornerBoundaryAssigner(PeriodicMap map)
        {
            this.map = map;
            visitedEdges = new LinkedList<EdgeID>();
        }

        public void SetEdges(MeshCell<T> cornerCell)
        {
            visitedEdges.Clear();
            IEnumerable<Edge<T>> firstEdges = PositiveRotationEdgesBeginningWith(cornerCell.Edges[0]);
            Edge<T> firstEdge = FindFirstExitEdge(firstEdges);
            firstEdgeID = firstEdge.Start.ID;

            FindWronglyAssignedBoundaries(firstEdge);
        }

        Edge<T> FindFirstExitEdge(IEnumerable<Edge<T>> boundaryEdges)
        {
            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(boundaryEdges);

            Edge<T> current = null;
            foreach (var edgePair in edges)
            {
                current = edgePair.Current;
                if (edgePair.Current.BoundaryEdgeNumber != edgePair.Previous.BoundaryEdgeNumber)
                {
                    map.PeriodicBoundaryCorrelation.TryGetValue(edgePair.Current.BoundaryEdgeNumber, out int pairedBoundary);
                    if (edgePair.Current.Twin.BoundaryEdgeNumber == pairedBoundary)
                    {
                        break;
                    }
                }
            }
            return current;
        }

        IEnumerable<Edge<T>> PositiveRotationEdgesBeginningWith(Edge<T> first)
        {
            MeshCell<T> cell = first.Cell;
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(cell.Edges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(cell.Edges, firstEdgeIndice);
            foreach (Edge<T> edge in edges)
            {
                if (edge.IsBoundary)
                {
                    yield return edge;
                }
            }
        }

        IEnumerable<Edge<T>> NegativeRotationEdgesBeginningWith(Edge<T> first)
        {
            MeshCell<T> cell = first.Cell;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(cell.Edges);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice);
            foreach (Edge<T> edge in edges)
            {
                if (edge.IsBoundary)
                {
                    yield return edge;
                }
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

        bool firstCircle;

        void FindWronglyAssignedBoundaries(Edge<T> edge)
        {
            firstCircle = true;
            while (firstCircle)
            {
                Edge<T> twin = edge.Twin;
                IEnumerable<Edge<T>> edges = NegativeRotationEdgesBeginningWith(twin);
                edge = AssignCornerTo(edges);
            }
        }

        Edge<T> AssignCornerTo(IEnumerable<Edge<T>> boundaryEdges)
        {
            Convolution<Edge<T>> edges = new Convolution<Edge<T>>(boundaryEdges);

            Edge<T> current = null;
            foreach (var edgePair in edges)
            {
                current = edgePair.Current;
                Edge<T> previous = edgePair.Previous;
                if (current.Start.ID != firstEdgeID)
                {
                    if (current.BoundaryEdgeNumber != previous.BoundaryEdgeNumber)
                    {
                        map.PeriodicBoundaryCorrelation.TryGetValue(current.BoundaryEdgeNumber, out int pairedBoundary);
                        if (current.Twin.BoundaryEdgeNumber == pairedBoundary)
                        {
                            break;
                        }
                    }
                    if (AlreadyVisited(current))
                    {
                        Corner periodicCorner = FromCell(current.Cell);
                        AssignEdge(current, periodicCorner);
                    };
                }
                else
                {
                    firstCircle = false;
                    break;
                }
            }
            return current;
        }

        Corner FromCell(MeshCell<T> cell)
        {
            int firstEdge = int.MaxValue;
            int secondEdge = int.MinValue;

            foreach(Edge<T> edge in cell.Edges)
            {
                int edgeNumber = edge.BoundaryEdgeNumber;
                firstEdge = Math.Min( Math.Max(edgeNumber, 0), firstEdge);
                secondEdge = Math.Max(edgeNumber, secondEdge);
            }

            Corner cellCorner = new Corner
            {
                FirstEdge = firstEdge,
                SecondEdge = secondEdge
            };
            return cellCorner;
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

        void AssignEdge(Edge<T> edge, Corner corner)
        {
            Corner twin = CreateTwinOf(corner);

            map.PeriodicCornerCorrelation.TryGetValue(corner, out int boundary);
            edge.BoundaryEdgeNumber = boundary;
            map.PeriodicCornerCorrelation.TryGetValue(twin, out int twinBoundary);
            edge.Twin.BoundaryEdgeNumber = twinBoundary;

            Debug.Assert(IsEdgeTransformation(map.PeriodicBoundaryTransformations[boundary], edge));
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
