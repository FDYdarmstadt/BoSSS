using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerBoundaryIdentifier<T>
    {
        struct Crossing
        {
            public Edge<T> First;

            public Edge<T> Second;

            public Edge<T> Third;
        }

        class CrossingFinder
        {
            public Stack<Edge<T>> Visited { get; private set; }

            readonly IEnumerator<Edge<T>> inside;

            readonly IEnumerator<Edge<T>> outside;

            public static PeriodicMap map;

            public Crossing Crossing {
                get {
                    return new Crossing
                    {
                        First = inside.Current,
                        Second = outside.Current.Twin,
                        Third = Visited.Peek().Twin
                    };
                }
            }

            public CrossingFinder(Edge<T> edge)
            {
                inside = PositiveRotationOfBoundaryEdgesBeginningWith(edge, 1).GetEnumerator();
                outside = NegativeRotationOfBoundaryEdgesBeginningWith(edge.Twin, 1).GetEnumerator();
                Visited = new Stack<Edge<T>>();
                Visited.Push(edge);
            }

            public bool MoveForwardAndCheckForCrossing()
            {
                if (inside.MoveNext() & outside.MoveNext())
                {
                    Edge<T> a = inside.Current;
                    Edge<T> b = outside.Current;
                    if (AreTwins(a, b))
                    {
                        Visited.Push(a);
                        return false;
                    }
                    else
                    {
                        return true;
                    }
                }
                else
                {
                    throw new Exception("Error in Boundary: must not end.");
                }
            }

            bool AreTwins(Edge<T> A, Edge<T> B)
            {
                if (map.PeriodicBoundaryCorrelation.TryGetValue(A.BoundaryEdgeNumber, out int pairedBoundary))
                {
                    return pairedBoundary == B.BoundaryEdgeNumber;
                }
                else
                {
                    throw new Exception("Boundary edge number is not registered in map.");
                }

            }
        }

        public PeriodicCornerBoundaryIdentifier(PeriodicMap map)
        {
            CrossingFinder.map = map;
        }

        public (Stack<Edge<T>>, Corner) FindEdges(MeshCell<T> cornerCell)
        {
            Edge<T> firstEdge = GetFirstEdgeOfBoundary(cornerCell);
            Crossing firstCrossing = FindFirstCrossing(firstEdge);

            Stack<Edge<T>> edges = default(Stack<Edge<T>>);
            Corner corner = default(Corner);
            bool stop = false;
            CrossingFinder first = new CrossingFinder(firstCrossing.First);
            CrossingFinder second = new CrossingFinder(firstCrossing.Second);
            CrossingFinder third = new CrossingFinder(firstCrossing.Third);
            while (!stop)
            {
                if (first.MoveForwardAndCheckForCrossing())
                {
                    stop = true;
                    edges = first.Visited;
                    corner = new Corner
                    {
                        FirstEdge = firstCrossing.Third.Twin.BoundaryEdgeNumber,
                        SecondEdge = first.Crossing.First.BoundaryEdgeNumber
                    };
                }
                else if (second.MoveForwardAndCheckForCrossing())
                {
                    stop = true;
                    edges = second.Visited;
                    corner = new Corner
                    {
                        FirstEdge = firstCrossing.First.Twin.BoundaryEdgeNumber,
                        SecondEdge = second.Crossing.First.BoundaryEdgeNumber
                    };
                }
                else if (third.MoveForwardAndCheckForCrossing())
                {
                    stop = true;
                    edges = third.Visited;
                    corner = new Corner
                    {
                        FirstEdge = firstCrossing.Second.Twin.BoundaryEdgeNumber,
                        SecondEdge = third.Crossing.First.BoundaryEdgeNumber
                    };
                }
            }
            return (edges, corner);
        }

        Crossing FindFirstCrossing(Edge<T> firstEdge)
        {
            CrossingFinder forward = new CrossingFinder(firstEdge);
            CrossingFinder backwards = new CrossingFinder(firstEdge.Twin);
            while (true)
            {
                if (forward.MoveForwardAndCheckForCrossing())
                {
                    return forward.Crossing;
                }
                else if (backwards.MoveForwardAndCheckForCrossing())
                {
                    return backwards.Crossing;
                };
            }
        }

        static Edge<T> GetFirstEdgeOfBoundary(MeshCell<T> cornerCell)
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

        static IEnumerable<Edge<T>> PositiveRotationOfBoundaryEdgesBeginningWith(Edge<T> first, int offset)
        {
            IEnumerable<Edge<T>> edges = PositiveEdgeRotationOfCellAfter(first, offset);
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
                        edges = PositiveEdgeRotationOfCellAfter(edge.Twin, 1);
                        break;
                    }
                }
            }
        }

        static CyclicArray<Edge<T>> PositiveEdgeRotationOfCellAfter(Edge<T> first, int offset)
        {
            MeshCell<T> cell = first.Cell;
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(cell.Edges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(cell.Edges, firstEdgeIndice + offset);
            return edges;
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

        static IEnumerable<Edge<T>> NegativeRotationOfBoundaryEdgesBeginningWith(Edge<T> first, int offset)
        {
            IList<Edge<T>> edges = NegativeEdgeRotationOfCellAfter(first, offset);
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

        static CyclicArray<Edge<T>> NegativeEdgeRotationOfCellAfter(Edge<T> first, int offset)
        {
            Edge<T>[] edgesOfCell = first.Cell.Edges;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(edgesOfCell);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice + offset);
            return edges;
        }
    }

    static class BoundaryEdgeEnumerator<T>
    {
        public static Edge<T> GetFirstEdgeOfBoundaryPositive(MeshCell<T> cornerCell)
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

        public static Edge<T> GetFirstEdgeOfBoundaryNegative(MeshCell<T> cornerCell)
        {
            foreach (var edge in new Convolution<Edge<T>>(cornerCell.Edges))
            {
                if (!edge.Current.IsBoundary && edge.Previous.IsBoundary)
                {
                    return edge.Previous;
                }
            }
            throw new Exception("Cell does not have a boundary.");
        }

        public static IEnumerable<Edge<T>> PositiveRotationOfBoundaryEdgesBeginningWith(Edge<T> first, int offset)
        {
            IList<Edge<T>> edges = PositiveEdgeRotationOfCellAfter(first, offset);
            while (edges !=  null)
            {
                foreach (Edge<T> edge in edges)
                {
                    if (edge.IsBoundary)
                    {
                        yield return edge;
                    }
                    else
                    {
                        edges = PositiveEdgeRotationOfCellAfter(edge.Twin, 1);
                        break;
                    }
                }
            }
        }

        public static CyclicArray<Edge<T>> PositiveEdgeRotationOfCellAfter(Edge<T> first, int offset)
        {
            MeshCell<T> cell = first.Cell;
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(cell.Edges, first);
            if(firstEdgeIndice != -1)
            {
                return new CyclicArray<Edge<T>>(cell.Edges, firstEdgeIndice + offset);
            }
            else
            {
                return null;
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
            return -1;
        }

        public static IEnumerable<Edge<T>> NegativeRotationOfBoundaryEdgesBeginningWith(Edge<T> first, int offset)
        {
            IList<Edge<T>> edges = NegativeEdgeRotationOfCellAfter(first, offset);
            while (edges != null)
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

        public static CyclicArray<Edge<T>> NegativeEdgeRotationOfCellAfter(Edge<T> first, int offset)
        {
            Edge<T>[] edgesOfCell = first.Cell.Edges;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(edgesOfCell);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            if(firstEdgeIndice != -1)
            {
                return new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice + offset);
            }
            else
            {
                return null;
            }
        }
    }
}
