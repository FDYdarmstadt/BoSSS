﻿using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerBoundaryIdentifier<T>
    {
        public PeriodicCornerBoundaryIdentifier(PeriodicMap map)
        {
            CrossingFinder.map = map;
        }

        struct Crossing
        {
            public Edge<T> First;

            public Edge<T> Second;

            public Edge<T> Third;
        }

        class CrossingFinder
        {
            public Stack<Edge<T>> Visited { get; private set; }

            IEnumerator<Edge<T>> inside;

            IEnumerator<Edge<T>> outside;

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

            public bool Increment()
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

        public (Stack<Edge<T>>, Corner) FindEdges(MeshCell<T> cornerCell)
        {
            Edge<T> firstEdge = GetFirstEdgeOfBoundary(cornerCell);
            Crossing firstCrossing = FindFirstCrossing(firstEdge);

            Stack<Edge<T>> edges = default;
            Corner corner = default;
            bool stop = false;
            CrossingFinder first = new CrossingFinder(firstCrossing.First);
            CrossingFinder second = new CrossingFinder(firstCrossing.Second);
            CrossingFinder third = new CrossingFinder(firstCrossing.Third);
            while (!stop)
            {
                if (first.Increment())
                {
                    stop = true;
                    edges = first.Visited;
                    corner = new Corner
                    {
                        FirstEdge = firstCrossing.Third.Twin.BoundaryEdgeNumber,
                        SecondEdge = first.Crossing.First.BoundaryEdgeNumber
                    };
                    
                }
                else if (second.Increment())
                {
                    stop = true;
                    edges = second.Visited;
                    corner = new Corner
                    {
                        FirstEdge = firstCrossing.First.Twin.BoundaryEdgeNumber,
                        SecondEdge = second.Crossing.First.BoundaryEdgeNumber
                    };
                }
                else if (third.Increment())
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
                if (forward.Increment())
                {
                    return forward.Crossing;
                }
                else if (backwards.Increment())
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
            IEnumerable<Edge<T>> edges = PositiveEdgeRotatationOfCellAfter(first, offset);
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

        static CyclicArray<Edge<T>> PositiveEdgeRotatationOfCellAfter(Edge<T> first, int offset)
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
            Edge<T>[] cell = first.Cell.Edges;
            Edge<T>[] reversedEdges = ArrayMethods.GetReverseOrderArray(cell);
            int firstEdgeIndice = FindIndiceOfEdgeInItsCell(reversedEdges, first);
            CyclicArray<Edge<T>> edges = new CyclicArray<Edge<T>>(reversedEdges, firstEdgeIndice + offset);
            return edges;
        }
    }
}
