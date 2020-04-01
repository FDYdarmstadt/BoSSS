using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryElementEnumerator<T>
    {
        InsideCellEnumerator<T> cells;

        public BoundaryElementEnumerator(Domain<T> mesh)
        {
            cells = new InsideCellEnumerator<T>(mesh);
        }

        public IEnumerable<MeshCell<T>> CycleCells()
        {
            IEnumerator<Edge<T>> edges = CycleEdges().GetEnumerator();
            if (edges.MoveNext())
            {
                Edge<T> edge = edges.Current;
                yield return edge.Cell;

                int previousCellID = edge.Cell.ID;
                while (edges.MoveNext())
                {
                    edge = edges.Current;
                    if (edge.Cell.ID != previousCellID)
                    {
                        yield return edge.Cell;
                        previousCellID = edge.Cell.ID;
                    }
                }
            }
        }

        public IEnumerable<Edge<T>> CycleEdges()
        {
            MeshCell<T> firstCell = cells.GetFirstCell();
            return FollowBoundaryEdges(firstCell);
        }

        public static IEnumerable<Edge<T>> FollowBoundaryEdges(MeshCell<T> cell)
        {
            Edge<T> firstEdge = FindFirstBoundaryEdge(cell);
            yield return firstEdge;
            foreach (Edge<T> edge in FollowBoundaryEdges(firstEdge))
            {
                yield return edge;
            }

           
        }

        public static IEnumerable<Edge<T>> FollowBoundaryEdges(Edge<T> edge)
        {
            Edge<T> currentEdge = edge;
            bool abort = false;
            int startID = currentEdge.Start.ID;
            do
            {
                CyclicArray<Edge<T>> edges = CycleThroughEdgesAfterEdge(currentEdge);
                for (int i = 0; i < edges.Length; ++i)
                {
                    currentEdge = edges[i];
                    if (currentEdge.IsBoundary)
                    {
                        if (currentEdge.Start.ID != startID)
                        {
                            yield return currentEdge;
                        }
                        else
                        {
                            abort = true;
                            break;
                        }
                    }
                    else
                    {
                        currentEdge = currentEdge.Twin;
                        break;
                    }
                }
            }
            while (!abort);
        }

        static Edge<T> FindFirstBoundaryEdge(MeshCell<T> cell)
        {
            foreach (Edge<T> edge in cell.Edges)
            {
                if (edge.IsBoundary)
                {
                    return edge;
                }
            }
            throw new Exception("Cell does not neighbor boundary.");
        }

        static CyclicArray<Edge<T>> CycleThroughEdgesAfterEdge(Edge<T> edge)
        {
            Edge<T>[] edges = edge.Cell.Edges;
            CyclicArray<Edge<T>> reorderedEdges = new CyclicArray<Edge<T>>(edges);
            for (int i = 0; i < edges.Length; ++i)
            {
                if (edge.End.ID == edges[i].Start.ID)
                {
                    reorderedEdges.Start = i;
                    return reorderedEdges;
                }
            }
            throw new Exception("Edge not found in edges.");
        }
    }
}
