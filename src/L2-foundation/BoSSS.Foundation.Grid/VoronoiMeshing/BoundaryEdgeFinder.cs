using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryEdgeFinder<T>
    {
        InsideCellEnumerator<T> cells;

        public BoundaryEdgeFinder(Mesh<T> mesh)
        {
            cells = new InsideCellEnumerator<T>(mesh);
        }

        public IEnumerable<Edge<T>> Edges(Vector start, int firstCellNodeIndice)
        {
            cells.SetFirstCell(start, firstCellNodeIndice);
            MeshCell<T> firstCell = cells.GetFirstCell();
            return IterativeYieldBoundaryCells(firstCell);
        }

        static IEnumerable<Edge<T>> IterativeYieldBoundaryCells(MeshCell<T> cell)
        {
            Edge<T> currentEdge = FindFirstBoundaryEdge(cell);
            yield return currentEdge;

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
                        if(currentEdge.Start.ID != startID)
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
