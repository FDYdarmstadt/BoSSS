using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System.Collections;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class MeshMethods
    {
        public static MeshCell<T> GetNeighbour<T>(Edge<T> ridge)
        {
            return ridge.Twin.Cell;
        }

        public static void InsertEdgesAndVertices<T>(params Edge<T>[] additionalEdges)
        {
            MeshCell<T> cell = additionalEdges[0].Cell;
            int countAdditionalRidges = additionalEdges.Length;
            Edge<T>[] newRidges = new Edge<T>[cell.Edges.Length + countAdditionalRidges];
            bool notInsertedYet = true;
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                if (notInsertedYet)
                {
                    newRidges[i] = cell.Edges[i];
                }
                else
                {
                    newRidges[i + countAdditionalRidges] = cell.Edges[i];
                }
                if (notInsertedYet && (additionalEdges[0].Start.ID == cell.Edges[i].End.ID))
                {
                    for (int k = 0; k < countAdditionalRidges; ++k)
                    {
                        newRidges[i + k + 1] = additionalEdges[k];
                    }

                    notInsertedYet = false;
                }
            }
            cell.Edges = newRidges;

            Vertex[] newVertices = new Vertex[cell.Vertices.Length + countAdditionalRidges];
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                newVertices[i] = newRidges[i].Start;
            }
            cell.Vertices = newVertices;
        }

        public static void InsertEdgesAndVertices<T>(Edge<T>[] newEdge, Edge<T>[] newNeighborEdges)
        {
            MeshCell<T> cell = newEdge[0].Cell;
            MeshCell<T> emptyNeighCell = newNeighborEdges[0].Cell;
            Edge<T>[] oldRidges = cell.Edges;

            List<Edge<T>> cellRidges = null;
            List<Edge<T>> emptyNeighCellRidges = null;
            List<Edge<T>> workerA = new List<Edge<T>>(newEdge.Length);
            List<Edge<T>> workerB = new List<Edge<T>>(newEdge.Length);
            bool workerAIsActive = true;
            List<Edge<T>> tmp = workerA;

            //Add new Ridges
            for (int i = 0; i < oldRidges.Length; ++i)
            {
                Edge<T> activeR = oldRidges[i];
                if (activeR.Start.ID == newEdge[0].Start.ID)
                {
                    cellRidges = tmp;
                    for (int j = 0; j < newEdge.Length; ++j)
                    {
                        cellRidges.Add(newEdge[j]);
                    }
                    tmp = workerAIsActive ? workerB : workerA;
                    workerAIsActive = !workerAIsActive;
                }
                if (activeR.Start.ID == newNeighborEdges[0].Start.ID)
                {
                    emptyNeighCellRidges = tmp;
                    for (int j = 0; j < newNeighborEdges.Length; ++j)
                    {
                        emptyNeighCellRidges.Add(newNeighborEdges[j]);
                    }

                    tmp = workerAIsActive ? workerB : workerA;
                    workerAIsActive = !workerAIsActive;
                }
                tmp.Add(activeR);
            }

            cell.Edges = cellRidges.ToArray();
            emptyNeighCell.Edges = emptyNeighCellRidges.ToArray();

            //Update AllRidges
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                cell.Edges[i].Cell = cell;
            }
            for (int i = 0; i < emptyNeighCell.Edges.Length; ++i)
            {
                emptyNeighCell.Edges[i].Cell = emptyNeighCell;
            }

            //Vertices
            Vertex[] newVertices = new Vertex[cell.Edges.Length];
            for (int i = 0; i < cell.Edges.Length; ++i)
            {
                newVertices[i] = cell.Edges[i].Start;
            }
            cell.Vertices = newVertices;

            Vertex[] newNeighVertices = new Vertex[emptyNeighCell.Edges.Length];
            for (int i = 0; i < emptyNeighCell.Edges.Length; ++i)
            {
                newNeighVertices[i] = emptyNeighCell.Edges[i].Start;
            }
            emptyNeighCell.Vertices = newNeighVertices;

        }

        public static void CreateBoundaryEdge<T>(
            Vertex[] vertices,
            MeshCell<T> cell,
            MeshCell<T> neighborCell,
            out Edge<T>[] ridges,
            out Edge<T>[] twinEdges,
            CyclicInterval boundaryCount)
        {
            int count = vertices.Length - 1;
            ridges = new Edge<T>[count];
            twinEdges = new Edge<T>[count];
            for (int i = 0; i < count; ++i)
            {
                Edge<T> ridge = new Edge<T>
                {
                    Start = vertices[i],
                    End = vertices[i + 1],
                    Cell = cell,
                    IsBoundary = true,
                    BoundaryEdgeNumber = boundaryCount.Current(),
                };
                Edge<T> twinRidge = new Edge<T>
                {
                    Start = vertices[i + 1],
                    End = vertices[i],
                    Twin = ridge,
                    Cell = neighborCell,
                    IsBoundary = true,
                    BoundaryEdgeNumber = boundaryCount.Current(),
                };
                ridge.Twin = twinRidge;
                ridges[i] = ridge;
                twinEdges[count - 1 - i] = twinRidge;
                boundaryCount.Previous();
            }
        }
    }
}
