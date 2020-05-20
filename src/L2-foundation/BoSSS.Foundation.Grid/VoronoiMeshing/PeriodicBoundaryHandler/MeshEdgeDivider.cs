using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class MeshEdgeDivider<T>
    {
        readonly IDMesh<T> mesh;

        static EdgeComparer<T> edgeComparer = new EdgeComparer<T>();

        struct IndiceEdge
        {
            public Edge<T> Edge;
            public int Indice;
        }

        public MeshEdgeDivider(IDMesh<T> mesh)
        {
            this.mesh = mesh;
        }

        public void DivideBoundary(IEnumerable<MeshCell<T>> cells)
        {
            var corners = new Convolution<IndiceEdge>(BoundaryOfEachCellClockwise(cells));
            foreach (Pair<IndiceEdge> corner in corners)
            {
                Divide(corner);
            }
        }

        public void DivideBoundary(MeshCell<T> first, MeshCell<T> last)
        {
            var corners = new Convolution<IndiceEdge>(BoundaryOfEachCellClockwise(first, last));
            foreach (Pair<IndiceEdge> corner in corners)
            {
                Divide(corner);
            }
        }

        static IEnumerable<IndiceEdge> BoundaryOfEachCellClockwise(MeshCell<T> first, MeshCell<T> last)
        {
            Edge<T> lastEdge = BoundaryEdgeEnumerator<T>.GetFirstEdgeOfBoundaryPositive(last);
            foreach (IndiceEdge edge in NegativeRotationOfBoundaryEdgesBeginningWith(first))
            {
                yield return edge; 
                if(edgeComparer.Equals(lastEdge, edge.Edge))
                {
                    break;
                }
            }
        }

        static IEnumerable<IndiceEdge> NegativeRotationOfBoundaryEdgesBeginningWith(MeshCell<T> cell)
        {
            bool boundaryEdgesLeft = true;
            while (boundaryEdgesLeft)
            {
                int first = FindFirstEdgeIndiceClockwise(cell);
                if(first != -1)
                {
                    for (int i = 0; i < cell.Edges.Length; ++i)
                    {
                        int indice = GoRound(first - i, cell.Edges.Length);
                        Edge<T> edge = cell.Edges[GoRound(first - i, cell.Edges.Length)];
                        if (edge.IsBoundary)
                        {
                            yield return new IndiceEdge
                            {
                                Edge = edge,
                                Indice = indice
                            };
                        }
                        else
                        {
                            cell = edge.Twin.Cell;
                            break;
                        }
                    }
                }
                else
                {
                    boundaryEdgesLeft = false;
                }
            }
        }

        static IEnumerable<IndiceEdge> BoundaryOfEachCellClockwise(IEnumerable<MeshCell<T>> cells)
        {
            foreach (MeshCell<T> cell in cells)
            {
                int first = FindFirstEdgeIndiceClockwise(cell);
                for (int i = 0; i < cell.Edges.Length; ++i)
                {
                    int indice = GoRound(first - i, cell.Edges.Length);
                    if(indice < 0)
                    {
                        break;
                    }
                    else
                    {
                        Edge<T> edge = cell.Edges[GoRound(first - i, cell.Edges.Length)];
                        if (edge.IsBoundary)
                        {
                            yield return new IndiceEdge
                            {
                                Edge = edge,
                                Indice = indice
                            };
                        }
                        else
                        {
                            break;
                        }
                    }
                }
            }
        }

        static int FindFirstEdgeIndiceClockwise(MeshCell<T> cell)
        {
            int indice = 0; 
            foreach(var followingEdges in new Convolution<Edge<T>>(cell.Edges))
            {
                if(followingEdges.Previous.IsBoundary && !followingEdges.Current.IsBoundary)
                {
                    return indice % cell.Edges.Length;
                }
                ++indice;
            }
            return -1;
        }

        void Divide(Pair<IndiceEdge> corner)
        {
            if(corner.Current.Edge.End == corner.Previous.Edge.Start)
            {
                Vertex clone = MeshElementCloner.Clone(corner.Current.Edge.End);
                mesh.AddVertex(clone);
                CircleAndInsert(clone, corner.Previous, corner.Current);
            }
            else
            {
                Vertex clone = MeshElementCloner.Clone(corner.Current.Edge.End);
                mesh.AddVertex(clone);
                InsertAtEnd(clone, corner.Current);

                Vertex clone2 = MeshElementCloner.Clone(corner.Previous.Edge.Start);
                mesh.AddVertex(clone2);
                InsertAtStart(clone2, corner.Previous);
            }
        }

        static void CircleAndInsert(Vertex clone, IndiceEdge from, IndiceEdge to)
        {
            IndiceEdge current = from;
            while (!edgeComparer.Equals(current.Edge, to.Edge.Twin))
            {
                current.Edge.Start = clone;
                current.Edge.Cell.Vertices[current.Indice] = clone;

                int nextEdgeIndice = GoRound(current.Indice - 1, current.Edge.Cell.Edges.Length);
                Edge<T> next = current.Edge.Cell.Edges[nextEdgeIndice];
                next.End = clone;

                current = new IndiceEdge
                {
                    Edge = next.Twin,
                    Indice = FindIndiceInCell(next.Twin)
                };
            }
        }

        static int FindIndiceInCell(Edge<T> edge)
        {
            for(int i = 0; i < edge.Cell.Vertices.Length; ++i)
            {
                if (edge.Start.ID == edge.Cell.Vertices[i].ID)
                    return i;
            }
            return -1;
        }

        static int GoRound(int i, int length)
        {
            return (i + length) % length;
        }
        
        static void InsertAtStart(Vertex clone, IndiceEdge edge)
        {
            edge.Edge.Start = clone;
            edge.Edge.Cell.Vertices[edge.Indice] = clone;
            edge.Edge.Cell.Edges[GoRound(edge.Indice - 1, edge.Edge.Cell.Edges.Length)]. End = clone;
        }

        static void InsertAtEnd(Vertex clone, IndiceEdge edge)
        {
            edge.Edge.End = clone;
            edge.Edge.Cell.Vertices[GoRound(edge.Indice + 1, edge.Edge.Cell.Vertices.Length)] = clone;
            edge.Edge.Cell.Edges[GoRound(edge.Indice + 1, edge.Edge.Cell.Edges.Length)].Start = clone;
        }
    }
}
