using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshIntersecter<T>
        where T: IMesherNode, new()
    {
        IDMesh<T> mesh;

        static EdgeComparer<T> edgeComparer = new EdgeComparer<T>();

        public class AfterCutEdgeEnumerator : IEnumerator<Edge<T>>
        {
            Edge<T>[] edges;
            Edge<T> current;
            int startIndex;
            int currentIndex;
            int counter;

            int length; //Changes after Reset

            public AfterCutEdgeEnumerator(Edge<T>[] Edges, Edge<T> StartEdge)
            {
                edges = Edges;
                counter = -1;
                length = Edges.Length - 2;
                for (int i = 0; i < Edges.Length; ++i)
                {
                    if (edgeComparer.Equals(StartEdge, Edges[i]))
                    {
                        startIndex = i;
                        currentIndex = i;
                        break;
                    }
                }
            }

            public Edge<T> Current => current;

            object IEnumerator.Current => Current;

            public void Dispose() { }

            public bool MoveNext()
            {

                if (++counter >= length)
                {
                    return false;
                }
                else
                {
                    currentIndex = (currentIndex + 1) % edges.Length;
                    current = edges[currentIndex];
                    return true;
                }
            }

            public void Reset()
            {
                counter = -1;
                length = edges.Length;
                currentIndex = startIndex;
            }
        }

        public MeshIntersecter(IDMesh<T> idMesh)
        {
            mesh = idMesh;
        }

        public IEnumerator<Edge<T>> GetFirstEnumerator(MeshCell<T> first)
        {
            AfterCutEdgeEnumerator enumerator = new AfterCutEdgeEnumerator(first.Edges, first.Edges[1]);
            enumerator.Reset();
            return enumerator;
        }

        public IEnumerator<Edge<T>> GetNeighborFromEdgeNeighbor(Edge<T> edge)
        {
            MeshCell<T> newCell = MeshMethods.GetNeighbour(edge);
            AfterCutEdgeEnumerator ridgeEnum = new AfterCutEdgeEnumerator(newCell.Edges, edge);

            return ridgeEnum;
        }

        public Edge<T> Subdivide(Edge<T> edge, List<BoundaryLine> lines, double alpha, CyclicInterval boundaryCount)
        {
            MeshCell<T> cell = edge.Cell;
            //Divide Ridge and update Ridge Arrays
            //-------------------------------------
            Vertex newVertex = DivideEdge(edge, alpha, out Edge<T> newRidge);
            edge.Twin.Cell.IntersectionVertex = newVertex.ID;
            //cell.IntersectionVertex = newVertex.ID;

            //Divide this cell
            //================================================================
            //NewVertices
            Vertex[] verticesOfNewRidgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewRidgeBoundary[0] = newVertex;
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1] = mesh.Vertices[cell.IntersectionVertex];
            //Add Vertices of lines
            for (int i = 1; i < verticesOfNewRidgeBoundary.Length - 1; ++i)
            {
                verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i] = lines[i - 1].End;
                int ID = mesh.AddVertex(verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newEdges;
            Edge<T>[] newNeighborEdges;
            MeshCell<T> newCell = new MeshCell<T> { Node = new T() };
            newCell.Node.Position = cell.Node.Position;
            mesh.AddCell(newCell);
            MeshMethods.CreateBoundaryEdge(verticesOfNewRidgeBoundary, cell, newCell, out newEdges, out newNeighborEdges, boundaryCount);
            //Link Ridges to old neighbors
            MeshMethods.InsertEdgesAndVertices(newEdges, newNeighborEdges);

            //dOnE, DoNe!
            return edge;
        }

        public Edge<T> FirstCut(Edge<T> edge, double alpha)
        {
            //Divide Ridge and update Ridge Arrays
            //-------------------------------------
            Vertex newVertex = DivideEdge(edge, alpha, out Edge<T> newRidge);
            edge.Twin.Cell.IntersectionVertex = newVertex.ID;

            //Find Intersection and insert Ridge
            edge.Cell.IntersectionVertex = newVertex.ID;

            return edge;
        }

        public Vertex DivideEdge(Edge<T> edge, double alpha, out Edge<T> newEdge)
        {
            Vertex newVertex = MeshMethods.DivideEdge(edge, alpha, out newEdge);
            mesh.AddVertex(newVertex);
            return newVertex;
        }

        public void CloseMesh(List<BoundaryLine> lines, Edge<T> firstCutEdge, CyclicInterval boundaryCount)
        {
            MeshCell<T> cell = firstCutEdge.Cell;
            //Divide this cell
            //================================================================
            //NewVertices
            Vertex[] verticesOfNewRidgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewRidgeBoundary[0] = firstCutEdge.End;
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1] = mesh.Vertices[cell.IntersectionVertex];
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 2] = lines[0].Start;
            int ID = mesh.AddVertex(lines[0].Start);

            //Add Vertices of lines
            for (int i = 2; i < verticesOfNewRidgeBoundary.Length - 1; ++i)
            {
                verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i] = lines[i - 1].End;
                ID = mesh.AddVertex(verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newEdges;
            Edge<T>[] newNeighborEdges;
            MeshCell<T> newCell = new MeshCell<T>();
            mesh.AddCell(newCell);
            MeshMethods.CreateBoundaryEdge(verticesOfNewRidgeBoundary, cell, newCell, out newEdges, out newNeighborEdges, boundaryCount);
            MeshMethods.InsertEdgesAndVertices(newEdges, newNeighborEdges);
        }

        public IEnumerator<Edge<T>> GetConnectedRidgeEnum(Edge<T> edge)
        {
            List<Edge<T>> outgoingEdges = new List<Edge<T>>();
            bool newNeighbor = true;
            Edge<T> worker = edge;
            while (newNeighbor)
            {
                worker = worker.Twin;
                Edge<T>[] workerRidges = worker.Cell.Edges;
                outgoingEdges.Add(worker);
                for (int i = 0; i < workerRidges.Length; ++i)
                {
                    if (edgeComparer.Equals(workerRidges[i], worker))
                    {
                        if (i == 0)
                            worker = workerRidges.Last();
                        else
                            worker = workerRidges[i - 1];
                        if (edgeComparer.Equals(worker, edge))
                        {
                            newNeighbor = false;
                        }
                        break;
                    }
                }
            }
            return new ArrayEnum<Edge<T>>(outgoingEdges);
        }

        public void VertexCut(Edge<T> edge, double alphaCut)
        {
            MeshCell<T> cell = edge.Cell;
            Vertex newOldVertex = edge.End;
            cell.IntersectionVertex = newOldVertex.ID;
            cell = edge.Twin.Cell;
            cell.IntersectionVertex = newOldVertex.ID;
        }

        public IEnumerator<Edge<T>> getNeighborFromLineDirection(Edge<T> edge, BoundaryLine line)
        {
            IEnumerator<Edge<T>> outgoingEdges = GetConnectedRidgeEnum(edge);
            Edge<T> A = null;
            Edge<T> B = null;
            if (outgoingEdges.MoveNext())
            {
                A = outgoingEdges.Current;
            }
            else
            {
                throw new Exception();
            }
            bool found = false;
            while (!found && outgoingEdges.MoveNext())
            {
                B = outgoingEdges.Current;
                if (IsBetween(A, line, B))
                {
                    found = true;
                }
                else
                {
                    A = B;
                }
            }

            return new AfterCutEdgeEnumerator(A.Cell.Edges, A);

            //Check if in positive rotation a, c, b order
            bool IsBetween(Edge<T> a, BoundaryLine b, Edge<T> c)
            {
                Vector A1 = a.End.Position - a.Start.Position;
                Vector C1 = b.End.Position - a.Start.Position;
                Vector B1 = c.End.Position - a.Start.Position;

                double crossAB = A1.CrossProduct2D(B1);
                double crossCB = C1.CrossProduct2D(B1);
                double crossAC = A1.CrossProduct2D(C1);
                if (crossAC > 0)
                {
                    if (crossCB > 0 || crossAB < 0)
                    {
                        return true;
                    }
                    return false;
                }
                else
                {
                    if (crossCB > 0 && crossAB < 0)
                    {
                        return true;
                    }
                    return false;
                }
            }
        }

        class MultiNeighRidgesAfterCutEnum : IEnumerator<Edge<T>>
        {
            IList<Edge<T>> enumEdges;
            IList<int> block;
            int pointer;
            int blockPointer;
            bool filter;
            public MultiNeighRidgesAfterCutEnum(IList<Edge<T>> edges, IList<int> blockInFirstRun)
            {
                enumEdges = edges;
                block = blockInFirstRun;
                pointer = -1;
                blockPointer = 0;
                filter = true;
            }
            public Edge<T> Current => enumEdges[pointer];

            object IEnumerator.Current => Current;

            public void Dispose()
            {
                throw new NotImplementedException();
            }

            public bool MoveNext()
            {
                ++pointer;
                if (filter)
                {
                    while (blockPointer < block.Count && pointer == block[blockPointer])
                    {
                        ++pointer;
                        ++blockPointer;
                    }
                }
                return (pointer < enumEdges.Count);
            }

            public void Reset()
            {
                pointer = -1;
                filter = false;
            }
        }

        public Edge<T> AddLineSegment(Edge<T> edge, double alpha, CyclicInterval boundaryCount)
        {
            DivideEdge(edge, alpha, out Edge<T> newRidge);
            edge.IsBoundary = true;
            edge.Twin.IsBoundary = true;
            edge.BoundaryEdgeNumber = boundaryCount.Current();
            edge.Twin.BoundaryEdgeNumber = boundaryCount.Current();
            edge.Cell.IntersectionVertex = edge.End.ID;
            edge.Twin.Cell.IntersectionVertex = edge.End.ID;
            return edge;
        }

        public Edge<T> SubdivideWithoutNewVertex(Edge<T> edge, List<BoundaryLine> lines, CyclicInterval boundaryCount)
        {
            MeshCell<T> cell = edge.Cell;
            Vertex cutVertex = edge.End;
            IEnumerator<Edge<T>> neighEdges = GetConnectedRidgeEnum(edge);
            while (neighEdges.MoveNext())
            {
                MeshCell<T> neighbor = neighEdges.Current.Cell;
                if (neighbor.ID != cell.ID)
                    neighbor.IntersectionVertex = cutVertex.ID;
            }
            //Divide this cell
            //================================================================
            Vertex[] verticesOfNewEdgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewEdgeBoundary[0] = cutVertex;
            verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1] = mesh.Vertices[cell.IntersectionVertex];
            //Add Vertices of lines
            for (int i = 1; i < verticesOfNewEdgeBoundary.Length - 1; ++i)
            {
                verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1 - i] = lines[i - 1].End;
                mesh.AddVertex(verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newEdges;
            Edge<T>[] newNeighborEdges;
            MeshCell<T> newCell = new MeshCell<T>();
            mesh.AddCell(newCell);
            MeshMethods.CreateBoundaryEdge(verticesOfNewEdgeBoundary, cell, newCell, out newEdges, out newNeighborEdges, boundaryCount);
            //Link Ridges to old neighbors
            MeshMethods.InsertEdgesAndVertices(newEdges, newNeighborEdges);


            //dOnE, DoNe!
            return edge;
        }

        public void AddEdge(Edge<T> edge, CyclicInterval boundaryCount)
        {
            edge.IsBoundary = true;
            edge.Twin.IsBoundary = true;
            edge.BoundaryEdgeNumber = boundaryCount.Current();
            edge.Twin.BoundaryEdgeNumber = boundaryCount.Current();
            edge.Cell.IntersectionVertex = edge.End.ID;
            edge.Twin.Cell.IntersectionVertex = edge.End.ID;
        }
    }
}
