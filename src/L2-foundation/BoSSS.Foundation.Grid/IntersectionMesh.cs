using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Foundation.Grid.Voronoi
{
    class ArrayEnum<T> : IEnumerator<T>
    {
        int pointer;
        IList<T> arr;
        public ArrayEnum(IList<T> Arr)
        {
            arr = Arr;
            pointer = -1;
        }

        public T Current => arr[pointer];

        object IEnumerator.Current => Current;

        public void Dispose()
        {
        }

        public bool MoveNext()
        {
            return (++pointer < arr.Count);
        }

        public void Reset()
        {
            pointer = -1;
        }
    }

    class Line
    {
        public static Line[] ToLines(Vector[] polygon)
        {
            Line[] lines = new Line[polygon.Length];
            Line line;
            for (int i = 0; i < polygon.Length - 1; ++i)
            {
                line = new Line { start = new Vertex(), end = new Vertex() };
                line.start.Position = polygon[i];
                line.end.Position = polygon[i + 1];
                lines[i] = line;
            }
            line = new Line { start = new Vertex(), end = new Vertex() };
            line.start.Position = polygon[polygon.Length - 1];
            line.end.Position = polygon[0];
            lines[lines.Length - 1] = line;

            return lines;
        }

        public static IEnumerator<Line> GetEnumerator(Vector[] polygon)
        {
            Line[] lines = ToLines(polygon);
            return new ArrayEnum<Line>(lines);
        }

        public Vertex start { get; set; }
        public Vertex end { get; set; }
    }

    class IntersectionMesh<T> : BoundaryMesh<T>, IIntersectableMesh<Cell<T>, Edge<T>, Line>
    {
        bool cutIsFresh = false;

        public bool CutIsFresh {
            set { cutIsFresh = value; } 
        }

        public IntersectionMesh(IIdMesh<T> Mesh, int firstCell_NodeIndice) 
            : base(Mesh)
        {
            insideCell = Cells[firstCell_NodeIndice];
        }

        public IntersectionMesh(IIdMesh<T> Mesh) 
            : base(Mesh)
        {
            insideCell = null;
        }

        static EdgeComparer<T> ridgeComparer = new EdgeComparer<T>();

        public override IReadOnlyList<Cell<T>> GetCells() {
            if (cutIsFresh)
            {
                DetermineInsideCells();
                cutIsFresh = false;
            }
            return base.GetCells();
        }

        static double accuracy = 1e-10;
        public (Cell<T>, IEnumerator<Edge<T>>) GetFirst(Line boundaryLine)
        {
            //Find cell that contains boundaryLine.Start;
            bool foundFirstCell = false;
            if (insideCell == null)
            {
                //SetFirst Cell: any random cell. Influences runtime, though
                insideCell = Cells[0];
                foundFirstCell = true;
            }

            //Check if boundaryLine.Start is still in cell, else search neighborhood
            foreach (Cell<T> cell in CellsOnSameSideOfBoundary_Iterative(insideCell))
            {
                Vector[] verts = Array.ConvertAll(cell.Vertices, item => (Vector)item);
                //At this point, every cell is convex!
                bool isInside = PolygonTesselation.PointInConvexPolygon(verts, (Vector)boundaryLine.start);
                if (isInside)
                {
                    foundFirstCell = true;
                    insideCell = cell;
                    break;
                }
            }
            if (foundFirstCell)
            {
                AfterCutRidgeEnumerator enumerator = new AfterCutRidgeEnumerator(insideCell.Edges, insideCell.Edges[1]);
                enumerator.Reset();
                return (insideCell, enumerator);
            }
            else
            {
                throw new Exception("First cell could not be found: boundaryLine.start not inside a cell");
            }
        }
        //Return end of ridge if parallel and overlapping.
        public bool intersect(Edge<T> edge, Line line, ref IntersectionCase intersectionCase, out double alpha)
        {
            double alpha2;
            Vector vector;
            bool notParallel = PolygonClipping.ComputeIntersection(
                edge.Start.Position,
                edge.End.Position,
                line.start.Position,
                line.end.Position,
                out alpha,
                out alpha2,
                out vector);
            if (notParallel)
            {
                if (accuracy <= alpha && 1 >= alpha)
                {
                    if (accuracy <= alpha2 && 1 >= alpha2)
                    {
                        if (alpha > 1 - accuracy)
                        {
                            alpha = 1;
                            intersectionCase = IntersectionCase.EndOfRidge;
                            if (alpha2 > 1 - accuracy)
                            {
                                intersectionCase = IntersectionCase.EndOfRidgeAndLine;
                            }
                            return true;
                        }
                        else
                        {
                            if (alpha2 > 1 - accuracy)
                            {
                                intersectionCase = IntersectionCase.EndOfLine;
                            }
                            else
                            {
                                intersectionCase = IntersectionCase.InMiddle;
                            }
                            return true;
                        }
                    }
                }
                return false;
            }
            else
            {
                if (Math.Abs(alpha) <= accuracy)
                {
                    double absRidge = (edge.Start.Position - edge.End.Position).L2Norm();
                    double absLine = (edge.Start.Position - line.end.Position).L2Norm();
                    //Only when in same direction!
                    double add = (edge.Start.Position - edge.End.Position + edge.Start.Position - line.end.Position).L2Norm();
                    if (add < absLine + absRidge)
                    {
                        return false;
                    }

                    alpha = absLine / absRidge;
                    if (alpha >= 1 - accuracy && alpha <= 1)
                    {
                        alpha = 1;
                        intersectionCase = IntersectionCase.EndOfRidgeAndLine;
                        return true;
                    }
                    if (absLine / absRidge < 1)
                    {
                        intersectionCase = IntersectionCase.EndOfLine;
                        return true;
                    }
                    alpha = 1;
                    intersectionCase = IntersectionCase.EndOfRidge;
                    return true;
                }
                return false;
            }
        }

        public (Cell<T>, IEnumerator<Edge<T>>) getNeighborFromEdgeNeighbor(Edge<T> edge)
        {
            Cell<T> newCell = getNeighbour(edge);
            AfterCutRidgeEnumerator ridgeEnum = new AfterCutRidgeEnumerator(newCell.Edges, edge);

            return (newCell, ridgeEnum);
        }

        public Edge<T> Subdivide(Edge<T> edge, List<Line> lines, double alpha)
        {
            Cell<T> cell = edge.Cell;
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
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1] = Vertices[cell.IntersectionVertex];
            //Add Vertices of lines
            for (int i = 1; i < verticesOfNewRidgeBoundary.Length - 1; ++i)
            {
                verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i] = lines[i - 1].end;
                int ID = AddVertex(verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newEdges;
            Edge<T>[] newNeighborEdges;
            Cell<T> newCell = new Cell<T> { Position = cell.Position };
            AddCell(newCell);
            CreateEdge(verticesOfNewRidgeBoundary, cell, newCell, out newEdges, out newNeighborEdges);
            //Link Ridges to old neighbors
            InsertEdgesAndVertices(newEdges, newNeighborEdges);


            //dOnE, DoNe!
            return edge;
        }

        class AfterCutRidgeEnumerator : IEnumerator<Edge<T>>
        {
            Edge<T>[] edges;
            Edge<T> current;
            int startIndex;
            int currentIndex;
            int counter;

            int length; //Changes after Reset

            public AfterCutRidgeEnumerator(Edge<T>[] Edges, Edge<T> StartEdge)
            {
                edges = Edges;
                counter = -1;
                length = Edges.Length - 2;
                for (int i = 0; i < Edges.Length; ++i)
                {
                    if (ridgeComparer.Equals(StartEdge, Edges[i]))
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

        public void CloseMesh(List<Line> lines, Edge<T> firstCutEdge)
        {
            Cell<T> cell = firstCutEdge.Cell;
            //Divide this cell
            //================================================================
            //NewVertices
            Vertex[] verticesOfNewRidgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewRidgeBoundary[0] = firstCutEdge.End;
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1] = Vertices[cell.IntersectionVertex];
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 2] = lines[0].start;
            int ID = AddVertex(lines[0].start);

            //Add Vertices of lines
            for (int i = 2; i < verticesOfNewRidgeBoundary.Length - 1; ++i)
            {
                verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i] = lines[i - 1].end;
                ID = AddVertex(verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newRidges;
            Edge<T>[] newNeighborRidges;
            Cell<T> newCell = new Cell<T>();
            AddCell(newCell);
            CreateEdge(verticesOfNewRidgeBoundary, cell, newCell, out newRidges, out newNeighborRidges);
            InsertEdgesAndVertices(newRidges, newNeighborRidges);
        }

        public IEnumerator<Edge<T>> getConnectedRidgeEnum(Edge<T> edge)
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
                    if (ridgeComparer.Equals(workerRidges[i], worker))
                    {
                        if (i == 0)
                            worker = workerRidges.Last();
                        else
                            worker = workerRidges[i - 1];
                        if (ridgeComparer.Equals(worker, edge))
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
            Cell<T> cell = edge.Cell;
            Vertex newOldVertex = edge.End;
            cell.IntersectionVertex = newOldVertex.ID;
            cell = edge.Twin.Cell;
            cell.IntersectionVertex = newOldVertex.ID;
        }

        public IEnumerator<Edge<T>> getNeighborFromLineDirection(Edge<T> edge, Line line)
        {
            IEnumerator<Edge<T>> outgoingEdges = getConnectedRidgeEnum(edge);
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

            return new AfterCutRidgeEnumerator(A.Cell.Edges, A);

            //Check if in positive rotation a, c, b order
            bool IsBetween(Edge<T> a, Line b, Edge<T> c)
            {
                Vector A1 = a.End.Position - a.Start.Position;
                Vector C1 = b.end.Position - a.Start.Position;
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

        public Edge<T> AddLineSegment(Edge<T> edge, double alpha)
        {
            DivideEdge(edge, alpha, out Edge<T> newRidge);
            edge.IsBoundary = true;
            edge.Twin.IsBoundary = true;
            edge.Cell.IntersectionVertex = edge.End.ID;
            edge.Twin.Cell.IntersectionVertex = edge.End.ID;
            return edge;
        }

        public Edge<T> SubdivideWithoutNewVertex(Edge<T> edge, List<Line> lines)
        {
            Cell<T> cell = edge.Cell;
            Vertex cutVertex = edge.End;
            IEnumerator<Edge<T>> neighEdges = getConnectedRidgeEnum(edge);
            while (neighEdges.MoveNext())
            {
                Cell<T> neighbor = neighEdges.Current.Cell;
                if (neighbor.ID != cell.ID)
                    neighbor.IntersectionVertex = cutVertex.ID;
            }
            //Divide this cell
            //================================================================
            Vertex[] verticesOfNewEdgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewEdgeBoundary[0] = cutVertex;
            verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1] = Vertices[cell.IntersectionVertex];
            //Add Vertices of lines
            for (int i = 1; i < verticesOfNewEdgeBoundary.Length - 1; ++i)
            {
                verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1 - i] = lines[i - 1].end;
                AddVertex(verticesOfNewEdgeBoundary[verticesOfNewEdgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Edge<T>[] newEdges;
            Edge<T>[] newNeighborEdges;
            Cell<T> newCell = new Cell<T>();
            AddCell(newCell);
            CreateEdge(verticesOfNewEdgeBoundary, cell, newCell, out newEdges, out newNeighborEdges);
            //Link Ridges to old neighbors
            InsertEdgesAndVertices(newEdges, newNeighborEdges);


            //dOnE, DoNe!
            return edge;
        }

        public void AddEdge(Edge<T> edge)
        {
            edge.IsBoundary = true;
            edge.Twin.IsBoundary = true;
            edge.Cell.IntersectionVertex = edge.End.ID;
            edge.Twin.Cell.IntersectionVertex = edge.End.ID;
        }
    }
}
