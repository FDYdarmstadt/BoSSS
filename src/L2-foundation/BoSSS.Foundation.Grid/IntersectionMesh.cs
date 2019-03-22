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

    class IntersectionMesh : Mesh, IIntersectableMesh<Cell, Ridge, Line>
    {
        public IntersectionMesh(IIdMesh Mesh, int firstCell_NodeIndice) : base(Mesh)
        {
            FirstCell = Cells[firstCell_NodeIndice];
        }

        public IntersectionMesh(IIdMesh Mesh) : base(Mesh)
        {
            FirstCell = null;
        }

        static RidgeComparer ridgeComparer = new RidgeComparer();

        Cell FirstCell;

        public (Cell, IEnumerator<Ridge>) getFirst(Line boundaryLine)
        {
            //Find cell that contains boundaryLine.Start;
            bool foundFirstCell = false;
            if(FirstCell == null)
            {
                //SetFirst Cell: any random cell. Influences runtime, though
                FirstCell = Cells[0];
            }
            else
            {
                //Check if boundaryLine.Start is still in cell, else search neighborhood
                foreach(Cell cell in ConnectedCells_Iterative(FirstCell))
                {
                    Vector[] verts = Array.ConvertAll(FirstCell.Vertices, item => (Vector)item);
                    //At this point, every cell is convex!
                    bool isInside = PolygonTesselation.PointInConvexPolygon(verts, (Vector)boundaryLine.start);
                    if (isInside)
                    {
                        foundFirstCell = true;
                        FirstCell = cell;
                        break;
                    }
                }
            }
            if (foundFirstCell)
            {
                AfterCutRidgeEnumerator enumerator = new AfterCutRidgeEnumerator(FirstCell.Ridges, FirstCell.Ridges[1]);
                enumerator.Reset();
                return (FirstCell, enumerator);
            }
            else
            {
                throw new Exception("First cell could not be found: boundaryLine.start not inside a cell");
            }
        }

        public IEnumerable<Cell> GetInsideCells()
        {
            return ConnectedCells_Iterative(FirstCell);
        }

        public AggregationGrid ToAggregationGrid()
        {
            return ToAggregationGrid(FirstCell);
        }

        static double accuracy = 1e-10;

        //Return end of ridge if parallel and overlapping.
        public bool intersect(Ridge ridge, Line line, ref IntersectionCase intersectionCase, out double alpha)
        {
            double alpha2;
            Vector vector;
            bool notParallel = PolygonClipping.ComputeIntersection(
                ridge.Start.Position,
                ridge.End.Position,
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
                            intersectionCase = IntersectionCase.IntersectionIsEndOfRidge;
                            if (alpha2 > 1 - accuracy)
                            {
                                intersectionCase = IntersectionCase.IntersectionIsEndOfRidgeAndLine;
                            }
                            return true;
                        }
                        else
                        {
                            if (alpha2 > 1 - accuracy)
                            {
                                intersectionCase = IntersectionCase.IntersectionIsEndOfLine;
                            }
                            else
                            {
                                intersectionCase = IntersectionCase.IntersectionInMiddle;
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
                    double absRidge = (ridge.Start.Position - ridge.End.Position).L2Norm();
                    double absLine = (ridge.Start.Position - line.end.Position).L2Norm();
                    //Only when in same direction!
                    double add = (ridge.Start.Position - ridge.End.Position + ridge.Start.Position - line.end.Position).L2Norm();
                    if (add < absLine + absRidge)
                    {
                        return false;
                    }

                    alpha = absLine / absRidge;
                    if (alpha >= 1 - accuracy && alpha <= 1)
                    {
                        alpha = 1;
                        intersectionCase = IntersectionCase.IntersectionIsEndOfRidgeAndLine;
                        return true;
                    }
                    if (absLine / absRidge < 1)
                    {
                        intersectionCase = IntersectionCase.IntersectionIsEndOfLine;
                        return true;
                    }
                    alpha = 1;
                    intersectionCase = IntersectionCase.IntersectionIsEndOfRidge;
                    return true;
                }
                return false;
            }
        }

        public (Cell, IEnumerator<Ridge>) getNeighborFromRidgeNeighbor(Ridge ridge)
        {
            Cell newCell = getNeighbour(ridge);
            AfterCutRidgeEnumerator ridgeEnum = new AfterCutRidgeEnumerator(newCell.Ridges, ridge);

            return (newCell, ridgeEnum);
        }

        public Ridge Subdivide(Ridge ridge, List<Line> lines, double alpha)
        {
            Cell cell = ridge.Cell;
            //Divide Ridge and update Ridge Arrays
            //-------------------------------------
            Vertex newVertex = DivideRidge(ridge, alpha, out Ridge newRidge);
            ridge.Twin.Cell.IntersectionVertex = newVertex.ID;
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
            Ridge[] newRidges;
            Ridge[] newNeighborRidges;
            Cell newCell = new Cell { VoronoiNode = cell.VoronoiNode };
            AddCell(newCell);
            CreateRidges(verticesOfNewRidgeBoundary, cell, newCell, out newRidges, out newNeighborRidges);
            //Link Ridges to old neighbors
            InsertRidgesAndVertices(newRidges, newNeighborRidges);


            //dOnE, DoNe!
            return ridge;
        }

        class AfterCutRidgeEnumerator : IEnumerator<Ridge>
        {
            Ridge[] ridges;
            Ridge current;
            int startIndex;
            int currentIndex;
            int counter;

            int length; //Changes after Reset

            public AfterCutRidgeEnumerator(Ridge[] Ridges, Ridge StartRidge)
            {
                ridges = Ridges;
                counter = -1;
                length = Ridges.Length - 2;
                for (int i = 0; i < Ridges.Length; ++i)
                {
                    if (ridgeComparer.Equals(StartRidge, Ridges[i]))
                    {
                        startIndex = i;
                        currentIndex = i;
                        break;
                    }
                }
            }

            public Ridge Current => current;

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
                    currentIndex = (currentIndex + 1) % ridges.Length;
                    current = ridges[currentIndex];
                    return true;
                }
            }

            public void Reset()
            {
                counter = -1;
                length = ridges.Length;
                currentIndex = startIndex;
            }
        }

        public Ridge FirstCut(Ridge ridge, double alpha)
        {
            //Divide Ridge and update Ridge Arrays
            //-------------------------------------
            Vertex newVertex = DivideRidge(ridge, alpha, out Ridge newRidge);
            ridge.Twin.Cell.IntersectionVertex = newVertex.ID;

            //Find Intersection and insert Ridge
            ridge.Cell.IntersectionVertex = newVertex.ID;

            return ridge;
        }

        public void CloseMesh(List<Line> lines, Ridge firstCutRidge)
        {
            Cell cell = firstCutRidge.Cell;
            //Divide this cell
            //================================================================
            //NewVertices
            Vertex[] verticesOfNewRidgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewRidgeBoundary[0] = firstCutRidge.End;
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
            Ridge[] newRidges;
            Ridge[] newNeighborRidges;
            Cell newCell = new Cell();
            AddCell(newCell);
            CreateRidges(verticesOfNewRidgeBoundary, cell, newCell, out newRidges, out newNeighborRidges);
            InsertRidgesAndVertices(newRidges, newNeighborRidges);
        }

        public IEnumerator<Ridge> getConnectedRidgeEnum(Ridge ridge)
        {
            List<Ridge> outgoingEdges = new List<Ridge>();
            bool newNeighbor = true;
            Ridge worker = ridge;
            while (newNeighbor)
            {
                worker = worker.Twin;
                Ridge[] workerRidges = worker.Cell.Ridges;
                outgoingEdges.Add(worker);
                for (int i = 0; i < workerRidges.Length; ++i)
                {
                    if (ridgeComparer.Equals(workerRidges[i], worker))
                    {
                        if (i == 0)
                            worker = workerRidges.Last();
                        else
                            worker = workerRidges[i - 1];
                        if (ridgeComparer.Equals(worker, ridge))
                        {
                            newNeighbor = false;
                        }
                        break;
                    }
                }
            }
            return new ArrayEnum<Ridge>(outgoingEdges);
        }

        public void VertexCut(Ridge ridge, double alphaCut)
        {
            Cell cell = ridge.Cell;
            Vertex newOldVertex = ridge.End;
            cell.IntersectionVertex = newOldVertex.ID;
            cell = ridge.Twin.Cell;
            cell.IntersectionVertex = newOldVertex.ID;
        }

        public IEnumerator<Ridge> getNeighborFromLineDirection(Ridge ridge, Line line)
        {
            IEnumerator<Ridge> outgoingRidges = getConnectedRidgeEnum(ridge);
            Ridge A = null;
            Ridge B = null;
            if (outgoingRidges.MoveNext())
            {
                A = outgoingRidges.Current;
            }
            else
            {
                throw new Exception();
            }
            bool found = false;
            while (!found && outgoingRidges.MoveNext())
            {
                B = outgoingRidges.Current;
                if (IsBetween(A, line, B))
                {
                    found = true;
                }
                else
                {
                    A = B;
                }
            }

            return new AfterCutRidgeEnumerator(A.Cell.Ridges, A);

            //Check if in positive rotation a, c, b order
            bool IsBetween(Ridge a, Line b, Ridge c)
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

        class MultiNeighRidgesAfterCutEnum : IEnumerator<Ridge>
        {
            IList<Ridge> enumRidges;
            IList<int> block;
            int pointer;
            int blockPointer;
            bool filter;
            public MultiNeighRidgesAfterCutEnum(IList<Ridge> ridges, IList<int> blockInFirstRun)
            {
                enumRidges = ridges;
                block = blockInFirstRun;
                pointer = -1;
                blockPointer = 0;
                filter = true;
            }
            public Ridge Current => enumRidges[pointer];

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
                return (pointer < enumRidges.Count);
            }

            public void Reset()
            {
                pointer = -1;
                filter = false;
            }
        }

        public Ridge AddLineSegment(Ridge ridge, double alpha)
        {
            DivideRidge(ridge, alpha, out Ridge newRidge);
            ridge.IsBoundary = true;
            ridge.Twin.IsBoundary = true;
            ridge.Cell.IntersectionVertex = ridge.End.ID;
            ridge.Twin.Cell.IntersectionVertex = ridge.End.ID;
            return ridge;
        }

        public Ridge SubdivideWithoutNewVertex(Ridge ridge, List<Line> lines)
        {
            Cell cell = ridge.Cell;
            Vertex cutVertex = ridge.End;
            IEnumerator<Ridge> neighRidges = getConnectedRidgeEnum(ridge);
            while (neighRidges.MoveNext())
            {
                Cell neighbor = neighRidges.Current.Cell;
                if (neighbor.ID != cell.ID)
                    neighbor.IntersectionVertex = cutVertex.ID;
            }
            //Divide this cell
            //================================================================
            Vertex[] verticesOfNewRidgeBoundary = new Vertex[lines.Count + 2];
            verticesOfNewRidgeBoundary[0] = cutVertex;
            verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1] = Vertices[cell.IntersectionVertex];
            //Add Vertices of lines
            for (int i = 1; i < verticesOfNewRidgeBoundary.Length - 1; ++i)
            {
                verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i] = lines[i - 1].end;
                AddVertex(verticesOfNewRidgeBoundary[verticesOfNewRidgeBoundary.Length - 1 - i]);
            }
            //New Ridges
            Ridge[] newRidges;
            Ridge[] newNeighborRidges;
            Cell newCell = new Cell();
            AddCell(newCell);
            CreateRidges(verticesOfNewRidgeBoundary, cell, newCell, out newRidges, out newNeighborRidges);
            //Link Ridges to old neighbors
            InsertRidgesAndVertices(newRidges, newNeighborRidges);


            //dOnE, DoNe!
            return ridge;
        }

        public void AddRidge(Ridge ridge)
        {
            ridge.IsBoundary = true;
            ridge.Twin.IsBoundary = true;
            ridge.Cell.IntersectionVertex = ridge.End.ID;
            ridge.Twin.Cell.IntersectionVertex = ridge.End.ID;
        }
    }
}
