using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi
{
    class Vertex
    {
        public Vector Position { get; set; }
        public int ID { get; set; }
        public static explicit operator Vector(Vertex vtx)
        {
            return vtx.Position;
        }
    }

    class Edge<T>
    {
        public bool IsBoundary = false;
        public Edge<T> Twin { get; set; }
        public Cell<T> Cell { get; set; }
        public Vertex Start { get; set; }
        public Vertex End { get; set; }
    }

    class EdgeComparer<T> : EqualityComparer<Edge<T>>
    {
        public override bool Equals(Edge<T> x, Edge<T> y)
        {
            return ((x.Start.ID == y.Start.ID && x.End.ID == y.End.ID) || (x.Start.ID == y.End.ID && x.End.ID == y.Start.ID));
        }

        public override int GetHashCode(Edge<T> obj)
        {
            //http://eternallyconfuzzled.com/tuts/algorithms/jsw_tut_hashing.aspx 
            //If x == y hash must be hash(x) = hash(y)
            int start = obj.End.ID > obj.Start.ID ? obj.Start.ID : obj.End.ID;
            int end = obj.End.ID < obj.Start.ID ? obj.Start.ID : obj.End.ID;

            int hash = 17;
            hash = hash * 23 + start.GetHashCode();
            hash = hash * 23 + end.GetHashCode();

            return hash;
        }
    }

    enum CType { Inside, NotDetermined }

    class Cell<T>
    {
        public CType type = CType.NotDetermined;

        public T Node;
        public int ID { get; set; }
        public Vertex[] Vertices { get; set; }
        public Edge<T>[] Edges { get; set; }
        public int IntersectionVertex { get; set; }

        public Vector Position { get; set; }
    }

    interface IIdMesh<T>
    {
        IReadOnlyList<Cell<T>> Cells { get; }
        IReadOnlyList<Vertex> Vertices { get; }
        int AddCell(Cell<T> cell);
        Cell<T> GetCell(int ID);
        int AddVertex(Vertex vert);
        Vertex GetVertex(int ID);
    }

    class SimpleIdMesh<T> : IIdMesh<T>
    {
        protected List<Cell<T>> cells;
        protected List<Vertex> vertices;

        public SimpleIdMesh()
        {
            cells = new List<Cell<T>>();
            vertices = new List<Vertex>();
        }

        public SimpleIdMesh(List<Cell<T>> Cells, List<Vertex> Vertices)
        {
            cells = Cells;
            vertices = Vertices;
        }

        public IReadOnlyList<Cell<T>> Cells {
            get {
                return cells;
            }
        }

        public IReadOnlyList<Vertex> Vertices {
            get {
                return vertices;
            }
        }

        public int AddCell(Cell<T> cell)
        {
            cell.ID = cells.Count;
            cells.Add(cell);
            return cell.ID;
        }

        public Cell<T> GetCell(int ID)
        {
            return cells[ID];
        }

        public int AddVertex(Vertex vert)
        {
            vert.ID = vertices.Count;
            vertices.Add(vert);
            return vert.ID;
        }

        public Vertex GetVertex(int ID)
        {
            return vertices[ID];
        }
    }

    class Mesh<T>
    {
        IIdMesh<T> mesh;

        public Mesh(IIdMesh<T> Mesh)
        {
            mesh = Mesh;
        }

        public Cell<T> getNeighbour(Edge<T> ridge)
        {
            return ridge.Twin.Cell;
        }

        public void InsertEdgesAndVertices(params Edge<T>[] additionalEdges)
        {
            Cell<T> cell = additionalEdges[0].Cell;
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

        public void InsertEdgesAndVertices(Edge<T>[] newEdge, Edge<T>[] newNeighborEdges)
        {
            Cell<T> cell = newEdge[0].Cell;
            Cell<T> emptyNeighCell = newNeighborEdges[0].Cell;
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

        public Vertex DivideEdge(Edge<T> edge, double alpha, out Edge<T> newEdge)
        {
            Vector start = edge.Start.Position;
            Vector end = edge.End.Position;

            Vector intersection = start * (1 - alpha) + end * alpha;
            Vertex newVertex = new Vertex
            {
                Position = intersection,
            };
            AddVertex(newVertex);

            newEdge = new Edge<T>
            {
                Start = newVertex,
                End = edge.End,
                Cell = edge.Cell
            };
            Edge<T> newRidgeTwin = new Edge<T>
            {
                End = newVertex,
                Start = edge.End,
                Cell = edge.Twin.Cell,
                Twin = newEdge
            };
            newEdge.Twin = newRidgeTwin;

            edge.End = newVertex;
            edge.Twin.Start = newVertex;

            InsertEdgesAndVertices(newEdge);
            InsertEdgesAndVertices(newRidgeTwin);

            return newVertex;
        }

        public void CreateEdge(Vertex[] vertices, Cell<T> cell, Cell<T> neighborCell, out Edge<T>[] ridges, out Edge<T>[] twinEdges)
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
                    Cell = cell
                };
                Edge<T> twinRidge = new Edge<T>
                {
                    Start = vertices[i + 1],
                    End = vertices[i],
                    Twin = ridge,
                    Cell = neighborCell,
                    IsBoundary = true
                };
                ridge.Twin = twinRidge;
                ridges[i] = ridge;
                twinEdges[count - 1 - i] = twinRidge;
            }
        }
        
        #region IIdMesh

        protected IReadOnlyList<Cell<T>> Cells => mesh.Cells;

        protected IReadOnlyList<Vertex> Vertices => mesh.Vertices;

        public int AddCell(Cell<T> cell)
        {
            return mesh.AddCell(cell);
        }

        public Cell<T> GetCell(int ID)
        {
            return mesh.GetCell(ID);
        }

        public int AddVertex(Vertex vert)
        {
            return mesh.AddVertex(vert);
        }

        public Vertex GetVertex(int ID)
        {
            return mesh.GetVertex(ID);
        }

        #endregion
    }

    class BoundaryMesh<T> : Mesh<T>
    {
        List<Cell<T>> insideCells;

        List<T> insideNodes;

        protected Cell<T> insideCell;

        public int Count {
            get {
                return insideCells.Count;
            }
        }

        public virtual IReadOnlyList<Cell<T>> GetCells()
        {
            if (insideCells == null)
            {
                DetermineInsideCells();
            }
            return insideCells;
            
        }

        public List<T> GetNodes()
        {
            if(insideNodes == null)
            {
                DetermineInsideNodes();
            }
            return insideNodes;
        }

        public BoundaryMesh(IIdMesh<T> mesh)
            : base(mesh)
        { }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        public IEnumerable<Cell<T>> CellsOnSameSideOfBoundary_Recursive(Cell<T> cell)
        {
            BitArray visited = new BitArray(Cells.Count);
            return YieldConnectedCells(cell, visited);
        }

        private IEnumerable<Cell<T>> YieldConnectedCells(Cell<T> cell, BitArray visited)
        {
            visited[cell.ID] = true;
            yield return cell;

            foreach (Edge<T> edge in cell.Edges)
            {
                Edge<T> newRidge = edge.Twin;
                if (!visited[newRidge.Cell.ID] && !newRidge.IsBoundary)
                {
                    foreach (Cell<T> neighbor in YieldConnectedCells(newRidge.Cell, visited))
                    {
                        yield return neighbor;
                    }
                }
            }
        }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        public IEnumerable<Cell<T>> CellsOnSameSideOfBoundary_Iterative(Cell<T> cell)
        {
            BitArray visited = new BitArray(Cells.Count);
            LinkedList<Cell<T>> cells = new LinkedList<Cell<T>>();
            cells.AddFirst(cell);
            visited[cell.ID] = true;
            while (cells.Count > 0)
            {
                Cell<T> current = cells.First.Value;
                yield return current;
                cells.RemoveFirst();
                foreach (Edge<T> edge in current.Edges)
                {
                    Edge<T> newEdge = edge.Twin;
                    if (!visited[newEdge.Cell.ID] && !newEdge.IsBoundary)
                    {
                        cells.AddLast(newEdge.Cell);
                        visited[newEdge.Cell.ID] = true;
                    }
                }
            }
        }

        protected void DetermineInsideCells()
        {
            insideCells = SetCellTypesOnSameSideOfBoundary(insideCell, CType.Inside, ref insideCells);
        }

        void DetermineInsideNodes()
        {
            IReadOnlyList<Cell<T>> cells = GetCells();
            insideNodes = new List<T>(cells.Count);
            foreach(Cell<T> cell in cells)
            {
                insideNodes.Add(cell.Node);
            }
        }

        protected List<Cell<T>> SetCellTypesOnSameSideOfBoundary(Cell<T> startingCell, CType type, ref List<Cell<T>> cells)
        {
            List<Cell<T>> sameSideCells = new List<Cell<T>>();
            foreach (Cell<T> cell in CellsOnSameSideOfBoundary_Iterative(startingCell))
            {
                cell.type = type;
                sameSideCells.Add(cell);
            }
            return sameSideCells;
        }
    }
}
