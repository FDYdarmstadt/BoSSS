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

    class Edge
    {
        public bool IsBoundary = false;
        public Edge Twin { get; set; }
        public Cell Cell { get; set; }
        public Vertex Start { get; set; }
        public Vertex End { get; set; }
    }

    class EdgeComparer : EqualityComparer<Edge>
    {
        public override bool Equals(Edge x, Edge y)
        {
            return ((x.Start.ID == y.Start.ID && x.End.ID == y.End.ID) || (x.Start.ID == y.End.ID && x.End.ID == y.Start.ID));
        }

        public override int GetHashCode(Edge obj)
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

    class Cell
    {
        public int ID { get; set; }
        public Vertex[] Vertices { get; set; }
        public Edge[] Edges { get; set; }
        public int IntersectionVertex { get; set; }
        public Vector VoronoiNode { get; set; }
    }

    interface IIdMesh
    {
        IReadOnlyList<Cell> Cells { get; }
        IReadOnlyList<Vertex> Vertices { get; }
        int AddCell(Cell cell);
        Cell GetCell(int ID);
        int AddVertex(Vertex vert);
        Vertex GetVertex(int ID);
    }

    class SimpleIdMesh : IIdMesh
    {
        protected List<Cell> cells;
        protected List<Vertex> vertices;

        public SimpleIdMesh()
        {
            cells = new List<Cell>();
            vertices = new List<Vertex>();
        }

        public SimpleIdMesh(List<Cell> Cells, List<Vertex> Vertices)
        {
            cells = Cells;
            vertices = Vertices;
        }

        public IReadOnlyList<Cell> Cells {
            get {
                return cells;
            }
        }

        public IReadOnlyList<Vertex> Vertices {
            get {
                return vertices;
            }
        }

        public int AddCell(Cell cell)
        {
            cell.ID = cells.Count;
            cells.Add(cell);
            return cell.ID;
        }

        public Cell GetCell(int ID)
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

    class Mesh : IIdMesh
    {
        IIdMesh mesh;

        public Mesh(IIdMesh Mesh)
        {
            mesh = Mesh;
        }

        public Cell getNeighbour(Edge ridge)
        {
            return ridge.Twin.Cell;
        }

        public void InsertEdgesAndVertices(params Edge[] additionalEdges)
        {
            Cell cell = additionalEdges[0].Cell;
            int countAdditionalRidges = additionalEdges.Length;
            Edge[] newRidges = new Edge[cell.Edges.Length + countAdditionalRidges];
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

        public void InsertEdgesAndVertices(Edge[] newEdge, Edge[] newNeighborEdges)
        {
            Cell cell = newEdge[0].Cell;
            Cell emptyNeighCell = newNeighborEdges[0].Cell;
            Edge[] oldRidges = cell.Edges;

            List<Edge> cellRidges = null;
            List<Edge> emptyNeighCellRidges = null;
            List<Edge> workerA = new List<Edge>(newEdge.Length);
            List<Edge> workerB = new List<Edge>(newEdge.Length);
            bool workerAIsActive = true;
            List<Edge> tmp = workerA;

            //Add new Ridges
            for (int i = 0; i < oldRidges.Length; ++i)
            {
                Edge activeR = oldRidges[i];
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

        public Vertex DivideEdge(Edge edge, double alpha, out Edge newEdge)
        {
            Vector start = edge.Start.Position;
            Vector end = edge.End.Position;

            Vector intersection = start * (1 - alpha) + end * alpha;
            Vertex newVertex = new Vertex
            {
                Position = intersection,
            };
            AddVertex(newVertex);

            newEdge = new Edge
            {
                Start = newVertex,
                End = edge.End,
                Cell = edge.Cell
            };
            Edge newRidgeTwin = new Edge
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

        public void CreateEdge(Vertex[] vertices, Cell cell, Cell neighborCell, out Edge[] ridges, out Edge[] twinEdges)
        {
            int count = vertices.Length - 1;
            ridges = new Edge[count];
            twinEdges = new Edge[count];
            for (int i = 0; i < count; ++i)
            {
                Edge ridge = new Edge
                {
                    Start = vertices[i],
                    End = vertices[i + 1],
                    Cell = cell
                };
                Edge twinRidge = new Edge
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

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        public IEnumerable<Cell> ConnectedCells_Recursive(Cell cell)
        {
            BitArray visited = new BitArray(Cells.Count);
            return YieldConnectedCells(cell, visited);
        }

        private IEnumerable<Cell> YieldConnectedCells(Cell cell, BitArray visited)
        {
            visited[cell.ID] = true;
            yield return cell;

            foreach (Edge edge in cell.Edges)
            {
                Edge newRidge = edge.Twin;
                if (!visited[newRidge.Cell.ID] && !newRidge.IsBoundary)
                {
                    foreach (Cell neighbor in YieldConnectedCells(newRidge.Cell, visited))
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
        public IEnumerable<Cell> ConnectedCells_Iterative(Cell cell)
        {
            BitArray visited = new BitArray(Cells.Count);
            LinkedList<Cell> cells = new LinkedList<Cell>();
            cells.AddFirst(cell);
            visited[cell.ID] = true;
            while (cells.Count > 0)
            {
                Cell current = cells.First.Value;
                yield return current;
                cells.RemoveFirst();
                foreach (Edge edge in current.Edges)
                {
                    Edge newEdge = edge.Twin;
                    if (!visited[newEdge.Cell.ID] && !newEdge.IsBoundary)
                    {
                        cells.AddLast(newEdge.Cell);
                        visited[newEdge.Cell.ID] = true;
                    }
                }
            }
        }

        (GridCommons grid, int[][] aggregation, MultidimensionalArray nodes) GetVoronoiData(Cell insideCell)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cellsBoSSS = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregation = new List<int[]>();
            MultidimensionalArray nodes = MultidimensionalArray.Create(Cells.Count, 2);
            foreach (Cell cell in ConnectedCells_Iterative(insideCell))
            {
                //Convert to BoSSSCell : Triangulate
                Vector[] VoronoiCell = cell.Vertices.Select(voVtx => voVtx.Position).ToArray();
                int[,] iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell);
                int[] Agg2Pt = new int[iVtxTri.GetLength(0)];

                for (int iTri = 0; iTri < iVtxTri.GetLength(0); iTri++)
                { // loop over triangles of voronoi cell
                    int iV0 = iVtxTri[iTri, 0];
                    int iV1 = iVtxTri[iTri, 1];
                    int iV2 = iVtxTri[iTri, 2];

                    Vector V0 = VoronoiCell[iV0];
                    Vector V1 = VoronoiCell[iV1];
                    Vector V2 = VoronoiCell[iV2];

                    Vector D1 = V1 - V0;
                    Vector D2 = V2 - V0;

                    if (D1.CrossProduct2D(D2) < 0)
                    {
                        int it = iV0;
                        iV0 = iV2;
                        iV2 = it;

                        Vector vt = V0;
                        V0 = V2;
                        V2 = vt;

                        D1 = V1 - V0;
                        D2 = V2 - V0;
                    }

                    Debug.Assert(D1.CrossProduct2D(D2) > 1.0e-8);


                    BoSSS.Foundation.Grid.Classic.Cell Cj = new BoSSS.Foundation.Grid.Classic.Cell();
                    Cj.GlobalID = cellsBoSSS.Count;
                    Cj.Type = CellType.Triangle_3;
                    Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                    Cj.NodeIndices = new int[] { cell.Vertices[iV0].ID, cell.Vertices[iV1].ID, cell.Vertices[iV2].ID };
                    Cj.TransformationParams.SetRowPt(0, V0);
                    Cj.TransformationParams.SetRowPt(1, V1);
                    Cj.TransformationParams.SetRowPt(2, V2);

                    Agg2Pt[iTri] = cellsBoSSS.Count;
                    cellsBoSSS.Add(Cj);
                }
                aggregation.Add(Agg2Pt);
            }

            GridCommons grd;
            grd = new Grid2D(Triangle.Instance);
            grd.Cells = cellsBoSSS.ToArray();
            return (grd, aggregation.ToArray(), nodes);
        } 

        public VoronoiGrid ToVoronoiGrid(Cell insideCell)
        {
            (GridCommons grid, int[][] aggregation, MultidimensionalArray nodes) = GetVoronoiData(insideCell);
            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes);
            return voronoiGrid;
        }

        #region IIdMesh

        public IReadOnlyList<Cell> Cells => mesh.Cells;

        public IReadOnlyList<Vertex> Vertices => mesh.Vertices;

        public int AddCell(Cell cell)
        {
            return mesh.AddCell(cell);
        }

        public Cell GetCell(int ID)
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

}
