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

    class Ridge
    {
        public bool IsBoundary = false;
        public Ridge Twin { get; set; }
        public Cell Cell { get; set; }
        public Vertex Start { get; set; }
        public Vertex End { get; set; }
    }

    class RidgeComparer : EqualityComparer<Ridge>
    {
        public override bool Equals(Ridge x, Ridge y)
        {
            return ((x.Start.ID == y.Start.ID && x.End.ID == y.End.ID) || (x.Start.ID == y.End.ID && x.End.ID == y.Start.ID));
        }

        public override int GetHashCode(Ridge obj)
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
        public Ridge[] Ridges { get; set; }
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

        public Cell getNeighbour(Ridge ridge)
        {
            return ridge.Twin.Cell;
        }

        public void InsertRidgesAndVertices(params Ridge[] additionalRidges)
        {
            Cell cell = additionalRidges[0].Cell;
            int countAdditionalRidges = additionalRidges.Length;
            Ridge[] newRidges = new Ridge[cell.Ridges.Length + countAdditionalRidges];
            bool notInsertedYet = true;
            for (int i = 0; i < cell.Ridges.Length; ++i)
            {
                if (notInsertedYet)
                {
                    newRidges[i] = cell.Ridges[i];
                }
                else
                {
                    newRidges[i + countAdditionalRidges] = cell.Ridges[i];
                }
                if (notInsertedYet && (additionalRidges[0].Start.ID == cell.Ridges[i].End.ID))
                {
                    for (int k = 0; k < countAdditionalRidges; ++k)
                    {
                        newRidges[i + k + 1] = additionalRidges[k];
                    }

                    notInsertedYet = false;
                }
            }
            cell.Ridges = newRidges;

            Vertex[] newVertices = new Vertex[cell.Vertices.Length + countAdditionalRidges];
            for (int i = 0; i < cell.Ridges.Length; ++i)
            {
                newVertices[i] = newRidges[i].Start;
            }
            cell.Vertices = newVertices;
        }

        public void InsertRidgesAndVertices(Ridge[] newRidges, Ridge[] newNeighborRidges)
        {
            Cell cell = newRidges[0].Cell;
            Cell emptyNeighCell = newNeighborRidges[0].Cell;
            Ridge[] oldRidges = cell.Ridges;

            List<Ridge> cellRidges = null;
            List<Ridge> emptyNeighCellRidges = null;
            List<Ridge> workerA = new List<Ridge>(newRidges.Length);
            List<Ridge> workerB = new List<Ridge>(newRidges.Length);
            bool workerAIsActive = true;
            List<Ridge> tmp = workerA;

            //Add new Ridges
            for (int i = 0; i < oldRidges.Length; ++i)
            {
                Ridge activeR = oldRidges[i];
                if (activeR.Start.ID == newRidges[0].Start.ID)
                {
                    cellRidges = tmp;
                    for (int j = 0; j < newRidges.Length; ++j)
                    {
                        cellRidges.Add(newRidges[j]);
                    }
                    tmp = workerAIsActive ? workerB : workerA;
                    workerAIsActive = !workerAIsActive;
                }
                if (activeR.Start.ID == newNeighborRidges[0].Start.ID)
                {
                    emptyNeighCellRidges = tmp;
                    for (int j = 0; j < newNeighborRidges.Length; ++j)
                    {
                        emptyNeighCellRidges.Add(newNeighborRidges[j]);
                    }

                    tmp = workerAIsActive ? workerB : workerA;
                    workerAIsActive = !workerAIsActive;
                }
                tmp.Add(activeR);
            }

            cell.Ridges = cellRidges.ToArray();
            emptyNeighCell.Ridges = emptyNeighCellRidges.ToArray();

            //Update AllRidges
            for (int i = 0; i < cell.Ridges.Length; ++i)
            {
                cell.Ridges[i].Cell = cell;
            }
            for (int i = 0; i < emptyNeighCell.Ridges.Length; ++i)
            {
                emptyNeighCell.Ridges[i].Cell = emptyNeighCell;
            }

            //Vertices
            Vertex[] newVertices = new Vertex[cell.Ridges.Length];
            for (int i = 0; i < cell.Ridges.Length; ++i)
            {
                newVertices[i] = cell.Ridges[i].Start;
            }
            cell.Vertices = newVertices;

            Vertex[] newNeighVertices = new Vertex[emptyNeighCell.Ridges.Length];
            for (int i = 0; i < emptyNeighCell.Ridges.Length; ++i)
            {
                newNeighVertices[i] = emptyNeighCell.Ridges[i].Start;
            }
            emptyNeighCell.Vertices = newNeighVertices;

        }

        public Vertex DivideRidge(Ridge ridge, double alpha, out Ridge newRidge)
        {
            Vector start = ridge.Start.Position;
            Vector end = ridge.End.Position;

            Vector intersection = start * (1 - alpha) + end * alpha;
            Vertex newVertex = new Vertex
            {
                Position = intersection,
            };
            AddVertex(newVertex);

            newRidge = new Ridge
            {
                Start = newVertex,
                End = ridge.End,
                Cell = ridge.Cell
            };
            Ridge newRidgeTwin = new Ridge
            {
                End = newVertex,
                Start = ridge.End,
                Cell = ridge.Twin.Cell,
                Twin = newRidge
            };
            newRidge.Twin = newRidgeTwin;

            ridge.End = newVertex;
            ridge.Twin.Start = newVertex;

            InsertRidgesAndVertices(newRidge);
            InsertRidgesAndVertices(newRidgeTwin);

            return newVertex;
        }

        public void CreateRidges(Vertex[] vertices, Cell cell, Cell neighborCell, out Ridge[] ridges, out Ridge[] twinRidges)
        {
            int count = vertices.Length - 1;
            ridges = new Ridge[count];
            twinRidges = new Ridge[count];
            for (int i = 0; i < count; ++i)
            {
                Ridge ridge = new Ridge
                {
                    Start = vertices[i],
                    End = vertices[i + 1],
                    Cell = cell
                };
                Ridge twinRidge = new Ridge
                {
                    Start = vertices[i + 1],
                    End = vertices[i],
                    Twin = ridge,
                    Cell = neighborCell,
                    IsBoundary = true
                };
                ridge.Twin = twinRidge;
                ridges[i] = ridge;
                twinRidges[count - 1 - i] = twinRidge;
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

            foreach (Ridge ridge in cell.Ridges)
            {
                Ridge newRidge = ridge.Twin;
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
                foreach (Ridge ridge in current.Ridges)
                {
                    Ridge newRidge = ridge.Twin;
                    if (!visited[newRidge.Cell.ID] && !newRidge.IsBoundary)
                    {
                        cells.AddLast(newRidge.Cell);
                        visited[newRidge.Cell.ID] = true;
                    }
                }
            }
        }

        public AggregationGrid ToAggregationGrid(Cell insideCell)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cellsBoSSS = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregationBoSSS = new List<int[]>();
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
                aggregationBoSSS.Add(Agg2Pt);
            }

            GridCommons grd;
            grd = new Grid2D(Triangle.Instance);
            grd.Cells = cellsBoSSS.ToArray();
            //grd.EdgeTagNames.Add(1, BoundaryType.Dirichlet.ToString());
            //grd.Plot2DGrid();
            //grd.DefineEdgeTags(X => (byte)1);

            // aggregation grid
            var agrd = new AggregationGrid(grd, aggregationBoSSS.ToArray());
            return agrd;
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
