using System.Collections;
using System.Collections.Generic;
using System;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class Divider<T>
    {
        IDMesh<T> mesh;

        MeshCell<T> insideCell;

        public Divider(IDMesh<T> mesh, int firstCell_NodeIndice)
            : this(mesh)
        {
            this.insideCell = mesh.Cells[firstCell_NodeIndice];
        }

        public Divider(IDMesh<T> mesh)
        {
            this.mesh = mesh;
        }

        public MeshCell<T> GetFirst(BoundaryLine boundaryLine)
        {
            //Find cell that contains boundaryLine.Start;
            bool foundFirstCell = false;
            if (insideCell == null)
            {
                //SetFirst Cell: any random cell. Influences runtime, though
                insideCell = mesh.Cells[0];
                foundFirstCell = true;
            }

            //Check if boundaryLine.Start is still in cell, else search neighborhood
            foreach (MeshCell<T> cell in CellsOnSameSideOfBoundary(insideCell))
            {
                Vector[] verts = Array.ConvertAll(cell.Vertices, item => (Vector)item);
                //At this point, every cell is convex!
                bool isInside = PolygonTesselation.PointInConvexPolygon(verts, (Vector)boundaryLine.Start);
                if (isInside)
                {
                    foundFirstCell = true;
                    insideCell = cell;
                    break;
                }
            }
            if (!foundFirstCell)
            {
                throw new Exception("First cell could not be found: boundaryLine.start not inside a cell");
            }
            return insideCell;
        }

        public void RemoveOutsideCells()
        {
            IdentifyInsideCells();
            List<MeshCell<T>> insideCells = CollectInsideCells();
            mesh.Cells = insideCells;
            List<T> nodes = CollectNodes(insideCells);
            mesh.Nodes = nodes;
        }

        List<MeshCell<T>> CollectInsideCells()
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach(MeshCell<T> cell in mesh.Cells)
            {
                if(cell.type == MeshCellType.Inside)
                {
                    cells.Add(cell);
                }
            }
            return cells;
        }

        static List<T> CollectNodes(List<MeshCell<T>> cells)
        {
            List<T> nodes = new List<T>(cells.Count);
            foreach(MeshCell<T> cell in cells)
            {
                nodes.Add(cell.Node);
            }
            return nodes;
        }

        void IdentifyInsideCells()
        {
            foreach (MeshCell<T> cell in CellsOnSameSideOfBoundary(insideCell))
            {
                cell.type = MeshCellType.Inside;
            }
        }

        IEnumerable<MeshCell<T>> CellsOnSameSideOfBoundary(MeshCell<T> cell)
        {
            BitArray visited = new BitArray(mesh.Cells.Count);
            return IterativeYieldConnectedCells(cell, visited);
        }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        private static IEnumerable<MeshCell<T>> RecursiveYieldConnectedCells(MeshCell<T> cell, BitArray visited)
        {
            visited[cell.ID] = true;
            yield return cell;

            foreach (Edge<T> edge in cell.Edges)
            {
                Edge<T> newRidge = edge.Twin;
                if (!visited[newRidge.Cell.ID] && !newRidge.IsBoundary)
                {
                    foreach (MeshCell<T> neighbor in RecursiveYieldConnectedCells(newRidge.Cell, visited))
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
        static IEnumerable<MeshCell<T>> IterativeYieldConnectedCells(MeshCell<T> cell, BitArray visited)
        {
            LinkedList<MeshCell<T>> cells = new LinkedList<MeshCell<T>>();
            cells.AddFirst(cell);
            visited[cell.ID] = true;
            while (cells.Count > 0)
            {
                MeshCell<T> current = cells.First.Value;
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

        
    }
}
