using System.Collections;
using System.Collections.Generic;
using System;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class MeshDivider<T>
    {
        IDMesh<T> mesh;

        List<MeshCell<T>> insideCells;

        List<T> insideNodes;

        public MeshCell<T> insideCell;

        public int Count {
            get {
                return insideCells.Count;
            }
        }

        public MeshDivider(IDMesh<T> mesh, int firstCell_NodeIndice)
            : this(mesh)
        {
            
            this.insideCell = mesh.Cells[firstCell_NodeIndice];
        }

        public MeshDivider(IDMesh<T> mesh)
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
            foreach (MeshCell<T> cell in CellsOnSameSideOfBoundary_Iterative(insideCell))
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

        public List<MeshCell<T>> GetInsideCells()
        {
            if (insideCells == null)
            {
                DetermineInsideCells();
            }
            return insideCells;

        }

        public void RemoveOutsideCells()
        {
            mesh.Cells = GetInsideCells();
            mesh.Nodes = GetInsideNodes();
        }

        public List<T> GetInsideNodes()
        {
            if (insideNodes == null)
            {
                DetermineInsideNodes();
            }
            return insideNodes;
        }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        public IEnumerable<MeshCell<T>> CellsOnSameSideOfBoundary_Recursive(MeshCell<T> cell)
        {
            BitArray visited = new BitArray(mesh.Cells.Count);
            return YieldConnectedCells(cell, visited);
        }

        private IEnumerable<MeshCell<T>> YieldConnectedCells(MeshCell<T> cell, BitArray visited)
        {
            visited[cell.ID] = true;
            yield return cell;

            foreach (Edge<T> edge in cell.Edges)
            {
                Edge<T> newRidge = edge.Twin;
                if (!visited[newRidge.Cell.ID] && !newRidge.IsBoundary)
                {
                    foreach (MeshCell<T> neighbor in YieldConnectedCells(newRidge.Cell, visited))
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
        public IEnumerable<MeshCell<T>> CellsOnSameSideOfBoundary_Iterative(MeshCell<T> cell)
        {
            BitArray visited = new BitArray(mesh.Cells.Count);
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

        public void DetermineInsideCells()
        {
            insideCells = SetCellTypesOnSameSideOfBoundary(insideCell, MeshCellType.Inside, ref insideCells);
        }

        void DetermineInsideNodes()
        {
            IReadOnlyList<MeshCell<T>> cells = GetInsideCells();
            insideNodes = new List<T>(cells.Count);
            foreach (MeshCell<T> cell in cells)
            {
                insideNodes.Add(cell.Node);
            }
        }

        protected List<MeshCell<T>> SetCellTypesOnSameSideOfBoundary(
            MeshCell<T> startingCell, 
            MeshCellType type, 
            ref List<MeshCell<T>> cells)
        {
            List<MeshCell<T>> sameSideCells = new List<MeshCell<T>>();
            foreach (MeshCell<T> cell in CellsOnSameSideOfBoundary_Iterative(startingCell))
            {
                cell.type = type;
                sameSideCells.Add(cell);
            }
            return sameSideCells;
        }
    }
}
