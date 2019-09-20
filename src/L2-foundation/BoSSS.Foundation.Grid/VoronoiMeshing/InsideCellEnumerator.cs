using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class InsideCellEnumerator<T>
    {
        protected Mesh<T> mesh;

        MeshCell<T> firstCell;

        public InsideCellEnumerator(Mesh<T> mesh)
        {
            firstCell = mesh.Cells[0];
            this.mesh = mesh;
        }

        public InsideCellEnumerator(Mesh<T> mesh, int firstCell_NodeIndice)
        {
            this.firstCell = mesh.Cells[firstCell_NodeIndice];
            this.mesh = mesh;
        }

        public void SetFirstCell(Vector start)
        {
            SetFirstCellToCellContaining(start);
        }

        public void SetFirstCell(Vector start, int firstCell_NodeIndice)
        {
            this.firstCell = mesh.Cells[firstCell_NodeIndice];
            SetFirstCellToCellContaining(start);
        }

        public MeshCell<T> GetFirstCell()
        {
            return firstCell;
        }

        void SetFirstCellToCellContaining(Vector start)
        {
            //Find cell that contains boundaryLine.Start;
            bool foundFirstCell = false;
            //Check if boundaryLine.Start is still in cell, else search neighborhood
            foreach (MeshCell<T> cell in Cells())
            {
                Vector[] verts = Array.ConvertAll(cell.Vertices, item => (Vector)item);
                //At this point, every cell is convex!
                bool isInside = PolygonTesselation.PointInConvexPolygon(verts, start);
                if (isInside)
                {
                    foundFirstCell = true;
                    firstCell = cell;
                    break;
                }
            }
            if (!foundFirstCell)
            {
                throw new Exception("First cell could not be found: boundaryLine.start not inside a cell");
            }
        }

        public virtual IEnumerable<MeshCell<T>> Cells()
        {
            Debug.Assert( firstCell != null, "Initialize before calling Cells()");
            HashSet<int> visited = new HashSet<int>();//mesh.Cells.Count);
            return IterativeYieldConnectedCells(firstCell, visited);
        }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <returns></returns>
        static IEnumerable<MeshCell<T>> RecursiveYieldConnectedCells(MeshCell<T> cell, HashSet<int> visited)
        {
            visited.Add(cell.ID);
            yield return cell;

            foreach (Edge<T> edge in cell.Edges)
            {
                Edge<T> newRidge = edge.Twin;
                if (!visited.Contains(newRidge.Cell.ID) && !newRidge.IsBoundary)
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
        static IEnumerable<MeshCell<T>> IterativeYieldConnectedCells(MeshCell<T> cell, HashSet<int> visited)
        {
            Queue<MeshCell<T>> cells = new Queue<MeshCell<T>>();
            cells.Enqueue(cell);
            visited.Add(cell.ID);
            while (cells.Count > 0)
            {
                MeshCell<T> current = cells.Dequeue();
                yield return current;
                foreach (Edge<T> edge in current.Edges)
                {
                    Edge<T> newEdge = edge.Twin;
                    if (!visited.Contains(newEdge.Cell.ID) && !newEdge.IsBoundary)
                    {
                        cells.Enqueue(newEdge.Cell);
                        visited.Add(newEdge.Cell.ID);
                    }
                }
            }
        }
    }
}
