using BoSSS.Platform;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class InsideCellEnumerator<T>
    {
        protected Mesh<T> mesh;

        readonly Boundary<T> boundary;

        public InsideCellEnumerator(Domain<T> mesh)
        {
            boundary = mesh.Boundary;
            this.mesh = mesh.Mesh;
        }

        public InsideCellEnumerator(Mesh<T> mesh, Boundary<T> boundary)
        {
            this.boundary = boundary;
            this.mesh = mesh;
        }

        public void SetFirstCell(Vector start)
        {
            SetFirstCellToCellContaining(start);
        }

        public MeshCell<T> GetFirstCell()
        {
            return boundary.FirstCorner;
        }

        void SetFirstCellToCellContaining(Vector start)
        {
            //Find cell that contains boundaryLine.Start;
            bool foundFirstCell = false;
            //Check if boundaryLine.Start is still in cell, else search neighborhood
            foreach (MeshCell<T> cell in EnumerateCellsInConcentricCircles())
            {
                Vector[] verts = Array.ConvertAll(cell.Vertices, item => (Vector)item);
                //At this point, every cell is convex!
                bool isInside = PolygonTesselation.PointInConvexPolygon(verts, start);
                if (isInside)
                {
                    foundFirstCell = true;
                    boundary.FirstCorner = cell;
                    break;
                }
            }
            if (!foundFirstCell)
            {
                throw new Exception("First cell could not be found: boundaryLine.start not inside a cell");
            }
        }

        public virtual IEnumerable<MeshCell<T>> EnumerateCellsInConcentricCircles()
        {
            Debug.Assert( boundary.FirstCorner != null, "Initialize before calling Cells()");
            return GetInsideCells(boundary.FirstCorner);
        }

        public static IEnumerable<MeshCell<T>> GetInsideCells(MeshCell<T> cell)
        {
            HashSet<int> visited = new HashSet<int>();
            return IterativeYieldConnectedCells(cell, visited);
        }

        /// <summary>
        /// Enumerates the cells inside the boundary of this mesh.
        /// Recursion produces Stack overflow, when to many cells in mesh.
        /// </summary>
        /// <param name="cell">
        /// Enumeration starts with this cell and then return its neighbors and so on.
        /// </param>
        /// <param name="visited"></param>
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
        /// <param name="visited"></param>
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
                    if (newEdge != null && !newEdge.IsBoundary && !visited.Contains(newEdge.Cell.ID) )
                    {
                        cells.Enqueue(newEdge.Cell);
                        visited.Add(newEdge.Cell.ID);
                    }
                }
            }
        }
    }
}
