using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerCellFinder<T>
    {
        readonly ICollection<int> visitedCorners;

        readonly Queue<MeshCell<T>> cornerCells;

        public PeriodicCornerCellFinder()
        {
            visitedCorners = new LinkedList<int>();
            cornerCells = new Queue<MeshCell<T>>();
        }

        public Queue<MeshCell<T>> FindInner(CellPairCollection<T> candidates)
        {
            foreach (CellPairCollection<T>.EdgeCombo edge in candidates.GetCollectedEdgeCombos())
            {
                if (edge.Inner.Count > 0)
                {
                    MeshCell<T> first = edge.Inner[0];
                    if (AlreadyVisited(first))
                    {
                        cornerCells.Enqueue(first);
                    }
                    //add last
                    if (edge.Inner.Count > 1)
                    {
                        MeshCell<T> last = edge.Inner[edge.Inner.Count - 1];
                        if (AlreadyVisited(last))
                        {
                            cornerCells.Enqueue(last);
                        }
                    }
                }
            }
            return cornerCells;
        }

        bool AlreadyVisited(MeshCell<T> cell)
        {
            if (visitedCorners.Contains(cell.ID))
            {
                visitedCorners.Remove(cell.ID);
                return true;
            }
            else
            {
                visitedCorners.Add(cell.ID);
                return false;
            }
        }
    }
}
