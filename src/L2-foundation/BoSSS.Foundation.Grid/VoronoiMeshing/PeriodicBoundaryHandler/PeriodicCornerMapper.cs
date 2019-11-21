using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerMapper<T>
        where T : ILocatable
    {
        PeriodicCornerBoundaryIdentifier<T> boundaryAssigner;

        Queue<MeshCell<T>> cornerCells;

        PeriodicCornerCellFinder<T> cornerFinder;

        public PeriodicCornerMapper(PeriodicMap map)
        {
            cornerFinder = new PeriodicCornerCellFinder<T>(map);
            boundaryAssigner = new PeriodicCornerBoundaryIdentifier<T>(map);
        }

        public void ConnectPeriodicCorners()
        {
            ConnectCorners(cornerCells);
        }

        public void FindPeriodicCorners(CellPairCollection<T> candidates)
        {
            cornerCells = cornerFinder.FindInner(candidates);
        }

        void ConnectCorners(Queue<MeshCell<T>> cornerCells)
        {
            while (cornerCells.Count > 0)
            {
                MeshCell<T> cornerCell = cornerCells.Dequeue();
                boundaryAssigner.SetEdges(cornerCell);
            }
        }
    }
}
