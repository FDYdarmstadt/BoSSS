using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerMapper<T>
    {
        PeriodicCornerBoundaryAssigner<T> boundaryAssigner;

        public PeriodicCornerMapper(PeriodicMap map)
        {
            boundaryAssigner = new PeriodicCornerBoundaryAssigner<T>(map);
        }

        public void ConnectPeriodicCorners(CellPairCollection<T> candidates)
        {
            var cornerFinder = new PeriodicCornerCellFinder<T>();
            Queue<MeshCell<T>> cornerCells = cornerFinder.FindInner(candidates);
            ConnectCorners(cornerCells);
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
