using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerMapper<T>
        where T : ILocatable
    {
        readonly PeriodicCornerBoundaryIdentifier<T> boundaryIdentifier;
        Queue<MeshCell<T>> cornerCells;
        readonly PeriodicCornerCellFinder<T> cornerFinder;
        readonly PeriodicCornerBoundaryAssigner<T> boundaryAssigner;

        public PeriodicCornerMapper(PeriodicMap map)
        {
            cornerFinder = new PeriodicCornerCellFinder<T>(map);
            boundaryIdentifier = new PeriodicCornerBoundaryIdentifier<T>(map);
            boundaryAssigner = new PeriodicCornerBoundaryAssigner<T>(map);
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
                (Stack<Edge<T>> wronglyAssignedEdges, Corner corner) = boundaryIdentifier.FindEdges(cornerCell);
                boundaryAssigner.AssignBoundariesOfPeriodicCorners(wronglyAssignedEdges, corner);
            }
        }
    }
}



