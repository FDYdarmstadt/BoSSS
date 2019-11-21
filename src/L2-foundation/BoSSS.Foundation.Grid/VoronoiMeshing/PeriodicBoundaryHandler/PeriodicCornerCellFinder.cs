using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerCellFinder<T>
        where T : ILocatable
    {
        readonly Queue<MeshCell<T>> cornerCells;

        readonly PeriodicMap map;

        public PeriodicCornerCellFinder(PeriodicMap map)
        {
            this.map = map;
            cornerCells = new Queue<MeshCell<T>>();
        }

        public Queue<MeshCell<T>> FindInner(CellPairCollection<T> candidates)
        {
            foreach (CellPairCollection<T>.EdgeCombo edge in candidates.GetCollectedEdgeCombos())
            {
                if (edge.Inner.Count > 0)
                {
                    MeshCell<T> first = edge.Inner[0];
                    if (IsPeriodicCorner(first))
                    {
                        cornerCells.Enqueue(first);
                    }
                }
            }
            return cornerCells;
        }

        bool IsPeriodicCorner(MeshCell<T> cell)
        {
            if (NodeIsInVoronoiCell(cell))
            {
                if (IsCorner(cell, out Corner corner))
                {
                    if (map.PeriodicCornerCorrelation.TryGetValue(corner, out int wayne))
                    {
                        return true;
                    }
                };
            }
            return false;
        }

        bool IsCorner(MeshCell<T> cell, out Corner corner)
        {
            foreach (Pair<Edge<T>> followingBoundaries in new Convolution<Edge<T>>(cell.Edges))
            {
                Edge<T> current = followingBoundaries.Current;
                Edge<T> following = followingBoundaries.Previous;
                if(current.IsBoundary && following.IsBoundary)
                {
                    if (Math.Abs(current.BoundaryEdgeNumber - following.BoundaryEdgeNumber) == 1)
                    {
                        corner = new Corner
                        {
                            FirstEdge = current.BoundaryEdgeNumber,
                            SecondEdge = following.BoundaryEdgeNumber
                        };
                        return true;
                    }
                }
            }
            corner = default;
            return false;
        }

        bool NodeIsInVoronoiCell(MeshCell<T> cell)
        {
            Vector[] vertices = Cast(cell.Vertices);
            return PolygonTesselation.PointInConvexPolygon(vertices, cell.Node.Position);
        }

        Vector[] Cast(IList<Vertex> children) 
        {
            Vector[] parents = new Vector[children.Count];
            for(int i = 0; i < children.Count; ++i)
            {
                parents[i] = (Vector)children[i];
            }
            return parents;
        }

    }
}
