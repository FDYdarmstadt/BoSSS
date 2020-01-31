using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerCellFinder<T>
        where T : ILocatable
    {
        readonly PeriodicMap map;

        public PeriodicCornerCellFinder(PeriodicMap map)
        {
            this.map = map;
        }

        public Queue<MeshCell<T>> FindInner(CellPairCollection<T> candidates)
        {
            Queue<MeshCell<T>> cornerCells = new Queue<MeshCell<T>>();
            foreach (CellPairCollection<T>.EdgeCombo edge in candidates.GetCollectedEdgeCombos())
            {
                if (edge.Inner.Count > 0)
                {
                    MeshCell<T> first = edge.Inner[0];
                    if (IsInnerPeriodicCorner(first))
                    {
                        cornerCells.Enqueue(first);
                    }
                }
            }
            return cornerCells;
        }

        bool IsInnerPeriodicCorner(MeshCell<T> cell)
        {
            if (NodeIsInVoronoiCell(cell))
            {
                if (IsRegisteredCorner(cell, out Corner corner))
                {
                    if (map.PeriodicCornerCorrelation.TryGetValue(corner, out int wayne))
                    {
                        return true;
                    }
                };
            }
            return false;
        }

        bool IsRegisteredCorner(MeshCell<T> cell, out Corner corner)
        {
            foreach (Pair<Edge<T>> edgePair in CornerPair(cell))
            {
                corner = new Corner
                {
                    SecondEdge = edgePair.Current.BoundaryEdgeNumber,
                    FirstEdge = edgePair.Previous.BoundaryEdgeNumber
                };
                if (map.PeriodicCornerCorrelation.ContainsKey(corner))
                {
                    return true;
                }
            }
            corner = default;
            return false;
        }

        IEnumerable<Pair<Edge<T>>> CornerPair(MeshCell<T> cell)
        {
            foreach (Pair<Edge<T>> followingBoundaries in new Convolution<Edge<T>>(cell.Edges))
            {
                Edge<T> current = followingBoundaries.Current;
                Edge<T> prevoius = followingBoundaries.Previous;
                if (current.IsBoundary && prevoius.IsBoundary)
                {
                    if (current.BoundaryEdgeNumber != prevoius.BoundaryEdgeNumber)
                    {
                        yield return followingBoundaries;
                    }
                }
            }
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

        public MeshCell<T> FindFirstCorner(Queue<MeshCell<T>> corners)
        {
            Corner second = new Corner
            {
                FirstEdge = 0,
                SecondEdge = 1
            };
            MeshCell<T> firstCorner = null;
            while (corners.Count > 0)
            {
                MeshCell<T> cell = corners.Dequeue();
                if(IsRegisteredCorner(cell, out Corner corner))
                {
                    if(corner.FirstEdge == 0 || corner.SecondEdge == 0)
                    {
                        if (!corner.Equals(second))
                        {
                            firstCorner = cell;
                        }
                    }
                }
            }
            return firstCorner;
        }
    }
}
