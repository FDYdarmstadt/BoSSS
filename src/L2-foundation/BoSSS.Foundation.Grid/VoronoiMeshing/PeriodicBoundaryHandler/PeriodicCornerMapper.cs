using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerMapper<T>
        where T : ILocatable
    {
        readonly PeriodicCornerBoundaryIdentifier<T> boundaryIdentifier;
        Queue<MeshCell<T>> innerCorners;
        MeshCell<T> firstCorner;
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
            ConnectCorners(innerCorners);
        }

        public MeshCell<T> GetFirstCornerCell()
        {
            if(firstCorner != null)
            {
                return firstCorner;
            }
            else
            {
                throw new Exception("First Cell not initialized.");
            }
        }

        public void FindPeriodicCorners(CellPairCollection<T> candidates)
        {
            innerCorners = cornerFinder.FindInner(candidates);
        }

        void ConnectCorners(Queue<MeshCell<T>> corners)
        {
            var needFixing = new List<(Stack<Edge<T>> wronglyAssignedEdges, Corner corner)>(corners.Count);

            //Find corner cells without modifying cells
            while (corners.Count > 0)
            {
                MeshCell<T> cornerCell = corners.Dequeue();
                needFixing.Add(boundaryIdentifier.FindEdges(cornerCell));
            }
            foreach (var fixMe in needFixing)
            {
                boundaryAssigner.AssignBoundariesOfPeriodicCorners(fixMe.wronglyAssignedEdges, fixMe.corner);
            }
        }
    }

    class CornerCleaner
    {
        readonly ICollection<int> cornerCandidates;

        public CornerCleaner(int corners)
        {
            this.cornerCandidates = new List<int>(corners * 2);
        }

        public void RemoveAlreadyDealtWithCornerCellMergePairsFrom<T>(CellPairCollection<T>.EdgeCombo edge)
        {
            if (edge.Outer.Count > 0)
            {
                MeshCell<T> firstCandidate = edge.Outer[0];
                if (edge.Inner.Count > 1)
                {
                    MeshCell<T> secondCandidate = edge.Outer[edge.Outer.Count - 1];
                    if (cornerCandidates.Contains(secondCandidate.ID))
                    {
                        edge.Inner.RemoveAt(edge.Outer.Count - 1);
                        edge.Outer.RemoveAt(edge.Outer.Count - 1);
                    }
                    else
                    {
                        cornerCandidates.Add(secondCandidate.ID);
                    }
                }
                if (cornerCandidates.Contains(firstCandidate.ID))
                {
                    edge.Inner.RemoveAt(0);
                    edge.Outer.RemoveAt(0);
                }
                else
                {
                    cornerCandidates.Add(firstCandidate.ID);
                }
            }
        }
    }
}



