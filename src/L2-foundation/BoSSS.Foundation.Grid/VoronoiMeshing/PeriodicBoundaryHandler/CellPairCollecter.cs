using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    static class CellPairCollecter<T>
        where T : ILocatable
    {
        static public CellPairCollection<T> FollowBoundaryAndCollectCandidates(IEnumerable<Edge<T>> periodicEdges)
        {
            LinkedList<MeshCell<T>> firstInnerCells = new LinkedList<MeshCell<T>>();
            LinkedList<MeshCell<T>> firstOuterCells = new LinkedList<MeshCell<T>>();
            CellPairCollection<T> candidates = new CellPairCollection<T>();
            int firstBoundaryEdgeNumber = -1;
            bool isFirstBoundaryLine = true;

            foreach (Edge<T> edge in periodicEdges)
            {
                if (isFirstBoundaryLine)
                {
                    isFirstBoundaryLine = IsStillFirstBoundary(edge.BoundaryEdgeNumber, ref firstBoundaryEdgeNumber);
                }
                if (isFirstBoundaryLine)
                {
                    if (NodeOfEdgeIsOutsideOfBoundary(edge))
                    {
                        firstOuterCells.AddLast(edge.Cell);
                    }
                    else
                    {
                        firstInnerCells.AddLast(edge.Cell);
                    }
                }
                else
                {
                    if (NodeOfEdgeIsOutsideOfBoundary(edge))
                    {
                        candidates.AddOuterCell(edge.Cell, edge.BoundaryEdgeNumber);
                    }
                    else
                    {
                        candidates.AddInnerCell(edge.Cell, edge.BoundaryEdgeNumber);
                    }
                }
            }
            candidates.AddOuterCells(firstOuterCells, firstBoundaryEdgeNumber);
            candidates.AddInnerCells(firstInnerCells, firstBoundaryEdgeNumber);
            return candidates;
        }

        static bool IsStillFirstBoundary(int boundaryEdgeNumber, ref int firstBoundaryEdgeNumber)
        {
            if (firstBoundaryEdgeNumber == -1)
            {
                firstBoundaryEdgeNumber = boundaryEdgeNumber;
            }
            return firstBoundaryEdgeNumber == boundaryEdgeNumber;
        }

        static bool NodeOfEdgeIsOutsideOfBoundary(Edge<T> edge)
        {
            Vector position = edge.Cell.Node.Position;
            return IsOnRightSide(position, edge);
        }

        static bool IsOnRightSide(Vector node, Line line)
        {
            Vector start = line.Start.Position;
            Vector end = line.End.Position;
            return ((end.x - start.x) * (node.y - start.y) - (end.y - start.y) * (node.x - start.x)) < 0;
        }
    }
}
