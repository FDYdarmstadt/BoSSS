using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler.NodeLocation;

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
            MeshCell<T> last;

            foreach (Edge<T> edge in periodicEdges)
            {
                if (isFirstBoundaryLine)
                {
                    isFirstBoundaryLine = IsStillFirstBoundary(edge.BoundaryEdgeNumber, ref firstBoundaryEdgeNumber);
                }
                if (isFirstBoundaryLine)
                {
                    if (NodeOfEdgeIsOnRightSideOfEdge(edge))
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
                    if (NodeOfEdgeIsOnRightSideOfEdge(edge))
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
    }
}
