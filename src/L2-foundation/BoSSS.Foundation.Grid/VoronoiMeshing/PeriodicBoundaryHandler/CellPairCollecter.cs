using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler.NodeLocation;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class CellPairCollecter<T>
        where T : ILocatable
    {
        readonly PeriodicMap map;

        public CellPairCollecter(PeriodicMap map)
        {
            this.map = map;
        }

        enum Position { left, right, leftUnsplit, rightUnsplit };

        public CellPairCollection<T> FollowBoundaryAndCollectCandidates(IEnumerable<Edge<T>> periodicEdges)
        {
            LinkedListDictionary<int, List<Edge<T>>> edges = DivideIntoBoundaries(periodicEdges);
            Debug.Assert(edges.Count % 2 == 0);
            CellPairCollection<T> candidates = new CellPairCollection<T>();


            while (edges.Count > 0)
            {
                List<Edge<T>> boundary = edges.First().Value;
                int boundaryNumber = edges.First().Key;
                Position[] positions = new Position[boundary.Count];
                for (int i = 0; i < boundary.Count; ++i)
                {
                    Edge<T> edge = boundary[i];
                    if (NodeIsOnRightSideOfEdge(edge.Cell.Node, edge))
                    {
                        if(edge.Twin.Cell.Node != null && NodeIsOnRightSideOfEdge(edge.Twin.Cell.Node, edge))
                        {
                            positions[i] = Position.right;
                            candidates.AddOuterSplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                        }
                        else
                        {
                            positions[i] = Position.rightUnsplit;
                            candidates.AddOuterUnsplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                        }
                    }
                    else
                    {
                        if (edge.Twin.Cell.Node!= null && !NodeIsOnRightSideOfEdge(edge.Twin.Cell.Node, edge))
                        {
                            positions[i] = Position.left;
                            candidates.AddInnerSplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                        }
                        else
                        {
                            positions[i] = Position.leftUnsplit;
                            candidates.AddInnerUnsplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                        }
                    }
                }
                edges.Remove(boundaryNumber);
                
                boundaryNumber = map.PeriodicBoundaryCorrelation[boundaryNumber];
                boundary = edges[boundaryNumber];
                for (int i = 0 ; i < boundary.Count; ++i)
                {
                    Edge<T> edge = boundary[i];
                    switch(positions[boundary.Count - 1 - i])
                    {
                        case Position.left:
                            candidates.AddOuterSplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                            break;
                        case Position.right:
                            candidates.AddInnerSplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                            break;
                        case Position.leftUnsplit:
                            candidates.AddOuterUnsplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                            break;
                        case Position.rightUnsplit:
                            candidates.AddInnerUnsplitCell(edge.Cell, edge.BoundaryEdgeNumber);
                            break;
                        default:
                            throw new Exception();
                    }
                }
                edges.Remove(boundaryNumber);
            }
            return candidates;
        }

        static LinkedListDictionary<int, List<Edge<T>>> DivideIntoBoundaries(IEnumerable<Edge<T>> periodicEdges)
        {
            LinkedListDictionary<int, List<Edge<T>>> edges = new LinkedListDictionary<int, List<Edge<T>>>();
            LinkedList<Edge<T>> firstEdges = new LinkedList<Edge<T>>();
            
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
                    firstEdges.AddLast(edge);
                }
                else
                {
                    Add(edges, edge);
                }
            }
            foreach(Edge<T> edge in firstEdges)
            {
                Add(edges, edge);
            }
            return edges;
        }

        static void Add(LinkedListDictionary<int, List<Edge<T>>> edges, Edge<T> edge)
        {
            if (edges.TryGetValue(edge.BoundaryEdgeNumber, out List<Edge<T>> boundary))
            {
                boundary.Add(edge);
            }
            else
            {
                List<Edge<T>> newBoundary = new List<Edge<T>>();
                newBoundary.Add(edge);
                edges.Add(edge.BoundaryEdgeNumber, newBoundary);
            }
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
