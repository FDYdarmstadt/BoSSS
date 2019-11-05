using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class CornerCell<T>
    {
        public MeshCell<T> Cell;

        public Transformation CornerTransformation;

        public int BoundaryEdgeNumber;
    }

    class CornerMapper<T>
    {
        struct Corner : IEquatable<Corner>
        {
            public int ID;

            public int boundaryEdgeNumber;

            public bool Equals(Corner other)
            {
                return other.ID == ID;
            }
        }

        readonly ICollection<Corner> visitedCorners;

        readonly PeriodicMap map;

        readonly Queue<CornerCell<T>> cornerCells;

        readonly CornerConnector<T> cornerConnector;

        public CornerMapper(PeriodicMap map)
        {
            visitedCorners = new LinkedList<Corner>();
            this.map = map;
            cornerCells = new Queue<CornerCell<T>>();
            cornerConnector = new CornerConnector<T>();
        }

        public void ConnectPeriodicCorners(CellPairCollection<T> candidates)
        {
            FindCorners(candidates);
            cornerConnector.ConnectCorners(cornerCells);
        }

        void FindCorners(CellPairCollection<T> candidates)
        {
            foreach (CellPairCollection<T>.EdgeCombo edge in candidates.GetCollectedEdgeCombos())
            {
                AddCorner(edge.Inner[0], edge.EdgeNumber);
                if (edge.Inner.Count > 1)
                {
                    AddCorner(edge.Inner[edge.Inner.Count - 1], edge.EdgeNumber);
                }
            }
        }

        void AddCorner(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            Corner corner = new Corner
            {
                ID = cell.ID,
                boundaryEdgeNumber = boundaryEdgeNumber
            };
            if (visitedCorners.Contains(corner))
            {
                visitedCorners.Remove(corner);
                Transformation first = map.PeriodicTransformationMap[corner.boundaryEdgeNumber];
                Transformation second = map.PeriodicTransformationMap[boundaryEdgeNumber];
                Transformation cornerTrafo = Transformation.Combine(first, second);

                int newBoundaryNumber = AddTransformationToPeriodicMap(cornerTrafo, map);
                CornerCell<T> cornerCell = new CornerCell<T>
                {
                    Cell = cell,
                    CornerTransformation = cornerTrafo,
                    BoundaryEdgeNumber = newBoundaryNumber
                };
                cornerCells.Enqueue(cornerCell);
            }
            else {
                visitedCorners.Add(corner);
            }
        }

        static int AddTransformationToPeriodicMap(Transformation transformation, PeriodicMap map)
        {
            int newEdgeNumber = map.PeriodicBoundaryMap.Count;
            while (map.PeriodicBoundaryMap.ContainsKey(newEdgeNumber))
            {
                ++newEdgeNumber;
            }
            int newEdgeNumberPair = newEdgeNumber + 1;
            while (map.PeriodicBoundaryMap.ContainsKey(newEdgeNumberPair))
            {
                ++newEdgeNumberPair;
            }
            map.PeriodicTransformationMap.Add(newEdgeNumber, transformation);
            map.PeriodicTransformationMap.Add(newEdgeNumberPair, Transformation.GetReverse(transformation));
            map.PeriodicBoundaryMap.Add(newEdgeNumber, newEdgeNumberPair);
            map.PeriodicBoundaryMap.Add(newEdgeNumberPair, newEdgeNumber);
            return newEdgeNumber;
        }
    }

    class CornerConnector<T>
    {
        public void  ConnectCorners(Queue<CornerCell<T>> cornerCells)
        {
            while (cornerCells.Count > 0)
            {
                CornerCell<T> cornerCell = cornerCells.Dequeue();
                SetCornerEdge(cornerCell);
            }
        }

        void SetCornerEdge(CornerCell<T> cornerCell)
        {

        }

        public IEnumerable<Edge<T>> FindEdgesThatArePartOfTwoBoundaries(CornerCell<T> cornerCell)
        {
            throw new NotSupportedException();
        }
    }

}
