using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryNodeCloner<T>
        where T : ILocatable, new()
    {
        IDictionary<int, Transformation> periodicTrafoMap;

        public BoundaryNodeCloner(IDictionary<int, Transformation> periodicTrafoMap)
        {
            this.periodicTrafoMap = periodicTrafoMap;
        }

        public List<T> CloneAndMirrorNodesOf(IEnumerable<Edge<T>> edges)
        {
            List<T> clones = new List<T>();
            foreach (Pair<Edge<T>> edgePair in new Convolution<Edge<T>>(edges))
            {
                Transformation transformation = GetBoundaryTransformationOf(edgePair.Current);
                T transformedClone = CloneAndTransFormAssociatedNodeOf(edgePair.Current, transformation);
                clones.Add(transformedClone);

                if (IsCorner(edgePair.Current, edgePair.Previous))
                {
                    Transformation previousTransformation = GetBoundaryTransformationOf(edgePair.Previous);
                    T doublyTransformedClone = CloneAndTransform(transformedClone, previousTransformation);
                    clones.Add(doublyTransformedClone);
                }
            }
            return clones;
        }

        Transformation GetBoundaryTransformationOf(Edge<T> edge)
        {
            Debug.Assert(periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out Transformation debug));

            periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out Transformation transformation);
            return transformation;
        }

        static T CloneAndTransFormAssociatedNodeOf(Edge<T> edge, Transformation transformation)
        {
            MeshCell<T> cell = edge.Cell;
            return CloneAndTransform(cell.Node, transformation);
        }

        static T CloneAndTransform(T node, Transformation transformation)
        {
            T transformedClone = new T
            {
                Position = transformation.Transform(node.Position)
            };
            return transformedClone;
        }

        static bool IsCorner(Edge<T> edge, Edge<T> preceedingEdge)
        {
            bool areDifferentBoundaries = (edge.BoundaryEdgeNumber != (preceedingEdge?.BoundaryEdgeNumber ?? int.MinValue));
            bool ofSameCell = edge.Cell.ID == preceedingEdge.Cell.ID;
            return areDifferentBoundaries & ofSameCell;
        }
    }
}
