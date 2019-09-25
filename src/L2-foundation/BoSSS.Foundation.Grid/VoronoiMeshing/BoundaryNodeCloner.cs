using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Voronoi;
using ilPSP;
using System.Diagnostics;
using System.Collections.Specialized;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryNodeCloner<T>
        where T : IMesherNode, new()
    {
        IDictionary<int, BoundaryTransformation> periodicTrafoMap;

        public BoundaryNodeCloner(IDictionary<int, BoundaryTransformation> periodicTrafoMap)
        {
            this.periodicTrafoMap = periodicTrafoMap;
        }

        public List<T> CloneAndMirrorNodesOf(IEnumerable<Edge<T>> edges)
        {
            List<T> clones = new List<T>();
            foreach (Pair<Edge<T>> edgePair in new Convolution<Edge<T>>(edges))
            {
                BoundaryTransformation transformation = GetBoundaryTransformationOf(edgePair.Current);
                T transformedClone = CloneAndTransFormAssociatedNodeOf(edgePair.Current, transformation);
                clones.Add(transformedClone);

                if (IsCorner(edgePair.Current, edgePair.Previous))
                {
                    BoundaryTransformation previousTransformation = GetBoundaryTransformationOf(edgePair.Previous);
                    T doublyTransformedClone = CloneAndTransform(transformedClone, previousTransformation);
                    clones.Add(doublyTransformedClone);
                }
            }
            return clones;
        }

        BoundaryTransformation GetBoundaryTransformationOf(Edge<T> edge)
        {
            Debug.Assert(periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out BoundaryTransformation debug));

            periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out BoundaryTransformation transformation);
            return transformation;
        }

        static T CloneAndTransFormAssociatedNodeOf(Edge<T> edge, BoundaryTransformation transformation)
        {
            
            MeshCell<T> cell = edge.Cell;
            return CloneAndTransform(cell.Node, transformation);
        }

        static T CloneAndTransform(T node, BoundaryTransformation transformation)
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
