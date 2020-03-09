using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryNodeMirrorer<T>
        where T : ICloneable<T>
    {
        readonly PeriodicMap map;

        public BoundaryNodeMirrorer(
            PeriodicMap map)
        {
            this.map = map;
        }

        ICollection<int> visited;

        List<T> clones;

        public List<T> CloneAndMirrorNodesOf(IEnumerable<Edge<T>> edges)
        {
            clones = new List<T>();
            visited = new LinkedList<int>();
            foreach (Pair<Edge<T>> edgePair in new Convolution<Edge<T>>(edges))
            {
                Transformation transformation = GetBoundaryTransformationOf(edgePair.Current);
                T transformedClone = CloneAndTransFormAssociatedNodeOf(edgePair.Current, transformation);
                clones.Add(transformedClone);

                if (IsCorner(edgePair.Current, edgePair.Previous))
                {
                    Transformation previousTransformation = GetBoundaryTransformationOf(edgePair.Previous);
                    Transformation edgeTransformation = Transformation.Combine(transformation, previousTransformation);

                    foreach(MeshCell<T> cell in GetSurroundingCells(edgePair.Previous))
                    {
                        TryToCloneAndTransform(cell, edgeTransformation);
                    }
                    foreach (MeshCell<T> cell in GetSurroundingCells(edgePair.Current))
                    {
                        TryToCloneAndTransform(cell, edgeTransformation);
                    }
                    visited.Clear();
                }
            }
            return clones;
        }

        Transformation GetBoundaryTransformationOf(Edge<T> edge)
        {
            Debug.Assert(map.PeriodicBoundaryTransformations.TryGetValue(edge.BoundaryEdgeNumber, out Transformation debug));

            map.PeriodicBoundaryTransformations.TryGetValue(edge.BoundaryEdgeNumber, out Transformation transformation);
            return transformation;
        }

        static T CloneAndTransFormAssociatedNodeOf(Edge<T> edge, Transformation transformation)
        {
            MeshCell<T> cell = edge.Cell;
            return CloneAndTransform(cell.Node, transformation);
        }

        static T CloneAndTransform(T node, Transformation transformation)
        {
            T transformedClone = node.Clone();
            transformedClone.Position = transformation.Transform(node.Position);
            
            return transformedClone;
        }

        static bool IsCorner(Edge<T> edge, Edge<T> preceedingEdge)
        {
            bool areDifferentBoundaries = (edge.BoundaryEdgeNumber != (preceedingEdge?.BoundaryEdgeNumber ?? int.MinValue));
            bool ofSameCell = edge.Cell.ID == preceedingEdge.Cell.ID;
            return areDifferentBoundaries & ofSameCell;
        }

        void TryToCloneAndTransform(MeshCell<T> cell, params Transformation[] trafos)
        {
            if (!visited.Contains(cell.ID))
            {
                foreach (Transformation trafo in trafos)
                {
                    CloneAndTransform(cell.Node, clones, trafo);
                }
                visited.Add(cell.ID);
            }
        }

        static void CloneAndTransform(T node, List<T> clones, Transformation trafo)
        {
            T transformedClone = CloneAndTransform(node, trafo);
            clones.Add(transformedClone);
        }

        static IEnumerable<MeshCell<T>> GetSurroundingCells(Edge<T> boundary)
        {
            yield return boundary.Cell;
            foreach(Edge<T> edge in boundary.Cell.Edges)
            {
                if (!edge.IsBoundary)
                {
                    yield return edge.Twin.Cell;
                }
            }
        }
    }
}
