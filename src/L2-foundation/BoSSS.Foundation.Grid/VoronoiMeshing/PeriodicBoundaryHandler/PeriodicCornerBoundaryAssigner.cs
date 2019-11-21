using BoSSS.Platform.LinAlg;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicCornerBoundaryAssigner<T>
    {
        readonly PeriodicMap map;

        public PeriodicCornerBoundaryAssigner(PeriodicMap map)
        {
            this.map = map;
        }

        public void AssignBoundariesOfPeriodicCorners(Queue< Edge<T>> periodicEdges, Corner periodicCorner)
        {
            while (periodicEdges.Count > 0)
            {
                Edge<T> current = periodicEdges.Dequeue();
                AssignEdge(current, periodicCorner);
            }
        }

        void AssignEdge(Edge<T> edge, Corner corner)
        {
            Corner twin = CreateTwinOf(corner);

            map.PeriodicCornerCorrelation.TryGetValue(corner, out int boundary);
            edge.BoundaryEdgeNumber = boundary;
            map.PeriodicCornerCorrelation.TryGetValue(twin, out int twinBoundary);
            edge.Twin.BoundaryEdgeNumber = twinBoundary;

            Debug.Assert(EdgeIsConnectedByTransformation(map.PeriodicBoundaryTransformations[boundary], edge));
        }

        Corner CreateTwinOf(Corner corner)
        {
            map.PeriodicBoundaryCorrelation.TryGetValue(corner.FirstEdge, out int cornerFirstTwin);
            map.PeriodicBoundaryCorrelation.TryGetValue(corner.SecondEdge, out int cornerSecondTwin);
            Corner currentTwinCorner = new Corner
            {
                FirstEdge = cornerFirstTwin,
                SecondEdge = cornerSecondTwin
            };

            Debug.Assert(AreRegisteredInMap(corner, currentTwinCorner));

            return currentTwinCorner;
        }

        bool AreRegisteredInMap(Corner From, Corner To)
        {
            bool areCorners = map.PeriodicCornerCorrelation.ContainsKey(From)
                && map.PeriodicCornerCorrelation.ContainsKey(To);
            return areCorners;
        }

        static bool EdgeIsConnectedByTransformation(Transformation transformation, Edge<T> edge)
        {
            Vector transformedStart = transformation.Transform(edge.Start.Position);
            double distance = Vector.Dist(transformedStart, edge.Twin.End.Position);
            return distance <= 1e-13;
        }
    }
}
