using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class PeriodicMapGenerator
    {
        public static PeriodicMap GeneratePeriodicMap(MeshingAlgorithm.Settings settings, VoronoiBoundary boundary)
        {
            PeriodicMap map = null;

            IDictionary<int, int> periodicBoundaryMap = ExtractPeriodicBoundaryMap(boundary.EdgeTags);
            if (periodicBoundaryMap.Count > 0)
            {
                IDictionary<int, Transformation> periodicBoundaryTransformationMap = CreatePeriodicTransformationsFrom(
                    settings.Boundary,
                    periodicBoundaryMap);
                map = new PeriodicMap
                {
                    PeriodicBoundaryCorrelation = periodicBoundaryMap,
                    PeriodicBoundaryTransformations = periodicBoundaryTransformationMap,
                };
                AddPeriodicCorners(map, boundary);
            }
            return map;
        }

        static IDictionary<int, int> ExtractPeriodicBoundaryMap(byte[] tags)
        {
            IDictionary<int, int> periodicBoundaryMap = new Dictionary<int, int>();
            IDictionary<byte, int> usedTags = new LinkedListDictionary<byte, int>();
            for (int i = 0; i < tags.Length; ++i)
            {
                byte tag = tags[i];
                if (tag >= GridCommons.FIRST_PERIODIC_BC_TAG && tag < 255)
                {
                    if (usedTags.TryGetValue(tag, out int index))
                    {
                        periodicBoundaryMap.Add(i, index);
                        periodicBoundaryMap.Add(index, i);
                    }
                    else
                    {
                        usedTags.Add(tag, i);
                    };
                }
            }
            return periodicBoundaryMap;
        }

        static IDictionary<int, Transformation> CreatePeriodicTransformationsFrom(
            Vector[] boundary,
            IDictionary<int, int> periodicBoundaryMap)
        {
            BoundaryLine[] boundaryLines = BoundaryLine.ToLines(boundary);
            return CreatePeriodicTransformationsFrom(boundaryLines, periodicBoundaryMap);
        }

        static IDictionary<int, Transformation> CreatePeriodicTransformationsFrom(
            BoundaryLine[] boundary,
            IDictionary<int, int> periodicBoundaryMap)
        {
            IDictionary<int, Transformation> periodicTrafoMap = new LinkedListDictionary<int, Transformation>();
            foreach (var boundaryPair in periodicBoundaryMap)
            {
                BoundaryLine source = boundary[boundaryPair.Key];
                BoundaryLine target = BoundaryLine.GetReverse(boundary[boundaryPair.Value]);

                Transformation transformation = new BoundaryTransformation(source, target);
                periodicTrafoMap.Add(boundaryPair.Key, transformation);
            }
            return periodicTrafoMap;
        }

        static void AddPeriodicCorners(PeriodicMap map, VoronoiBoundary boundary)
        {
            IDictionary<Corner, int> periodicCornerCorrelation = ExtractPeriodicCornerMap(
                boundary.EdgeTags.Length,
                map.PeriodicBoundaryCorrelation);
            map.PeriodicCornerCorrelation = periodicCornerCorrelation;
            AddPeriodicCornerTransformations(periodicCornerCorrelation, map.PeriodicBoundaryTransformations);
            AddPeriodicCornerBoundaryCorrelation(periodicCornerCorrelation, map.PeriodicBoundaryCorrelation);
        }

        static IDictionary<Corner, int> ExtractPeriodicCornerMap(
            int totalEdges,
            IDictionary<int, int> periodicBoundaryMap)
        {
            LinkedListDictionary<Corner, int> periodicCornerMap = new LinkedListDictionary<Corner, int>();
            int cornerEdge = totalEdges;
            for (int edge = 0; edge < totalEdges; ++edge)
            {
                int followingEdge = (edge + 1) % totalEdges;
                if (periodicBoundaryMap.TryGetValue(edge, out int pairedEdge)
                    && periodicBoundaryMap.TryGetValue(followingEdge, out int pairedFollowingEdge))
                {
                    Corner periodicCorner = new Corner
                    {
                        FirstEdge = edge,
                        SecondEdge = followingEdge
                    };
                    periodicCornerMap.Add(periodicCorner, cornerEdge);
                    ++cornerEdge;
                }
            }
            return periodicCornerMap;
        }

        static void AddPeriodicCornerTransformations(
            IDictionary<Corner, int> periodicCornerMap,
            IDictionary<int, Transformation> periodicBoundaryTransformations)
        {
            foreach (var cornerPair in periodicCornerMap)
            {
                Corner corner = cornerPair.Key;
                int boundary = cornerPair.Value;
                periodicBoundaryTransformations.TryGetValue(corner.FirstEdge, out Transformation first);
                periodicBoundaryTransformations.TryGetValue(corner.SecondEdge, out Transformation second);

                Transformation cornerTrafo = Transformation.Combine(first, second);
                periodicBoundaryTransformations.Add(boundary, cornerTrafo);
            }
        }

        static void AddPeriodicCornerBoundaryCorrelation(
            IDictionary<Corner, int> periodicCornerMap,
            IDictionary<int, int> periodicBoundaryCorrelation)
        {
            foreach (var cornerPair in periodicCornerMap)
            {
                Corner corner = cornerPair.Key;
                int boundary = cornerPair.Value;

                periodicBoundaryCorrelation.TryGetValue(corner.SecondEdge, out int firstEdgePaired);
                periodicBoundaryCorrelation.TryGetValue(corner.FirstEdge, out int secondEdgePaired);
                Corner pairedCorner = new Corner
                {
                    FirstEdge = firstEdgePaired,
                    SecondEdge = secondEdgePaired
                };

                if(periodicCornerMap.TryGetValue(pairedCorner, out int pairedBoundary))
                {
                    periodicBoundaryCorrelation.Add(boundary, pairedBoundary);
                }
                else
                {
                    throw new System.Exception("Periodic corner missing.");
                }
            }
        }

        
    }
}
