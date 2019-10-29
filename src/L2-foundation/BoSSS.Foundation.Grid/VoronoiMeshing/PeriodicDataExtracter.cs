using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class PeriodicDataExtracter
    {
        VoronoiBoundary boundary;

        public PeriodicDataExtracter(VoronoiBoundary boundary)
        {
            this.boundary = boundary;
        }

        public void InitializePeriodic(MeshingAlgorithm.Settings settings)
        {
            IDictionary<int, int> periodicBoundaryMap = ExtractPeriodicBoundaryMap(boundary);
            if (periodicBoundaryMap.Count > 0)
            {
                settings.PeriodicBoundaryMap = periodicBoundaryMap;
                settings.PeriodicTransformations = CreatePeriodicTransformationMapFrom(settings.Boundary, periodicBoundaryMap);
            }
        }

        static IDictionary<int, int> ExtractPeriodicBoundaryMap(VoronoiBoundary boundary)
        {
            byte[] tags = boundary.EdgeTags;
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

        static IDictionary<int, Transformation> CreatePeriodicTransformationMapFrom(
            Vector[] boundary,
            IDictionary<int, int> periodicBoundaryMap)
        {
            BoundaryLine[] boundaryLines = BoundaryLine.ToLines(boundary);
            return CreatePeriodicTransformationMapFrom(boundaryLines, periodicBoundaryMap);
        }

        static IDictionary<int, Transformation> CreatePeriodicTransformationMapFrom(
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
    }
}
