using BoSSS.Foundation.Voronoi;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class Mesher
    {
        public class Settings
        {
            public VoronoiInfo GridInfo;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        MeshingAlgorithm.Settings ConvertToMesherSettings(Settings settings)
        {
            MeshingAlgorithm.Settings mesherSettings = new MeshingAlgorithm.Settings
            {
                Boundary = settings.GridInfo.Boundary,
                BoundingBox = settings.GridInfo.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
                FirstCellNode_indice = settings.FirstCellNode_indice
            };
            return mesherSettings;
        }

        internal BoundaryMesh<T> CreateMesh<T>(List<T> nodes, Settings settings)
            where T : IMesherNode, new()
        {
            MeshingAlgorithm.Settings meshingSettings = ConvertToMesherSettings(settings);
            BoundaryMesh<T> mesh = MeshingAlgorithm.ComputeMesh(nodes, meshingSettings);
            return mesh;
        }
    }
}
