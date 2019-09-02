using BoSSS.Foundation.Voronoi;
using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class VoronoiMesher<T>
        where T : IMesherNode, IVoronoiNodeCastable, new()
    {
        internal BoundaryMesh<T> mesh;

        public class Settings
        {
            public VoronoiBoundary Boundary;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        MeshingAlgorithm.Settings ConvertToMeshingAlgoSettings(Settings settings)
        {
            MeshingAlgorithm.Settings mesherSettings = new MeshingAlgorithm.Settings
            {
                Boundary = settings.Boundary.Polygon,
                BoundingBox = settings.Boundary.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
                FirstCellNode_indice = settings.FirstCellNode_indice
            };
            return mesherSettings;
        }

        void CreateMesh(List<T> nodes, Settings settings)
        {
            MeshingAlgorithm.Settings meshingSettings = ConvertToMeshingAlgoSettings(settings);
            mesh = MeshingAlgorithm.ComputeMesh(nodes, meshingSettings);
        }

        public VoronoiGrid CreateGrid(List<T> nodes, Settings settings)
        {
            CreateMesh(nodes, settings); 
            VoronoiGrid grid = GridConverter.Convert2VoronoiGrid(mesh, settings.Boundary);
            return grid;
        }
    }
}
