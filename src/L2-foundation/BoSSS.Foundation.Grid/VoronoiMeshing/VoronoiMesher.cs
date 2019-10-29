using System.Collections.Generic;
using System;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class VoronoiMesher<T>
        where T : ILocatable, IVoronoiNodeCastable, new()
    {
        public class Settings
        {
            public VoronoiBoundary Boundary;

            public int NumberOfLloydIterations = 10;

            public int FirstCellNode_indice = 0;
        }

        internal Mesh<T> mesh;

        MeshingAlgorithm.Settings mesherSettings;

        Settings settings;

        GridConverter<T> gridConverter; 

        public VoronoiMesher(VoronoiBoundary boundary)
        {
            settings = new Settings()
            {
                Boundary = boundary
            };
            Initialize();
        }

        public VoronoiMesher(Settings settings)
        {
            this.settings = settings;
            Initialize();
        }

        void Initialize()
        {
            CreateSetupForMeshingAlgorithm();
            CreateGridConverter();
        }

        void CreateSetupForMeshingAlgorithm()
        {
            mesherSettings = new MeshingAlgorithm.Settings
            {
                Boundary = settings.Boundary.Polygon,
                BoundingBox = settings.Boundary.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
                FirstCellNode_indice = settings.FirstCellNode_indice,
            };
            PeriodicDataExtracter periodicDataExtracter = new PeriodicDataExtracter(settings.Boundary);
            periodicDataExtracter.InitializePeriodic(mesherSettings);
        }

        void CreateGridConverter()
        {
            gridConverter = new GridConverter<T>(settings.Boundary);
        }

        protected void CreateMesh(List<T> nodes)
        {
            mesh = MeshingAlgorithm.ComputeMesh(nodes, mesherSettings);
        }

        public VoronoiGrid CreateGrid(List<T> nodes)
        {
            CreateMesh(nodes); 
            VoronoiGrid grid = gridConverter.Convert(mesh);
            return grid;
        }
    }
}
