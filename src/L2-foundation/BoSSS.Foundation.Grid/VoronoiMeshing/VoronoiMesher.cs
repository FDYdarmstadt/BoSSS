using System.Collections.Generic;
using System;
using BoSSS.Foundation.Grid.Voronoi.Meshing.Converter;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class VoronoiMesher<T>
        where T : ICloneable<T>, IVoronoiNodeCastable, new()
    {
        public class Settings
        {
            public VoronoiBoundary Boundary;

            public int NumberOfLloydIterations = 10;
        }

        internal IMesh<T> mesh;

        MeshingAlgorithm.State mesherSettings;

        readonly Settings settings;

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
            mesherSettings = new MeshingAlgorithm.State
            {
                Boundary = settings.Boundary.Polygon,
                BoundingBox = settings.Boundary.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
            };
            mesherSettings.PeriodicMap = PeriodicMapGenerator.GeneratePeriodicMap(
                mesherSettings, 
                settings.Boundary);
        }

        void CreateGridConverter()
        {
            if(mesherSettings.PeriodicMap != null)
            {
                gridConverter = new GridConverter<T>(settings.Boundary, mesherSettings.PeriodicMap);
            }
            else
            {
                gridConverter = new GridConverter<T>(settings.Boundary);
            }
        }

        protected void CreateMesh(List<T> nodes, int firstCornerNodeIndice = 0)
        {
            mesherSettings.FirstCellNodeIndice = firstCornerNodeIndice;
            mesh = MeshingAlgorithm.ComputeMesh(nodes, mesherSettings);
        }

        public VoronoiGrid CreateGrid(List<T> nodes, int firstCornerNodeIndice)
        {
            CreateMesh(nodes, firstCornerNodeIndice); 
            VoronoiGrid grid = gridConverter.ConvertToVoronoiGrid(mesh);
            return grid;
        }
    }
}
