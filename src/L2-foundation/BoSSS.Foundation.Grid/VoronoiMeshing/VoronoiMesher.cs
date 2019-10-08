using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class VoronoiMesher<T>
        where T : ILocatable, IVoronoiNodeCastable, new()
    {
        internal Mesh<T> mesh;

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
                FirstCellNode_indice = settings.FirstCellNode_indice,
                PeriodicBoundaryMap = ExtractPeriodicBoundaryMap(settings.Boundary)
            };
            return mesherSettings;
        }

        IDictionary<int,int> ExtractPeriodicBoundaryMap(VoronoiBoundary boundary)
        {
            byte[] tags = boundary.EdgeTags;
            IDictionary<int, int> periodicBoundaryMap = new Dictionary<int, int>();
            IDictionary<byte, int> usedTags = new LinkedListDictionary<byte, int>(); 
            for(int i = 0; i < tags.Length; ++i)
            {
                byte tag = tags[i];
                if(tag >= GridCommons.FIRST_PERIODIC_BC_TAG && tag < 255)
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

        protected void CreateMesh(List<T> nodes, Settings settings)
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
