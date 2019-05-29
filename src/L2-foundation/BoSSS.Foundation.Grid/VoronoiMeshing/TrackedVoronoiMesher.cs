using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Voronoi;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class TrackedVoronoiGrid
    {
        public VoronoiGrid Result;
        public OneWayArrayMap InputNodesToResultNodes;
        public OneWayArrayMap ResultNodesToInputNodes;
    }

    class TrackableNode : IMesherNode, IVoronoiNodeCastable
    {
        VoronoiNode node;

        public ArrayConnection Type { get; set; }

        public TrackableNode(VoronoiNode node, int j, Connection type)
        {
            this.node = node;
            Type = new ArrayConnection
            {
                J = j,
                Type = type
            };
        }

        public TrackableNode()
        {
            node = new VoronoiNode();
            Type = new ArrayConnection
            {
                J = -1,
                Type = Connection.Created
            };
        }

        public Vector Position {
            get { return node.Position; }
            set { node.Position = value; }
        }

        public VoronoiNode AsVoronoiNode()
        {
            return node;
        }
    }

    public class TrackedVoronoiMesher
    {
        public class Settings
        {
            public VoronoiInfo GridInfo;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        public static TrackedVoronoiGrid CreateGrid(VoronoiNodes nodes, Settings settings)
        {
            BoundaryMesh<TrackableNode> mesh = CreateBoundaryMesh(nodes, settings);

            VoronoiGrid grid = GridConverter.Convert2VoronoiGrid(mesh, settings.GridInfo);
            OneWayArrayMap resultMap = ExtractMap(mesh.GetNodes());
            OneWayArrayMap inputMap = GetInputMap(resultMap, nodes.Count);
            TrackedVoronoiGrid movingGrid = new TrackedVoronoiGrid
            {
                Result = grid,
                InputNodesToResultNodes = inputMap,
                ResultNodesToInputNodes = resultMap
            };
            return movingGrid;
        }

        static BoundaryMesh<TrackableNode> CreateBoundaryMesh(VoronoiNodes nodes, Settings settings)
        {
            Mesher.Settings mesherSettings = ConvertToMesherSettings(settings);
            List<TrackableNode> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            
            BoundaryMesh<TrackableNode> mesh = Mesher.Create(mesherNodes, mesherSettings);
            return mesh;
        }

        static Mesher.Settings ConvertToMesherSettings(Settings settings)
        {
            Mesher.Settings mesherSettings = new Mesher.Settings
            {
                Boundary = settings.GridInfo.Boundary,
                BoundingBox = settings.GridInfo.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
                FirstCellNode_indice = settings.FirstCellNode_indice
            };
            return mesherSettings;
        }

        static List<TrackableNode> WrapInMesherNodes(IList<VoronoiNode> voronoiNodes)
        {
            List<TrackableNode> wrappedNodes = new List<TrackableNode>(voronoiNodes.Count);
            for (int i = 0; i < voronoiNodes.Count; ++i)
            {
                TrackableNode wrappedNode = new TrackableNode(voronoiNodes[i], i, Connection.Remained);
                wrappedNodes.Add(wrappedNode);
            }
            return wrappedNodes;
        }

        static OneWayArrayMap ExtractMap(IList<TrackableNode> processed)
        {
            ArrayConnection[] connections = new ArrayConnection[processed.Count];
            for (int i = 0; i < processed.Count; ++i)
            {
                connections[i] = processed[i].Type;
            }
            OneWayArrayMap backwardsMap = new OneWayArrayMap(connections);
            return backwardsMap;
        }

        static OneWayArrayMap GetInputMap(OneWayArrayMap result, int noOfInitialNodes)
        {
            OneWayArrayMap input = OneWayArrayMap.CreateEmpty(Connection.Removed, noOfInitialNodes);
            input.AddReverse(result);
            return input;
        }
    }
}
