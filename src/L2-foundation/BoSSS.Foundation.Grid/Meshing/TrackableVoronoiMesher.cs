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

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class TrackableVoronoiGrid
    {
        public VoronoiGrid Result;
        public ArrayMap IncomingNodesToResultNodes;
        public ArrayMap ResultNodesToIncomingNodes;
    }

    class TrackableNode : IMesherNode, IVoronoiNodeCastable
    {
        VoronoiNode node;

        public ArrayConnection Type { get; set; }

        public TrackableNode(VoronoiNode node, int j)
        {
            this.node = node;
            Type = new ArrayConnection
            {
                J = j,
                Type = Connection.Remained
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

    public class TrackableVoronoiMesher
    {
        public class Settings
        {
            public VoronoiInfo GridInfo;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        public static TrackableVoronoiGrid CreateGrid(VoronoiNodes nodes, Settings settings)
        {
            BoundaryMesh<TrackableNode> mesh = CreateBoundaryMesh(nodes, settings);

            ArrayMap node2nodeMap = ExtractMapping(mesh.GetNodes());
            VoronoiGrid grid = GridConverter.Convert2VoronoiGrid(mesh, settings.GridInfo);
            TrackableVoronoiGrid movingGrid = new TrackableVoronoiGrid
            {
                Result = grid,
                ResultNodesToIncomingNodes = node2nodeMap,
                IncomingNodesToResultNodes = 
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
                TrackableNode wrappedNode = new TrackableNode(voronoiNodes[i], i);
                wrappedNodes.Add(wrappedNode);
            }
            return wrappedNodes;
        }

        static ArrayMap ExtractMapping(IList<TrackableNode> processed)
        {
            ArrayConnection[] connections = new ArrayConnection[processed.Count];
            for (int i = 0; i < processed.Count; ++i)
            {
                connections[i] = processed[i].Type;
            }
            ArrayMap backwardsMap = new ArrayMap(connections);
            return backwardsMap;
        }
    }
}
