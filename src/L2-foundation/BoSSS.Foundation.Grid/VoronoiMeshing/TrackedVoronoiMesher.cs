using BoSSS.Foundation.Voronoi;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class MappedVoronoiGrid
    {
        public VoronoiGrid Result;
        public ConnectionMap InputNodesToResultNodes;
        public ConnectionMap ResultNodesToInputNodes;
    }

    public struct NodeConnection
    {
        public Connection Type;

        public int J;
    }

    public class TrackableNode : Node
    {
        public NodeConnection Type { get; set; }

        public TrackableNode(VoronoiNode node, int j, Connection type)
            : base(node)
        {
            Type = new NodeConnection
            {
                J = j,
                Type = type
            };
        }

        public TrackableNode()
            :base()
        {
            Type = new NodeConnection
            {
                J = -1,
                Type = Connection.Created
            };
        }
    }

    public class NodeTrackingVoronoiMesher : VoronoiMesher<TrackableNode>
    {
        public MappedVoronoiGrid CreateGrid(VoronoiNodes nodes, Settings settings)
        {
            List<TrackableNode> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            VoronoiGrid grid = CreateGrid(mesherNodes, settings);

            ConnectionMap resultMap = ExtractMap(base.mesh.Nodes);
            ConnectionMap inputMap = GetInputMap(resultMap, nodes.Count);
            MappedVoronoiGrid movingGrid = new MappedVoronoiGrid
            {
                Result = grid,
                InputNodesToResultNodes = inputMap,
                ResultNodesToInputNodes = resultMap
            };
            return movingGrid;
        }

        List<TrackableNode> WrapInMesherNodes(IList<VoronoiNode> voronoiNodes)
        {
            List<TrackableNode> wrappedNodes = new List<TrackableNode>(voronoiNodes.Count);
            for (int i = 0; i < voronoiNodes.Count; ++i)
            {
                TrackableNode wrappedNode = new TrackableNode(voronoiNodes[i], i, Connection.Remained);
                wrappedNodes.Add(wrappedNode);
            }
            return wrappedNodes;
        }

        static ConnectionMap ExtractMap(IList<TrackableNode> processed)
        {
            Connection[] connections = new Connection[processed.Count];
            int[] map = new int[processed.Count];

            for (int i = 0; i < processed.Count; ++i)
            {
                connections[i] = processed[i].Type.Type;
                map[i] = processed[i].Type.J;
            }

            ConnectionMap backwardsMap = new ConnectionMap(connections, map);
            return backwardsMap;
        }

        static ConnectionMap GetInputMap(ConnectionMap result, int noOfInitialNodes)
        {
            ConnectionMap input = ConnectionMap.CreateEmpty(Connection.Removed, noOfInitialNodes);
            input.AddReverse(result);
            return input;
        }
    }
}
