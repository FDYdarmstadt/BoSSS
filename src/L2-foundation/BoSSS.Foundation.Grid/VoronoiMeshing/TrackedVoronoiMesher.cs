using BoSSS.Foundation.Voronoi;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class TrackedVoronoiGrid
    {
        public VoronoiGrid Result;
        public OneWayArrayMap InputNodesToResultNodes;
        public OneWayArrayMap ResultNodesToInputNodes;
    }

    public class TrackableNode : Node
    {
        public ArrayConnection Type { get; set; }

        public TrackableNode(VoronoiNode node, int j, Connection type)
            : base(node)
        {
            Type = new ArrayConnection
            {
                J = j,
                Type = type
            };
        }

        public TrackableNode()
            :base()
        {
            Type = new ArrayConnection
            {
                J = -1,
                Type = Connection.Created
            };
        }
    }

    public class TrackedVoronoiMesher : Mesher<TrackableNode>
    {
        public TrackedVoronoiGrid CreateGrid(VoronoiNodes nodes, Settings settings)
        {
            List<TrackableNode> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            BoundaryMesh<TrackableNode> mesh = CreateMesh(mesherNodes, settings);
            VoronoiGrid grid = Convert2VoronoiGrid(mesh, settings);

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
