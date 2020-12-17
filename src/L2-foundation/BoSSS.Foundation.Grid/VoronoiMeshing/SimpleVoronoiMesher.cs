using BoSSS.Foundation.Grid.Voronoi.Meshing.Converter;
using ilPSP;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class Node : ICloneable<Node>, IVoronoiNodeCastable
    {
        readonly VoronoiNode node;

        public Node()
        {
            node = new VoronoiNode();
        }

        public Node(VoronoiNode node)
        {
            this.node = node;
        }

        public Vector Position {
            get { return node.Position; }
            set {
                node.Position = value;
            }
        }

        public VoronoiNode AsVoronoiNode()
        {
            return node;
        }

        Node ICloneable<Node>.Clone()
        {
            return new Node();
        }
    } 

    public class VoronoiMesher : VoronoiMesher<Node>
    {
        public VoronoiMesher(Settings settings) : base(settings) { }

        public VoronoiMesher(VoronoiBoundary boundary) : base(boundary) { }

        List<Node> WrapInMesherNodes(IList<VoronoiNode> voronoiNodes)
        {
            List<Node> wrappedNodes = new List<Node>(voronoiNodes.Count);
            for (int i = 0; i < voronoiNodes.Count; ++i)
            {
                Node wrappedNode = new Node(voronoiNodes[i]);
                wrappedNodes.Add(wrappedNode);
            }
            return wrappedNodes;
        }

        public VoronoiGrid CreateGrid(VoronoiNodes nodes, int firstCornerNodeIndice)
        {
            List<Node> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            VoronoiGrid grid = CreateGrid(mesherNodes, firstCornerNodeIndice);
            return grid;
        }
    }
}
