﻿using BoSSS.Foundation.Voronoi;
using BoSSS.Platform.LinAlg;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class Node : ILocatable, IVoronoiNodeCastable
    {
        VoronoiNode node;

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
            set { node.Position = value; }
        }

        public VoronoiNode AsVoronoiNode()
        {
            return node;
        }
    } 

    public class VoronoiMesher : VoronoiMesher<Node>
    {
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

        public VoronoiGrid CreateGrid(VoronoiNodes nodes, Settings settings)
        {
            List<Node> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            VoronoiGrid grid = CreateGrid(mesherNodes, settings);
            return grid;
        }
    }
}
