using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    static class NodeLocation
    {
        public static bool NodeIsOnRightSideOfEdge<T>(T Node, Edge<T> edge)
            where T : ILocatable
        {
            return IsOnRightSide(Node.Position, edge);
        }

        public static bool NodeOfEdgeIsOnRightSideOfEdge<T>(Edge<T> edge)
            where T : ILocatable
        {
            Vector position = edge.Cell.Node.Position;
            return IsOnRightSide(position, edge);
        }

        static bool IsOnRightSide(Vector node, Line line)
        {
            Vector start = line.Start.Position;
            Vector end = line.End.Position;
            return ((end.x - start.x) * (node.y - start.y) - (end.y - start.y) * (node.x - start.x)) < 0;
        }
    }
}
