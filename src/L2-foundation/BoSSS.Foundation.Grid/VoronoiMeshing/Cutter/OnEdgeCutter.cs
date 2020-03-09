using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Cutter
{
    class OnEdgeCutter<T>
        where T : ILocatable, new()
    {
        MeshIntersecter<T> meshIntersecter;

        BoundaryLineEnumerator boundary;

        public OnEdgeCutter(
            MeshIntersecter<T> meshIntersecter,
            BoundaryLineEnumerator boundary)
        {
            this.meshIntersecter = meshIntersecter;
            this.boundary = boundary;
        }

        public (IEnumerator<Edge<T>>, IntersectionCase ) CutOut(
            IEnumerator<Edge<T>> edgeEnum,
            Edge<T> outerEdge)
        {
            bool keepOnRunning = true;
            Edge<T> activeEdge = default(Edge<T>);
            BoundaryLine activeLine = default(BoundaryLine);
            IntersectionCase intersectionCase = IntersectionCase.NotIntersecting;
            while (keepOnRunning)
            {
                keepOnRunning = false;
                double alphaCut = 0;
                activeLine = boundary.Current;
                while (edgeEnum.MoveNext() && !keepOnRunning)
                {
                    activeEdge = edgeEnum.Current;
                    keepOnRunning = LineIntersect.Find(activeEdge, activeLine, ref intersectionCase, out alphaCut);
                }
                switch (intersectionCase)
                {
                    case IntersectionCase.NotIntersecting:
                        break;
                    case IntersectionCase.EndOfLine:
                        activeEdge = meshIntersecter.AddLineSegment(activeEdge, alphaCut, boundary.LineIndex());
                        edgeEnum = meshIntersecter.GetConnectedEdgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.EndOfEdge:
                        meshIntersecter.AddEdge(activeEdge, boundary.LineIndex());
                        edgeEnum = meshIntersecter.GetConnectedEdgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        break;
                    case IntersectionCase.EndOfEdgeAndLine:
                        meshIntersecter.AddEdge(activeEdge, boundary.LineIndex());
                        edgeEnum = meshIntersecter.GetConnectedEdgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.StartOfLine:
                    case IntersectionCase.InMiddle:
                    default:
                        throw new InvalidOperationException();
                }
            }
            edgeEnum = meshIntersecter.GetNeighborFromLineDirection(outerEdge, activeLine);
            return (edgeEnum, intersectionCase);
        }
    }
}
