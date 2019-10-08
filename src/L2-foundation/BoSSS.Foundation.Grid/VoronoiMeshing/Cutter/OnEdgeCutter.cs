using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
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

        public IEnumerator<Edge<T>> CutOut(
            IEnumerator<Edge<T>> edgeEnum,
            Edge<T> outerEdge)
        {
            bool keepOnRunning = true;
            Edge<T> activeEdge = default(Edge<T>);
            BoundaryLine activeLine = default(BoundaryLine);
            while (keepOnRunning)
            {
                keepOnRunning = false;
                IntersectionCase intersectionCase = IntersectionCase.NotIntersecting;
                double alphaCut = 0;
                activeLine = boundary.Current;
                while (edgeEnum.MoveNext() && !keepOnRunning)
                {
                    activeEdge = edgeEnum.Current;
                    keepOnRunning = LineIntersection.Intersect(activeEdge, activeLine, ref intersectionCase, out alphaCut);
                }
                switch (intersectionCase)
                {
                    case IntersectionCase.NotIntersecting:
                        break;
                    case IntersectionCase.EndOfLine:
                        activeEdge = meshIntersecter.AddLineSegment(activeEdge, alphaCut, boundary.LineIndex);
                        edgeEnum = meshIntersecter.GetConnectedRidgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.EndOfRidge:
                        meshIntersecter.AddEdge(activeEdge, boundary.LineIndex);
                        edgeEnum = meshIntersecter.GetConnectedRidgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        break;
                    case IntersectionCase.EndOfRidgeAndLine:
                        meshIntersecter.AddEdge(activeEdge, boundary.LineIndex);
                        edgeEnum = meshIntersecter.GetConnectedRidgeEnum(activeEdge);
                        outerEdge = activeEdge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.InMiddle:
                    default:
                        throw new InvalidOperationException();
                }
            }
            edgeEnum = meshIntersecter.getNeighborFromLineDirection(outerEdge, activeLine);
            return edgeEnum;
        }
    }
}
