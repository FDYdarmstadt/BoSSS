using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.Linq;
using System.Resources;
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
            Edge<T> outerEdge,
            CutterState<Edge<T>> state)
        {
            bool keepOnRunning = true;
            Edge<T> activeEdge = default(Edge<T>);
            IntersectionCase intersectionCase = IntersectionCase.NotIntersecting;
            while (keepOnRunning)
            {
                keepOnRunning = false;
                double alphaCut = 0;
                state.ActiveLine = boundary.Current;
                while (edgeEnum.MoveNext() && !keepOnRunning)
                {
                    activeEdge = edgeEnum.Current;
                    keepOnRunning = LineIntersect.Find(activeEdge, state.ActiveLine, ref intersectionCase, out alphaCut);
                    if(alphaCut < LineIntersect.accuracy)
                    {
                        intersectionCase = IntersectionCase.NotIntersecting;
                        keepOnRunning = false;
                    }
                }
                switch (intersectionCase)
                {
                    case IntersectionCase.NotIntersecting:
                    case IntersectionCase.InMiddle:
                        //keepOnRunning = false;
                        //IEnumerator<Edge<T>> cellEnum = meshIntersecter.GetAfterCutEdgeEnumerator(state.ActiveEdge.Cell.Edges, state.ActiveEdge);
                        //return (cellEnum, intersectionCase);
                        keepOnRunning = false;
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
                        //edgeEnum.MoveNext();
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
                    default:
                        throw new InvalidOperationException();
                }
            }
            MeshIntersecter<T>.AfterCutEdgeEnumerator cellEnumerator = meshIntersecter.GetNeighborFromLineDirection(outerEdge, state.ActiveLine);
            cellEnumerator.Cell.IntersectionVertex = outerEdge.End.ID;
            return (cellEnumerator, intersectionCase);
        }
    }
}
