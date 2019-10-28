using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class LineIntersect
    {
        const double accuracy = 1e-10;

        //Return end of edge if parallel and overlapping.
        public static bool FindFirst(Line edge, Line line, ref IntersectionCase intersectionCase, out double alpha)
        {
            bool notParallel = PolygonClipping.ComputeIntersection(
                edge.Start.Position,
                edge.End.Position,
                line.Start.Position,
                line.End.Position,
                out alpha,
                out double alpha2,
                out Vector vector);
            if (notParallel)
            {
                if (1 >= alpha2)
                {
                    if (accuracy <= alpha2)
                    {
                        if (accuracy <= alpha && 1 >= alpha)
                        {
                            if (alpha > 1 - accuracy)
                            {
                                alpha = 1;
                                intersectionCase = IntersectionCase.EndOfEdge;
                                if (alpha2 > 1 - accuracy)
                                {
                                    intersectionCase = IntersectionCase.EndOfEdgeAndLine;
                                }
                                return true;
                            }
                            else
                            {
                                if (alpha2 > 1 - accuracy)
                                {
                                    intersectionCase = IntersectionCase.EndOfLine;
                                }
                                else
                                {
                                    intersectionCase = IntersectionCase.InMiddle;
                                }
                                return true;
                            }
                        }
                    }
                    else
                    {
                        if (alpha2 > -accuracy)
                        {
                            intersectionCase = IntersectionCase.StartOfLine;
                            return true;
                        }
                    }
                }
                intersectionCase = IntersectionCase.NotIntersecting;
                return false;
            }
            else
            {
                if (Math.Abs(alpha) <= accuracy)
                {
                    double absRidge = (edge.Start.Position - edge.End.Position).L2Norm();
                    double absLine = (edge.Start.Position - line.End.Position).L2Norm();
                    //Only when in same direction!
                    double add = (edge.Start.Position - edge.End.Position + edge.Start.Position - line.End.Position).L2Norm();
                    if (add < absLine + absRidge)
                    {
                        return false;
                    }

                    alpha = absLine / absRidge;
                    if (alpha >= 1 - accuracy && alpha <= 1)
                    {
                        alpha = 1;
                        intersectionCase = IntersectionCase.EndOfEdgeAndLine;
                        return true;
                    }
                    if (absLine / absRidge < 1)
                    {
                        intersectionCase = IntersectionCase.EndOfLine;
                        return true;
                    }
                    alpha = 1;
                    intersectionCase = IntersectionCase.EndOfEdge;
                    return true;
                }
                intersectionCase = IntersectionCase.NotIntersecting;
                return false;
            }
        }

        public static bool Find(Line edge, Line line, ref IntersectionCase intersectionCase, out double alpha)
        {
            bool notParallel = PolygonClipping.ComputeIntersection(
                edge.Start.Position,
                edge.End.Position,
                line.Start.Position,
                line.End.Position,
                out alpha,
                out double alpha2,
                out Vector vector);
            if (notParallel)
            {
                if (1 + accuracy >= alpha2)
                {
                    if (accuracy <= alpha2)
                    {
                        if (accuracy <= alpha && 1 + accuracy >= alpha)
                        {
                            if (alpha > 1 - accuracy)
                            {
                                alpha = 1;
                                intersectionCase = IntersectionCase.EndOfEdge;
                                if (alpha2 > 1 - accuracy)
                                {
                                    intersectionCase = IntersectionCase.EndOfEdgeAndLine;
                                }
                                return true;
                            }
                            else
                            {
                                if (alpha2 > 1 - accuracy)
                                {
                                    intersectionCase = IntersectionCase.EndOfLine;
                                }
                                else
                                {
                                    intersectionCase = IntersectionCase.InMiddle;
                                }
                                return true;
                            }
                        }
                    }
                }
                intersectionCase = IntersectionCase.NotIntersecting;
                return false;
            }
            else
            {
                if (Math.Abs(alpha) <= accuracy)
                {
                    double absRidge = (edge.Start.Position - edge.End.Position).L2Norm();
                    double absLine = (edge.Start.Position - line.End.Position).L2Norm();
                    //Only when in same direction!
                    double add = (edge.Start.Position - edge.End.Position + edge.Start.Position - line.End.Position).L2Norm();
                    if (add < absLine + absRidge)
                    {
                        return false;
                    }

                    alpha = absLine / absRidge;
                    if (alpha >= 1 - accuracy && alpha <= 1)
                    {
                        alpha = 1;
                        intersectionCase = IntersectionCase.EndOfEdgeAndLine;
                        return true;
                    }
                    if (absLine / absRidge < 1)
                    {
                        intersectionCase = IntersectionCase.EndOfLine;
                        return true;
                    }
                    alpha = 1;
                    intersectionCase = IntersectionCase.EndOfEdge;
                    return true;
                }
                intersectionCase = IntersectionCase.NotIntersecting;
                return false;
            }
        }
    }
}
