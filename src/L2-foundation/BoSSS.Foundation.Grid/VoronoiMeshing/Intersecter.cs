using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Platform;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    enum IntersectionCase
    {
        NotIntersecting,
        InMiddle,
        EndOfLine,
        EndOfRidge,
        EndOfRidgeAndLine
    }

    interface IBoundaryEnumerator<T> : IEnumerator<T>
    {
        CyclicInterval LineIndex { get; }    
    }

    interface IIntersectableMesh<TCell, TEdge, TLine>
    {
        (TCell, IEnumerator<TEdge>) GetFirst(TLine boundaryLine);

        bool Intersect(TEdge ridge, TLine line, ref IntersectionCase intersectionCase, out double alphaCut);

        (TCell, IEnumerator<TEdge>) GetNeighborFromEdgeNeighbor(TEdge ridge);

        TEdge Subdivide(TEdge ridge, List<TLine> lines, double alphaCut, CyclicInterval boundaryCount);

        TEdge SubdivideWithoutNewVertex(TEdge ridge, List<TLine> lines, CyclicInterval boundaryCount);

        TEdge FirstCut(TEdge ridge, double alphaCut);

        void VertexCut(TEdge ridge, double alphaCut);

        void CloseMesh(List<TLine> lines, TEdge ridgeOfFirstCut, CyclicInterval boundaryCount);

        IEnumerator<TEdge> GetConnectedRidgeEnum(TEdge ridge);

        IEnumerator<TEdge> getNeighborFromLineDirection(TEdge ridge, TLine line);

        TEdge AddLineSegment(TEdge ridge, double alpha, CyclicInterval boundaryCount);

        void AddEdge(TEdge ridge, CyclicInterval boundaryCount);

        bool CutIsFresh{ set; } 
    }

    static class Intersecter
    {
        public static void Intersect<TCell, TEdge, TLine>(
            IIntersectableMesh<TCell, TEdge, TLine> vMesh,
            IBoundaryEnumerator<TLine> boundary)
        {
            //Setup
            TCell activeCell;
            IEnumerator<TEdge> ridgeEnum;
            IEnumerator<TEdge> runningEnum;
            boundary.MoveNext();
            (activeCell, ridgeEnum) = vMesh.GetFirst(boundary.Current);
            boundary.Reset();
            List<TLine> lines = new List<TLine>(10);

            List<TLine> linesFirstCell = null;
            TEdge firstCellCutRidge = default(TEdge);
            IntersectionCase intersectionCase = default(IntersectionCase); ;

            bool firstCut = true;
            TLine activeLine = default(TLine);
            while (boundary.MoveNext())
            {
                bool found = true;
                while (found)
                {
                    found = false;
                    activeLine = boundary.Current;
                    TEdge activeRidge = default(TEdge);
                    double alphaCut = default(double);
                    intersectionCase = IntersectionCase.NotIntersecting;

                    //Check for intersection of Cell and Line
                    while (ridgeEnum.MoveNext() && !found)
                    {
                        activeRidge = ridgeEnum.Current;
                        found = vMesh.Intersect(activeRidge, activeLine, ref intersectionCase, out alphaCut);
                    }
                    //Handle intersection
                    if (firstCut)
                    {
                        //First cut ever
                        //-----------------------------------------------------------
                        switch (intersectionCase)
                        {
                            case IntersectionCase.NotIntersecting:
                                ridgeEnum.Reset();
                                break;
                            case IntersectionCase.InMiddle:
                                //if intersection was successfull, select next cell
                                activeRidge = vMesh.FirstCut(activeRidge, alphaCut);
                                (activeCell, ridgeEnum) = vMesh.GetNeighborFromEdgeNeighbor(activeRidge);
                                break;
                            case IntersectionCase.EndOfLine:
                                activeRidge = vMesh.FirstCut(activeRidge, alphaCut);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            case IntersectionCase.EndOfRidge:
                                vMesh.VertexCut(activeRidge, alphaCut);
                                runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                break;
                            case IntersectionCase.EndOfRidgeAndLine:
                                vMesh.VertexCut(activeRidge, alphaCut);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            default:
                                throw new Exception();
                        }
                        //Start with empty line set
                        if (intersectionCase != IntersectionCase.NotIntersecting)
                        {
                            firstCut = false;
                            lines.Clear();
                            firstCellCutRidge = activeRidge;
                            //Information needed for the last cell, when boundary is closed
                            linesFirstCell = new List<TLine>(lines);
                            if (linesFirstCell.Count == 0)
                            {
                                linesFirstCell.Add(boundary.Current);
                            }
                        }
                    }
                    else
                    {
                        //All other cuts and subdivisions etc.
                        //-----------------------------------------------------------
                        switch (intersectionCase)
                        {
                            case IntersectionCase.NotIntersecting:
                                ridgeEnum.Reset();
                                break;
                            case IntersectionCase.InMiddle:
                                //if intersection was successfull, select next cell
                                activeRidge = vMesh.Subdivide(activeRidge, lines, alphaCut, boundary.LineIndex);
                                (activeCell, ridgeEnum) = vMesh.GetNeighborFromEdgeNeighbor(activeRidge);
                                break;
                            case IntersectionCase.EndOfLine:
                                activeRidge = vMesh.Subdivide(activeRidge, lines, alphaCut, boundary.LineIndex);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            case IntersectionCase.EndOfRidge:
                                activeRidge = vMesh.SubdivideWithoutNewVertex(activeRidge, lines, boundary.LineIndex);
                                runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                break;
                            case IntersectionCase.EndOfRidgeAndLine:
                                activeRidge = vMesh.SubdivideWithoutNewVertex(activeRidge, lines, boundary.LineIndex);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            default:
                                throw new NotImplementedException();
                        }
                        if (intersectionCase != IntersectionCase.NotIntersecting)
                        {
                            lines.Clear();
                        }
                    }
                }
                lines.Add(activeLine);
            }
            //Handle last cell
            //-----------------------------------------------------------
            switch (intersectionCase)
            {
                case IntersectionCase.NotIntersecting:
                    vMesh.CloseMesh(linesFirstCell, firstCellCutRidge, boundary.LineIndex);
                    break;
                case IntersectionCase.EndOfLine:
                case IntersectionCase.EndOfRidgeAndLine:
                    break;
                case IntersectionCase.EndOfRidge:
                case IntersectionCase.InMiddle:
                default:
                    throw new Exception();
            }
            boundary.Reset();
            vMesh.CutIsFresh = true;
        }

        static IEnumerator<TRidge> RunningOnRidges<TCell, TRidge, TLine>(IIntersectableMesh<TCell, TRidge, TLine> vMesh,
            IBoundaryEnumerator<TLine> boundary,
            IEnumerator<TRidge> ridgeEnum,
            TRidge outerRidge)
        {
            bool keepOnRunning = true;
            TRidge activeRidge = default(TRidge);
            TLine activeLine = default(TLine);
            while (keepOnRunning)
            {
                keepOnRunning = false;
                IntersectionCase intersectionCase = IntersectionCase.NotIntersecting;
                double alphaCut = 0;
                activeLine = boundary.Current;
                while (ridgeEnum.MoveNext() && !keepOnRunning)
                {
                    activeRidge = ridgeEnum.Current;
                    keepOnRunning = vMesh.Intersect(activeRidge, activeLine, ref intersectionCase, out alphaCut);
                }
                switch (intersectionCase)
                {
                    case IntersectionCase.NotIntersecting:
                        break;
                    case IntersectionCase.EndOfLine:
                        activeRidge = vMesh.AddLineSegment(activeRidge, alphaCut, boundary.LineIndex);
                        ridgeEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                        outerRidge = activeRidge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.EndOfRidge:
                        vMesh.AddEdge(activeRidge, boundary.LineIndex);
                        ridgeEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                        outerRidge = activeRidge;
                        break;
                    case IntersectionCase.EndOfRidgeAndLine:
                        vMesh.AddEdge(activeRidge, boundary.LineIndex);
                        ridgeEnum = vMesh.GetConnectedRidgeEnum(activeRidge);
                        outerRidge = activeRidge;
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
            ridgeEnum = vMesh.getNeighborFromLineDirection(outerRidge, activeLine);
            return ridgeEnum;
        }
    }
}
