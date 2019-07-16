using System;
using System.Collections.Generic;
using System.Linq;
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

    interface IIntersectableMesh<TCell, TRidge, TLine>
    {
        (TCell, IEnumerator<TRidge>) GetFirst(TLine boundaryLine);

        bool intersect(TRidge ridge, TLine line, ref IntersectionCase intersectionCase, out double alphaCut);

        (TCell, IEnumerator<TRidge>) getNeighborFromEdgeNeighbor(TRidge ridge);

        TRidge Subdivide(TRidge ridge, List<TLine> lines, double alphaCut);

        TRidge SubdivideWithoutNewVertex(TRidge ridge, List<TLine> lines);

        TRidge FirstCut(TRidge ridge, double alphaCut);

        void VertexCut(TRidge ridge, double alphaCut);

        void CloseMesh(List<TLine> lines, TRidge ridgeOfFirstCut);

        IEnumerator<TRidge> getConnectedRidgeEnum(TRidge ridge);

        IEnumerator<TRidge> getNeighborFromLineDirection(TRidge ridge, TLine line);

        TRidge AddLineSegment(TRidge ridge, double alpha);

        void AddEdge(TRidge ridge);

        bool CutIsFresh{ set; } 
    }

    static class Intersecter
    {
        public static void Intersect<TCell, TRidge, TLine>(IIntersectableMesh<TCell, TRidge, TLine> vMesh,
            IEnumerator<TLine> boundary)
        {
            //Setup
            TCell activeCell;
            IEnumerator<TRidge> ridgeEnum;
            IEnumerator<TRidge> runningEnum;
            boundary.MoveNext();
            (activeCell, ridgeEnum) = vMesh.GetFirst(boundary.Current);
            boundary.Reset();
            List<TLine> lines = new List<TLine>(10);

            List<TLine> linesFirstCell = null;
            TRidge firstCellCutRidge = default(TRidge);
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
                    TRidge activeRidge = default(TRidge);
                    double alphaCut = default(double);
                    intersectionCase = IntersectionCase.NotIntersecting;

                    //Check for intersection of Cell and Line
                    while (ridgeEnum.MoveNext() && !found)
                    {
                        activeRidge = ridgeEnum.Current;
                        found = vMesh.intersect(activeRidge, activeLine, ref intersectionCase, out alphaCut);
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
                                (activeCell, ridgeEnum) = vMesh.getNeighborFromEdgeNeighbor(activeRidge);
                                break;
                            case IntersectionCase.EndOfLine:
                                activeRidge = vMesh.FirstCut(activeRidge, alphaCut);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            case IntersectionCase.EndOfRidge:
                                vMesh.VertexCut(activeRidge, alphaCut);
                                runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                                ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                break;
                            case IntersectionCase.EndOfRidgeAndLine:
                                vMesh.VertexCut(activeRidge, alphaCut);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
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
                                activeRidge = vMesh.Subdivide(activeRidge, lines, alphaCut);
                                (activeCell, ridgeEnum) = vMesh.getNeighborFromEdgeNeighbor(activeRidge);
                                break;
                            case IntersectionCase.EndOfLine:
                                activeRidge = vMesh.Subdivide(activeRidge, lines, alphaCut);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                                    ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                }
                                else
                                {
                                    found = false;
                                }
                                break;
                            case IntersectionCase.EndOfRidge:
                                activeRidge = vMesh.SubdivideWithoutNewVertex(activeRidge, lines);
                                runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                                ridgeEnum = RunningOnRidges(vMesh, boundary, runningEnum, activeRidge);
                                break;
                            case IntersectionCase.EndOfRidgeAndLine:
                                activeRidge = vMesh.SubdivideWithoutNewVertex(activeRidge, lines);
                                if (boundary.MoveNext())
                                {
                                    runningEnum = vMesh.getConnectedRidgeEnum(activeRidge);
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
                    vMesh.CloseMesh(linesFirstCell, firstCellCutRidge);
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
            IEnumerator<TLine> boundary,
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
                    keepOnRunning = vMesh.intersect(activeRidge, activeLine, ref intersectionCase, out alphaCut);
                }
                switch (intersectionCase)
                {
                    case IntersectionCase.NotIntersecting:
                        break;
                    case IntersectionCase.EndOfLine:
                        activeRidge = vMesh.AddLineSegment(activeRidge, alphaCut);
                        ridgeEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                        outerRidge = activeRidge;
                        if (!boundary.MoveNext())
                        {
                            keepOnRunning = false;
                        }
                        break;
                    case IntersectionCase.EndOfRidge:
                        vMesh.AddEdge(activeRidge);
                        ridgeEnum = vMesh.getConnectedRidgeEnum(activeRidge);
                        outerRidge = activeRidge;
                        break;
                    case IntersectionCase.EndOfRidgeAndLine:
                        vMesh.AddEdge(activeRidge);
                        ridgeEnum = vMesh.getConnectedRidgeEnum(activeRidge);
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
