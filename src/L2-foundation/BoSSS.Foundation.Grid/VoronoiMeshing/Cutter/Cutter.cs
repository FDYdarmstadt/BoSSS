using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Cutter
{
    class Cutter<T>
        where T :ILocatable, new()
    {
        BoundaryLineEnumerator boundaryLines;

        MeshIntersecter<T> meshIntersecter;

        Divider<T> greatDivider;

        CutterState<Edge<T>> state;

        FirstCutInfo firstCell;

        OnEdgeCutter<T> edgeCutter;

        bool isFirstIntersection;

        class FirstCutInfo
        {
            public List<BoundaryLine> linesFirstCell = null;

            public Edge<T> firstCellCutRidge = default(Edge<T>);

            public IntersectionCase CutCase = default(IntersectionCase);
        }

        public void CutOut(IDMesh<T> mesh, Boundary<T> boundary)
        {
            this.boundaryLines = BoundaryLineEnumerator.GetEnumerator(BoundaryLine.Copy(boundary.BoundaryLines));
            Initialize(mesh, boundary);
            CutOut();
        }

        void Initialize(
            IDMesh<T> mesh,
            Boundary<T> boundary)
        {
            this.meshIntersecter = new MeshIntersecter<T>(mesh);
            this.greatDivider = new Divider<T>(mesh, boundary);
            
            state = new CutterState<Edge<T>>();
            edgeCutter = new OnEdgeCutter<T>(this.meshIntersecter, this.boundaryLines);
            firstCell = null;
        }

        void CutOut()
        {
            IEnumerator<Edge<T>> edgeEnumerator = FirstCellEdgeEnumerator();
            List<BoundaryLine> lines = new List<BoundaryLine>(10);
            BoundaryLine activeLine = default(BoundaryLine);

            bool firstCut = true;
            isFirstIntersection = true;

            while (boundaryLines.MoveNext())
            {
                bool circleCells = true;
                while (circleCells)
                {
                    circleCells = false;
                    activeLine = boundaryLines.Current;
                    ResetState();

                    //Check for intersection of Cell and Line
                    while (edgeEnumerator.MoveNext() && !circleCells)
                    {
                        circleCells = FindIntersection(edgeEnumerator.Current, activeLine);
                    }
                    isFirstIntersection = false;

                    //Handle intersection
                    if (firstCut)
                    {
                        edgeEnumerator = TryFirstCut(ref circleCells, edgeEnumerator);
                        //Start with empty line set
                        if (state.Case != IntersectionCase.NotIntersecting)
                        {
                            firstCut = false;
                            FinishFirstCut(lines);
                        }
                    }
                    else
                    {
                        edgeEnumerator = TryCut(ref circleCells, edgeEnumerator, lines);
                    }
                    if (state.Case != IntersectionCase.NotIntersecting)
                    {
                        lines.Clear();
                    }
                }
                lines.Add(activeLine);
            }
            HandleLastCell();
            boundaryLines.Reset();
            greatDivider.RemoveOutsideCells();
        }

        IEnumerator<Edge<T>> FirstCellEdgeEnumerator()
        {
            IEnumerator<Edge<T>> ridgeEnum;
            boundaryLines.MoveNext();
            MeshCell<T> first = greatDivider.GetFirst(boundaryLines.Current);
            ridgeEnum = meshIntersecter.GetFirstEnumerator(first);
            boundaryLines.Reset();
            return ridgeEnum;
        }

        void ResetState()
        {
            state.Case = IntersectionCase.NotIntersecting;
            state.ActiveEdge = default(Edge<T>);
            state.AlphaCut = default(double);
        }

        IEnumerator<Edge<T>> TryFirstCut(ref bool circleCells, IEnumerator<Edge<T>> ridgeEnum)
        {
            IEnumerator<Edge<T>> runningEnum;
            Edge<T> activeRidge = state.ActiveEdge;

            //First cut ever
            //-----------------------------------------------------------
            switch (state.Case)
            {
                case IntersectionCase.NotIntersecting:
                    ridgeEnum.Reset();
                    break;
                case IntersectionCase.InMiddle:
                    //if intersection was successfull, select next cell
                    activeRidge = meshIntersecter.FirstCut(activeRidge, state.AlphaCut);
                    ridgeEnum = meshIntersecter.GetNeighborFromEdgeNeighbor(activeRidge);
                    break;
                case IntersectionCase.StartOfLine:
                    activeRidge = meshIntersecter.FirstCut(activeRidge, state.AlphaCut);
                    runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                    ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    break;
                case IntersectionCase.EndOfLine:
                    activeRidge = meshIntersecter.FirstCut(activeRidge, state.AlphaCut);
                    if (boundaryLines.MoveNext())
                    {
                        runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                        ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    }
                    else
                    {
                        circleCells = false;
                    }
                    break;
                case IntersectionCase.EndOfEdge:
                    meshIntersecter.VertexCut(activeRidge, state.AlphaCut);
                    runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                    ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    break;
                case IntersectionCase.EndOfEdgeAndLine:
                    meshIntersecter.VertexCut(activeRidge, state.AlphaCut);
                    if (boundaryLines.MoveNext())
                    {
                        runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                        ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    }
                    else
                    {
                        circleCells = false;
                    }
                    break;
                default:
                    throw new Exception();
            }
            return ridgeEnum;
        }

        void FinishFirstCut(List<BoundaryLine> lines)
        {
            //Information needed for the last cell, when boundary is closed
            firstCell = new FirstCutInfo
            {
                firstCellCutRidge = state.ActiveEdge,
                linesFirstCell = new List<BoundaryLine>(lines),
                CutCase = state.Case
            };
            if (firstCell.linesFirstCell.Count == 0)
            {
                firstCell.linesFirstCell.Add(boundaryLines.Current);
            }
        }

        IEnumerator<Edge<T>> TryCut(
            ref bool circleCells, 
            IEnumerator<Edge<T>> ridgeEnum, 
            List<BoundaryLine> lines)
        {
            IEnumerator<Edge<T>> runningEnum;
            Edge<T> activeRidge = state.ActiveEdge;
            //All other cuts and subdivisions etc.
            //-----------------------------------------------------------
            switch (state.Case)
            {
                case IntersectionCase.NotIntersecting:
                    ridgeEnum.Reset();
                    break;
                case IntersectionCase.InMiddle:
                    //if intersection was successfull, select next cell
                    activeRidge = meshIntersecter.Subdivide(activeRidge, lines, state.AlphaCut, boundaryLines.LineIndex());
                    ridgeEnum = meshIntersecter.GetNeighborFromEdgeNeighbor(activeRidge);
                    break;
                case IntersectionCase.EndOfLine:
                    activeRidge = meshIntersecter.Subdivide(activeRidge, lines, state.AlphaCut, boundaryLines.LineIndex());
                    if (boundaryLines.MoveNext())
                    {
                        runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                        ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    }
                    else
                    {
                        circleCells = false;
                    }
                    break;
                case IntersectionCase.EndOfEdge:
                    activeRidge = meshIntersecter.SubdivideWithoutNewVertex(activeRidge, lines, boundaryLines.LineIndex());
                    runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                    ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    break;
                case IntersectionCase.EndOfEdgeAndLine:
                    activeRidge = meshIntersecter.SubdivideWithoutNewVertex(activeRidge, lines, boundaryLines.LineIndex());
                    if (boundaryLines.MoveNext())
                    {
                        runningEnum = meshIntersecter.GetConnectedEdgeEnum(activeRidge);
                        ridgeEnum = edgeCutter.CutOut(runningEnum, activeRidge);
                    }
                    else
                    {
                        circleCells = false;
                    }
                    break;
                default:
                    throw new NotImplementedException();
            }
            return ridgeEnum;
            
        }

        bool FindIntersection(Edge<T> edge, BoundaryLine line)
        {
            state.ActiveEdge = edge;
            bool foundIntersection;
            if (isFirstIntersection)
            {
                foundIntersection = LineIntersect.FindFirst(edge, line, ref state.Case, out state.AlphaCut);
            }
            else
            {
                foundIntersection = LineIntersect.Find(edge, line, ref state.Case, out state.AlphaCut);
            }
            return foundIntersection;
        }

        void HandleLastCell()
        {
            //Handle last cell
            //-----------------------------------------------------------
            switch (state.Case)
            {
                case IntersectionCase.NotIntersecting:
                    meshIntersecter.CloseMesh(firstCell.linesFirstCell, firstCell.firstCellCutRidge, boundaryLines.LineIndex());
                    break;
                case IntersectionCase.EndOfLine:
                case IntersectionCase.EndOfEdgeAndLine:
                    if(firstCell.CutCase != IntersectionCase.StartOfLine)
                    {
                        CyclicInterval lineIndex = boundaryLines.LineIndex();
                        lineIndex.Previous();
                        meshIntersecter.CloseMesh(firstCell.firstCellCutRidge, lineIndex);
                    }
                    break;
                case IntersectionCase.EndOfEdge:
                case IntersectionCase.InMiddle:
                default:
                    throw new Exception();
            }
        }
    }
}
