using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class BoundaryRecomposer<T>
        where T : ILocatable
    {
        readonly PeriodicMap map;

        CellDetacher<T> cellDetacher;

        CornerCleaner cleaner;

        public BoundaryRecomposer(PeriodicMap map)
        {
            this.map = map;
        }

        public void RecomposePeriodicEdges(Domain<T> mesh, IEnumerable<Edge<T>> periodicEdges)
        {
            CellPairCollecter<T> collecter = new CellPairCollecter<T>(map);
            CellPairCollection<T> candidates = collecter.FollowBoundaryAndCollectCandidates(periodicEdges);
            RecomposeCutCells(mesh, candidates);
        }

        void RecomposeCutCells(Domain<T> mesh, CellPairCollection<T> candidates)
        {
            cellDetacher = new CellDetacher<T>(mesh, map);
            cleaner = new CornerCleaner(map.PeriodicCornerCorrelation.Count);

            int i = 0;
            //MatlabPlotter.Plot(mesh, "Anfang");
            foreach (CellPairCollection<T>.EdgeCombo mergePair in MergePairsOfEachEdge(candidates))
            {
                //MatlabPlotter.Plot(mesh, i + "aMerge");
                Debug.Assert(CellNodePositionsMatch(mergePair));
                MergeAtBoundary(mergePair);
                //MatlabPlotter.Plot(mesh, i + "Merge");
                ++i;
            }

            //MatlabPlotter.Plot(mesh,"2aRemove");
            RemoveOuterCellsFromMesh(mesh);
            //MatlabPlotter.Plot(mesh, "final");
        }

        IEnumerable<CellPairCollection<T>.EdgeCombo> MergePairsOfEachEdge(
            CellPairCollection<T> candidates)
        {
            foreach (Pair<CellPairCollection<T>.EdgeCombo> opposingEdges in candidates.GetCollectedEdgeComboPairs(map))
            {
                CellPairCollection<T>.EdgeCombo mergePair = ExtractMergePairs(opposingEdges.Previous, candidates);
                CellPairCollection<T>.EdgeCombo opposedMergePair = ExtractMergePairs(opposingEdges.Current, candidates);

                MoveBoundary(mergePair);
                MoveBoundary(opposedMergePair);

                Transform(mergePair);
                Transform(opposedMergePair);

                yield return mergePair;
                yield return opposedMergePair;
            }
        }

        CellPairCollection<T>.EdgeCombo ExtractMergePairs(
            CellPairCollection<T>.EdgeCombo edgePair,
            CellPairCollection<T> candidates)
        {
            var mergePair = new CellPairCollection<T>.EdgeCombo(edgePair.EdgeNumber);
            int pairedBoundary = map.PeriodicBoundaryCorrelation[edgePair.EdgeNumber];
            candidates.TryGetOuterCells(pairedBoundary, out List<(MeshCell<T>, bool)> pairedOuterCells);
            mergePair.Outer = pairedOuterCells;
            mergePair.Inner = new List<(MeshCell<T>, bool)>(ArrayMethods.GetReverseOrderArray(edgePair.Inner));
            Debug.Assert(mergePair.Outer.Count == mergePair.Inner.Count);
            cleaner.RemoveAlreadyDealtWithCornerCellMergePairsFrom(mergePair);
            InitializeGlueMapOf(mergePair);

            return mergePair;
        }

        void InitializeGlueMapOf(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            (int outerEdge, int innerEdge, bool glue)[] glueMap= ExtractEdgeGlueMap(cellsOfABoundary);
            cellsOfABoundary.GlueMap = new List<(int, int, bool)>(glueMap);
        }

        (int outerEdge, int innerEdge, bool glue)[] ExtractEdgeGlueMap(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            int outerBoundaryNumber = map.PeriodicBoundaryCorrelation[cellsOfABoundary.EdgeNumber];
            int innerBoundaryNumber = cellsOfABoundary.EdgeNumber;
            (int outerEdge, int innerEdge, bool glue)[] glueMap = new (int, int, bool)[cellsOfABoundary.Outer.Count];
            for(int i = 0; i < cellsOfABoundary.Outer.Count; ++i)
            {
                int outerIndice = default(int);
                Edge<T>[] outer = cellsOfABoundary.Outer[i].cell.Edges;
                for (int j = 0; j < outer.Length; ++j)
                {
                    if(outer[j].BoundaryEdgeNumber == outerBoundaryNumber)
                    {
                        outerIndice = j;
                        break;
                    }
                }

                int innerIndice = default(int);
                Edge<T>[] inner = cellsOfABoundary.Inner[i].cell.Edges;
                for (int j = 0; j < inner.Length; ++j)
                {
                    if (inner[j].BoundaryEdgeNumber == innerBoundaryNumber)
                    {
                        innerIndice = j;
                        break;
                    }
                }

                Debug.Assert(cellsOfABoundary.Outer[i].isSplit == cellsOfABoundary.Inner[i].isSplit);
                bool glue = cellsOfABoundary.Outer[i].isSplit;
                glueMap[i] = (outerIndice, innerIndice, glue);
            }
            return glueMap;
        }

        void MoveBoundary(CellPairCollection<T>.EdgeCombo mergePair)
        {
            int pairedBoundary = map.PeriodicBoundaryCorrelation[mergePair.EdgeNumber];
            cellDetacher.DetachCells(mergePair.Outer, pairedBoundary, mergePair.EdgeNumber);
        }

        void Transform(CellPairCollection<T>.EdgeCombo mergePair)
        {
            int pairedBoundary = map.PeriodicBoundaryCorrelation[mergePair.EdgeNumber];
            Transform(mergePair.Outer, pairedBoundary);
        }

        void Transform(
            List<(MeshCell<T>, bool)> cells,
            int boundaryNumber)
        {
            map.PeriodicBoundaryTransformations.TryGetValue(boundaryNumber, out Transformation trafo);
            TransformCells(cells, trafo);
        }

        static void TransformCells(IList<(MeshCell<T>, bool)> cells, Transformation transformation)
        {
            HashSet<int> transformed = new HashSet<int>();
            foreach ((MeshCell<T> cell, bool isSplit) item in cells)
            {
                if (item.isSplit)
                {
                    TransformCell(item.cell, transformation, transformed);
                }
            }
        }

        static void TransformCell(
            MeshCell<T> cell, 
            Transformation transformation, 
            HashSet<int> transformed)
        {
            for (int i = 0; i < cell.Vertices.Length; ++i)
            {
                Vertex vertex = cell.Vertices[i];
                if (!transformed.Contains(vertex.ID))
                {
                    transformed.Add(vertex.ID);
                    cell.Vertices[i].Position = transformation.Transform(cell.Vertices[i].Position);
                }
            }
            cell.Node.Position = transformation.Transform(cell.Node.Position);
        }

        static bool CellNodePositionsMatch(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            IList<(MeshCell<T> cell, bool isSplit)> cellsA = cellsOfABoundary.Inner;
            IList<(MeshCell<T> cell, bool isSplit)> cellsB = cellsOfABoundary.Outer;
            if (cellsA.Count != cellsB.Count)
            {
                throw new Exception("Array dimension mismatch.");
            }
            for (int i = 0; i < cellsA.Count; ++i)
            {
                if (cellsA[i].isSplit)
                {
                    if (!CellNodePositionMatches(cellsA[i].cell, cellsB[i].cell))
                    {
                        return false;
                    }
                }
                else if (cellsA[i].isSplit != cellsB[i].isSplit)
                {
                    return false;
                }
            }
            return true;
        }

        static bool CellNodePositionMatches(MeshCell<T> cellA, MeshCell<T> cellB)
        {
            if ((cellA.Node.Position - cellB.Node.Position) * (cellA.Node.Position - cellB.Node.Position) > 1e-6)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        static void MergeAtBoundary(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            BoundaryCellMerger<T>.MergeAtBoundary(
                cellsOfABoundary.Outer,
                cellsOfABoundary.Inner,
                cellsOfABoundary.GlueMap
                );
        }

        static void RemoveOuterCellsFromMesh(Domain<T> mesh)
        {
            InsideCellEnumerator<T> insideCells = new InsideCellEnumerator<T>(mesh);
            RemoveOuterCellsFromMesh(mesh.Mesh, insideCells);
        }

        static void RemoveOuterCellsFromMesh(Mesh<T> mesh, InsideCellEnumerator<T> insideCells)
        {
            List<MeshCell<T>> cells = new List<MeshCell<T>>(mesh.Cells.Count);
            foreach (MeshCell<T> cell in insideCells.EnumerateCellsInConcentricCircles())
            {
                cells.Add(cell);
            }
            for(int i = 0; i < cells.Count; ++i)
            {
                cells[i].ID = i;
            }
            mesh.Cells = cells;
            mesh.Nodes.Clear();
            foreach (MeshCell<T> cell in mesh.Cells)
            {
                mesh.Nodes.Add(cell.Node);
            }
        }
    }
}
