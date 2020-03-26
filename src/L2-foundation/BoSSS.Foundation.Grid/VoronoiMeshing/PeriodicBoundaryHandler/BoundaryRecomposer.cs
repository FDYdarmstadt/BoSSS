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

        readonly PeriodicCornerMapper<T> cornerMapper;

        CellDetacher<T> cellDetacher;

        CornerCleaner cleaner;

        public BoundaryRecomposer(PeriodicMap map)
        {
            this.map = map;
            cornerMapper = new PeriodicCornerMapper<T>(map);
        }

        public void RecomposePeriodicEdges(Domain<T> mesh, IEnumerable<Edge<T>> periodicEdges)
        {
            CellPairCollection<T> candidates = CellPairCollecter<T>.FollowBoundaryAndCollectCandidates(periodicEdges);
            cornerMapper.FindPeriodicCorners(candidates);
            RecomposeCutCells(mesh, candidates);
            cornerMapper.ConnectPeriodicCorners();
        }

        void RecomposeCutCells(Domain<T> mesh, CellPairCollection<T> candidates)
        {
            cellDetacher = new CellDetacher<T>(mesh);
            cleaner = new CornerCleaner(map.PeriodicCornerCorrelation.Count);

            int i = 0;
            //MatlabPlotter.Plot(mesh, "Anfang");
            foreach (CellPairCollection<T>.EdgeCombo mergePair in CreateMergePairsOfEachEdge(candidates))
            {
                //MatlabPlotter.Plot(mesh, i + "aMerge");
                Debug.Assert(CellNodePositionsMatch(mergePair));
                MergeAtBoundary(mergePair);
                //MatlabPlotter.Plot(mesh, i + "bMerge");
                ++i;
            }
            //MatlabPlotter.Plot(mesh,"2aRemove");
            RemoveOuterCellsFromMesh(mesh);
            //MatlabPlotter.Plot(mesh, "2bRemove");
        }

        IEnumerable<CellPairCollection<T>.EdgeCombo> CreateMergePairsOfEachEdge(
            CellPairCollection<T> candidates)
        {
            foreach (Pair<CellPairCollection<T>.EdgeCombo> opposingEdges in candidates.GetCollectedEdgeComboPairs(map))
            {
                CellPairCollection<T>.EdgeCombo mergePair = ExtractMergePair(opposingEdges.Previous, candidates);
                CellPairCollection<T>.EdgeCombo opposedMergePair = ExtractMergePair(opposingEdges.Current, candidates);

                MoveBoundary(mergePair);
                MoveBoundary(opposedMergePair);

                Transform(mergePair);
                Transform(opposedMergePair);

                yield return mergePair;
                yield return opposedMergePair;
            }
        }

        CellPairCollection<T>.EdgeCombo ExtractMergePair(
            CellPairCollection<T>.EdgeCombo edgePair,
            CellPairCollection<T> candidates)
        {
            var mergePair = new CellPairCollection<T>.EdgeCombo(edgePair.EdgeNumber);
            int pairedBoundary = map.PeriodicBoundaryCorrelation[edgePair.EdgeNumber];
            candidates.TryGetOuterCells(pairedBoundary, out List<MeshCell<T>> pairedOuterCells);
            mergePair.Outer = pairedOuterCells;
            mergePair.Inner = new List<MeshCell<T>>(ArrayMethods.GetReverseOrderArray(edgePair.Inner));
            Debug.Assert(mergePair.Outer.Count == mergePair.Inner.Count);

            cleaner.RemoveAlreadyDealtWithCornerCellMergePairsFrom(mergePair);
            InitializeGlueMapOf(mergePair);
            return mergePair;
        }

        void ExtractEdgeGlueMap(CellPairCollection<T> pairsForMerging)
        {
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairsForMerging.GetCollectedEdgeCombos())
            {
                IList<(int outerEdge, int innerEdge)> glueMap = ExtractEdgeGlueMap(cellsOfABoundary);
                pairsForMerging.AddGlueMap(glueMap, cellsOfABoundary.EdgeNumber);
            }
        }

        void InitializeGlueMapOf(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            (int outerEdge, int innerEdge)[] glueMap= ExtractEdgeGlueMap(cellsOfABoundary);
            cellsOfABoundary.GlueMap = new List<(int outerEdge, int innerEdge)>(glueMap);
        }

        (int outerEdge, int innerEdge)[] ExtractEdgeGlueMap(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            int outerBoundaryNumber = map.PeriodicBoundaryCorrelation[cellsOfABoundary.EdgeNumber];
            int innerBoundaryNumber = cellsOfABoundary.EdgeNumber;
            (int outerEdge, int innerEdge)[] glueMap = new (int outerEdge, int innerEdge)[cellsOfABoundary.Outer.Count];
            for(int i = 0; i < cellsOfABoundary.Outer.Count; ++i)
            {
                int outerIndice = default(int);
                Edge<T>[] outer = cellsOfABoundary.Outer[i].Edges;
                for (int j = 0; j < outer.Length; ++j)
                {
                    if(outer[j].BoundaryEdgeNumber == outerBoundaryNumber)
                    {
                        outerIndice = j;
                        break;
                    }
                }

                int innerIndice = default(int);
                Edge<T>[] inner = cellsOfABoundary.Inner[i].Edges;
                for (int j = 0; j < inner.Length; ++j)
                {
                    if (inner[j].BoundaryEdgeNumber == innerBoundaryNumber)
                    {
                        innerIndice = j;
                        break;
                    }
                }
                
                glueMap[i] = (outerIndice, innerIndice);
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
            List<MeshCell<T>> cells,
            int boundaryNumber)
        {
            map.PeriodicBoundaryTransformations.TryGetValue(boundaryNumber, out Transformation trafo);
            TransformCells(cells, trafo);
        }

        static void TransformCells(IList<MeshCell<T>> cells, Transformation transformation)
        {
            HashSet<int> transformed = new HashSet<int>();
            foreach (MeshCell<T> cell in cells)
            {
                TransformCell(cell, transformation, transformed);
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

        static bool CellNodePositionsMatch(CellPairCollection<T> pairs)
        {
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairs.GetCollectedEdgeCombos())
            {
                return CellNodePositionsMatch(cellsOfABoundary);
            }
            return true;
        }

        static bool CellNodePositionsMatch(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            if (!CellNodePositionsMatch(cellsOfABoundary.Inner, cellsOfABoundary.Outer))
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        static bool CellNodePositionsMatch(IList<MeshCell<T>> cellsA, IList<MeshCell<T>> cellsB)
        {
            if (cellsA.Count != cellsB.Count)
            {
                throw new Exception("Array dimension mismatch.");
            }
            for (int i = 0; i < cellsA.Count; ++i)
            {
                MeshCell<T> cellA = cellsA[i];
                MeshCell<T> cellB = cellsB[i];
                if ((cellA.Node.Position - cellB.Node.Position) * (cellA.Node.Position - cellB.Node.Position) > 1e-6)
                {
                    return false;
                }
            }
            return true;
        }

        static void MergeAtBoundary(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
        {
            // 4) Merge cloned and transformed outer cells to corresponding innner cells
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
