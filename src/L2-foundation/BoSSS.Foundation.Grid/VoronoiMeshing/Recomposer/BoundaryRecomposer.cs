using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Voronoi;
using ilPSP;
using System.Diagnostics;
using System.Collections.Specialized;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryRecomposer<T>
        where T : ILocatable, new()
    {
        readonly IDictionary<int, Transformation> periodicTrafoMap;

        readonly IDictionary<int, int> periodicBoundaryMap;

        readonly int firstCellNodeIndice;

        public BoundaryRecomposer(
            IDictionary<int, int> periodicBoundaryMap, 
            IDictionary<int, Transformation> periodicTrafoMap,
            int firstCellNodeIndice)
        {
            this.periodicBoundaryMap = periodicBoundaryMap;
            this.periodicTrafoMap = periodicTrafoMap;
            this.firstCellNodeIndice = firstCellNodeIndice;
        }

        public void RecomposePeriodicEdges(IDMesh<T> mesh, IEnumerable<Edge<T>> periodicEdges)
        {
            CellPairCollection<T> candidates = CellPairCollecter<T>.FollowBoundaryAndCollectCandidates(periodicEdges);
            MeshCellCopier<T> cellCopier = new MeshCellCopier<T>(mesh);
            BoundaryCellRemover<T> remover = new BoundaryCellRemover<T>(mesh, firstCellNodeIndice);
            //MatlabPlotter plotter = new MatlabPlotter();
            //int i = 0;

            foreach (CellPairCollection<T>.EdgeCombo mergePair in CreateMergePairsOfEachEdge(candidates, cellCopier, remover))
            {
                Debug.Assert(CellNodePositionsMatch(mergePair));
                InitializeGlueMapOf(mergePair);
                MergeAtBoundary(mergePair);
                //plotter.Plot(mesh, "intermediate" + i);
                //++i;
            }
            //plotter.Plot(mesh, "final");
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
            int outerBoundaryNumber = periodicBoundaryMap[cellsOfABoundary.EdgeNumber];
            int innerBoundaryNumber = cellsOfABoundary.EdgeNumber;
            (int outerEdge, int innerEdge)[] map = new (int outerEdge, int innerEdge)[cellsOfABoundary.Outer.Count];
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
                
                map[i] = (outerIndice, innerIndice);
            }
            return map;
        }

        IEnumerable<CellPairCollection<T>.EdgeCombo> CreateMergePairsOfEachEdge(
            CellPairCollection<T> candidates,
            MeshCellCopier<T> cellCopier,
            BoundaryCellRemover<T> remover)
        {
            foreach (CellPairCollection<T>.EdgeCombo edgePair in candidates.GetCollectedEdgeCombos())
            {
                CellPairCollection<T>.EdgeCombo mergePair = new CellPairCollection<T>.EdgeCombo(edgePair.EdgeNumber);
                int pairedBoundary = periodicBoundaryMap[edgePair.EdgeNumber];
                candidates.TryGetOuterCells(pairedBoundary, out List<MeshCell<T>> pairedOuterCells);
                mergePair.Outer = TransformedCopyOfOuter(pairedOuterCells, pairedBoundary, cellCopier);
                mergePair.Inner = new List<MeshCell<T>>(ArrayMethods.GetReverseOrderArray(edgePair.Inner));
                remover.EnqueueForRemoval(pairedOuterCells, pairedBoundary);

                yield return mergePair;
            }
            remover.RemoveQueuedCells();
        }

        List<MeshCell<T>> TransformedCopyOfOuter(
            List<MeshCell<T>> candidates,
            int boundaryNumber,
            MeshCellCopier<T> cellCopier)
        {
            //CloneAndTransformOuterCells
            periodicTrafoMap.TryGetValue(boundaryNumber, out Transformation trafo);
            List<MeshCell<T>> clones = cellCopier.Copy(candidates);
            TransformCells(clones, trafo);
            return clones;
        }

        static void TransformCells(IEnumerable<MeshCell<T>> cells, Transformation transformation)
        {
            foreach (MeshCell<T> cell in cells)
            {
                TransformCell(cell, transformation);
            }
        }

        static void TransformCell(MeshCell<T> cell, Transformation transformation)
        {
            for (int i = 0; i < cell.Vertices.Length; ++i)
            {
                cell.Vertices[i].Position = transformation.Transform(cell.Vertices[i].Position);
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
    }
}
