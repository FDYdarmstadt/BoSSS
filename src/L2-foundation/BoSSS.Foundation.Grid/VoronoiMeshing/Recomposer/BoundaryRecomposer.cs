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
        readonly IDictionary<int, BoundaryTransformation> periodicTrafoMap;

        readonly IDictionary<int, int> periodicBoundaryMap;

        readonly int firstCellNodeIndice;

        public BoundaryRecomposer(
            IDictionary<int, int> periodicBoundaryMap, 
            IDictionary<int, BoundaryTransformation> periodicTrafoMap,
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
            CellPairCollection<T> pairsForMerging = CopyAndTransformOuterCells(candidates, cellCopier);
            Debug.Assert(CellNodePositionsMatch(pairsForMerging));
            ExtractEdgeGlueMap(pairsForMerging);

            RemoveCollectedOuterCellsFromMesh(mesh, candidates);
            MergeAtBoundary(pairsForMerging);
        }

        void ExtractEdgeGlueMap(CellPairCollection<T> pairsForMerging)
        {
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairsForMerging.GetCollectedEdgeCombos())
            {
                IList<(int outerEdge, int innerEdge)> glueMap = ExtractEdgeGlueMap(cellsOfABoundary);
                pairsForMerging.AddGlueMap(glueMap, cellsOfABoundary.EdgeNumber);
            }
        }

        IList<(int outerEdge, int innerEdge)> ExtractEdgeGlueMap(CellPairCollection<T>.EdgeCombo cellsOfABoundary)
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

        CellPairCollection<T> CopyAndTransformOuterCells(
            CellPairCollection<T> candidates, 
            MeshCellCopier<T> cellCopier)
        {
            CellPairCollection<T> pairCollecter = new CellPairCollection<T>();
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in candidates.GetCollectedEdgeCombos())
            {
                //CloneAndTransformOuterCells
                periodicTrafoMap.TryGetValue(cellsOfABoundary.EdgeNumber, out BoundaryTransformation trafo);
                MeshCell<T>[] clones = cellCopier.Copy(cellsOfABoundary.Outer);
                TransformCells(clones, trafo);

                //Collect pairs
                pairCollecter.AddOuterCells(clones, periodicBoundaryMap[cellsOfABoundary.EdgeNumber]);
                MeshCell<T>[] reversedInner = ArrayMethods.GetReverseOrderArray(cellsOfABoundary.Inner);
                pairCollecter.AddInnerCells(reversedInner, cellsOfABoundary.EdgeNumber);

            }
            return pairCollecter;
        }

        static void TransformCells(IEnumerable<MeshCell<T>> cells, BoundaryTransformation transformation)
        {
            foreach (MeshCell<T> cell in cells)
            {
                TransformCell(cell, transformation);
            }
        }

        static void TransformCell(MeshCell<T> cell, BoundaryTransformation transformation)
        {
            for (int i = 0; i < cell.Vertices.Length; ++i)
            {
                cell.Vertices[i].Position = transformation.Transform(cell.Vertices[i].Position);
            }
            foreach (Edge<T> edge in cell.Edges)
            {
                edge.Start.Position = transformation.Transform(edge.Start.Position);
                edge.End.Position = transformation.Transform(edge.End.Position);
            }
            cell.Node.Position = transformation.Transform(cell.Node.Position);
        }

        static bool CellNodePositionsMatch(CellPairCollection<T> pairs)
        {
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairs.GetCollectedEdgeCombos())
            {
                if (!CellNodePositionsMatch(cellsOfABoundary.Inner, cellsOfABoundary.Outer))
                {
                    return false;
                }
            }
            return true;
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

        void RemoveCollectedOuterCellsFromMesh(Mesh<T> mesh, CellPairCollection<T> candidates)
        {
            // Remove collected outer cells from mesh
            BoundaryCellRemover<T> remover = new BoundaryCellRemover<T>(mesh, firstCellNodeIndice);
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in candidates.GetCollectedEdgeCombos())
            {
                remover.SetAsBoundary(cellsOfABoundary.Outer, cellsOfABoundary.EdgeNumber);
            }
            remover.RemoveOuterCellsFromMesh();
        }

        static void MergeAtBoundary(CellPairCollection<T> pairsForMerging)
        {
            // 4) Merge cloned and transformed outer cells to corresponding innner cells
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairsForMerging.GetCollectedEdgeCombos())
            {

                BoundaryCellMerger<T>.MergeAtBoundary(
                    cellsOfABoundary.Outer,
                    cellsOfABoundary.Inner,
                    cellsOfABoundary.GlueMap
                    );
            }
        }
    }
}
