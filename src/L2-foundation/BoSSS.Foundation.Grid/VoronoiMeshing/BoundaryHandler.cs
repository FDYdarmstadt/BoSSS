using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Voronoi;
using ilPSP;
using System.Diagnostics;
using System.Collections.Specialized;
using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryHandler<T>
        where T : IMesherNode, new()
    {
        readonly IDictionary<int, BoundaryTransformation> periodicTrafoMap;

        readonly IDictionary<int, int> periodicBoundaryMap;

        readonly BoundaryLine[] boundary;

        readonly int firstCellNodeIndice;

        public bool ContainsPeriodicBoundaries { get; private set; }

        public BoundaryHandler(
            BoundaryLine[] boundary,
            IDictionary<int, int> periodicBoundaryMap,
            int firstCellNodeIndice)
        {
            this.boundary = boundary;
            this.firstCellNodeIndice = firstCellNodeIndice;
            if (periodicBoundaryMap != null)
            {
                this.periodicBoundaryMap = periodicBoundaryMap;
                ContainsPeriodicBoundaries = true;
                periodicTrafoMap = CreateMappingFrom(boundary, periodicBoundaryMap);
            }
            else
            {
                ContainsPeriodicBoundaries = false;
            }
        }

        static IDictionary<int, BoundaryTransformation> CreateMappingFrom(
            BoundaryLine[] boundary, 
            IDictionary<int, int> periodicBoundaryMap)
        {
            IDictionary<int, BoundaryTransformation> periodicTrafoMap = new LinkedListDictionary< int, BoundaryTransformation>();
            foreach(var boundaryPair in periodicBoundaryMap)
            {
                BoundaryLine source = boundary[boundaryPair.Key];
                BoundaryLine target = boundary[boundaryPair.Value];
                BoundaryTransformation transformation = new BoundaryTransformation(source, target);
                periodicTrafoMap.Add(boundaryPair.Key, transformation);
            }
            return periodicTrafoMap;
        }

        public void CloneNodesAlongPeriodicBoundaries(Mesh<T> mesh, IList<T> nodes) 
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);

            List<T> clones = new List<T>();
            Edge<T> preceedingEdge = null;
            //Follow boundary and
            foreach (Edge<T> edge in PeriodicEdges(mesh))
            {
                // 1) Collect nodes that should be mirrored
                // 2) Clone and teleport/transform to paired boundary
                periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out BoundaryTransformation transformation);
                MeshCell<T> cell = edge.Cell;
                T clone = new T
                {
                    Position = transformation.Transform(cell.Node.Position)
                };
                clones.Add(clone);

                if (IsCorner(edge, preceedingEdge))
                {
                    periodicTrafoMap.TryGetValue(preceedingEdge.BoundaryEdgeNumber, out BoundaryTransformation preceedingTransformation);
                    T cornerClone = new T
                    {
                        Position = preceedingTransformation.Transform(transformation.Transform(cell.Node.Position))
                    };
                    clones.Add(cornerClone);
                }
                preceedingEdge = edge;
            }
            nodes.AddRange(clones);
        }

        IEnumerable<Edge<T>> PeriodicEdges(Mesh<T> mesh)
        {
            BoundaryEdgeEnumerator<T> edgeCells = new BoundaryEdgeEnumerator<T>(mesh);
            Vector startFromEnclosingCell = boundary[0].Start.Position;
            foreach (Edge<T> edge in edgeCells.Edges(startFromEnclosingCell, firstCellNodeIndice))
            {
                if (periodicTrafoMap.ContainsKey(edge.BoundaryEdgeNumber))
                {
                    yield return edge;
                }
            }
        }

        static bool IsCorner(Edge<T> edge, Edge<T> preceedingEdge)
        {
            bool areDifferentBoundaries = (edge.BoundaryEdgeNumber != (preceedingEdge?.BoundaryEdgeNumber ?? int.MinValue));
            bool ofSameCell = edge.Cell.ID == preceedingEdge.Cell.ID;
            return areDifferentBoundaries & ofSameCell;
        }

        bool NodeOfEdgeIsOutsideOfBoundary(Edge<T> edge)
        {
            Vector position = edge.Cell.Node.Position;
            Line boundaryEdge = boundary[edge.BoundaryEdgeNumber];
            return IsLeft(position, boundaryEdge);
        }

        static bool IsLeft(Vector node, Line line)
        {
            Vector start = line.Start.Position;
            Vector end = line.End.Position;
            return ((end.x - start.x) * (node.y - start.y) - (end.y - start.y) * (node.x - start.x)) > 0;
        }

        public void RecomposePeriodicEdges(Mesh<T> mesh)
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);

            CellPairCollection<T> candidates = new CellPairCollection<T>();
            //Merge/Shift: Follow boundary and
            // 1) Collect
            foreach (Edge<T> edge in PeriodicEdges(mesh))
            {
                if (NodeOfEdgeIsOutsideOfBoundary(edge))
                {
                    candidates.AddOuterCell(edge.Cell, edge.BoundaryEdgeNumber);
                }
                else
                {
                    candidates.AddInnerCell(edge.Cell, edge.BoundaryEdgeNumber);
                }
            }
            // 2) CloneAndTransformOuterCells
            CellPairCollection<T> pairsForMerging = CloneAndTransformOuterCells(candidates);
            Debug.Assert(CellNodePositionsMatch( pairsForMerging));

            BoundaryCellRemover<T> remover = new BoundaryCellRemover<T>(mesh, firstCellNodeIndice);
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in candidates.GetCollectedEdgeCombos())
            {
                // 3) Remove collected outer cells from mesh
                remover.Remove(cellsOfABoundary.Outer);
            }
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in pairsForMerging.GetCollectedEdgeCombos())
            {
                // 4) Merge collected outer cells to corresponding innner cells
                BoundaryCellMerger<T>.MergeAtBoundary(
                    cellsOfABoundary.Inner, 
                    cellsOfABoundary.EdgeNumber, 
                    cellsOfABoundary.Outer, 
                    periodicBoundaryMap[cellsOfABoundary.EdgeNumber]);
            }
        }

        CellPairCollection<T> CloneAndTransformOuterCells(CellPairCollection<T> candidateCollecter)
        {
            CellPairCollection<T> pairCollecter = new CellPairCollection<T>();
            foreach (CellPairCollection<T>.EdgeCombo cellsOfABoundary in candidateCollecter.GetCollectedEdgeCombos())
            {
                //CloneAndTransformOuterCells
                periodicTrafoMap.TryGetValue(cellsOfABoundary.EdgeNumber, out BoundaryTransformation trafo);
                MeshCell<T>[] clones = MeshCloner.Clone(cellsOfABoundary.Outer);
                TransformCells(clones, trafo);
                
                //Collect pairs
                pairCollecter.AddOuterCells(clones, periodicBoundaryMap[cellsOfABoundary.EdgeNumber]);
                pairCollecter.AddInnerCells(cellsOfABoundary.Inner, cellsOfABoundary.EdgeNumber);
                
            }
            return pairCollecter;
        }

        bool CellNodePositionsMatch(CellPairCollection<T> pairs)
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

        bool CellNodePositionsMatch(IList<MeshCell<T>> cellsA, IList<MeshCell<T>> cellsB)
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

        static void TransformCells(IEnumerable<MeshCell<T>> cells, BoundaryTransformation transformation)
        {
            foreach(MeshCell<T> cell in cells)
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
            foreach(Edge<T> edge in cell.Edges)
            {
                edge.Start.Position = transformation.Transform(edge.Start.Position);
                edge.End.Position = transformation.Transform(edge.End.Position);
            }
            cell.Node.Position = transformation.Transform(cell.Node.Position);
        }
    }
}
