using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class CellPairCollection<T>
    {
        readonly IDictionary<int, EdgeCombo> periodicEdges;

        public class EdgeCombo
        {
            public int EdgeNumber;

            public EdgeCombo(int edgeNumber)
            {
                EdgeNumber = edgeNumber;
                Inner = new List<(MeshCell<T>, bool)>();
                Outer = new List<(MeshCell<T>, bool)>();
                GlueMap = new List<(int outerEdge, int innerEdge, bool glue)>();
            }

            public List<(MeshCell<T> cell, bool isSplit)> Inner;

            public List<(MeshCell<T> cell, bool isSplit)> Outer;

            public List<(int outerEdge, int innerEdge, bool glue)> GlueMap;
        }

        public CellPairCollection()
        {
            periodicEdges = new LinkedListDictionary<int, EdgeCombo>();
        }

        public void AddInnerSplitCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Inner.Add((cell, true));
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Inner.Add((cell, true));
            }
        }
        
        public void AddInnerUnsplitCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Inner.Add((cell, false));
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Inner.Add((cell, false));
            }
        }

        public void AddOuterSplitCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Outer.Add((cell, true));
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Outer.Add((cell, true));
            }
        }

        public void AddOuterUnsplitCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Outer.Add((cell, false));
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Outer.Add((cell, false));
            }
        }

        public bool TryGetOuterCells(int boundaryEdgeNumber, out List<(MeshCell<T>, bool)> cells)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                cells = edgeCells.Outer;
                return true;
            }
            else
            {
                cells = null;
                return false;
            }
        }

        public EdgeCombo[] GetCollectedEdgeCombos()
        {
            EdgeCombo[] combos = new EdgeCombo[periodicEdges.Count];
            foreach((int i, KeyValuePair<int, EdgeCombo> combo) in CountingEnumerable.Wrap(periodicEdges))
            {
                combos[i] = combo.Value;
            }
            return combos;
        }

        public bool TryGetEdgeComboOf(int boundaryNumber, out EdgeCombo edgeCombo)
        {
            edgeCombo = null;
            if (periodicEdges.TryGetValue(boundaryNumber, out edgeCombo))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public IEnumerable<Pair<EdgeCombo>> GetCollectedEdgeComboPairs(PeriodicMap map)
        {
            IDictionary<int, int> boundaryMap = map.PeriodicBoundaryCorrelation;
            EdgeCombo[] combos = EdgeComboAsPairs(boundaryMap);

            Pair<EdgeCombo>[] pairedEdgeCombos = new Pair<EdgeCombo>[combos.Length / 2];
            for ( int i = 0; i < pairedEdgeCombos.Length; ++i)
            {
                pairedEdgeCombos[i] = new Pair<EdgeCombo>
                {
                    Previous = combos[2 * i],
                    Current = combos[2 * i + 1],
                };
            }
            return pairedEdgeCombos;
        }

        EdgeCombo[] EdgeComboAsPairs( IDictionary<int, int> boundaryMap)
        {
            EdgeCombo[] combos = GetCollectedEdgeCombos();
            Debug.Assert(combos.Length % 2 == 0);
            for (int i = 0; i < combos.Length - 3; i += 2)
            {
                EdgeCombo active = combos[i];
                boundaryMap.TryGetValue(active.EdgeNumber, out int pairedBoundary);
                for (int j = i + 1; j < combos.Length; ++j)
                {
                    EdgeCombo matchCandidate = combos[j];
                    if (matchCandidate.EdgeNumber == pairedBoundary)
                    {
                        EdgeCombo temp = combos[i + 1];
                        combos[i + 1] = matchCandidate;
                        combos[j] = temp;
                    }
                }
            }
            return combos;
        }
    }
}
