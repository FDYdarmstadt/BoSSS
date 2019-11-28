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
                Inner = new List<MeshCell<T>>();
                Outer = new List<MeshCell<T>>();
                GlueMap = new List<(int outerEdge, int innerEdge)>();
            }

            public List<MeshCell<T>> Inner;

            public List<MeshCell<T>> Outer;

            public List<(int outerEdge, int innerEdge)> GlueMap;
        }

        public CellPairCollection()
        {
            periodicEdges = new LinkedListDictionary<int, EdgeCombo>();
        }

        public void AddInnerCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Inner.Add(cell);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Inner.Add(cell);
            }
        }

        public void AddInnerCells(IEnumerable<MeshCell<T>> cells, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Inner.AddRange(cells);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Inner.AddRange(cells);
            }
        }

        public void AddOuterCell(MeshCell<T> cell, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Outer.Add(cell);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Outer.Add(cell);
            }
        }

        public void AddOuterCells(IEnumerable<MeshCell<T>> cells, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.Outer.AddRange(cells);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.Outer.AddRange(cells);
            }
        }

        public void AddGlueMap(IList<(int outerEdge, int innerEdge)> glueMap, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.GlueMap.AddRange(glueMap);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.GlueMap.AddRange(glueMap);
            }
        }

        public void AddGlueMapEntry((int outerEdge, int innerEdge) edge2EdgeMap, int boundaryEdgeNumber)
        {
            EdgeCombo edgeCells;
            if (periodicEdges.TryGetValue(boundaryEdgeNumber, out edgeCells))
            {
                edgeCells.GlueMap.Add(edge2EdgeMap);
            }
            else
            {
                edgeCells = new EdgeCombo(boundaryEdgeNumber);
                periodicEdges.Add(boundaryEdgeNumber, edgeCells);
                edgeCells.GlueMap.Add(edge2EdgeMap);
            }
        }

        public bool TryGetOuterCells(int boundaryEdgeNumber, out List<MeshCell<T>> cells)
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
