using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class CellPairCollection<T>
    {
        IDictionary<int, EdgeCombo> periodicEdges;

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

        public void AddInnerCells(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
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

        public void AddOuterCells(IList<MeshCell<T>> cells, int boundaryEdgeNumber)
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

        public IEnumerable<EdgeCombo> GetCollectedEdgeCombos()
        {
            EdgeCombo[] combos = new EdgeCombo[periodicEdges.Count];
            foreach((int i, KeyValuePair<int, EdgeCombo> combo) in CountingEnumerable.Wrap(periodicEdges))
            {
                combos[i] = combo.Value;
            }
            return combos;
        }
    }
}
