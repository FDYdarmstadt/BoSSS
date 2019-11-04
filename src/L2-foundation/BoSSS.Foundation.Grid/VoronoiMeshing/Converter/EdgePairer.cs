using BoSSS.Foundation.Grid.Classic;
using System.Collections.Generic;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    class EdgePairer
    {
        readonly Dictionary<int, Edge> periodicEdges;

        class Edge
        {
            readonly public CellFaceTag[] tags;

            readonly int tagIndice;

            readonly long cellID;

            public Edge(int tagIndice, CellFaceTag[] tags, long cellID)
            {
                this.cellID = cellID;
                this.tags = tags;
                this.tagIndice = tagIndice;
            }

            public long CellID {
                get {
                    return cellID;
                }
            }

            public long NeighCell_GlobalID {
                set {
                    tags[tagIndice].NeighCell_GlobalID = value;
                }
            }
        }

        public EdgePairer()
        {
            periodicEdges = new Dictionary<int, Edge>();
        }

        public void SetNeighborCell(int faceID, int tagIndice, CellFaceTag[] tags, long cellID)
        {
            if (periodicEdges.TryGetValue(faceID, out Edge edge))
            {
                tags[tagIndice].NeighCell_GlobalID = edge.CellID;
                edge.NeighCell_GlobalID = cellID;
            }
            else
            {
                Edge newEdge = new Edge(tagIndice, tags, cellID);
                periodicEdges.Add(faceID, newEdge);
            }
        }
    }
}
