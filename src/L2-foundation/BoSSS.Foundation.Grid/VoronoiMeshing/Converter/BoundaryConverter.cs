using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    class BoundaryConverter
    {
        VoronoiBoundary boundary;

        readonly PeriodicBoundaryConverter periodicBoundary;

        public BoundaryConverter(VoronoiBoundary boundary, PeriodicMap periodicMap = null)
        {
            this.periodicBoundary = new PeriodicBoundaryConverter(boundary, periodicMap);
            this.boundary = boundary;
        }

        public void RegisterBoundaries(Cell cell, List<BoundaryFace> tags)
        {
            CreateCellFaceTags(cell, tags);
            AddPeriodicNeighbor(cell, tags);
        }

        void CreateCellFaceTags(Cell cell, List<BoundaryFace> tags)
        {
            foreach (BoundaryFace tag in tags)
            {
                byte edgeTag = boundary.GetEdgeTagOfPolygonEdge(tag.BoundaryEdgeNumber);
                int faceIndice = tag.Face;
                CellFaceTag CFT = new CellFaceTag()
                {
                    EdgeTag = edgeTag,
                    FaceIndex = faceIndice,
                    NeighCell_GlobalID = long.MinValue
                };
                CFT.AddToArray(ref cell.CellFaceTags);
            }
        }

        void AddPeriodicNeighbor(Cell cell, List<BoundaryFace> tags)
        {
            for (int i = 0; i < tags.Count; ++i)
            {
                CellFaceTag CFT = cell.CellFaceTags[i];
                BoundaryFace face2EdgeMap = tags[i];
                if (CFT.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG)
                {
                    Debug.Assert(periodicBoundary != null);
                    periodicBoundary.SetPeriodicData(
                        i,
                        cell.CellFaceTags,
                        face2EdgeMap.BoundaryEdgeNumber,
                        cell.GlobalID,
                        face2EdgeMap.ID,
                        face2EdgeMap.NeighborID);
                }
            }
        }

        public void RegisterEdgesTo(GridCommons grid)
        {
            if (boundary.EdgeTagNames != null)
            {
                foreach (KeyValuePair<byte, string> tagName in boundary.EdgeTagNames)
                {
                    grid.EdgeTagNames.Add(tagName);
                }
            }
            if(periodicBoundary != null)
            {
                foreach (var pair in periodicBoundary.GetPeriodicTrafos())
                {
                    AffineTrafo trafo = pair.Value;
                    byte edgeTag = pair.Key;
                    byte edgeTagInGrid = grid.AddPeriodicEdgeTrafo(trafo);
                    Debug.Assert(edgeTag == edgeTagInGrid);
                }
            }
        }
    }
}
