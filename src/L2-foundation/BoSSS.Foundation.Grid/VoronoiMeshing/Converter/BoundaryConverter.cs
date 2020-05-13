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
        readonly VoronoiBoundary boundary;

        readonly PeriodicBoundaryConverter periodicBoundaryConverter;

        public BoundaryConverter(VoronoiBoundary boundary, PeriodicMap periodicMap = null)
        {
            if(periodicMap != null)
            {
                this.periodicBoundaryConverter = new PeriodicBoundaryConverter(boundary, periodicMap);
            }
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
                byte edgeTag = GetEdgeTagOf(tag.BoundaryEdgeNumber);
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

        byte GetEdgeTagOf(int boundaryEdgeNumber)
        {
            byte edgeTag = default(byte);
            if(periodicBoundaryConverter != null)
            {
                edgeTag = periodicBoundaryConverter.GetEdgeTagOf(boundaryEdgeNumber);
            }
            else
            {
                edgeTag = boundary.GetEdgeTagOfPolygonEdge(boundaryEdgeNumber);
            }
            return edgeTag;
        }

        void AddPeriodicNeighbor(Cell cell, List<BoundaryFace> tags)
        {
            for (int i = 0; i < tags.Count; ++i)
            {
                CellFaceTag CFT = cell.CellFaceTags[i];
                BoundaryFace face2EdgeMap = tags[i];
                if (CFT.EdgeTag >= GridCommons.FIRST_PERIODIC_BC_TAG)
                {
                    Debug.Assert(periodicBoundaryConverter != null);
                    periodicBoundaryConverter.SetPeriodicData(
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
                RegisterTagNamesTo(grid);
            }
            if(periodicBoundaryConverter != null)
            {
                periodicBoundaryConverter.RegisterPeriodicBoundariesTo(grid);
            }
        }

        void RegisterTagNamesTo(GridCommons grid)
        {
            foreach (KeyValuePair<byte, string> tagName in boundary.EdgeTagNames)
            {
                grid.EdgeTagNames.Add(tagName);
            }
        }

        public void Clear()
        {
            if(periodicBoundaryConverter != null)
            {
                periodicBoundaryConverter.Clear();
            }
        }
    }
}
