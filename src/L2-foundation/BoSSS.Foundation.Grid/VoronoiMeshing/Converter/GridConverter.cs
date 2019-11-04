using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    public interface IVoronoiNodeCastable
    {
        VoronoiNode AsVoronoiNode();
    }

    class GridConverter<T>
        where T : IVoronoiNodeCastable
    {
        readonly VoronoiBoundary boundary;

        readonly BoundaryConverter boundaryConverter;

        public GridConverter(VoronoiBoundary boundary, PeriodicBoundary periodicBoundary = null)
        {
            this.boundary = boundary;
            boundaryConverter = new BoundaryConverter(boundary, periodicBoundary);
        }

        public VoronoiGrid ConvertToVoronoiGrid(
            Mesh<T> mesh)
        {
            (GridCommons grid, int[][] aggregation) = ExtractGridCommonsAndCellAggregation(mesh.Cells, boundaryConverter);
            VoronoiNodes nodes = ExtractVoronoiNodes(mesh);
            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, boundary);

            return voronoiGrid;
        }

        static VoronoiNodes ExtractVoronoiNodes(Mesh<T> mesh)
        {
            IList<T> nodeList = mesh.Nodes;
            IList<VoronoiNode> voronoiNodeList = CastAsVoronoiNodes(nodeList);
            VoronoiNodes nodes = new VoronoiNodes(voronoiNodeList);
            return nodes;
        }

        static (GridCommons, int[][]) ExtractGridCommonsAndCellAggregation(
            IEnumerable<MeshCell<T>> cells,
            BoundaryConverter boundaryConverter)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cells_GridCommons = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregation = new List<int[]>();

            foreach (MeshCell<T> cell in cells)
            {
                //Convert to BoSSSCell : Triangulate
                Vector[] VoronoiCell = cell.Vertices.Select(voVtx => voVtx.Position).ToArray();
                int[,] iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell);
                int[] Agg2Pt = new int[iVtxTri.GetLength(0)];

                bool isBoundaryCell = IsBoundary(cell);

                for (int iTri = 0; iTri < iVtxTri.GetLength(0); iTri++)
                { // loop over triangles of voronoi cell
                    int iV0 = iVtxTri[iTri, 0];
                    int iV1 = iVtxTri[iTri, 1];
                    int iV2 = iVtxTri[iTri, 2];

                    Vector V0 = VoronoiCell[iV0];
                    Vector V1 = VoronoiCell[iV1];
                    Vector V2 = VoronoiCell[iV2];

                    Vector D1 = V1 - V0;
                    Vector D2 = V2 - V0;

                    if (D1.CrossProduct2D(D2) < 0)
                    {
                        int it = iV0;
                        iV0 = iV2;
                        iV2 = it;

                        Vector vt = V0;
                        V0 = V2;
                        V2 = vt;

                        D1 = V1 - V0;
                        D2 = V2 - V0;
                    }

                    Debug.Assert(D1.CrossProduct2D(D2) > 1.0e-8);

                    Cell Cj = new Cell()
                    {
                        GlobalID = cells_GridCommons.Count,
                        Type = CellType.Triangle_3,
                    };
                    Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                    Cj.TransformationParams.SetRowPt(0, V0);
                    Cj.TransformationParams.SetRowPt(1, V1);
                    Cj.TransformationParams.SetRowPt(2, V2);

                    Agg2Pt[iTri] = cells_GridCommons.Count;
                    cells_GridCommons.Add(Cj);

                    //Save BoundaryInformation
                    if (isBoundaryCell)
                    {
                        List<BoundaryFace> tags = GetBoundaryFacesOf(cell, iV0, iV1, iV2);
                        boundaryConverter.RegisterBoundaries(Cj, tags);
                    }
                    Cj.NodeIndices = new int[]
                    {
                        cell.Vertices[iV0].ID,
                        cell.Vertices[iV1].ID,
                        cell.Vertices[iV2].ID
                    };
                }
                aggregation.Add(Agg2Pt);
            }
            GridCommons grid = new Grid2D(Triangle.Instance)
            {
                Cells = cells_GridCommons.ToArray()
            };
            boundaryConverter.RegisterEdgesTo(grid);
            return (grid, aggregation.ToArray());
        }

        static List<BoundaryFace> GetBoundaryFacesOf(MeshCell<T> cell, int iV0, int iV1, int iV2)
        {
            //Indices are debug magic. FML
            List<BoundaryFace> tags = new List<BoundaryFace>(3);
            int max = cell.Edges.Length;
            if (iV0 + 1 == iV1 || iV0 - max + 1 == iV1)
            {
                IfIsBoundaryAddEdge2Tags(iV0, 0);
            }
            if (iV1 + 1 == iV2 || iV1 - max + 1 == iV2)
            {
                IfIsBoundaryAddEdge2Tags(iV1, 1);
            }
            if (iV2 + 1 == iV0 || iV2 - max + 1 == iV0)
            {
                IfIsBoundaryAddEdge2Tags(iV2, 2);
            }
            return tags;

            void IfIsBoundaryAddEdge2Tags(int iV, int iFace)
            {
                Edge<T> edge = cell.Edges[iV];
                if (edge.IsBoundary)
                {
                    BoundaryFace tag = new BoundaryFace
                    {
                        Face = iFace,
                        BoundaryEdgeNumber = edge.BoundaryEdgeNumber,
                        ID = edge.Start.ID,
                        NeighborID = edge.Twin.Start.ID
                    };
                    tags.Add(tag);
                }
            }
        }

        static bool IsBoundary(MeshCell<T> cell)
        {
            foreach (Edge<T> edge in cell.Edges)
            {
                if (edge.IsBoundary)
                    return true;
            }
            return false;
        }
        
        static IList<VoronoiNode> CastAsVoronoiNodes(IList<T> nodes)
        {
            IList<VoronoiNode> voronoiNodes = new List<VoronoiNode>(nodes.Count);
            for (int i = 0; i < nodes.Count; ++i)
            {
                voronoiNodes.Add(nodes[i].AsVoronoiNode());
            }
            return voronoiNodes;
        }
    }
}
