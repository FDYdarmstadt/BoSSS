using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    interface IVoronoiNodeCastable
    {
        VoronoiNode AsVoronoiNode();
    }

    class GridConverter
    {
        public static VoronoiGrid Convert2VoronoiGrid<T>(BoundaryMesh<T> mesh, VoronoiInfo info)
            where T : IVoronoiNodeCastable
        {
            IReadOnlyList<MeshCell<T>> cells = mesh.GetCells();
            (GridCommons grid, int[][] aggregation) = ExtractGridCommonsAndCellAggregation(cells);

            IList<T> nodeList = mesh.GetNodes();
            IList<VoronoiNode> voronoiNodeList = Cast(nodeList);
            VoronoiNodes nodes = new VoronoiNodes(voronoiNodeList);

            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, info);
            return voronoiGrid;
        }

        static (GridCommons grid, int[][] aggregation) ExtractGridCommonsAndCellAggregation<T>(IEnumerable<MeshCell<T>> cells)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cells_GridCommons = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregation = new List<int[]>();
            foreach (MeshCell<T> cell in cells)
            {
                //Convert to BoSSSCell : Triangulate
                Vector[] VoronoiCell = cell.Vertices.Select(voVtx => voVtx.Position).ToArray();
                int[,] iVtxTri = PolygonTesselation.TesselatePolygon(VoronoiCell);
                int[] Agg2Pt = new int[iVtxTri.GetLength(0)];

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


                    BoSSS.Foundation.Grid.Classic.Cell Cj = new BoSSS.Foundation.Grid.Classic.Cell();
                    Cj.GlobalID = cells_GridCommons.Count;
                    Cj.Type = CellType.Triangle_3;
                    Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                    Cj.NodeIndices = new int[] { cell.Vertices[iV0].ID, cell.Vertices[iV1].ID, cell.Vertices[iV2].ID };
                    Cj.TransformationParams.SetRowPt(0, V0);
                    Cj.TransformationParams.SetRowPt(1, V1);
                    Cj.TransformationParams.SetRowPt(2, V2);

                    Agg2Pt[iTri] = cells_GridCommons.Count;
                    cells_GridCommons.Add(Cj);
                }
                aggregation.Add(Agg2Pt);
            }

            GridCommons grid;
            grid = new Grid2D(Triangle.Instance);
            grid.Cells = cells_GridCommons.ToArray();
            return (grid, aggregation.ToArray());
        }

        static IList<VoronoiNode> Cast<T>(IList<T> nodes)
            where T : IVoronoiNodeCastable
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
