using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using log4net.Core;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    class GridConverter<T>
        where T : IVoronoiNodeCastable
    {
        readonly VoronoiBoundary boundary;

        readonly BoundaryConverter boundaryConverter;

        public GridConverter(VoronoiBoundary boundary, PeriodicMap periodicMap = null)
        {
            this.boundary = boundary;
            boundaryConverter = new BoundaryConverter(boundary, periodicMap);
        }

        public VoronoiGrid ConvertToVoronoiGrid(
            IMesh<T> mesh)
        {
            boundaryConverter.Clear();
            (GridCommons grid, int[][] aggregation) = ExtractGridCommonsAndCellAggregation(mesh.Cells, boundaryConverter);
            VoronoiNodes nodes = ExtractVoronoiNodes(mesh);
            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, boundary);
            voronoiGrid.FirstCornerNodeIndice = mesh.CornerNodeIndice;

            return voronoiGrid;
        }

        static VoronoiNodes ExtractVoronoiNodes(IMesh<T> mesh)
        {
            IList<T> nodeList = mesh.Nodes;
            IList<VoronoiNode> voronoiNodeList = CastAsVoronoiNodes(nodeList);
            VoronoiNodes nodes = new VoronoiNodes(voronoiNodeList);
            return nodes;
        }

        static void SortDescendingArea(int[,] triangles, Vector[] vertices)
        {
            //Quicksort for fun, lazy https://en.wikipedia.org/wiki/Quicksort implementation
            double[] areas = new double[triangles.GetLength(0)];
            for(int i = 0; i < triangles.GetLength(0); ++i)
            {
                int iV0 = triangles[i, 0];
                int iV1 = triangles[i, 1];
                int iV2 = triangles[i, 2];

                Vector V0 = vertices[iV0];
                Vector V1 = vertices[iV1];
                Vector V2 = vertices[iV2];

                Vector D1 = V1 - V0;
                Vector D2 = V2 - V0;
                areas[i] = D1.CrossProduct2D(D2);
            }

            Quicksort(0, areas.Length -1);

            void Quicksort(int low, int high)
            {
                if(low < high)
                {
                    int p = Partition(low, high);
                    Quicksort(low, p -1);
                    Quicksort(p + 1, high);
                }
            }

            int Partition(int low, int high)
            {
                double pivot = areas[high];
                int i = low;
                for(int j = low; j < high; ++j)
                {
                    if(areas[j] > pivot)
                    {
                        double temp = areas[j];
                        areas[j] = areas[i];
                        areas[i] = temp;
                        SwitchTriangles(i, j);
                        ++i;
                    }
                }
                double temp1 = areas[high];
                areas[high] = areas[i];
                areas[i] = temp1;
                SwitchTriangles(i, high);
                return i;
            }

            void SwitchTriangles(int i, int j)
            {
                for(int k = 0; k < 3; ++k)
                {
                    int temp = triangles[j, k];
                    triangles[j, k] = triangles[i, k];
                    triangles[i, k] = temp;
                }
            }
        }

        static (int, int, int) Rearrange(int iV0, int iV1, int iV2, Vector V0, Vector V1, Vector V2)
        {
            double D01 = (V1 - V0).AbsSquare();  
            double D12 = (V2 - V1).AbsSquare();  
            double D20 = (V0 - V2).AbsSquare();
            if(D01 > D20 && D12 > D20)
            {
                return (iV1, iV2, iV0);
            }
            else if ( D12 > D01 && D20 > D01)
            {
                return (iV2, iV0, iV1);
            }
            else
            {
                return ( iV0, iV1, iV2);
            }
        }

        static (GridCommons, int[][]) ExtractGridCommonsAndCellAggregation(
            IEnumerable<MeshCell<T>> cells,
            BoundaryConverter boundaryConverter)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cellsGridCommons = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregation = new List<int[]>();

            foreach (MeshCell<T> cell in cells)
            {
                //Convert to BoSSSCell : Triangulate
                Vector[] voronoiCellVertices = cell.Vertices.Select(voVtx => voVtx.Position).ToArray();
                int[,] iVtxTri = PolygonTesselation.TesselatePolygon(voronoiCellVertices);
                //SortDescendingArea(iVtxTri, voronoiCellVertices);
                int[] Agg2Pt = new int[iVtxTri.GetLength(0)];

                bool isBoundaryCell = IsBoundary(cell);

                for (int iTri = 0; iTri < iVtxTri.GetLength(0); iTri++)
                { // loop over triangles of voronoi cell
                    int iV0 = iVtxTri[iTri, 0];
                    int iV1 = iVtxTri[iTri, 1];
                    int iV2 = iVtxTri[iTri, 2];
                    
                    Vector V0 = voronoiCellVertices[iV0];
                    Vector V1 = voronoiCellVertices[iV1];
                    Vector V2 = voronoiCellVertices[iV2];

                    (iV0, iV1, iV2) = Rearrange(iV0, iV1, iV2, V0, V1, V2);

                    V0 = voronoiCellVertices[iV0];
                    V1 = voronoiCellVertices[iV1];
                    V2 = voronoiCellVertices[iV2];

                    Vector D1 = V1 - V0;
                    Vector D2 = V2 - V0;

                    bool positive = true;
                    if (D1.CrossProduct2D(D2) < 0)
                    {
                        if (isBoundaryCell)
                        {
                            Console.WriteLine($"Cell{cell.Node.AsVoronoiNode().Position} not positive. It is probably not convex.");
                            positive = false;
                            int it = iV0;
                            iV0 = iV2;
                            iV2 = it;

                            Vector vt = V0;
                            V0 = V2;
                            V2 = vt;

                            D1 = V1 - V0;
                            D2 = V2 - V0;
                        }
                        else
                        {
                            throw new Exception($"Created negatively oriented cell. Node = {cell.Node.AsVoronoiNode().Position}");
                        }

                    }

                    //if(D1.CrossProduct2D(D2) < 1.0e-7)
                    //{
                    //    Console.WriteLine($"Created very small cell. Node = {cell.Node.AsVoronoiNode().Position}, Area= {D1.CrossProduct2D(D2)}, positive = {positive}");
                    //}

                    Cell Cj = new Cell()
                    {
                        GlobalID = cellsGridCommons.Count,
                        Type = CellType.Triangle_3,
                    };
                    Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                    Cj.TransformationParams.SetRowPt(0, V0);
                    Cj.TransformationParams.SetRowPt(1, V1);
                    Cj.TransformationParams.SetRowPt(2, V2);

                    Agg2Pt[iTri] = cellsGridCommons.Count;
                    cellsGridCommons.Add(Cj);

                    //Save BoundaryInformation
                    if (isBoundaryCell)
                    {
                        List<BoundaryFace> tags;
                        if (positive)
                        {
                            tags = GetBoundaryFacesOfTriangle(cell, iV0, iV1, iV2);
                        }
                        else
                        {
                            tags = GetBoundaryFacesOfNegativeTriangle(cell, iV0, iV1, iV2);
                        }
                        boundaryConverter.RegisterBoundaries(Cj, tags);
                    }
                    Cj.NodeIndices = new long[]
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
                Cells = cellsGridCommons.ToArray()
            };
            PrintEdgeTags(cellsGridCommons);
            boundaryConverter.RegisterEdgesTo(grid);
            return (grid, aggregation.ToArray());
        }

        static void PrintEdgeTags(List<BoSSS.Foundation.Grid.Classic.Cell> cellsGridCommons)
        {
            using (StreamWriter sw = new StreamWriter("EdgeTags.txt"))
            {
                foreach (Cell cell in cellsGridCommons)
                {
                    for (int i = 0; i < (cell.CellFaceTags?.Length ?? 0); ++i)
                    {
                        CellFaceTag tag = cell.CellFaceTags[i];
                        double x = cell.TransformationParams[(tag.FaceIndex) % 3, 0]
                            + cell.TransformationParams[(tag.FaceIndex + 1) % 3, 0];
                        x /= 2;
                        double y = cell.TransformationParams[(tag.FaceIndex) % 3, 1]
                            + cell.TransformationParams[(tag.FaceIndex + 1) % 3, 1];
                        y /= 2;
                        sw.WriteLine($"{x}, {y}, {tag.EdgeTag}");
                    }
                }
            }
        }

        static List<BoundaryFace> GetBoundaryFacesOfTriangle(MeshCell<T> cell, int iV0, int iV1, int iV2)
        {
            //Indices are debug magic. FML
            List<BoundaryFace> tags = new List<BoundaryFace>(3);
            int max = cell.Edges.Length;
            if (iV0 + 1 == iV1 || iV0 - max + 1 == iV1)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV0], 0, tags);
            }
            if (iV1 + 1 == iV2 || iV1 - max + 1 == iV2)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV1], 1, tags);
            }
            if (iV2 + 1 == iV0 || iV2 - max + 1 == iV0)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV2], 2, tags);
            }
            return tags;
        }

        static void IfIsBoundaryAddEdge2Tags(Edge<T> edge, int iFace, List<BoundaryFace> tags)
        {
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

        static List<BoundaryFace> GetBoundaryFacesOfNegativeTriangle(MeshCell<T> cell, int iV0, int iV1, int iV2)
        {
            //Indices are debug magic. FML
            List<BoundaryFace> tags = new List<BoundaryFace>(3);
            int max = cell.Edges.Length;
            if (iV0 - 1 == iV1 || iV1 - max + 1 == iV0)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV1], 0, tags);
            }
            if (iV1 - 1 == iV2 || iV2 - max + 1 == iV1)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV2], 1, tags);
            }
            if (iV2 - 1 == iV0 || iV0 - max + 1 == iV2)
            {
                IfIsBoundaryAddEdge2Tags(cell.Edges[iV0], 2, tags);
            }
            return tags;
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
