using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace BoSSS.Foundation.Grid.Voronoi
{
    

    static class VoronoiMesher
    {
        public class Settings
        {
            public VoronoiInfo GridInfo;
            public int NumberOfLloydIterations = 10;
            public int FirstCellNode_indice = 0;
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="nodes">
        /// Voronoi nodes: Center of each agglomerated cell. Will not be considered if outside of PolygonBoundary.
        /// </param>
        /// <param name="settings">
        /// Contains information of Voronoi grid,e.g. a boundary polygon.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Create(
            VoronoiNodes nodes,
            Settings settings)
        {
            BoundaryMesh<VoronoiNode> mesh = Create(nodes.Nodes, settings);
            VoronoiGrid grid = Convert2VoronoiGrid(mesh, settings.GridInfo);

            return grid;
        }

        static BoundaryMesh<T> Create<T>(IList<T> nodeList, Settings settings)
            where T : INode, new()
        {
            // Create Voronoi mesh
            // =================================
            IEnumerator<Line> boundaryLines = Line.GetEnumerator(settings.GridInfo.Boundary);
            IntersectionMesh<T> mesh = null;

            for (int iLloyd = 0; iLloyd <= settings.NumberOfLloydIterations; ++iLloyd)
            {
                // Mesh generation
                //-------------------------------------
                AddDistantBoundingNodes(nodeList, settings.GridInfo.BoundingBox);
                mesh = IntersectionMeshGenerator.CreateMesh(nodeList, settings.FirstCellNode_indice);
                //Clip
                //-------------------------------------
                Intersecter.Intersect(mesh, boundaryLines);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                IReadOnlyList<Cell<T>> cells = mesh.GetCells(); //GetInsideCell : Return Cells in order of MeshArray.
                if (iLloyd != settings.NumberOfLloydIterations)
                {
                    MoveNodesTowardsCellCenter(cells, ref settings.FirstCellNode_indice); 
                }
                nodeList = mesh.GetNodes();
            }
            return mesh;
        }

        static void AddDistantBoundingNodes<T>(IList<T> nodes, Vector[] boundingBox)
            where T : INode, new()
        {
            Debug.Assert(nodes.Count > 0);
            foreach(Vector corner in boundingBox)
            {
                T cornerNode = new T();
                cornerNode.Position = corner * 10;
                nodes.Add(cornerNode);
            }
        }

        static void MoveNodesTowardsCellCenter<T>(IReadOnlyList<Cell<T>> Cells, ref int FirstCellNode_indice)
            where T : INode, new()
        {
            //Mark inside nodes
            //Use Original Nodes List and update. Use LinkedList?! Only iterate and cut some nodes, or insert
            //Let's give it a try!
            for (int i = 0; i < Cells.Count; ++i)
            {
                Cell<T> cell = Cells[i];
                double relaxValue = 0.1;
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                CenterOfGravity = CenterOfGravity * relaxValue + new Vector(cell.Node.Position) * (1 - relaxValue);

                cell.Node.Position = CenterOfGravity;
                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = i;
                }
            }
        }

        static VoronoiGrid Convert2VoronoiGrid(BoundaryMesh<VoronoiNode> mesh, VoronoiInfo info)
        {
            IReadOnlyList<Cell<VoronoiNode>> cells = mesh.GetCells();
            (GridCommons grid, int[][] aggregation) = ExtracGridCommonsAndCellAggregation(cells);

            List<VoronoiNode> nodeList = mesh.GetNodes();
            VoronoiNodes nodes = new VoronoiNodes(nodeList);

            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, info);
            return voronoiGrid;
        }

        static (GridCommons grid, int[][] aggregation) ExtracGridCommonsAndCellAggregation(IEnumerable<Cell<VoronoiNode>> cells)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cells_GridCommons = new List<BoSSS.Foundation.Grid.Classic.Cell>();
            List<int[]> aggregation = new List<int[]>();
            foreach (Cell<VoronoiNode> cell in cells)
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
    }
}
