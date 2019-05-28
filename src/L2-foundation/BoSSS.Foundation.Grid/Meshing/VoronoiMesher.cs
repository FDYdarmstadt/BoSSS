using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public static class VoronoiMesher
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
        public static VoronoiGrid CreateGrid( VoronoiNodes nodes, Settings settings)
        {
            //Also possible with BoundaryMesh<VoronoiNode> and Ducktyping, maybe some other time.
            BoundaryMesh<Node> mesh = CreateBoundaryMesh(nodes, settings);
            VoronoiGrid grid = Convert2VoronoiGrid(mesh, settings.GridInfo);

            return grid;
        }

        public static MovingGrid CreateMovingGrid(VoronoiNodes nodes, Settings settings)
        {
            BoundaryMesh<Node> mesh = CreateBoundaryMesh(nodes, settings);

            ArrayMap node2nodeMap = ExtractMapping(mesh.GetNodes());
            VoronoiGrid grid = Convert2VoronoiGrid(mesh, settings.GridInfo);
            MovingGrid movingGrid = new MovingGrid
            {
                Grid = grid,
                Grid2PredecessorGrid = node2nodeMap
            };
            return movingGrid;
        }

        static BoundaryMesh<Node> CreateBoundaryMesh(VoronoiNodes nodes, Settings settings)
        {
            Mesher.Settings mesherSettings = ConvertToMesherSettings(settings);
            List<Node> mesherNodes = WrapInMesherNodes(nodes.Nodes);
            BoundaryMesh<Node> mesh = Mesher.Create(mesherNodes, mesherSettings);
            return mesh;
        }

        static Mesher.Settings ConvertToMesherSettings(Settings settings)
        {
            Mesher.Settings mesherSettings = new Mesher.Settings
            {
                Boundary = settings.GridInfo.Boundary,
                BoundingBox = settings.GridInfo.BoundingBox,
                NumberOfLloydIterations = settings.NumberOfLloydIterations,
                FirstCellNode_indice = settings.FirstCellNode_indice
            };
            return mesherSettings;
        }

        class Node : IMesherNode
        {
            VoronoiNode node;

            public ArrayConnection Type { get; set; }

            public Node(VoronoiNode node, int j)
            {
                this.node = node;
                Type = new ArrayConnection
                {
                    J = j,
                    Type = Connection.Remained
                };
            }

            public Node()
            {
                node = new VoronoiNode();
                Type = new ArrayConnection
                {
                    J = -1,
                    Type = Connection.Created
                };
            }
            
            public Vector Position {
                get { return node.Position; }
                set { node.Position = value; }
            }

            public static explicit operator VoronoiNode(Node node)
            {
                return node.node;
            }
        }

        static List<Node> WrapInMesherNodes(IList<VoronoiNode> voronoiNodes)
        {
            List<Node> wrappedNodes = new List<Node>(voronoiNodes.Count);
            for(int i = 0; i < voronoiNodes.Count; ++i)
            {
                Node wrappedNode = new Node(voronoiNodes[i], i);
                wrappedNodes.Add(wrappedNode);
            }
            return wrappedNodes;
        }

        static VoronoiGrid Convert2VoronoiGrid(BoundaryMesh<Node> mesh, VoronoiInfo info)
        {
            IReadOnlyList<MeshCell<Node>> cells = mesh.GetCells();
            (GridCommons grid, int[][] aggregation) = ExtractGridCommonsAndCellAggregation(cells);

            IList<Node> nodeList = mesh.GetNodes();
            IList<VoronoiNode> voronoiNodeList = Cast(nodeList);
            VoronoiNodes nodes = new VoronoiNodes(voronoiNodeList);

            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, info);
            return voronoiGrid;
        }

        static IList<VoronoiNode> Cast(IList<Node> nodes)
        {
            IList<VoronoiNode> voronoiNodes = new List<VoronoiNode>(nodes.Count);
            for(int i = 0; i < nodes.Count; ++i)
            {
                voronoiNodes.Add((VoronoiNode)nodes[i]);
            }
            return voronoiNodes;
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

        static ArrayMap ExtractMapping(IList<Node> processed)
        {
            ArrayConnection[] connections = new ArrayConnection[processed.Count];
            for(int i = 0; i < processed.Count; ++i)
            {
                connections[i] = processed[i].Type;
            }
            ArrayMap backwardsMap = new ArrayMap(connections);
            return backwardsMap;
        }
    }
}
