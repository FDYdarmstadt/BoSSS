using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform.LinAlg;
using BoSSS.Platform;

namespace BoSSS.Foundation.Grid.Voronoi
{
    /// <summary>
    /// Static methods to create Voronoi Meshes
    /// </summary>
    public static class VoronoiGrid2D
    {
        /// <summary>
        /// Create a random voronoi mesh inside a rectangle
        /// </summary>
        /// <param name="width"></param>
        /// <param name="height"></param>
        /// <param name="numberOfNodes"></param>
        /// <returns></returns>
        public static VoronoiGrid Rectangle(double width, double height, int numberOfNodes) {
            Vector[] polygonBoundary = new Vector[]
            {
                new Vector(-width / 2, height / 2),
                new Vector(width / 2, height / 2),
                new Vector(width / 2, -height / 2),
                new Vector(-width / 2, -height / 2)
            };
            return Polygonal(polygonBoundary, 10, numberOfNodes);
        }

        /// <summary>
        /// Creates a random voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="PolygonBoundary">
        /// Outer boundary of mesh. Is expected to be closed and must be non-overlapping.
        /// </param>
        /// <param name="NoOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="noOfNodeSeed">
        /// Number of random nodes that are placed in the bounding box of the PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter,
            int noOfNodeSeed)
        {
            //creates random nodes in bounding box of PolygonBoundary, first node ist first entry of polygon boundary
            Vector[] boundingBox = BoundingBox(PolygonBoundary);
            MultidimensionalArray nodePositions = RandomVoronoiNodesInBoundingBox(boundingBox, noOfNodeSeed);
            nodePositions.SetRowPt(0, PolygonBoundary[0]);
            VoronoiNodes nodes = new VoronoiNodes(nodePositions);

            VoronoiInfo info = new VoronoiInfo {
                MesherInfo = new VoronoiMesherInfo
                {
                    BoundingBox = boundingBox,
                    Boundary = PolygonBoundary,
                    NumberOfLloydIterations = NoOfLyyodsIter
                },
            };
            return Polygonal(nodes, info);
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="Nodes">
        /// Voronoi nodes: Center of each agglomerated cell. Will not be considered if outside of PolygonBoundary.
        /// </param>
        /// <param name="PolygonBoundary">
        /// Outer boundary of mesh. Is expected to be closed and must be non-overlapping.
        /// </param>
        /// <param name="NoOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="FirstCellNode_Indice">
        /// Indice of node where the algorithm will start looking for the first Vector of PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            MultidimensionalArray positions,
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter,
            int FirstCellNode_Indice)
        {
            //Short hack
            VoronoiNodes nodes = new VoronoiNodes(positions);
            Vector[] boundingBox = BoundingBox(PolygonBoundary);

            VoronoiInfo info = new VoronoiInfo
            {
                MesherInfo = new VoronoiMesherInfo
                {
                    BoundingBox = boundingBox,
                    Boundary = PolygonBoundary,
                    NumberOfLloydIterations = NoOfLyyodsIter,
                    FirstCellNode_indice = FirstCellNode_Indice
                },
            };
            return Polygonal(nodes, info);
        }

        static Vector[] BoundingBox(Vector[] polygon)
        {
            double[,] intervals = FindMaxAndMinInEachDimension(polygon);

            int dim = polygon[0].Dim;
            int numberOfCorners = (int)Math.Pow(2, dim);
            Vector[] corners = new Vector[numberOfCorners];
            for(int i = 0; i < numberOfCorners; ++i)
            {
                corners[i] = new Vector(dim);
            }
            SetBoxCorners(corners, intervals);
            return corners;
        }

        static void SetBoxCorners(Vector[] corners, double[,] intervals)
        {
            int dim = corners[0].Dim;
            int repeats = 1;
            for (int i_dim = 0; i_dim < dim; ++i_dim)
            {
                for (int i = 0, counter = 0; i < corners.Length; i += repeats)
                {
                    for (int j = 0; j < repeats; ++j)
                    {
                        corners[i + j][i_dim] = intervals[i_dim, counter];
                    }
                    counter = (counter + 1) % 2;
                }
                repeats *= 2;
            }
        }

        static double[,] FindMaxAndMinInEachDimension(Vector[] PolygonBoundary)
        {
            Debug.Assert(PolygonBoundary.Length > 0);
            int dim = PolygonBoundary[0].Dim;
            double[,] max_min = new double[dim, 2];
            for (int i = 0; i < dim; ++i)
            {
                max_min[i, 0] = -double.MaxValue; //Look for maximum
                max_min[i, 1] = double.MaxValue; //Look for minimum
            }
            foreach (Vector vec in PolygonBoundary)
            {
                for (int i = 0; i < dim; ++i)
                {
                    if (max_min[i, 0] < vec[i])
                    {
                        max_min[i, 0] = vec[i];
                    }
                    if (max_min[i, 1] > vec[i])
                    {
                        max_min[i, 1] = vec[i];
                    }
                }
            }
            return max_min;
        }

        static MultidimensionalArray RandomVoronoiNodesInBoundingBox(Vector[] boundingBox, int nSeedVoronois)
        {
            Debug.Assert(boundingBox.Length > 0);
            int dim = boundingBox[0].Dim;
            MultidimensionalArray positions = MultidimensionalArray.Create(nSeedVoronois, dim);
            Random rnd = new Random(0);

            double[] scales = ScalesFromRandomIntervalToBoundingBox(boundingBox);
            Vector center = CenterOfBoundingBox(boundingBox);

            for (int i = 0; i < nSeedVoronois; ++i)
            {
                Vector randomNodePosition = RandomNodePositioninBoundingBox(rnd, scales, center);
                positions.SetRowPt(i, randomNodePosition);
            }
            
            return positions;
        }

        static double[] ScalesFromRandomIntervalToBoundingBox(Vector[] boundingBox) {
            int dim = boundingBox[0].Dim;
            double[,] max_min = FindMaxAndMinInEachDimension(boundingBox);

            double[] scales = new double[dim];
            for (int j = 0; j < dim; ++j)
            {
                scales[j] = (max_min[j, 0] - max_min[j, 1]);
            }
            return scales;
        }

        static Vector CenterOfBoundingBox(Vector[] boundingBox) {
            int dim = boundingBox[0].Dim;
            int numberOfCorners = boundingBox.Length;
            Vector center = new Vector(dim);

            for (int j = 0; j < numberOfCorners; ++j)
            {
                center += boundingBox[j];
            }
            center /= numberOfCorners;

            return center;
        }

        static Vector RandomNodePositioninBoundingBox(Random rnd, double[] scales, Vector center)
        {
            int dim = center.Dim;
            Vector randomPosition = new Vector(dim);
            for (int j = 0; j < dim; ++j)
            {
                randomPosition[j] = center[j] - 0.5 * scales[j] + rnd.NextDouble() * scales[j];
            }
            return randomPosition;
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="nodes">
        /// Voronoi nodes: Center of each agglomerated cell. Will not be considered if outside of PolygonBoundary.
        /// </param>
        /// <param name="info">
        /// Contains information of Voronoi grid,e.g. a boundary polygon.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            VoronoiNodes nodes,
            VoronoiInfo info)
        {
            BoundaryMesh<VoronoiNode> mesh = VoronoiMesher<VoronoiNode>.Create(nodes.Nodes, info.MesherInfo);
            VoronoiGrid grid = Convert2VoronoiGrid(mesh, info);
            return grid;
        }

        static VoronoiGrid Convert2VoronoiGrid(BoundaryMesh<VoronoiNode> mesh, VoronoiInfo info)
        {
            IReadOnlyList<Cell<VoronoiNode>>cells = mesh.GetCells();
            (GridCommons grid, int[][] aggregation) = GetVoronoiData(cells);

            List<VoronoiNode> nodeList = mesh.GetNodes();
            VoronoiNodes nodes = new VoronoiNodes(nodeList);

            VoronoiGrid voronoiGrid = new VoronoiGrid(grid, aggregation, nodes, info);
            return voronoiGrid;
        }

        static (GridCommons grid, int[][] aggregation) GetVoronoiData(IReadOnlyList<Cell<VoronoiNode>> cells)
        {
            List<BoSSS.Foundation.Grid.Classic.Cell> cellsBoSSS = new List<BoSSS.Foundation.Grid.Classic.Cell>();
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
                    Cj.GlobalID = cellsBoSSS.Count;
                    Cj.Type = CellType.Triangle_3;
                    Cj.TransformationParams = MultidimensionalArray.Create(3, 2);
                    Cj.NodeIndices = new int[] { cell.Vertices[iV0].ID, cell.Vertices[iV1].ID, cell.Vertices[iV2].ID };
                    Cj.TransformationParams.SetRowPt(0, V0);
                    Cj.TransformationParams.SetRowPt(1, V1);
                    Cj.TransformationParams.SetRowPt(2, V2);

                    Agg2Pt[iTri] = cellsBoSSS.Count;
                    cellsBoSSS.Add(Cj);
                }
                aggregation.Add(Agg2Pt);
            }

            GridCommons grd;
            grd = new Grid2D(Triangle.Instance);
            grd.Cells = cellsBoSSS.ToArray();
            return (grd, aggregation.ToArray());
        }
    }
}

