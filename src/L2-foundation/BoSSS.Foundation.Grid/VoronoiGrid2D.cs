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
            MultidimensionalArray nodePositions = CreateVoronoiNodesFromPolygon(PolygonBoundary, noOfNodeSeed);
            VoronoiNodes nodes = new VoronoiNodes(nodePositions);

            VoronoiInfo info = new VoronoiInfo {
                Boundary = new VoronoiBoundary {
                    Polygon = PolygonBoundary
                },
                NumberOfLloydIterations = NoOfLyyodsIter
            };
            return Polygonal(nodes, info);
        }

        static MultidimensionalArray CreateVoronoiNodesFromPolygon(Vector[] PolygonBoundary, int nSeedVoronois)
        {
            Debug.Assert(PolygonBoundary.Length > 0);
            int dim = PolygonBoundary[0].Dim;
            MultidimensionalArray positions = MultidimensionalArray.Create(nSeedVoronois, dim);
            Random rnd = new Random(0);

            double[,] max_min = FindMaxAndMin(PolygonBoundary);
            double[] scales = ScalesFromRandomIntervalToBoundingBox(dim, max_min);
            Vector center = CenterOfBoundingBox(dim, max_min);

            for (int i = 0; i < nSeedVoronois; ++i)
            {
                Vector randomNodePosition = RandomNodePositioninBoundingBox(rnd, scales, center);
                positions.SetRowPt(i, randomNodePosition);
            }
            positions.SetRowPt(0, PolygonBoundary[0]); 
            return positions;
        }

        static double[] ScalesFromRandomIntervalToBoundingBox(int dim, double[,] max_min) {
            double[] scales = new double[dim];
            for (int j = 0; j < dim; ++j)
            {
                scales[j] = (max_min[j, 0] - max_min[j, 1]);
            }
            return scales;
        }

        static Vector CenterOfBoundingBox(int dim, double[,] max_min) {
            Vector center = new Vector(dim);
            for (int j = 0; j < dim; ++j)
            {
                center[j] = (max_min[j, 0] + max_min[j, 1]) / 2;
            }
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

            VoronoiInfo info = new VoronoiInfo
            {
                Boundary = new VoronoiBoundary
                {
                    Polygon = PolygonBoundary
                },
                NumberOfLloydIterations = NoOfLyyodsIter,
                FirstCellNode_indice = FirstCellNode_Indice
            };
            return Polygonal(nodes, info);
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="nodes">
        /// Voronoi nodes: Center of each agglomerated cell. Will not be considered if outside of PolygonBoundary.
        /// </param>
        /// <param name="PolygonBoundary">
        /// Outer boundary of mesh. Is expected to be closed and must be non-overlapping.
        /// </param>
        /// <param name="NoOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="FirstCellNode_indice">
        /// Indice of node where the algorithm will start looking for the first Vector of PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            VoronoiNodes nodes,
            VoronoiInfo info)
        {
            List<VoronoiNode> nodeList = nodes.Nodes;
            Vector[] PolygonBoundary = info.Boundary.Polygon;
            int NoOfLyyodsIter = info.NumberOfLloydIterations;
            int FirstCellNode_indice = info.FirstCellNode_indice;

            // Create Voronoi mesh
            // =================================

            IEnumerator<Line> boundaryLines = Line.GetEnumerator(PolygonBoundary);
            IntersectionMesh voronoiMesh = null;
            Func<List<VoronoiNode>, int, IntersectionMesh> CreateMesh = MIConvexHullMeshGenerator.CreateMesh;

            for (int iLloyd = 0; iLloyd <= NoOfLyyodsIter; ++iLloyd)
            {
                // Voronoi generation
                //-------------------------------------
                AddFarNodes(nodeList, PolygonBoundary);
                voronoiMesh = CreateMesh(nodeList, FirstCellNode_indice);
                //Clip
                //-------------------------------------
                Intersecter.Intersect(voronoiMesh, boundaryLines);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                if (iLloyd != NoOfLyyodsIter)
                {
                    nodeList = RelaxNodes(voronoiMesh.GetInsideCells(), ref FirstCellNode_indice);
                }
            }

            return voronoiMesh.ToVoronoiGrid();
        }
        
        static List<VoronoiNode> RelaxNodes(IEnumerable<Cell> Cells, ref int FirstCellNode_indice)
        {
            List<VoronoiNode> nodes = new List<VoronoiNode>();
            foreach (Cell cell in Cells)
            {
                double relaxValue = 0.1;
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                CenterOfGravity = CenterOfGravity * relaxValue + new Vector(cell.Position) * (1 - relaxValue);
                nodes.Add(new VoronoiNode { Position = CenterOfGravity, GlobalID = cell.GlobalID});
                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = nodes.Count - 1;
                }
            }
            return nodes;
        }

        //creates random nodes in bounding box of PolygonBoundary, first node ist first entry of polygon boundary

        static void AddFarNodes(List<VoronoiNode> nodes, Vector[] PolygonBoundary)
        {
            Debug.Assert(nodes.Count > 0);
            double[,] max_min = FindMaxAndMin(PolygonBoundary);
            int dim = nodes[0].Position.Dim;
            for(int i = 0; i < dim; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    Vector far = new Vector(dim);
                    far[i] = max_min[i,j] * 10;
                    nodes.Add(new VoronoiNode { Position = far});
                }
            }
        }

        static double[,] FindMaxAndMin(Vector[] PolygonBoundary)
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
    }
}

