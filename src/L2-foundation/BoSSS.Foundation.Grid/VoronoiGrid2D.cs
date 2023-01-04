using BoSSS.Foundation.Grid.Voronoi.Meshing;
using BoSSS.Platform.LinAlg;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Diagnostics;

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
        /// Create a regular checkerboard mesh
        /// </summary>
        /// <param name="xTics">
        /// Ascending tics of grid in x direction.
        /// </param>
        /// <param name="yTics">
        /// Ascending tics of grid in y direction.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Rectangle(double[] xTics, double[] yTics)
        {
            Vector[] rectangle = new Vector[]
            {
                new Vector(xTics[0], yTics[0]),
                new Vector(xTics[0], yTics[yTics.Length - 1]),
                new Vector(xTics[xTics.Length - 1], yTics[yTics.Length - 1]),
                new Vector(xTics[xTics.Length - 1], yTics[0])
            };
            MultidimensionalArray nodes = Checkerize(xTics, yTics);

            return Polygonal(nodes, rectangle, 0, 0);
        }

        static MultidimensionalArray Checkerize(double[] a, double[] b)
        {
            MultidimensionalArray abT = MultidimensionalArray.Create((a.Length - 1) * (b.Length - 1), 2);
            for(int i = 0; i < a.Length - 1; ++i)
            {
                for (int j = 0; j < b.Length - 1; ++j)
                {
                    abT[i * (b.Length - 1) + j, 0] = (a[i] + a[i +1]) / 2.0;
                    abT[i * (b.Length - 1) + j, 1] = (b[j] + b[j + 1]) / 2.0;
                }
            }
            return abT;
        }

        /// <summary>
        /// Creates a random voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="polygonBoundary">
        /// Outer boundary of mesh. Is expected to be closed and must be non-overlapping.
        /// </param>
        /// <param name="noOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="noOfNodeSeed">
        /// Number of random nodes that are placed in the bounding box of the PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            Vector[] polygonBoundary,
            int noOfLyyodsIter,
            int noOfNodeSeed)
        {
            Vector[] boundingBox = BoundingBox(polygonBoundary);
            VoronoiMesher.Settings mesherSettings = new VoronoiMesher.Settings
            {
                Boundary = new VoronoiBoundary()
                {
                    BoundingBox = boundingBox,
                    Polygon = polygonBoundary,
                    EdgeTags = DefaultEdgeTags(polygonBoundary.Length),
                },
                NumberOfLloydIterations = noOfLyyodsIter
            };

            VoronoiNodes nodes = GetVoronoiNodesIn(mesherSettings.Boundary, noOfNodeSeed);
            VoronoiMesher mesher = new VoronoiMesher(mesherSettings);

            return mesher.CreateGrid(nodes, 0);
        }

        //creates random nodes in bounding box of PolygonBoundary, 
        //first node coincides with first corner of polygon boundary
        static VoronoiNodes GetVoronoiNodesIn(VoronoiBoundary boundary, int amount)
        {
            MultidimensionalArray nodePositions = RandomVoronoiNodesInBoundingBox(boundary.BoundingBox, amount);
            nodePositions.SetRowPt(0, boundary.Polygon[0] + new Vector(0.01, - 0.01));
            VoronoiNodes nodes = new VoronoiNodes(nodePositions);
            return nodes;
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="boundary">
        /// Specifies polygonal boundary, edgetags, bounding box, ...
        /// </param>
        /// <param name="noOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="noOfNodeSeed">
        /// Number of random nodes that are placed in the bounding box of the PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            VoronoiBoundary boundary, 
            int noOfLyyodsIter,
            int noOfNodeSeed)
        {
            if (boundary.BoundingBox == null)
            {
                boundary.BoundingBox = BoundingBox(boundary.Polygon);
            }
            VoronoiMesher.Settings mesherSettings = new VoronoiMesher.Settings
            {
                Boundary = boundary,
                NumberOfLloydIterations = noOfLyyodsIter,
            };

            VoronoiNodes nodes = GetVoronoiNodesIn(mesherSettings.Boundary, noOfNodeSeed);
            VoronoiMesher mesher = new VoronoiMesher(mesherSettings);
            return mesher.CreateGrid(nodes, 0);
        }

        static byte[] DefaultEdgeTags(int count)
        {
            byte[] edgetag = new byte[count];
            edgetag.SetAll((byte)1);
            return edgetag;
        }

        /// <summary>
        /// Creates a voronoi mesh inside a polygon.
        /// </summary>
        /// <param name="nodePositions">
        /// Voronoi nodes: Center of each agglomerated cell. Will not be considered if outside of PolygonBoundary.
        /// </param>
        /// <param name="polygonBoundary">
        /// Outer boundary of mesh. Is expected to be closed and must be non-overlapping.
        /// </param>
        /// <param name="noOfLyyodsIter">
        /// Number of smoothing iterations.
        /// </param>
        /// <param name="firstCellNodeIndice">
        /// Indices of node where the algorithm will start looking for the first Vector of PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid Polygonal(
            MultidimensionalArray nodePositions,
            Vector[] polygonBoundary,
            int noOfLyyodsIter,
            int firstCellNodeIndice)
        {
            //Short hack
            VoronoiNodes nodes = new VoronoiNodes(nodePositions);
            Vector[] boundingBox = BoundingBox(polygonBoundary);

            VoronoiMesher.Settings settings = new VoronoiMesher.Settings
            {
                Boundary = new VoronoiBoundary()
                {
                    BoundingBox = boundingBox,
                    Polygon = polygonBoundary,
                    EdgeTags = DefaultEdgeTags(polygonBoundary.Length),
                },
                NumberOfLloydIterations = noOfLyyodsIter,
            };

            VoronoiMesher mesher = new VoronoiMesher(settings);
            return mesher.CreateGrid(nodes, firstCellNodeIndice);
        }

        public static VoronoiGrid Polygonal(
            MultidimensionalArray nodePositions, 
            VoronoiBoundary boundary, 
            int noOfLyyodsIter,
            int firstCellNodeIndice)
        {
            VoronoiNodes nodes = new VoronoiNodes(nodePositions);
            if(boundary.BoundingBox == null)
            {
                boundary.BoundingBox = BoundingBox(boundary.Polygon);
            }
            VoronoiMesher.Settings settings = new VoronoiMesher.Settings
            {
                Boundary = boundary,
                NumberOfLloydIterations = noOfLyyodsIter,
            };
            VoronoiMesher mesher = new VoronoiMesher(settings);
            return mesher.CreateGrid(nodes, firstCellNodeIndice);
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
    }
}