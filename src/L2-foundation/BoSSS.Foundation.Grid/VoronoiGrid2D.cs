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
        public static VoronoiGrid FromPolygonalDomain(
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter,
            int noOfNodeSeed)
        {
            //creates random nodes in bounding box of PolygonBoundary, first node ist first entry of polygon boundary
            List<Vector> nodes = CreateVoronoiNodesFromPolygon(PolygonBoundary, noOfNodeSeed);
            return FromPolygonalDomain(nodes, PolygonBoundary, NoOfLyyodsIter, 0);
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
        /// <param name="FirstCellNode_indice">
        /// Indice of node where the algorithm will start looking for the first Vector of PolygonBoundary.
        /// </param>
        /// <returns></returns>
        public static VoronoiGrid FromPolygonalDomain(
            MultidimensionalArray Nodes,
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter,
            int FirstCellNode_indice)
        {
            //Short hack
            List<Vector> nodes = new List<Vector>(Nodes.NoOfRows);
            for (int i = 0; i < Nodes.NoOfRows; ++i)
            {
                nodes.Add(new Vector(Nodes.GetRow(i)));
            }
            return FromPolygonalDomain(nodes, PolygonBoundary, NoOfLyyodsIter, FirstCellNode_indice);
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
        public static VoronoiGrid FromPolygonalDomain(
                List<Vector> nodes,
                Vector[] PolygonBoundary,
                int NoOfLyyodsIter,
                int FirstCellNode_indice
                )
        {
            // Create Voronoi mesh
            // =================================

            IEnumerator<Line> lines = Line.GetEnumerator(PolygonBoundary);
            IntersectionMesh voronoiMesh = null;
            Func<List<Vector>, int, IntersectionMesh> CreateMesh = MIConvexHullMeshGenerator.CreateMesh;

            for (int iLloyd = 0; iLloyd <= NoOfLyyodsIter; ++iLloyd)
            {
                // Voronoi generation
                //-------------------------------------
                AddFarNodes(nodes, PolygonBoundary);
                voronoiMesh = CreateMesh(nodes, FirstCellNode_indice);
                //Clip
                //-------------------------------------
                Intersecter.Intersect(voronoiMesh, lines);

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                if (iLloyd != NoOfLyyodsIter)
                {
                    nodes = RelaxVoronois(voronoiMesh.GetInsideCells(), ref FirstCellNode_indice);
                }
            }

            return voronoiMesh.ToVoronoiGrid();
        }

        static List<Vector> RelaxVoronois(IEnumerable<Cell> Cells, ref int FirstCellNode_indice)
        {
            List<Vector> nodes = new List<Vector>();
            foreach (Cell cell in Cells)
            {
                double relaxValue = 0.1;
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                nodes.Add(CenterOfGravity * relaxValue + new Vector(cell.VoronoiNode) * (1 - relaxValue));
                if (cell.ID == FirstCellNode_indice)
                {
                    FirstCellNode_indice = nodes.Count - 1;
                }
            }
            return nodes;
        }

        //creates random nodes in bounding box of PolygonBoundary, first node ist first entry of polygon boundary
        static List<Vector> CreateVoronoiNodesFromPolygon(Vector[] PolygonBoundary, int nSeedVoronois)
        {
            Debug.Assert(PolygonBoundary.Length > 0);
            int dim = PolygonBoundary[0].Dim;
            Vector[] nodes = new Vector[nSeedVoronois];
            Random rnd = new Random(0);

            double[,] max_min = FindMaxMin(PolygonBoundary);
            double[] scl = new double[dim];
            Vector center = new Vector(dim);
            for (int j = 0; j < dim; ++j)
            {
                scl[j] = (max_min[j, 0] - max_min[j, 1]);
                center[j] = (max_min[j, 0] + max_min[j, 1]) / 2;
            }

            for (int i = 0; i < nSeedVoronois; ++i)
            {
                nodes[i] = new Vector(dim);
                for(int j = 0; j < dim; ++j)
                {
                    nodes[i][j] = center[j] - 0.5 * scl[j] + rnd.NextDouble() * scl[j];
                }
            }
            nodes[0] = PolygonBoundary[0];
            return nodes.ToList();
        }

        static double[,] FindMaxMin(Vector[] PolygonBoundary)
        {
            Debug.Assert(PolygonBoundary.Length > 0);
            int dim = PolygonBoundary[0].Dim;
            double[,] max_min = new double[dim,2];
            for(int i = 0; i < dim; ++i)
            {
                max_min[i,0] = -double.MaxValue; //Look for maximum
                max_min[i, 1] = double.MaxValue; //Look for minimum
            }
            foreach(Vector vec in PolygonBoundary)
            {
                for (int i = 0; i < dim; ++i)
                {
                    if(max_min[i,0] < vec[i])
                    {
                        max_min[i,0] = vec[i];
                    }
                    if (max_min[i, 1] > vec[i])
                    {
                        max_min[i, 1] = vec[i];
                    }
                }
            }
            return max_min;
        }

        static void AddFarNodes(List<Vector> nodes, Vector[] PolygonBoundary)
        {
            Debug.Assert(nodes.Count > 0);
            double[,] max_min = FindMaxMin(PolygonBoundary);
            int dim = nodes[0].Dim;
            for(int i = 0; i < dim; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    Vector far = new Vector(dim);
                    far[i] = max_min[i,j] * 10;
                    nodes.Add(far);
                }
            }
        }
    }
}

