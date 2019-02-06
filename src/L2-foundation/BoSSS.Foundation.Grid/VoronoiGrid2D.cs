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
    public static class VoronoiGrid2D
    {
        public static AggregationGrid FromPolygonalDomain(
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter,
            int noOfNodeSeed)
        {
            List<Vector> nodes = CreateVoronoiNodesFromPolygon(PolygonBoundary, noOfNodeSeed);
            return FromPolygonalDomain(nodes, PolygonBoundary, NoOfLyyodsIter);
        }

        public static AggregationGrid FromPolygonalDomain(
            MultidimensionalArray Nodes,
            Vector[] PolygonBoundary,
            int NoOfLyyodsIter)
        {
            //Short hack
            List<Vector> nodes = new List<Vector>(Nodes.NoOfRows);
            for (int i = 0; i < Nodes.NoOfRows; ++i)
            {
                nodes.Add(new Vector(Nodes.GetRow(i)));
            }
            return FromPolygonalDomain(nodes, PolygonBoundary, NoOfLyyodsIter);
        }

        public static AggregationGrid FromPolygonalDomain(
                List<Vector> nodes,
                Vector[] PolygonBoundary,
                int NoOfLyyodsIter
                )
        {
            // Create Voronoi mesh
            // =================================

            IEnumerator<Line> lines = Line.GetEnumerator(PolygonBoundary);
            IntersectionMesh voronoiMesh = null;
            Func<List<Vector>, IntersectionMesh> CreateMesh = MIConvexHullMeshGenerator.CreateMesh;

            for (int iLloyd = 0; iLloyd <= NoOfLyyodsIter; ++iLloyd)
            {
                // Voronoi generation
                //-------------------------------------
                Stopwatch stopwatch = new Stopwatch();
                stopwatch.Start();
                AddFarNodes(nodes, PolygonBoundary);
                voronoiMesh = CreateMesh(nodes);
                stopwatch.Stop();
                Console.WriteLine(stopwatch.ElapsedMilliseconds);
                //Clip
                //-------------------------------------
                stopwatch.Restart();
                Intersecter.Intersect(voronoiMesh, lines);
                stopwatch.Stop();
                Console.WriteLine(stopwatch.ElapsedMilliseconds);
                //Console.ReadKey();

                // Lloyds algorithm (Voronoi relaxation)
                // -------------------------------------
                if (iLloyd != NoOfLyyodsIter)
                {
                    nodes = RelaxVoronois(voronoiMesh.GetInnerCells());
                }
            }
            return voronoiMesh.ToAggregationGrid();
        }

        static List<Vector> RelaxVoronois(IEnumerable<Cell> Cells)
        {
            int hack = 0;
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
                if (cell.ID == 3)
                {
                    hack = nodes.Count - 1;
                }
            }
            //------------------------------Achtung Hack!--------------------------------------------------------------
            Vector tmp = nodes[3];
            nodes[3] = new Vector(new double[] { -1, 1 });
            nodes[hack] = tmp;
            return nodes;
        }

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
            nodes[3] = PolygonBoundary[0];
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

