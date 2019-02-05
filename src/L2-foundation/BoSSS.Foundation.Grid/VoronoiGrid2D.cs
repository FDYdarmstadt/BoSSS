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
    class VoronoiGrid2D //: AggregationGrid
    {
        VoronoiGrid2D(MultidimensionalArray Nodes, Vector[] PolygonBoundary, double NoOfLyyodsIter)
        {
            //bla bla bla 
        }
    }

    class ArrayEnum<T> : IEnumerator<T>
        {
            int pointer;
            IList<T> arr;
            public ArrayEnum(IList<T> Arr)
            {
                arr = Arr;
                pointer = -1;
            }

            public T Current => arr[pointer];

            object IEnumerator.Current => Current;

            public void Dispose()
            {
            }

            public bool MoveNext()
            {
                return (++pointer < arr.Count);
            }

            public void Reset()
            {
                pointer = -1;
            }
        }

    class Line
    {
        public static Line[] ToLines(Vector[] polygon)
        {
            Line[] lines = new Line[polygon.Length];
            Line line;
            for (int i = 0; i < polygon.Length - 1; ++i)
            {
                line = new Line { start = new Vertex(), end = new Vertex() };
                line.start.Position = polygon[i];
                line.end.Position = polygon[i + 1];
                lines[i] = line;
            }
            line = new Line { start = new Vertex(), end = new Vertex() };
            line.start.Position = polygon[polygon.Length - 1];
            line.end.Position = polygon[0];
            lines[lines.Length - 1] = line;

            return lines;
        }
        public static IEnumerator<Line> GetEnumerator(Vector[] polygon)
        {
            Line[] lines = ToLines(polygon);
            return new ArrayEnum<Line>(lines);
        }

        public Vertex start { get; set; }
        public Vertex end { get; set; }
    }

    static class VoronoiMeshGen
    {
        static void RelaxVoronois(IReadOnlyList<Cell> Cells, MultidimensionalArray Nodes)
        {
            for (int i = 8; i < Nodes.NoOfRows; ++i)
            {
                double relaxValue = 0.1;
                Cell cell = Cells[i];
                Vector CenterOfGravity = new Vector(0, 0);
                foreach (Vertex vertex in cell.Vertices)
                {
                    CenterOfGravity += vertex.Position;
                }
                CenterOfGravity.Scale(1.0 / cell.Vertices.Length);
                Nodes.SetRowPt(i, CenterOfGravity * relaxValue + Nodes.GetRowPt(i) * (1 - relaxValue));
            }
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

        /* boundaries of domain as 
            * 
            * 
            * 
            */
        static public AggregationGrid FromPolygonalDomain(  MultidimensionalArray Nodes,
                                                            Vector[] PolygonBoundary,
                                                            double NoOfLyyodsIter)
        {
            //Short hack
            List<Vector> nodes = new List<Vector>(Nodes.NoOfRows);
            for (int i = 0; i < Nodes.NoOfRows; ++i)
            {
                nodes.Add(new Vector(Nodes.GetRow(i)));
            }
            // Create Voronoi mesh
            // =================================

            IEnumerator<Line> lines = Line.GetEnumerator(PolygonBoundary);
            IntersectionMesh voronoiMesh = null;
            Func<List<Vector>, IntersectionMesh> CreateMesh = MIConvexHullMeshGenerator.CreateMesh;

            int LloydIterations = (int)Math.Ceiling(NoOfLyyodsIter);
            for (int iLloyd = 0; iLloyd <= LloydIterations; ++iLloyd)
            {
                // Voronoi generation
                //-------------------------------------
                Stopwatch stopwatch = new Stopwatch();
                stopwatch.Start();
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
                if (iLloyd != LloydIterations)
                {
                    nodes = RelaxVoronois(voronoiMesh.GetInnerCells());
                    //---------------------------------------------------------------Achtung Hack!-----------------------
                    nodes.Add(new Vector(new double[] { -10, 0 }));
                    nodes.Add(new Vector(new double[] { 10, 0 }));
                    nodes.Add(new Vector(new double[] { 0, 10 }));
                    nodes.Add(new Vector(new double[] { 0, -10 }));
                }
            }
            return voronoiMesh.ToAggregationGrid();

        }
    }
}

