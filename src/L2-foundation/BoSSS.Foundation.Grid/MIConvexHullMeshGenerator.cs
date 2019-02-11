using System;
using System.Collections.Generic;
using ilPSP;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi
{
    static class MIConvexHullMeshGenerator
    {
        class MICHVertex : MIConvexHull.IVertex
        {
            static int ID_Counter = 0;
            public static void Clear()
            {
                ID_Counter = 0;
            }
            public readonly int ID;
            public MICHVertex()
            {
                ID = ID_Counter;
                ++ID_Counter;
            }

            public double[] Position { get; set; }
        }

        class MICHDelaunayCell : MIConvexHull.TriangulationCell<MICHVertex, MICHDelaunayCell>
        {
            static int ID_Counter = 0;
            public readonly int ID;
            public static void Clear()
            {
                ID_Counter = 0;
            }

            public MICHDelaunayCell()
            {
                ID = ID_Counter;
                ++ID_Counter;
            }

            public bool done = false;

            Vertex circumCenter;

            public Vertex Circumcenter 
            {
                get 
                {
                    circumCenter = circumCenter ?? GetCircumCenter();
                    return circumCenter; 
                }
                set 
                {
                    circumCenter = value;
                }
            }

            double Det(double[,] m)
            {
                return m[0, 0] * ((m[1, 1] * m[2, 2]) - (m[2, 1] * m[1, 2])) - m[0, 1] * (m[1, 0] * m[2, 2] - m[2, 0] * m[1, 2]) + m[0, 2] * (m[1, 0] * m[2, 1] - m[2, 0] * m[1, 1]);
            }

            double LengthSquared(double[] v)
            {
                double norm = 0;
                for (int i = 0; i < v.Length; i++)
                {
                    var t = v[i];
                    norm += t * t;
                }
                return norm;
            }

            Vertex GetCircumCenter()
            {
                // From MathWorld: http://mathworld.wolfram.com/Circumcircle.html

                var points = Vertices;

                double[,] m = new double[3, 3];

                // x, y, 1
                for (int i = 0; i < 3; i++)
                {
                    m[i, 0] = points[i].Position[0];
                    m[i, 1] = points[i].Position[1];
                    m[i, 2] = 1;
                }
                var a = Det(m);

                // size, y, 1
                for (int i = 0; i < 3; i++)
                {
                    m[i, 0] = LengthSquared(points[i].Position);
                }
                var dx = -Det(m);

                // size, x, 1
                for (int i = 0; i < 3; i++)
                {
                    m[i, 1] = points[i].Position[0];
                }
                var dy = Det(m);

                // size, x, y
                for (int i = 0; i < 3; i++)
                {
                    m[i, 2] = points[i].Position[1];
                }
                var c = -Det(m);

                var s = -1.0 / (2.0 * a);
                var r = System.Math.Abs(s) * System.Math.Sqrt(dx * dx + dy * dy - 4 * a * c);

                return new Vertex { Position = new Vector(s * dx, s * dy),
                                    ID = this.ID };
            }
        }

        class VariableCell : Cell
        {
            public VariableCell()
            {
                ridgeList = new List<Ridge>();
            }

            public bool boundary = false;

            List<Ridge> ridgeList;

            public void Insert(Ridge a)
            {
                ridgeList.Add(a);
            }

            public void Init()
            {
                if (!boundary)
                {
                    Ridges = GetRidges();
                    Vertices = GetVertices();
                }
            }

            //This should be possible in nlogn time!!111
            Ridge[] GetRidges()
            {
                //Remove zero Ridges
                for (int i = 0; i < ridgeList.Count(); ++i)
                {
                    Ridge c_i = ridgeList[i];
                    if (c_i.End.ID == c_i.Start.ID)
                    {
                        ridgeList.RemoveAt(i);
                    }
                }
                for (int i = 0; i < ridgeList.Count() - 2; ++i)
                {
                    Ridge c_i = ridgeList[i];
                    for ( int j = i + 1; j < ridgeList.Count(); ++j)
                    {
                        Ridge c_j = ridgeList[j];
                        if (c_i.End.ID == c_j.Start.ID)
                        {
                            //Ridge is in order
                            if(j != i + 1)
                            {
                                Ridge tmp = ridgeList[i + 1];
                                ridgeList[i + 1] = c_j;
                                ridgeList[j] = tmp;
                            }
                            break;
                        }
                    }
                }
                return ridgeList.ToArray();
            }

            Vertex[] GetVertices()
            {
                Vertex[] verts = new Vertex[Ridges.Length]; 
                for(int i = 0; i < Ridges.Length ; ++i)
                {
                    verts[i] = Ridges[i].Start;
                }
                return verts;
            }

        }

        class MICHIdMesh : SimpleIdMesh
        {
            static readonly double accuracy = 1e-20;

            public MICHIdMesh(  
                IEnumerable<MICHDelaunayCell> delaCells, 
                IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex, MICHDelaunayCell>> delaEdges, 
                int numberOfVoronois)
            {
                (VariableCell[] vCells, Vertex[] arrVertices) = 
                    CreateMeshLists(delaCells, delaEdges, numberOfVoronois);
                for(int i = 0; i < vCells.Length; ++i)
                {
                    vCells[i].Init();
                }
                cells = new List<Cell>(vCells);
                vertices = new List<Vertex>(arrVertices);

            }

            //Helper function
            static (int, int) OpposingIndices(int k)
            {
                if (k == 0)
                {
                    return (1, 2);
                }
                else if (k == 1)
                {
                    return (2, 0);
                }
                else
                {
                    return (0, 1);
                }
            }

            static (VariableCell[] cells, Vertex[] vertices) CreateMeshLists(   
                IEnumerable<MICHDelaunayCell> delaCells,
                IEnumerable<MIConvexHull.VoronoiEdge
                    <MICHVertex, MICHDelaunayCell>> edges,
                int numberOfVoronois)
            {
                //remove zero-Edges
                foreach(var edge in edges)
                {
                    Vertex start = edge.Source.Circumcenter;
                    Vertex end = edge.Target.Circumcenter;
                    double length = (start.Position - end.Position).AbsSquare();
                    if(length < accuracy)
                    {
                        edge.Source.Circumcenter = edge.Target.Circumcenter;
                    }
                }

                Vertex[] vertices = new Vertex[delaCells.Count()];
                VariableCell[] cells = new VariableCell[numberOfVoronois];

                foreach (MICHDelaunayCell delaCell in delaCells)
                {
                    delaCell.done = true;
                    MICHDelaunayCell[] neighbors = delaCell.Adjacency;
                    VariableCell[] voronoiCells = new VariableCell[3];
                    vertices[delaCell.ID] = delaCell.Circumcenter; 

                    //Create Voronoi Cells
                    //--------------------------------------------------------------
                    for(int i = 0; i < 3; ++i)
                    {
                        MICHVertex vert = delaCell.Vertices[i];
                        //Create or pick cell
                        VariableCell voronoiCell = null;
                        if (cells[vert.ID] == null)
                        {
                            voronoiCell = new VariableCell { ID = vert.ID, VoronoiNode = new Vector(vert.Position)};
                            cells[vert.ID] = voronoiCell;
                        }
                        else
                        {
                            voronoiCell = cells[vert.ID];
                        }
                        voronoiCells[i] = voronoiCell;
                    }
                    //Create Ridges for each neighbor
                    //--------------------------------------------------------------
                    for (int i = 0; i < 3; ++i)
                    {
                        MICHDelaunayCell neighbor = neighbors[i];
                        if (neighbor != null)
                        {
                            if (!neighbor.done)
                            {
                                (int k, int j) = OpposingIndices(i);
                                //Create Ridges
                                Ridge ridgeOutwards = new Ridge
                                {
                                    Start = delaCell.Circumcenter,
                                    End = neighbor.Circumcenter,
                                    Cell = voronoiCells[k]
                                };
                                Ridge ridgeInwards = new Ridge
                                {
                                    Start = neighbor.Circumcenter,
                                    End = delaCell.Circumcenter,
                                    Cell = voronoiCells[j]
                                };
                                ridgeInwards.Twin = ridgeOutwards;
                                ridgeOutwards.Twin = ridgeInwards;
                                //add to the two cells
                                
                                voronoiCells[k].Insert(ridgeOutwards);
                                voronoiCells[j].Insert(ridgeInwards);
                            }
                        }
                        else
                        {
                            (int k, int j) = OpposingIndices(i);
                            voronoiCells[k].boundary = true;
                            voronoiCells[j].boundary = true;
                        }
                    }
                }
                return (cells, vertices);
            }

        }

        static IntersectionMesh CreateMesh(IList<MICHVertex> startNodes)
        {
            var mICHMesh = MIConvexHull.VoronoiMesh.Create<
                MICHVertex,
                MICHDelaunayCell,
                MIConvexHull.VoronoiEdge<MICHVertex, MICHDelaunayCell>
                >(startNodes);
            var delaunayCells = mICHMesh.Vertices;
            var delaunayEdges = mICHMesh.Edges;
            MICHVertex.Clear();
            MICHDelaunayCell.Clear();

            return new IntersectionMesh(new MICHIdMesh(delaunayCells, delaunayEdges, startNodes.Count));
        }

        public static IntersectionMesh CreateMesh(MultidimensionalArray nodes)
        {
            MICHVertex[] startNodes = new MICHVertex[nodes.NoOfRows];
            for(int i = 0; i < nodes.NoOfRows; ++i)
            {
                startNodes[i] = new MICHVertex { Position = nodes.GetRowPt(i) };
            }
            return CreateMesh(startNodes); 
        }

        public static IntersectionMesh CreateMesh(IList<Vector> nodes)
        {
            MICHVertex[] startNodes = new MICHVertex[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                startNodes[i] = new MICHVertex { Position = nodes[i] };
            }
            return CreateMesh(startNodes);
        }

    }
}
