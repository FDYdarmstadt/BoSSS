using System;
using System.Collections.Generic;
using ilPSP;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Platform.LinAlg;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class IntersectionMeshGenerator
    {
        public static IntersectionMesh<T> CreateMesh<T>(IList<T> nodes)
            where T : IMesherNode, new()
        {
            ResetDataIDCounters<T>();
            MICHVertex<T>[] startNodes = new MICHVertex<T>[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                startNodes[i] = new MICHVertex<T>(nodes[i].Position, nodes[i]);
            }
            IntersectionMesh<T> mesh = CreateMesh(startNodes);
            return mesh;
        }

        public static IntersectionMesh<T> CreateMesh<T>(IList<T> nodes, int startCell_NodeIndice)
            where T : IMesherNode, new()
        {
            ResetDataIDCounters<T>();
            MICHVertex<T>[] startNodes = new MICHVertex<T>[nodes.Count];
            for (int i = 0; i < nodes.Count; ++i)
            {
                startNodes[i] = new MICHVertex<T>(nodes[i].Position, nodes[i]);
            }
            IntersectionMesh<T> mesh = CreateMesh(startNodes, startCell_NodeIndice);
            return mesh;
        }

        static void ResetDataIDCounters<T>()
        {
            MICHVertex<T>.Clear();
            MICHDelaunayCell<T>.Clear();
        }

        static IntersectionMesh<T> CreateMesh<T>(IList<MICHVertex<T>> startNodes)
            where T : IMesherNode, new()
        {
            SimpleIdMesh<T> mesh = MICMesher<T>.Create(startNodes);
            IntersectionMesh<T> intersectionMesh = new IntersectionMesh<T>(mesh);

            return intersectionMesh;
        }

        static IntersectionMesh<T> CreateMesh<T>(IList<MICHVertex<T>> startNodes, int startCell_NodeIndice)
            where T : IMesherNode, new()
        {
            SimpleIdMesh<T> mesh = MICMesher<T>.Create(startNodes);
            IntersectionMesh<T> intersectionMesh = new IntersectionMesh<T>(mesh, startCell_NodeIndice);
            return intersectionMesh;
        }
    }
    
    class MICHVertex<T> : MIConvexHull.IVertex
    {
        public T Node;

        static int ID_Counter = 0;
        public static void Clear()
        {
            ID_Counter = 0;
        }
        public readonly int ID;

        public MICHVertex(Vector position, T node)
        {
            Position = position;
            Node = node;
            ID = ID_Counter;
            ++ID_Counter;
        }

        public double[] Position { get; set; }
    }

    class MICHDelaunayCell<T> : MIConvexHull.TriangulationCell<MICHVertex<T>, MICHDelaunayCell<T>>
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

    class VariableCell<T> : MeshCell<T>
    {
        public VariableCell()
        {
            ridgeList = new List<Edge<T>>();
        }

        public bool boundary = false;

        List<Edge<T>> ridgeList;

        public void Insert(Edge<T> a)
        {
            ridgeList.Add(a);
        }

        public void Init()
        {
            if (!boundary)
            {
                Edges = GetRidges();
                Vertices = GetVertices();
            }
            else
            {
                //Empty Ridges
                Edges = ridgeList.ToArray();
                for(int i = 0; i < Edges.Length; ++i)
                {
                    Edges[i].IsBoundary = true;
                }
            }
                
        }

        //This should be possible in nlogn time!!111
        Edge<T>[] GetRidges()
        {
            //Remove zero Ridges
            for (int i = 0; i < ridgeList.Count(); ++i)
            {
                Edge<T> c_i = ridgeList[i];
                if (c_i.End.ID == c_i.Start.ID)
                {
                    ridgeList.RemoveAt(i);
                }
            }
            for (int i = 0; i < ridgeList.Count() - 2; ++i)
            {
                Edge<T> c_i = ridgeList[i];
                for ( int j = i + 1; j < ridgeList.Count(); ++j)
                {
                    Edge<T> c_j = ridgeList[j];
                    if (c_i.End.ID == c_j.Start.ID)
                    {
                        //Ridge is in order
                        if(j != i + 1)
                        {
                            Edge<T> tmp = ridgeList[i + 1];
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
            Vertex[] verts = new Vertex[Edges.Length]; 
            for(int i = 0; i < Edges.Length ; ++i)
            {
                verts[i] = Edges[i].Start;
            }
            return verts;
        }

    }

    static class MICMesher<T>
    {
        static readonly double accuracy = 1e-20;

        public static SimpleIdMesh<T> Create(IList<MICHVertex<T>> startNodes)
        {
            var mICHMesh = CreateVoronoiMeshFromDLL(startNodes);
            var delaunayCells = mICHMesh.Vertices;
            var delaunayEdges = mICHMesh.Edges;
           
            SimpleIdMesh<T> mesh = MeshFromCellsAndEdges(delaunayCells, delaunayEdges, startNodes.Count);
            return mesh;
        }

        static MIConvexHull.VoronoiMesh<MICHVertex<T>, MICHDelaunayCell<T>,
            MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>>
            CreateVoronoiMeshFromDLL(IList<MICHVertex<T>> nodes)
        {
            var mICHMesh = MIConvexHull.VoronoiMesh.Create<
                    MICHVertex<T>,
                    MICHDelaunayCell<T>,
                    MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>
                >(nodes);
            
            return mICHMesh;
        }

        static SimpleIdMesh<T> MeshFromCellsAndEdges(  
            IEnumerable<MICHDelaunayCell<T>> delaCells, 
            IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> delaEdges, 
            int numberOfVoronois)
        {
            (VariableCell<T>[] vCells, Vertex[] arrVertices) = 
                CreateMeshLists(delaCells, delaEdges, numberOfVoronois);

            for(int i = 0; i < vCells.Length; ++i)
            {
                vCells[i].Init();
            }
            List<MeshCell<T>> cells = new List<MeshCell<T>>(vCells);
            List<Vertex> vertices = new List<Vertex>(arrVertices);
            SimpleIdMesh<T> mesh = new SimpleIdMesh<T>(cells, vertices);
            return mesh;
        }

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

        static (VariableCell<T>[] cells, Vertex[] vertices) CreateMeshLists(   
            IEnumerable<MICHDelaunayCell<T>> delaCells,
            IEnumerable<MIConvexHull.VoronoiEdge
                <MICHVertex<T>, MICHDelaunayCell<T>>> edges,
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
            VariableCell<T>[] cells = new VariableCell<T>[numberOfVoronois];

            foreach (MICHDelaunayCell<T> delaCell in delaCells)
            {
                delaCell.done = true;
                MICHDelaunayCell<T>[] neighbors = delaCell.Adjacency;
                VariableCell<T>[] voronoiCells = new VariableCell<T>[3];
                vertices[delaCell.ID] = delaCell.Circumcenter; 

                //Create Voronoi Cells
                //--------------------------------------------------------------
                for(int i = 0; i < 3; ++i)
                {
                    MICHVertex<T> vert = delaCell.Vertices[i];
                    //Create or pick cell
                    VariableCell<T> voronoiCell = null;
                    if (cells[vert.ID] == null)
                    {
                        voronoiCell = new VariableCell<T> { ID = vert.ID, Node = vert.Node};
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
                    MICHDelaunayCell<T> neighbor = neighbors[i];
                    if (neighbor != null)
                    {
                        if (!neighbor.done)
                        {
                            (int k, int j) = OpposingIndices(i);
                            //Create Ridges
                            Edge<T> ridgeOutwards = new Edge<T>
                            {
                                Start = delaCell.Circumcenter,
                                End = neighbor.Circumcenter,
                                Cell = voronoiCells[k]
                            };
                            Edge<T> ridgeInwards = new Edge<T>
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
}
