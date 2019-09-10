using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    static class MICMesher<T>
    {
        static readonly double accuracy = 1e-20;

        public static IDMesh<T> Create(IList<MICHVertex<T>> startNodes)
        {
            var mICHMesh = CreateVoronoiMeshFromDLL(startNodes);
            var delaunayCells = mICHMesh.Vertices;
            var delaunayEdges = mICHMesh.Edges;

            IDMesh<T> mesh = MeshFromCellsAndEdges(delaunayCells, delaunayEdges, startNodes.Count);
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

        static IDMesh<T> MeshFromCellsAndEdges(
            IEnumerable<MICHDelaunayCell<T>> delaCells,
            IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> delaEdges,
            int numberOfVoronois)
        {
            (VariableCell<T>[] vCells, Vertex[] arrVertices) =
                CreateMeshLists(delaCells, delaEdges, numberOfVoronois);

            for (int i = 0; i < vCells.Length; ++i)
            {
                vCells[i].Init();
            }
            List<MeshCell<T>> cells = new List<MeshCell<T>>(vCells);
            List<Vertex> vertices = new List<Vertex>(arrVertices);
            IDMesh<T> mesh = new IDMesh<T>
            {
                Cells = cells,
                Vertices = vertices
            };
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
            foreach (var edge in edges)
            {
                Vertex start = edge.Source.Circumcenter;
                Vertex end = edge.Target.Circumcenter;
                double length = (start.Position - end.Position).AbsSquare();
                if (length < accuracy)
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
                for (int i = 0; i < 3; ++i)
                {
                    MICHVertex<T> vert = delaCell.Vertices[i];
                    //Create or pick cell
                    VariableCell<T> voronoiCell = null;
                    if (cells[vert.ID] == null)
                    {
                        voronoiCell = new VariableCell<T> { ID = vert.ID, Node = vert.Node };
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
