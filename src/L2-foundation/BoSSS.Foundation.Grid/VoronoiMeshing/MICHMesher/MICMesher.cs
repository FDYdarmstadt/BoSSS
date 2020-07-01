using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.MICHMesher
{
    static class MICMesher<T>
    {
        static readonly double accuracy = 1e-6;
        static readonly double zeroAccuracy = 1e-7;

        public static IDMesh<T> Create(IList<MICHVertex<T>> startNodes)
        {
            MIConvexHull.VoronoiMesh<MICHVertex<T>, MICHDelaunayCell<T>, MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> mICHMesh = CreateDealaunayTriangulation(startNodes);
            IEnumerable<MICHDelaunayCell<T>> delaunayVertices = mICHMesh.Vertices;
            IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> delaunayEdges = mICHMesh.Edges;

            IDMesh<T> mesh = MeshFromCellsAndEdges(delaunayVertices, delaunayEdges, startNodes.Count);
            return mesh;
        }

        static MIConvexHull.VoronoiMesh<MICHVertex<T>, MICHDelaunayCell<T>,
            MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>>
            CreateDealaunayTriangulation(IList<MICHVertex<T>> nodes)
        {
            var mICHMesh = MIConvexHull.VoronoiMesh.Create<
                    MICHVertex<T>,
                    MICHDelaunayCell<T>,
                    MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>
                >(nodes);

            return mICHMesh;
        }

        static IDMesh<T> MeshFromCellsAndEdges(
            IEnumerable<MICHDelaunayCell<T>> delaunayCells,
            IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> delaunayEdges,
            int numberOfVoronois)
        {
            (VariableCell<T>[] vCells, Vertex[] arrVertices) =
                CreateMeshLists(delaunayCells, delaunayEdges, numberOfVoronois);

            for (int i = 0; i < vCells.Length; ++i)
            {
                vCells[i].Init();
            }
            Debug.Assert(VerticesMatch(vCells));

            List<MeshCell<T>> cells = new List<MeshCell<T>>(vCells);
            List<Vertex> vertices = new List<Vertex>(arrVertices);
            IDMesh<T> mesh = new IDMesh<T>
            {
                Cells = cells,
                Vertices = vertices
            };
            return mesh;
        }

        static bool VerticesMatch(VariableCell<T>[] vCells)
        {
            foreach(VariableCell<T> cell in vCells)
            {
                foreach(Edge<T> edge in cell.Edges)
                {
                    if (!edge.Start.Equals(edge.Twin.End))
                    {
                        Console.WriteLine($"Found unmatching Vertices");
                    }
                    if (!edge.End.Equals(edge.Twin.Start))
                    {
                        Console.WriteLine($"Found unmatching Vertices");
                    }
                    if (!BoSSS.Platform.FloatingPointArithmetic.IsEqual(edge.Start.Position, edge.Twin.End.Position, accuracy, zeroAccuracy) )
                    {
                        Console.WriteLine($"Vertices of an edge between cell {cell.ID} " +
                            $"and cell {edge.Twin.Cell.ID} do not match geometrically");
                        return false;
                    };
                    if(BoSSS.Platform.FloatingPointArithmetic.IsEqual(edge.Start.Position, edge.End.Position, accuracy, zeroAccuracy))
                    {
                        Console.WriteLine("Detected zero edge");
                        return false;
                    }
                    if (edge.Start.ID == edge.End.ID)
                    {
                        Console.WriteLine("Detected zero edge");
                        return false;
                    }
                }
                if (!cell.boundary)
                {
                    Edge<T>[] edges = cell.Edges;
                    for (int i = 0; i < edges.Length; ++i)
                    {
                        if (!BoSSS.Platform.FloatingPointArithmetic.IsEqual(edges[i].End.Position, edges[(i + 1) % edges.Length].Start.Position, accuracy, zeroAccuracy))
                        {
                            Console.WriteLine($"Vertices of an edge of cell {cell.ID} do not match geometrically");
                            return false;
                        };
                        
                        if (edges[i].End.ID != edges[ (i+ 1) % edges.Length].Start.ID)
                        {
                            Console.WriteLine($"Detected error in Vertex Ids of edge in cell {cell.ID}");
                            return false;
                        }
                    }
                }
            }
            return true;
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
            IEnumerable<MICHDelaunayCell<T>> delaunayCells,
            IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> edges,
            int numberOfVoronois)
        {
            //Merge Position and ID of Vertices of zero-Edges
            MergeVerticesOfZeroEdges(edges);

            Vertex[] vertices = new Vertex[delaunayCells.Count()];
            VariableCell<T>[] cells = new VariableCell<T>[numberOfVoronois];

            foreach (MICHDelaunayCell<T> delaCell in delaunayCells)
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
                            if(delaCell.Circumcenter.ID != neighbor.Circumcenter.ID)
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

        static void MergeVerticesOfZeroEdges(IEnumerable<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> edges)
        {
            Dictionary<int, LinkedList<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>>> mergedVertices 
                = new Dictionary<int, LinkedList<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>>>();
            foreach (var edge in edges)
            {
                Vertex start = edge.Source.Circumcenter;
                Vertex end = edge.Target.Circumcenter;

                if (BoSSS.Platform.FloatingPointArithmetic.IsEqual(start.Position, end.Position, accuracy, zeroAccuracy))
                {
                    if (!mergedVertices.TryGetValue(start.ID, out LinkedList<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> mergeMe))
                    {
                        mergeMe = new LinkedList<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>>();
                        mergedVertices.Add(start.ID, mergeMe);
                    }
                    mergeMe.AddLast(edge);

                    if (mergedVertices.TryGetValue(end.ID, out LinkedList<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> mergedEnd))
                    {
                        mergeMe.AddRange(mergedEnd);
                    }
                    
                    MergeVertices(mergeMe, start.ID);
                }
            }
        }

        static void MergeVertices(ICollection<MIConvexHull.VoronoiEdge<MICHVertex<T>, MICHDelaunayCell<T>>> edges, int ID) 
        {
            Vector mergedPosition = new Vector(0,0); 
            foreach(var edge in edges)
            {
                mergedPosition += (edge.Source.Circumcenter.Position + edge.Target.Circumcenter.Position) / 2;
            };
            mergedPosition /= edges.Count;

            Vertex merged = edges.First().Source.Circumcenter;
            merged.Position = mergedPosition;
            merged.ID = ID;

            foreach (var edge in edges)
            {
                edge.Source.Circumcenter = merged;
                edge.Target.Circumcenter = merged;
            };
        }
    }
}
