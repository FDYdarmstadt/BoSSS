using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.PeriodicBoundaryHandler
{
    class PeriodicBoundaryHandler<T>
        where T : ILocatable, new()
    {
        readonly IDictionary<int, Transformation> periodicTrafoMap;

        readonly BoundaryLine[] boundary;

        readonly int firstCellNodeIndice;

        readonly BoundaryNodeCloner<T> nodeCloner;

        readonly BoundaryRecomposer<T> recomposer;

        public bool ContainsPeriodicBoundaries { get; private set; }

        public PeriodicBoundaryHandler(
            BoundaryLine[] boundary,
            int firstCellNodeIndice,
            PeriodicMap map = null)
        {
            if (map != null)
            {
                ContainsPeriodicBoundaries = true;

                this.boundary = boundary;
                this.firstCellNodeIndice = firstCellNodeIndice;
                periodicTrafoMap = map.PeriodicBoundaryTransformations;
                nodeCloner = new BoundaryNodeCloner<T>(map);
                recomposer = new BoundaryRecomposer<T>(map, firstCellNodeIndice); 
            }
            else
            {
                ContainsPeriodicBoundaries = false;
            }
        }

        IEnumerable<Edge<T>> PeriodicEdgesOf(Mesh<T> mesh)
        {
            Vector startFromEnclosingCell = boundary[0].Start.Position;
            return PeriodicEdgesOf(mesh, startFromEnclosingCell);
        }

        IEnumerable<Edge<T>> PeriodicEdgesOf(Mesh<T> mesh, Vector start)
        {
            foreach (Edge<T> edge in BoundaryEdgesOf(mesh, start))
            {
                if (periodicTrafoMap.ContainsKey(edge.BoundaryEdgeNumber))
                {
                    yield return edge;
                }
            }
        }

        IEnumerable<Edge<T>> BoundaryEdgesOf(Mesh<T> mesh, Vector start)
        {
            BoundaryEdgeFinder<T> edgeCells = new BoundaryEdgeFinder<T>(mesh);
            
            foreach (Edge<T> edge in edgeCells.Edges(start, firstCellNodeIndice))
            {
                yield return edge;
            }
        }

        public IList<T> CloneNodesAlongPeriodicBoundaries(Mesh<T> mesh)
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);

            IEnumerable<Edge<T>> periodicEdges = PeriodicEdgesOf(mesh);
            List<T> clones = nodeCloner.CloneAndMirrorNodesOf(periodicEdges);
            IList<T> nodes = mesh.Nodes;
            nodes.AddRange(clones);
            return nodes;
        }

        public void RecomposePeriodicEdges(IDMesh<T> mesh)
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);

            IEnumerable<Edge<T>> periodicEdges = PeriodicEdgesOf(mesh);
            recomposer.RecomposePeriodicEdges(mesh, periodicEdges);

            Debug.Assert(PeriodicEdgeIDsLineUp(mesh));
        }

        bool PeriodicEdgeIDsLineUp(Mesh<T> mesh)
        {
            IEnumerable<Edge<T>> periodicEdges = PeriodicEdgesOf(mesh, mesh.Cells[firstCellNodeIndice].Node.Position);
            List<Vertex> innerIds = new List<Vertex>();
            List<Vertex> outerIds = new List<Vertex>();
            foreach(var edge in periodicEdges)
            {
                innerIds.Add(edge.Start);
                outerIds.Add(edge.Twin.Start);
            }
            HashSet< int> innerEdges = new HashSet<int>();
            int missCounter = 0;
            int duplicateCounter = 0;
            foreach(Vertex i in innerIds)
            {
                if (innerEdges.Contains(i.ID))
                {
                    Console.WriteLine($"({i.Position.x}, {i.Position.y} ) is a duplicate.");
                    ++duplicateCounter;
                }
                innerEdges.Add(i.ID);
            }
            HashSet<int> outerEdges = new HashSet<int>();
            foreach (Vertex i in outerIds)
            {
                if (outerEdges.Contains(i.ID))
                {
                    Console.WriteLine($"({i.Position.x}, {i.Position.y} ) is a duplicate.");
                    ++duplicateCounter;
                }
                outerEdges.Add(i.ID);

                if (!innerEdges.Contains(i.ID))
                {
                    Console.WriteLine($"({i.Position.x}, {i.Position.y} ) is a miss.");
                    ++missCounter;
                }
            }
            foreach (Vertex i in innerIds)
            {
                if (!outerEdges.Contains(i.ID))
                {
                    Console.WriteLine($"({i.Position.x}, {i.Position.y} ) is a miss.");
                    ++missCounter;
                }
            }
            if (missCounter > 0 || duplicateCounter > 0)
                return false;
            else 
                return true;
        }
    }
}
