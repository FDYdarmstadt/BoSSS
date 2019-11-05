using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using ilPSP;
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
            PeriodicMap map = null
            )
        {
            if (map != null)
            {
                ContainsPeriodicBoundaries = true;

                this.boundary = boundary;
                this.firstCellNodeIndice = firstCellNodeIndice;
                periodicTrafoMap = map.PeriodicTransformationMap;
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
            foreach (Edge<T> edge in BoundaryEdgesOf(mesh))
            {
                if (periodicTrafoMap.ContainsKey(edge.BoundaryEdgeNumber))
                {
                    yield return edge;
                }
            }
        }

        IEnumerable<Edge<T>> BoundaryEdgesOf(Mesh<T> mesh)
        {
            BoundaryEdgeFinder<T> edgeCells = new BoundaryEdgeFinder<T>(mesh);
            Vector startFromEnclosingCell = boundary[0].Start.Position;
            foreach (Edge<T> edge in edgeCells.Edges(startFromEnclosingCell, firstCellNodeIndice))
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
        }
    }
}
