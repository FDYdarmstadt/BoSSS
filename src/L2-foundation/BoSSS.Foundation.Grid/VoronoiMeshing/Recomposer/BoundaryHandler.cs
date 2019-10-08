using BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Recomposer
{
    class BoundaryHandler<T>
        where T : ILocatable, new()
    {
        readonly IDictionary<int, BoundaryTransformation> periodicTrafoMap;

        readonly BoundaryLine[] boundary;

        readonly int firstCellNodeIndice;

        readonly BoundaryNodeCloner<T> nodeCloner;

        readonly BoundaryRecomposer<T> recomposer;

        public bool ContainsPeriodicBoundaries { get; private set; }

        public BoundaryHandler(
            BoundaryLine[] boundary,
            IDictionary<int, int> periodicBoundaryMap,
            int firstCellNodeIndice)
        {
            if (periodicBoundaryMap != null)
            {
                this.boundary = boundary;
                this.firstCellNodeIndice = firstCellNodeIndice;
                ContainsPeriodicBoundaries = true;
                periodicTrafoMap = CreatePeriodicTransformationMapFrom(boundary, periodicBoundaryMap);
                nodeCloner = new BoundaryNodeCloner<T>(periodicTrafoMap);
                recomposer = new BoundaryRecomposer<T>(periodicBoundaryMap, periodicTrafoMap, firstCellNodeIndice); 
            }
            else
            {
                ContainsPeriodicBoundaries = false;
            }
        }

        static IDictionary<int, BoundaryTransformation> CreatePeriodicTransformationMapFrom(
            BoundaryLine[] boundary, 
            IDictionary<int, int> periodicBoundaryMap)
        {
            IDictionary<int, BoundaryTransformation> periodicTrafoMap = new LinkedListDictionary< int, BoundaryTransformation>();
            foreach(var boundaryPair in periodicBoundaryMap)
            {
                BoundaryLine source = boundary[boundaryPair.Key];
                BoundaryLine target = BoundaryLine.GetReverse(boundary[boundaryPair.Value]);
                
                BoundaryTransformation transformation = new BoundaryTransformation(source, target);
                periodicTrafoMap.Add(boundaryPair.Key, transformation);
            }
            return periodicTrafoMap;
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

        public void RecomposePeriodicEdges(Mesh<T> mesh)
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);

            IEnumerable<Edge<T>> periodicEdges = PeriodicEdgesOf(mesh);
            recomposer.RecomposePeriodicEdges(mesh, periodicEdges); 
        }
    }
}
