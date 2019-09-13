using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Voronoi;
using ilPSP;
using System.Diagnostics;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryHandler<T>
        where T : IMesherNode, new()
    {
        Dictionary<int, BoundaryTransformation> periodicTrafoMap;

        BoundaryLine[] boundary;

        int firstCellNodeIndice;

        public bool ContainsPeriodicBoundaries { get; private set; }

        public BoundaryHandler(
            BoundaryLine[] boundary,
            Dictionary<int, int> periodicBoundaryMap,
            int firstCellNodeIndice)
        {
            this.boundary = boundary;
            this.firstCellNodeIndice = firstCellNodeIndice;
            if (periodicBoundaryMap != null)
            {
                ContainsPeriodicBoundaries = true;
                periodicTrafoMap = CreateMappingFrom(boundary, periodicBoundaryMap);
            }
            else
            {
                ContainsPeriodicBoundaries = false;
            }
        }

        static Dictionary<int, BoundaryTransformation> CreateMappingFrom(BoundaryLine[] boundary, IDictionary<int, int> map)
        {
            Dictionary<int, BoundaryTransformation> periodicTrafoMap = new Dictionary<int, BoundaryTransformation>(map.Count);
            foreach(var boundaryPair in map)
            {
                BoundaryLine source = boundary[boundaryPair.Key];
                BoundaryLine target = boundary[boundaryPair.Value];
                BoundaryTransformation transformation = new BoundaryTransformation(source, target);
                periodicTrafoMap.Add(boundaryPair.Key, transformation);
            }
            return periodicTrafoMap;
        }

        public void CloneNodesAlongPeriodicBoundaries(Mesh<T> mesh, IList<T> nodes) 
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);
            List<T> clones = new List<T>();

            //Follow boundary and
            BoundaryEdgeEnumerator<T> edgeCells = new BoundaryEdgeEnumerator<T>(mesh);
            Vector start = (Vector)boundary[0].Start;
            foreach (Edge<T> edge in edgeCells.Edges(start, firstCellNodeIndice))
            {
                // 1) Collect nodes that should be mirrored
                // 2) Clone and teleport/transform to paired boundary
                if (periodicTrafoMap.TryGetValue(edge.BoundaryEdgeNumber, out BoundaryTransformation transformation))
                {
                    MeshCell<T> cell = edge.Cell;
                    T clone = new T
                    {
                        Position = transformation.Transform(cell.Node.Position)
                    };
                    clones.Add(clone); 
                }
            }
            nodes.AddRange(clones);
        }

        public void RecomposePeriodicEdges(Mesh<T> mesh)
        {
            Debug.Assert(ContainsPeriodicBoundaries == true);
            throw new NotImplementedException();
            //Merge/Shift: Follow boundary 2 times and
            // 1) Collect
            // 2) Remove & Add
        }
    }
}
