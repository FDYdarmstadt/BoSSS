using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Voronoi;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryHandler
    {
        bool boundaryHasPeriodicEdges; 

        public BoundaryHandler(Vector[] boundary, ArrayMap periodicEdgeMapping)
        {
            
            if (periodicEdgeMapping != null)
            {
                boundaryHasPeriodicEdges = true;
                CheckIfValidInput(boundary, periodicEdgeMapping);
            }
            else
            {
                boundaryHasPeriodicEdges = false;
            }
        }

        static void CheckIfValidInput(Vector[] boundary, ArrayMap periodicEdgeMapping)
        {
            if(!periodicEdgeMapping.IsBijective())
            {
                throw new Exception();
            }
            if(boundary.Length - 1 > periodicEdgeMapping.MaxIndice())
            {
                throw new Exception();
            }
        }

        public IntersectionMesh<T> ImposePeriodicity<T>(IntersectionMesh<T> mesh, IList<T> nodeList) 
            where T : IMesherNode, new()
        {
            if (boundaryHasPeriodicEdges)
            {
                //Mirror : Follow boundary 2 times and
                // 1) Collect nodes that should be copied
                // 2) Copy

                //Cut

                //Merge/Shift: Follow boundary 2 times and
                // 1) Collect
                // 2) Remove & Add
            }
            return mesh;
        }
    }
}
