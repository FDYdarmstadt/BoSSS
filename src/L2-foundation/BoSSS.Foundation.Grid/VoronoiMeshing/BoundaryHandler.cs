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

        Line[] boundary;

        Map periodicEdgeMapping;

        public BoundaryHandler(Line[] boundary, Map periodicEdgeMapping)
        {
            if (periodicEdgeMapping != null)
            {
                boundaryHasPeriodicEdges = true;
                this.boundary = boundary;
                this.periodicEdgeMapping = periodicEdgeMapping;
                CheckIfValid();
            }
            else
            {
                boundaryHasPeriodicEdges = false;
            }
        }

        void CheckIfValid()
        {
            if(!periodicEdgeMapping.IsBijective())
            {
                throw new Exception();
            }
            if(periodicEdgeMapping.Max() > boundary.Length - 1 || periodicEdgeMapping.Min() < 0)
            {
                throw new Exception();
            }
            if ( !PeriodicBoundariesAreSameSize())
            {
                throw new Exception();
            }
        }

        bool PeriodicBoundariesAreSameSize()
        {
            for(int i = 0; i < periodicEdgeMapping.Length; ++i)
            {
                Line boundaryLineA = boundary[i];
                Line boundaryLineB = periodicEdgeMapping.GetCorrespondingEntry(i, boundary);

                double lengthA = boundaryLineA.Length();
                double lengthB = boundaryLineB.Length();
                
                if( Math.Abs(lengthA - lengthB) > 1e-14 )
                {
                    return false;
                }
            }
            return true;
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
