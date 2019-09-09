using BoSSS.Platform.LinAlg;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Voronoi;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BoundaryHandler<T>
        where T : IMesherNode, new()
    {
        public bool boundaryHasPeriodicEdges { get; private set; }

        BoundaryLine[] boundary;

        Map periodicEdgeMapping;

        public BoundaryHandler(BoundaryLine[] boundary, Map periodicEdgeMapping)
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
                BoundaryLine boundaryLineA = boundary[i];
                BoundaryLine boundaryLineB = periodicEdgeMapping.GetCorrespondingEntry(i, boundary);

                double lengthA = boundaryLineA.Length();
                double lengthB = boundaryLineB.Length();
                
                if( Math.Abs(lengthA - lengthB) > 1e-14 )
                {
                    return false;
                }
            }
            return true;
        }

        public void MirrorNodes(MeshIntersecter<T> mesh, IList<T> nodeList) 
        {
            //Mirror : Follow boundary and
            // 1) Collect nodes that should be mirrored
            /*while (boundaryEdges.MoveNext())
            {
                Edge<T> edge = boundaryEdges.Current;
            }
            for ()
            {

            }
                
            // 2) Mirror

            //Cut

            //Merge/Shift: Follow boundary 2 times and
            // 1) Collect
            // 2) Remove & Add
            */
        }

        public void UnifyMirroredCells(IList<MeshCell<T>> cells)
        {

        }

    }

    class Collecter<T>
    {
        IEnumerator<Edge<T>> boundaryEdges;

    }

    class Mirror<T>
    {

    }
}
