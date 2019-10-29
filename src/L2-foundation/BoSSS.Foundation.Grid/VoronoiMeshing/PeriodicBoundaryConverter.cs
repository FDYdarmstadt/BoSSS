using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform.LinAlg;
using System;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class PeriodicBoundaryConverter
    {
        AffineTrafo[] PeriodicTrafos;

        public PeriodicBoundaryConverter()
        {

        }

        public void RegisterPeriodicBoundariesTo(GridCommons grid)
        {
            byte tag = grid.AddPeriodicEdgeTrafo(PeriodicTrafos[0]);
        }

        public bool IsPeriodicInverse(int boundaryEdgeNumber)
        {
            throw new NotImplementedException();
        }
    }
}
