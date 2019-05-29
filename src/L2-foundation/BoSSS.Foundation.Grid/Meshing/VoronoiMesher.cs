using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public static class VoronoiMesher
    {
        public static VoronoiGrid CreateGrid(VoronoiNodes nodes, TrackableVoronoiMesher.Settings settings)
        {
            TrackableVoronoiGrid grid = TrackableVoronoiMesher.CreateGrid(nodes, settings);
            return grid.Result;
        }
    }

}
