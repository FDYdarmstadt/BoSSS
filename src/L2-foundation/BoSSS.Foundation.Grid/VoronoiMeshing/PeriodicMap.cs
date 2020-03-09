using BoSSS.Foundation.Grid.Classic;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class PeriodicMap
    {
        public IDictionary<int, int> PeriodicBoundaryCorrelation;

        public IDictionary<int, Transformation> PeriodicBoundaryTransformations;

        public IDictionary<Corner, int> PeriodicCornerCorrelation;
    }
}
