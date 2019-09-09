using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class CutterState<TEdge>
    {
        public TEdge ActiveEdge;

        public IntersectionCase Case;

        public double AlphaCut;
    }
}
