using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Cutter
{
    class CutterState<TEdge>
    {
        public TEdge ActiveEdge;

        public IntersectionCase Case;

        public double AlphaCut;

        public BoundaryLine ActiveLine;

        public void Reset()
        {
            AlphaCut = default(double);
            ActiveEdge = default(TEdge);
            Case = IntersectionCase.NotIntersecting;
            ActiveLine = default(BoundaryLine);
        }
    }
}
