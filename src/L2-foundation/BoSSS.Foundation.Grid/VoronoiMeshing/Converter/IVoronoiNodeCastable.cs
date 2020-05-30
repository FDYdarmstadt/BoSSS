using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.Converter
{
    public interface IVoronoiNodeCastable
    {
        VoronoiNode AsVoronoiNode();
    }
}
