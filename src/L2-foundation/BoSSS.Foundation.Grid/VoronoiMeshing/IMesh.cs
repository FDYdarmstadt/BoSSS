using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    interface IMesh<T>
    {
        T CornerNode { get; }

        int CornerNodeIndice { get; }

        List<MeshCell<T>> Cells { get; }

        List<T> Nodes { get; }
    }
}
