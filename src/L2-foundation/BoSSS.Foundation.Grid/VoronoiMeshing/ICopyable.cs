using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public interface ICloneable<T> : ILocatable
    {
        T Clone();
    }
}
