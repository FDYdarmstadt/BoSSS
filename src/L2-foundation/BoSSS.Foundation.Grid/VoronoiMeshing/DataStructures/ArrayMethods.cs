using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    class ArrayMethods
    {
        public static S[] GetCopyInReverseOrder<S>(IList<S> list)
        {
            S[] reverse = new S[list.Count];
            for (int i = 0, j = list.Count - 1; i < list.Count; ++i, --j)
            {
                reverse[j] = list[i];
            }
            return reverse;
        }
    }
}
