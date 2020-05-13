using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    static class CountingEnumerable
    { 
        public static IEnumerable<(int i, T value)>Wrap<T>(IEnumerable<T> enumerable)
        {
            int i = 0; 
            foreach(T value in enumerable)
            {
                yield return (i, value);
                ++i;
            }
        } 

    }
}
