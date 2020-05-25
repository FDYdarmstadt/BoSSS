using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    class Convolution<T> : IEnumerable<Pair<T>>
    {
        readonly ConvolutionEnumerator<T> enumerator;

        public Convolution(IEnumerable<T> enumerable)
        {
            enumerator = new ConvolutionEnumerator<T>(enumerable.GetEnumerator()); 
        }

        public IEnumerator<Pair<T>> GetEnumerator()
        {
            return enumerator;
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return enumerator;
        }
    }
}
