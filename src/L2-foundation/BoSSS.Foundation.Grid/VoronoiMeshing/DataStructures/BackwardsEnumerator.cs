using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    class BackwardsEnumerator<T> : IEnumerator<T>
    {
        int pointer; 

        readonly IList<T> array;

        public BackwardsEnumerator(IList<T> array)
        {
            this.array = array;
        }

        public T Current => array[pointer];

        object IEnumerator.Current => Current;

        public void Dispose()
        {
        }

        public bool MoveNext()
        {
            return (--pointer > -1);
        }

        public void Reset()
        {
            pointer = array.Count + 1;
        }
    }
}
