using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    public class ArrayEnumerator<T> : IEnumerator<T>
    {
        int pointer;

        readonly IList<T> array;

        public ArrayEnumerator(IList<T> array)
        {
            this.array = array;
            pointer = -1;
        }

        public T Current => array[pointer];

        object IEnumerator.Current => Current;

        public void Dispose()
        {
        }

        public bool MoveNext()
        {
            return (++pointer < array.Count);
        }

        public void Reset()
        {
            pointer = -1;
        }
    }
}
