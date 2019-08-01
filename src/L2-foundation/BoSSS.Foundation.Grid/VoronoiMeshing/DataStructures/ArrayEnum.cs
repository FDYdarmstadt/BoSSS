using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    public class ArrayEnum<T> : IEnumerator<T>
    {
        int pointer;

        IList<T> arr;

        public ArrayEnum(IList<T> Arr)
        {
            arr = Arr;
            pointer = -1;
        }

        public T Current => arr[pointer];

        object IEnumerator.Current => Current;

        public void Dispose()
        {
        }

        public bool MoveNext()
        {
            return (++pointer < arr.Count);
        }

        public void Reset()
        {
            pointer = -1;
        }
    }
}
