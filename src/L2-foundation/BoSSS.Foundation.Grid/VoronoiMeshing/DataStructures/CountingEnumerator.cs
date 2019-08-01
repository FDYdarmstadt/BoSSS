using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing
{
    /// <summary>
    /// Wrapper arround IEnumerator object that counts the enumerated objects
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class CountingEnumerator<T> : IEnumerator<T>
    {
        IEnumerator<T> enumerator;

        T current;

        public CountingEnumerator(IEnumerator<T> Enumerator)
        {
            Counter = -1;
            this.enumerator = Enumerator;
        }

        public int Counter { get; private set; }

        public T Current => current;

        object IEnumerator.Current => current;

        public void Dispose()
        {
            enumerator.Dispose();
        }

        public bool MoveNext()
        {
            ++Counter;
            bool moved = enumerator.MoveNext();
            if (moved)
            {
                current = enumerator.Current;
            }
            return moved;
        }

        public void Reset()
        {
            Counter = -1;
            enumerator.Reset();
        }
    }
}
