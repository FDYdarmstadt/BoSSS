using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Grid.Voronoi.Meshing.DataStructures
{
    class ConvolutionEnumerator<T> : IEnumerator<Pair<T>>
    {
        bool moveNext;

        T first;

        T previous;

        T current;

        IEnumerator<T> enumerator;

        public ConvolutionEnumerator(IEnumerator<T> enumerator)
        {
            this.enumerator = enumerator;
            Initialize();
        }

        Pair<T> IEnumerator<Pair<T>>.Current => currentPair;

        Pair<T> currentPair {
            get 
            {
                return new Pair<T>
                {
                    Current = current,
                    Previous = previous,
                };
            }
        }

        void Initialize()
        {
            moveNext = true;
            if (this.enumerator.MoveNext())
            {
                current = this.enumerator.Current;
                first = this.enumerator.Current;
            }
            else
            {
                moveNext = false;
            }
        }

        object IEnumerator.Current => currentPair;

        public void Dispose()
        {
            enumerator.Dispose();
        }

        public bool MoveNext()
        {
            if (moveNext)
            {
                previous = current;
                if (enumerator.MoveNext())
                {
                    current = enumerator.Current;
                }
                else
                {
                    moveNext = false;
                    current = first; 
                }
                return true;
            }
            else
            {
                return false;
            };
        }

        public void Reset()
        {
            enumerator.Reset();
            Initialize();
        }
    }
}
