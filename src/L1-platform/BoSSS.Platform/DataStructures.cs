using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Platform
{
    public class CyclicArray<T> : IList<T>
    {
        protected readonly T[] data;

        public readonly int Length;

        protected int start = 0;

        public CyclicArray(int Length)
        {
            data = new T[Length];
            this.Length = Length;
        }

        public CyclicArray(T[] array)
        {
            data = array;
            Length = array.Length;

        }

        public CyclicArray(T[] array, int start)
            : this(array)
        {
            this.start = start; 
        }

        public int Start {
            get {
                return start;
            }
            set {
                this.start = value % Length;
            }
        }

        public int Count => Length;

        public bool IsReadOnly => false;

        public T this[int i] {
            get {
                if (i >= Length || i < 0)
                {
                    throw new IndexOutOfRangeException();
                }
                else
                {
                    int index = goesRoundLikeARecordBaby(i);
                    return data[index];
                }
            }
            set {
                if (i >= Length || i < 0)
                {
                    throw new IndexOutOfRangeException();
                }
                else
                {
                    int index = goesRoundLikeARecordBaby(i);
                    data[index] = value;
                }
            }

        }

        private int goesRoundLikeARecordBaby(int indice)
        {
            indice = (indice + start) % Length;
            indice = Math.Abs(indice);
            return indice;
        }

        public int IndexOf(T item)
        {
            throw new NotImplementedException();
        }

        public void Insert(int index, T item)
        {
            throw new NotImplementedException();
        }

        public void RemoveAt(int index)
        {
            throw new NotImplementedException();
        }

        public void Add(T item)
        {
            throw new NotImplementedException();
        }

        public void Clear()
        {
            throw new NotImplementedException();
        }

        public bool Contains(T item)
        {
            throw new NotImplementedException();
        }

        public void CopyTo(T[] array, int arrayIndex)
        {
            throw new NotImplementedException();
        }

        public bool Remove(T item)
        {
            throw new NotImplementedException();
        }

        public IEnumerator<T> GetEnumerator()
        {
            return new CyclicArrayEnumerator<T>(this);
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    public class CyclicArrayEnumerator<T> : IEnumerator<T>
    {
        int pointer;

        readonly CyclicArray<T> array;

        public CyclicArrayEnumerator(CyclicArray<T> array)
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
            return (++pointer < array.Length);
        }

        public void Reset()
        {
            pointer = -1;
        }
    }

}
