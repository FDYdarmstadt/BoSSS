using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Platform
{
    public class CyclicArray<T>
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

        public T this[int i] {
            get {
                if (i >= Length)
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
                if (i >= Length)
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
    }
}
