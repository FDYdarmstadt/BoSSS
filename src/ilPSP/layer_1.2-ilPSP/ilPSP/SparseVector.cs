/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace ilPSP {

    /// <summary>
    /// a sparse vector (for some arbitrary type); looks like a dense vector, i.e. implements <see cref="IList{T}"/>;
    /// If one wants to access only nonzeros, use <see cref="SparseStruct"/>.
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public interface ISparseVector<T> : IList<T> {

        /// <summary>
        /// returns the number of nonzeros
        /// </summary>
        int NonZeros { get; }

        /// <summary>
        /// <see cref="NonZeros"/> over vector length (Count)
        /// </summary>
        double Sparsity { get; }

        /// <summary>
        /// gets a dictionary of all nonzero-entries (the keys are the indices);
        /// Manipulation of the returned object affects the sparse vector itself.
        /// </summary>
        IDictionary<int, T> SparseStruct { get; }
    }
    
    
    /// <summary>
    /// class to store sparse vectors, e.g. for some right-hand-sides that are only populated in boundary cells, ...
    /// </summary>
    /// <typeparam name="T">
    /// </typeparam>
    /// <remarks>
    /// The length of the vector (i.e. <see cref="Count"/>) is set during the construction, and remains constant afterwards;
    /// For entries which are not set, the default value of <typeparamref name="T"/> is returned (usually 0 or null);
    /// </remarks>
    public class SparseVector<T> : ISparseVector<T> {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="Count">
        /// initial overall length of the vector
        /// </param>
        public SparseVector(int Count) {
            m_Count = Count;
        }

        /// <summary>
        /// ctor; initializes the sparse structure from a dense array (only entries unequal to the default
        /// value of <typeparamref name="T"/> will be stored);
        /// </summary>
        /// <param name="Full"></param>
        public SparseVector(T[] Full) : this(Full.Length) {
            for (int i = 0; i < Full.Length; i++) {
                if (!Full[i].Equals(default(T)))
                    m_Content.Add(i, Full[i]);
            }
        }

        /// <summary>
        /// converts this vector to a dense array;
        /// </summary>
        /// <returns></returns>
        public T[] ToArray() {
            T[] ret = new T[this.Count];
            foreach (KeyValuePair<int, T> kv in this.m_Content) {
                ret[kv.Key] = kv.Value;
            }
            return ret;
        }

        
        int m_Count;

        SortedList<int,T> m_Content = new SortedList<int,T>();
        
        #region IList<T> Members

        /// <summary>
        /// as defined by interface <see cref="IList{T}"/>
        /// </summary>
        public int IndexOf(T item) {
            int i = m_Content.Values.IndexOf(item);
            return m_Content.Keys[i];
        }

        /// <summary>
        /// not supported: vector is not resizeable
        /// </summary>
        public void Insert(int index, T item) {
            throw new NotSupportedException("Count/Size of Vector is constant");
        }

        /// <summary>
        /// not supported: vector is not resizeable
        /// </summary>
        public void RemoveAt(int index) {
            throw new NotSupportedException("Count/Size of Vector is constant");
        }

        /// <summary>
        /// as defined by interface <see cref="IList{T}"/>
        /// </summary>
        public T this[int index] {
            get {
                if (index < 0 || index >= m_Count)
                    throw new IndexOutOfRangeException();
                T ret; 
                if(m_Content.TryGetValue(index,out ret)) {
                    return ret;
                } else {
                    return default(T);
                }
            }
            set {
                if (index < 0 || index >= m_Count)
                    throw new IndexOutOfRangeException();
                if (value.Equals(default(T))) {
                    m_Content.Remove(index);
                } else {
                    m_Content[index] = value;
                }
            }
        }

        #endregion

        #region ICollection<T> Members

        /// <summary>
        /// not supported: vector is not resizeable
        /// </summary>
        public void Add(T item) {
            throw new NotSupportedException("Count/Size of Vector is constant");
        }

        /// <summary>
        /// sets all entries to default value (usually 0 or 0.0 or null)
        /// </summary>
        public void Clear() {
            m_Content.Clear();
        }

        /// <summary>
        /// as defined by interface <see cref="IList{T}"/>
        /// </summary>
        public bool Contains(T item) {
            //return m_Content.Values.Contains<T>(item);
			throw new NotImplementedException();
        }

        /// <summary>
        /// copies content to array
        /// </summary>
        /// <param name="array">
        /// length of array must be at least <see cref="Count"/>+<paramref name="arrayIndex"/>
        /// </param>
        /// <param name="arrayIndex">
        /// </param>
        public void CopyTo(T[] array, int arrayIndex) {
            Array.Clear(array, arrayIndex, this.Count);
            foreach (var entry in m_Content) {
                array[entry.Key] = entry.Value;
            }
        }

        /// <summary>
        /// length of vector
        /// </summary>
        public int Count {
            get { return m_Count; }
        }

        /// <summary>
        /// false;
        /// </summary>
        public bool IsReadOnly {
            get { return false; }
        }

        /// <summary>
        /// not supported: vector is not resizeable
        /// </summary>
        public bool Remove(T item) {
            throw new NotSupportedException("Count/Size of Vector is constant");
        }

        #endregion

        #region IEnumerable<T> Members
        
        class MyEnum : IEnumerator<T>, System.Collections.IEnumerator {


            int m_Count;
            IEnumerator<KeyValuePair<int, T>> m_en;

            public MyEnum(int cnt, IEnumerator<KeyValuePair<int,T>> en) {
                m_Count = cnt;
                m_en = en;
            }

            int ind = -1;

            #region IEnumerator<T> Members

            public T Current {
                get {
                    if (ind < 0)
                        throw new InvalidOperationException();

                    if (!StillContVals)
                        return default(T);

                    if (ind == next.Key) {
                        return next.Value;
                    }

                    return default(T);
                }
            }

            #endregion

            #region IDisposable Members

            public void Dispose() {
                m_en.Dispose();
            }

            #endregion

            #region IEnumerator Members

            object System.Collections.IEnumerator.Current {
                get { return this.Current; }
            }

            KeyValuePair<int, T> next;
            bool StillContVals = true;
            bool next_worked = true;

            public bool MoveNext() {
                next_worked = (ind==next.Key);

                ind++;
                if (StillContVals) {
                    if (next_worked) {
                        StillContVals = m_en.MoveNext();

                        if (!StillContVals) {
                            next = m_en.Current;
                            next_worked = false;
                        }
                    }
                }

                return (ind < m_Count);
            }

            public void Reset() {
                StillContVals = true;
                next_worked = true;
                ind = -1;
                m_en.Reset();
            }

            #endregion
        }


        /// <summary>
        /// gets an enumerator
        /// </summary>
        public IEnumerator<T> GetEnumerator() {
            return new MyEnum(m_Count, this.m_Content.GetEnumerator());
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// gets an enumerator
        /// </summary>
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return new MyEnum(m_Count, this.m_Content.GetEnumerator());
        }

        #endregion

        #region ISparseVector<T> Members

        /// <summary>
        /// as defined by <see cref="ISparseVector{T}.NonZeros"/>
        /// </summary>
        public int NonZeros {
            get { return m_Content.Count; }
        }

        /// <summary>
        /// as defined by <see cref="ISparseVector{T}.Sparsity"/>
        /// </summary>
        public double Sparsity {
            get { return (double)this.NonZeros / (double)this.Count; }
        }

        /// <summary>
        /// as defined by <see cref="ISparseVector{T}.SparseStruct"/>
        /// </summary>
        public IDictionary<int, T> SparseStruct {
            get { return m_Content; }
        }

        #endregion
    }
}
