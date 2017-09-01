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
using System.Diagnostics;
using System.Text;

namespace ilPSP.LinSolvers.monkey {
    
    /// <summary>
    /// baseclass for distribuetd vectors.
    /// </summary>
    public abstract partial class VectorBase : LockAbleObject, IList<double> {

        /// <summary>
        /// .
        /// </summary>
        /// <param name="P">
        /// <see cref="Part"/>
        /// </param>
        protected VectorBase(IPartitioning P) {
            if (P.IsMutable) {
                throw new NotSupportedException(P.GetType().Name + " is marked as mutable partitioning.");
            }
            m_Part = P;
        }

        /// <summary>
        /// <see cref="Part"/>
        /// </summary>
        protected IPartitioning m_Part;

        /// <summary>
        /// distribution of this vector among MPI processes
        /// </summary>
        public IPartitioning Part {
            get { return m_Part; }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vals"></param>
        /// <param name="arrayIndex">
        /// The zero-based index into <paramref name="vals"/> at which copying (reading) begins.
        /// </param>
        /// <param name="insertAt">
        /// The zero-based index into tihs vector at which copying (writing to this array) begins.
        /// </param>
        /// <param name="Length">
        /// No of elemtns to copy;
        /// </param>
        abstract public void SetValues<T>(T vals, int arrayIndex, int insertAt, int Length) 
            where T : IList<double>;

        /// <summary>
        /// Copies from this array to <paramref name="vals"/>
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vals"></param>
        /// <param name="arrayIndex">
        /// The zero-based index into <paramref name="vals"/> in array at which copying (writeing) begins.
        /// </param>
        /// <param name="readAt">
        /// The zero-based index into tihs vector at which copying (reading from this array) begins.
        /// </param>
        /// <param name="Length">
        /// No of elemtns to copy;
        /// </param>
        abstract public void GetValues<T>(T vals, int arrayIndex, int readAt, int Length)
            where T : IList<double>;

        /// <summary>
        /// not supported
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public int IndexOf(double item) {
            throw new NotSupportedException("not appopriate.");
        }

        /// <summary>
        /// not supported
        /// </summary>
        /// <param name="index"></param>
        /// <param name="item"></param>
        public void Insert(int index, double item) {
            throw new NotSupportedException("resizeing is not supported.");
        }

        /// <summary>
        /// not supported
        /// </summary>
        /// <param name="index"></param>
        public void RemoveAt(int index) {
            throw new NotSupportedException("resizeing is not supported.");
        }

        /// <summary>
        /// gets/sets an entry of this vector
        /// </summary>
        /// <param name="index">a local index (i.e. only valid on the current MPI processor)</param>
        /// <returns></returns>
        abstract public double this[int index] { get; set; }

        /// <summary>
        /// not supported;
        /// </summary>
        /// <param name="item"></param>
        public void Add(double item) {
            throw new NotSupportedException("resizeing is not supported.");
        }

        /// <summary>
        /// sets all entries to 0.0; works in locke and unlocked mode;
        /// </summary>
        abstract public void Clear();

        /// <summary>
        /// not supported
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Contains(double item) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="array"></param>
        /// <param name="arrayIndex"></param>
        public void CopyTo(double[] array, int arrayIndex) {
            GetValues(array, arrayIndex, 0, this.Count);
        }

        /// <summary>
        /// number of items on the current processor
        /// </summary>
        virtual public int Count {
            get { return m_Part.LocalLength; }
        }

        /// <summary>
        /// allwasy false;
        /// </summary>
        public bool IsReadOnly {
            get { return false; }
        }

        /// <summary>
        /// not supported
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Remove(double item) {
            throw new NotSupportedException("resizeing is not supported.");
        }
        

        #region IEnumerable<double> Members

        class MyEnum : IEnumerator<double> {

            VectorBase m_owner;

            public MyEnum(VectorBase owner) {
                m_owner = owner;
            }

            #region IEnumerator<double> Members

            public double Current {
                get { return m_owner[ind]; }
            }

            #endregion

            #region IDisposable Members

            public void Dispose() {
                // nop
            }

            #endregion

            #region IEnumerator Members

            object System.Collections.IEnumerator.Current {
                get { return Current; }
            }

            int ind = -1;

            public bool MoveNext() {
                ind++;
                return (ind < m_owner.Count);
            }

            public void Reset() {
                ind = -1;
            }

            #endregion
        }

        /// <summary>
        /// like defined in interface
        /// </summary>
        /// <returns></returns>
        public IEnumerator<double> GetEnumerator() {
            return new MyEnum(this);
        }

        #endregion

        #region IEnumerable Members

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }

        #endregion

        /// <summary>
        /// multiplies this vector by <paramref name="alpha"/>. Works in locked mode only;
        /// </summary>
        /// <param name="alpha"></param>
        public abstract void Scale(double alpha);


        /// <summary>
        /// this = this + <paramref name="alpha"/>*<paramref name="other"/>;
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="other"></param>
        abstract public void Acc(double alpha, VectorBase other);


        /// <summary>
        /// For each <em>j</em>, <br/>
        /// this[j] = this[j]*<paramref name="other"/>[j]
        /// </summary>
        /// <param name="other"></param>
        abstract public void MultiplyElementWise(VectorBase other);
        
        ///// <summary>
        ///// this = this*<paramref name="beta"/> + <paramref name="other"/>
        ///// </summary>
        ///// <param name="alpha"></param>
        ///// <param name="other"></param>
        ///// <param name="beta"></param>
        //abstract public void daxpy(VectorBase other, double beta);

        /// <summary>
        /// the square of the two-norm of this vector;
        /// </summary>
        /// <returns></returns>
        /// <remarks>
        /// this is a MPI-collective operation;
        /// works only in locked mode;
        /// </remarks>
        abstract public double TwoNormSquare();

        /// <summary>
        /// computes the inner product of this and another vector;
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        /// <remarks>
        /// this is a MPI-collective operation;
        /// works only in locked mode
        /// </remarks>
        abstract public double InnerProd(VectorBase other);

        /// <summary>
        /// initilizes this vector to be a copy of <paramref name="other"/>
        /// </summary>
        /// <remarks>
        /// works only in locked mode
        /// </remarks>
        abstract public void CopyFrom(VectorBase other);

        /// <summary>
        /// swaps the content of this vector with <paramref name="other"/>
        /// </summary>
        /// <remarks>
        /// works only in locked mode
        /// </remarks>
        abstract public void Swap(VectorBase other);


        /// <summary>
        /// 
        /// </summary>
        /// <param name="_src">
        /// source vector to copy data from
        /// </param>
        /// <param name="IdxThis">
        /// indices into this vector, where to copy to
        /// </param>
        /// <param name="PerThis">
        /// must be a divider of the local length of this vector
        /// </param>
        /// <param name="IdxSrc">
        /// indices into vector <paramref name="_src"/>, where to copy from;
        /// length of this array must be equal to length of <paramref name="IdxThis"/>
        /// </param>
        /// <param name="PerSrc">
        /// must be a divider of the local length of <paramref name="_src"/>
        /// </param>
        public abstract void CopyPartlyFrom(VectorBase _src, int[] IdxThis, int PerThis, int[] IdxSrc, int PerSrc);

    }
}
