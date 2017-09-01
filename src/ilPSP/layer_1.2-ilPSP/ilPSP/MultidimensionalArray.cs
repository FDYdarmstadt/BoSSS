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
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Serialization;
using ilPSP.Utils;

namespace ilPSP {

    /// <summary>
    /// The main difference to standard multidimensional arrays in .Net
    /// is the ability to extract sub-arrays
    /// (<see cref="ExtractSubArrayShallow(int[])"/>) or to reshape the array
    /// without reallocation of memory.
    /// </summary>
    /// <remarks>
    /// If this array is 2D (i.e. <see cref="Dimension"/>==2)
    /// it can be represented as a matrix <see cref="ilPSP.Utils.IMatrix"/>
    /// (see <see cref="NoOfRows"/>, see <see cref="NoOfCols"/>);
    /// </remarks>
    [Serializable]
    [DataContract]
    public partial class MultidimensionalArray : IMatrix, ICloneable, IEquatable<MultidimensionalArray> {

        /// <summary>
        /// The actual data.
        /// </summary>
        [DataMember]
        private double[] m_Storage;

        /// <summary>
        /// Indicates whether this instance is a shallow copy of another
        /// instance (i.e., whether <see cref="m_Storage"/> is shared)
        /// </summary>
        [DataMember]
        private bool m_ShallowCopy;

        /// <summary>
        /// Internal information about the layout of this array, i.e. the
        /// proper indexing of <see cref="m_Storage"/>
        /// </summary>
        [Serializable]
        [StructLayout(LayoutKind.Sequential)]
        private struct StorageLayout {

            /// <summary>
            /// Number of entries in the i-th dimension
            /// </summary>
            public int m_Length0;
            public int m_Length1;
            public int m_Length2;
            public int m_Length3;
            public int m_Length4;
            public int m_Length5;
            public int m_Length6;
            public int m_Length7;
            public int m_Length8;

            /// <summary>
            /// Length of the cycle corresponding to the i-th dimension. Since
            /// non-continuous sub-arrays can be extracted, this number may be
            /// larger than the corresponding length.
            /// </summary>
            public int m_Cycle0;
            public int m_Cycle1;
            public int m_Cycle2;
            public int m_Cycle3;
            public int m_Cycle4;
            public int m_Cycle5;
            public int m_Cycle6;
            public int m_Cycle7;
            public int m_Cycle8;
        }

        /// <summary>
        /// <see cref="StorageLayout"/>
        /// </summary>
        [DataMember]
        private StorageLayout m_StorageLayout;

        /// <summary>
        /// See <see cref="Dimension"/>
        /// </summary>
        [DataMember]
        private int m_Dimension;

        /// <summary>
        /// Indicates whether this instance may be edited
        /// </summary>
        [DataMember]
        private bool m_LockedForever;

        /// <summary>
        /// An optional offset into <see cref="m_Storage"/> of this instance is
        /// a shallow sub-array.
        /// </summary>
        [DataMember]
        private int m_Offset;

        /// <summary>
        /// Retrieves the value of the cycle (see
        /// <see cref="StorageLayout.m_Cycle0"/>,
        /// <see cref="StorageLayout.m_Cycle1"/>, ...) indexed by
        /// <paramref name="cycleIndex"/>
        /// </summary>
        /// <param name="cycleIndex">
        /// The index of the selected cycle.
        /// </param>
        /// <returns>
        /// The length of the cycle indexed by <paramref name="cycleIndex"/>.
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int GetCycle(int cycleIndex) {
            Debug.Assert(cycleIndex >= 0 && cycleIndex < m_Dimension);
            unsafe {
                fixed (int* pCycle0 = &this.m_StorageLayout.m_Cycle0) {
                    return *(pCycle0 + cycleIndex);
                }
            }
        }

        /// <summary>
        /// Alters the value of the cycle (see
        /// <see cref="StorageLayout.m_Cycle0"/>,
        /// <see cref="StorageLayout.m_Cycle1"/>, ...) indexed by
        /// <paramref name="cycleIndex"/>
        /// </summary>
        /// <param name="cycleIndex">
        /// The index of the selected cycle.
        /// </param>
        /// <param name="value">
        /// The value to be set.
        /// </param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void SetCycle(int cycleIndex, int value) {
            Debug.Assert(cycleIndex >= 0 && cycleIndex < m_Dimension);
            unsafe {
                fixed (int* pCycle0 = &this.m_StorageLayout.m_Cycle0) {
                    *(pCycle0 + cycleIndex) = value;
                }
            }
        }

        /// <summary>
        /// The length of the <paramref name="i"/>-th dimension.
        /// </summary>
        /// <param name="i">
        /// The selected dimension
        /// </param>
        /// <returns></returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int GetLength(int i) {
            Debug.Assert(i >= 0 && i < m_Dimension);
            unsafe {
                fixed (int* pCycle0 = &this.m_StorageLayout.m_Length0) {
                    return *(pCycle0 + i);
                }
            }
        }

        /// <summary>
        /// Alters the value of the cycle length (see
        /// <see cref="StorageLayout.m_Length0"/>,
        /// <see cref="StorageLayout.m_Length1"/>, ...) indexed by
        /// <paramref name="cycleIndex"/>
        /// </summary>
        /// <param name="cycleIndex">
        /// The index of the selected cycle length.
        /// </param>
        /// <param name="value">
        /// The value to be set.
        /// </param>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void SetLength(int cycleIndex, int value) {
            Debug.Assert(cycleIndex >= 0 && cycleIndex < m_Dimension);
            unsafe {
                fixed (int* pCycle0 = &this.m_StorageLayout.m_Length0) {
                    *(pCycle0 + cycleIndex) = value;
                }
            }
        }

        /// <summary>
        /// Dimension ("number of lengths") of this multidimensional array;
        /// </summary>
        public int Dimension {
            get {
                return m_Dimension;
            }
        }

        /// <summary>
        /// The lengths of all dimensions
        /// </summary>
        public int[] Lengths {
            get {
                var ret = new int[this.Dimension];
                for (int i = 0; i < ret.Length; i++)
                    ret[i] = this.GetLength(i);
                return ret;
            }
        }

        /// <summary>
        /// the array that stores all entries;<br/>
        /// Internally, the array is organized in C-Order, e. g. for a 3D-array 
        /// indexed by [i,j,k], k varies fastest and i varies slowest.
        /// </summary>
        public double[] Storage {
            get {
                return m_Storage;
            }
        }

        /// <summary>
        /// true if all entries of the array lie in a continuous index region of the 
        /// memory <see cref="Storage"/>.
        /// </summary>
        public bool IsContinious {
            get {
                int Dim = this.Dimension;
                int Cyc_d_cont = 1;
                for (int d = Dimension - 1; d >= 0; d--) {
                    int L_d = this.GetLength(d);
                    if (L_d > 1 && this.GetCycle(d) != Cyc_d_cont) {
                        return false;
                    }

                    Cyc_d_cont *= L_d;
                }
                return true;
            }
        }

        /// <summary>
        /// set/get one element; 2D - version;
        /// </summary>
        /// <param name="i"></param>
        /// <param name="j"></param>
        /// <returns></returns>
        public double this[int i, int j] {
            get {
                return m_Storage[Index(i, j)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 1D - version;
        /// </summary>
        public double this[int i] {
            get {
                return m_Storage[Index(i)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 3D - version;
        /// </summary>
        public double this[int i, int j, int k] {
            get {
                return m_Storage[Index(i, j, k)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j, k)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 4D - version;
        /// </summary>
        public double this[int i, int j, int k, int l] {
            get {
                return m_Storage[Index(i, j, k, l)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j, k, l)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 5D - version;
        /// </summary>
        public double this[int i, int j, int k, int l, int m] {
            get {
                return m_Storage[Index(i, j, k, l, m)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j, k, l, m)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 6D - version;
        /// </summary>
        public double this[int i, int j, int k, int l, int m, int n] {
            get {
                return m_Storage[Index(i, j, k, l, m, n)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j, k, l, m, n)] = value;
            }
        }

        /// <summary>
        /// set/get one element; 7D - version;
        /// </summary>
        public double this[int i, int j, int k, int l, int m, int n, int o] {
            get {
                return m_Storage[Index(i, j, k, l, m, n, o)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(i, j, k, l, m, n, o)] = value;
            }
        }

        /// <summary>
        /// set/get one element; version for general indices;
        /// </summary>
        /// <param name="indices"></param>
        /// <returns></returns>
        public double this[params int[] indices] {
            get {
                return m_Storage[Index(indices)];
            }
            set {
                Debug.Assert(!m_LockedForever, "illegal call - object is locked.");
                m_Storage[Index(indices)] = value;
            }
        }

        /// <summary>
        /// calculate the index of a 1D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i) {
            Debug.Assert(m_Dimension == 1, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "Index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            return ind;
        }

        /// <summary>
        /// calculate the index of a 2D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j) {
            Debug.Assert(m_Dimension == 2, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            return ind;
        }

        /// <summary>
        /// calculate the index of a 3D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j, int k) {
            Debug.Assert(m_Dimension == 3, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");
            Debug.Assert(k >= 0 && k < GetLength(2), "Third index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            ind += m_StorageLayout.m_Cycle2 * k;
            return ind;

        }

        /// <summary>
        /// calculate the index of a 4D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j, int k, int l) {
            Debug.Assert(m_Dimension == 4, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");
            Debug.Assert(k >= 0 && k < GetLength(2), "Third index out of range");
            Debug.Assert(l >= 0 && l < GetLength(3), "Fourth index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            ind += m_StorageLayout.m_Cycle2 * k;
            ind += m_StorageLayout.m_Cycle3 * l;
            return ind;
        }

        /// <summary>
        /// calculate the index of a 5D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j, int k, int l, int m) {
            Debug.Assert(m_Dimension == 5, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");
            Debug.Assert(k >= 0 && k < GetLength(2), "Third index out of range");
            Debug.Assert(l >= 0 && l < GetLength(3), "Fourth index out of range");
            Debug.Assert(m >= 0 && m < GetLength(4), "Fifth index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            ind += m_StorageLayout.m_Cycle2 * k;
            ind += m_StorageLayout.m_Cycle3 * l;
            ind += m_StorageLayout.m_Cycle4 * m;
            return ind;
        }

        /// <summary>
        /// calculate the index of a 6D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j, int k, int l, int m, int n) {
            Debug.Assert(m_Dimension == 6, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");
            Debug.Assert(k >= 0 && k < GetLength(2), "Third index out of range");
            Debug.Assert(l >= 0 && l < GetLength(3), "Fourth index out of range");
            Debug.Assert(m >= 0 && m < GetLength(4), "Fifth index out of range");
            Debug.Assert(n >= 0 && n < GetLength(5), "Sixth index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            ind += m_StorageLayout.m_Cycle2 * k;
            ind += m_StorageLayout.m_Cycle3 * l;
            ind += m_StorageLayout.m_Cycle4 * m;
            ind += m_StorageLayout.m_Cycle5 * n;
            return ind;
        }

        /// <summary>
        /// calculate the index of a 7D-array into <see cref="Storage"/>;
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(int i, int j, int k, int l, int m, int n, int o) {
            Debug.Assert(m_Dimension == 7, "Wrong array dimension");
            Debug.Assert(i >= 0 && i < GetLength(0), "First index out of range");
            Debug.Assert(j >= 0 && j < GetLength(1), "Second index out of range");
            Debug.Assert(k >= 0 && k < GetLength(2), "Third index out of range");
            Debug.Assert(l >= 0 && l < GetLength(3), "Fourth index out of range");
            Debug.Assert(m >= 0 && m < GetLength(4), "Fifth index out of range");
            Debug.Assert(n >= 0 && n < GetLength(5), "Sixth index out of range");
            Debug.Assert(o >= 0 && o < GetLength(6), "Seventh index out of range");

            int ind = m_Offset;
            ind += m_StorageLayout.m_Cycle0 * i;
            ind += m_StorageLayout.m_Cycle1 * j;
            ind += m_StorageLayout.m_Cycle2 * k;
            ind += m_StorageLayout.m_Cycle3 * l;
            ind += m_StorageLayout.m_Cycle4 * m;
            ind += m_StorageLayout.m_Cycle5 * n;
            ind += m_StorageLayout.m_Cycle6 * o;
            return ind;
        }

        /// <summary>
        /// Calculate the index of a nD-array into <see cref="Storage"/>;
        /// </summary>
        /// <param name="indices">
        /// Indices into each dimension of this instance.
        /// </param>
        /// <returns>
        /// The corresponding index into <see cref="Storage"/>
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Index(params int[] indices) {
            Debug.Assert(
                m_Dimension == indices.Length,
                "Wrong array dimension");
            Debug.Assert(
                indices.All(i => i >= 0),
                "Index out of range");
            Debug.Assert(indices.Where(
                (index, i) => index >= this.Lengths[i]).Count() == 0,
                "Index out of range");

            return IndexUnchecked(indices);
        }

        /// <summary>
        /// Calculates the index of a nD-array into <see cref="Storage"/>
        /// while omitting all kinds of checks on the
        /// <paramref name="indices"/>. This is required for the 
        /// extraction of shallow sub-arrays.
        /// </summary>
        /// <param name="indices">
        /// Indices into each dimension of this instance.
        /// </param>
        /// <returns>
        /// The corresponding index into <see cref="Storage"/>
        /// </returns>
        /// <remarks>
        /// Bjoern to himself: One does not simply outsmart Florian on these
        /// things...
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private int IndexUnchecked(int[] indices) {
            int ind = m_Offset;
            for (int i = this.m_Dimension - 1; i >= 0; i--) {
                ind += GetCycle(i) * indices[i];
            }
            return ind;
        }

        /// <summary>
        /// Calculates the index of a nD-array into <see cref="Storage"/>
        /// while omitting all kinds of checks on the
        /// <paramref name="indices"/>. This is required for the 
        /// extraction of shallow sub-arrays.
        /// </summary>
        /// <param name="indices">
        /// Indices into each dimension of this instance.
        /// </param>
        /// <returns>
        /// The corresponding index into <see cref="Storage"/>
        /// </returns>
        /// <remarks>
        /// Bjoern to himself: One does not simply outsmart Florian on these
        /// things...
        /// </remarks>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        unsafe private int IndexUnchecked(int* indices) {
            int ind = m_Offset;
            for(int i = this.m_Dimension - 1; i >= 0; i--) {
                ind += GetCycle(i) * indices[i];
            }
            return ind;
        }

        /// <summary>
        /// total number of elements
        /// </summary>
        public int Length {
            get {
                int l = 1;
                for (int i = 0; i < m_Dimension; i++) {
                    l *= GetLength(i);
                }
                return l;
            }
        }

        /// <summary>
        /// creates a resized, shallow copy of this array; 
        /// the copy uses the same storage as this array (shallow copy)
        /// so changes in the copy also change this array;
        /// this operation is relatively light-weighted;
        /// </summary>
        /// <param name="NewLengths">
        /// The new length of each individual dimension. The total number of
        /// entries must stay constant
        /// </param>
        /// <returns>
        /// A shallow, resized copy.
        /// </returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public MultidimensionalArray ResizeShallow(params int[] NewLengths) {
            if (this.IsContinious == false)
                throw new NotSupportedException("currently, resizing is only available for continuous arrays.");

            // check
            int newLen = 1;
            int dimNew = NewLengths.Length;
            for (int i = dimNew - 1; i >= 0; i--)
                newLen *= NewLengths[i];
            if (newLen != Length)
                throw new ArgumentOutOfRangeException("total number of elements must stay constant.");

            // create new array
            MultidimensionalArray ret = new MultidimensionalArray(NewLengths.Length);
            ret.m_Offset = this.m_Offset;
            ret.m_Storage = this.m_Storage;
            ret.m_ShallowCopy = true;
            for (int i = dimNew - 1; i >= 0; i--) {
                ret.SetLength(i, NewLengths[i]);
            }
            ret.CreateCycles();

            // return
            return ret;
        }

        /// <summary>
        /// extracts a sub-array from this array.
        /// The returned array is a shallow copy, i.e. changing its entries 
        /// will also change the corresponding entry in this array.
        /// </summary>
        /// <param name="Istart">start index (including)</param>
        /// <param name="Iend">end index (including)</param>
        /// <returns></returns>
        /// <remarks>
        /// For the i-th dimension, the range from <paramref name="Istart"/>[i]
        /// (including) to <paramref name="Iend"/>[i] (including) will be
        /// extracted. If end index is smaller than starting index, the
        /// dimension will be omitted in the returned array. If start and end
        /// index are equal, the dimension will be present in the returned
        /// array.
        /// </remarks>
        public MultidimensionalArray ExtractSubArrayShallow(int[] Istart, int[] Iend) {
            { //unsafe {
                int DIM = this.m_Dimension;

                // check Arguments
                if(Istart.Length != DIM)
                    throw new ArgumentOutOfRangeException("Mismatch between array dimension and length of 1st argument.");
                if(Iend.Length != DIM)
                    throw new ArgumentOutOfRangeException("Mismatch between array dimension and length of 1st argument.");

                bool isThis = true;
                for(int i = 0; i < DIM; i++) {
                    if(Istart[i] != 0) {
                        isThis = false;
                        break;
                    }
                    if(Iend[i] != GetLength(i) - 1) {
                        isThis = false;
                        break;
                    }
                }
                if(isThis)
                    return this;


                // dimension of SubArray
                int[] rStart = (int[])Istart.Clone(); // we want to overwrite some entries
                int[] rEnd = (int[])Iend.Clone();

                //int* rStart = stackalloc int[DIM];
                //int* rEnd = stackalloc int[DIM];
                

                int dims = 0;
                for(int i = 0; i < m_Dimension; i++) {
                    rStart[i] = Istart[i];
                    rEnd[i] = Iend[i];

                    int li = rEnd[i] - rStart[i] + 1;

                    if(rStart[i] < 0 || rStart[i] >= GetLength(i))
                        throw new IndexOutOfRangeException("subarray start range out of range for dimension " + i);
                    if(rEnd[i] >= GetLength(i))
                        throw new IndexOutOfRangeException("subarray end range out of range for dimension " + i);

                    if(li > 0)
                        dims++;
                    else {
                        rEnd[i] = -2;
                    }
                }
                if(dims <= 0)
                    throw new ArgumentOutOfRangeException("0-dimensional arrays are not supported.");

                // initialize SubArray
                MultidimensionalArray ret = new MultidimensionalArray(dims);
                ret.m_Storage = this.m_Storage;
                ret.m_ShallowCopy = true;
                ret.m_LockedForever = this.m_LockedForever;

                // calculate offset 
                ret.m_Offset = this.IndexUnchecked(rStart);

                // calculate cycles and lengths
                int ii = 0;
                for(int i = 0; i < m_Dimension; i++) {
                    int li = rEnd[i] - rStart[i] + 1;

                    int Cyc = int.MinValue;
                    if(GetLength(i) >= 1) {
                        int i0 = this.IndexUnchecked(rStart);
                        rStart[i]++;
                        //if (li > 1)
                        //    Cyc = this.IndexUnchecked(rStart) - i0;
                        //else
                        //    Cyc = 1;
                        Cyc = this.IndexUnchecked(rStart) - i0;
                        rStart[i]--;
                    }

                    if(li > 0) {
                        ret.SetLength(ii, li);
                        ret.SetCycle(ii, Cyc);
                        ii++;
                    }
                }

                // continuous piece of memory?
                //ret.m_IsContinious = true;
                //if (ret.m_Cycles[0] != 1) {
                //    ret.m_IsContinious = false;
                //} else {
                //    for (int i = 1; i < ret.m_Lengths.Length; i++) {
                //        if (ret.m_Cycles[i] != ret.m_Lengths[i-1]) {
                //            ret.m_IsContinious = false;
                //            break;
                //        }
                //    }
                //}


                // return
                return ret;
            }
        }

        /// <summary>
        /// copies one row or column or ... to a vector
        /// </summary>
        /// <param name="vec">output</param>
        /// <param name="d">rank to extract</param>
        /// <param name="st">first element (in the d-th index) to extract</param>
        /// <param name="len">number of elements to copy</param>
        /// <param name="index">
        /// indices for all dimensions
        /// </param>
        public void ExtractVector<T>(T vec, int d, int st, int len, params int[] index) where T : IList<double> {
            if (d < 0 || d >= this.Dimension)
                throw new ArgumentOutOfRangeException("d");

            int en = st + len;
            for (index[d] = st; index[d] < en; index[d]++) {
                vec[index[d] - st] = this[index];
            }
        }

        /// <summary>
        /// Extract a sub-array; the sub-array is a so-called "shallow copy", that means, changes on the
        /// sub-array also affect this array; Because there is no internal storage allocation,
        /// this is a relatively light-weighted operation;
        /// </summary>
        /// <param name="subArrayIndices">a negative entry means that the associated dimension
        /// is in the sub-array; Example: let be A[*,*,*] the 3D mother array; than {-1,2,-1} 
        /// extracts the 2D - array A[*,2,*];</param>
        /// <returns></returns>
        public MultidimensionalArray ExtractSubArrayShallow(params int[] subArrayIndices) {
            // check Arguments
            if (subArrayIndices.Length != this.Dimension) {
                throw new ArgumentOutOfRangeException("Mismatch between array dimension and length of argument.");
            }

            // dimension of SubArray
            int dims = 0;
            for (int i = 0; i < subArrayIndices.Length; i++) {
                if (subArrayIndices[i] < 0)
                    dims++;
            }
            if (dims <= 0)
                throw new ArgumentOutOfRangeException("0-dimensional arrays are not supported.");

            // initialize SubArray
            MultidimensionalArray ret = new MultidimensionalArray(dims);
            ret.m_Storage = this.m_Storage;
            ret.m_ShallowCopy = true;
            ret.m_LockedForever = this.m_LockedForever;

            // calculate offset ...
            int[] indices = (int[])subArrayIndices.Clone();
            for (int j = 0; j < indices.Length; j++) if (indices[j] < 0)
                    indices[j] = 0;
            ret.m_Offset = this.Index(indices);

            // ... and cycles:
            int cnt = 0;
            for (int i = 0; i < subArrayIndices.Length; i++) {
                indices = (int[])subArrayIndices.Clone();

                if (indices[i] >= 0)
                    continue; // this dimension is not included in the 
                // returned array
                if (this.GetLength(i) <= 1) {
                    // the only valid index is 0, 
                    // so the cycle is of no interrest
                    // (and cannot be calculated);

                    ret.SetCycle(cnt, Int32.MinValue);
                    ret.SetLength(cnt, this.GetLength(i));
                } else {

                    for (int j = 0; j < indices.Length; j++) if (indices[j] < 0)
                            indices[j] = 0;
                    indices[i] = 1;

                    int i1 = this.Index(indices);

                    ret.SetCycle(cnt, i1 - ret.m_Offset);
                    ret.SetLength(cnt, this.GetLength(i));
                }
                cnt++;
            }

            // return
            return ret;
        }

        /// <summary>
        /// clears all entries
        /// </summary>
        public void Clear() {
            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");

            //int[] ind = new int[this.Dimension];
            if (this.IsContinious) {
                int L = this.Length;
                if (L > 0)
                    Array.Clear(m_Storage, this.m_Offset, L);
            } else {
                ApplyAll(d => 0);
            }
        }

        /// <summary>
        /// checks all entries for infinity or NAN - values, and
        /// throws an <see cref="ArithmeticException"/> if found;
        /// </summary>
        public int CheckForNanOrInf(bool CheckForInf = true, bool CheckForNan = true, bool ExceptionIfFound = true) {

            int Problems = 0;

            this.ApplyAll(delegate(int[] index, double a) {

                
                if (CheckForNan)
                    if (double.IsNaN(a)) {
                        if(ExceptionIfFound) {
                            using(StringWriter stw = new StringWriter()) {
                                stw.Write("[");
                                for(int i = 0; i < index.Length - 1; i++) {
                                    stw.Write(index[i]);
                                    stw.Write(", ");
                                }
                                stw.Write(index.Last());
                                stw.Write("]");

                                throw new ArithmeticException("NaN found at " + stw.ToString() + "-th entry.");
                            }
                        } else {
                            Problems++;
                        }
                    }

                if (CheckForInf)
                    if (double.IsInfinity(a)) {

                        if(ExceptionIfFound) {
                            using(StringWriter stw = new StringWriter()) {
                                stw.Write("[");
                                for(int i = 0; i < index.Length - 1; i++) {
                                    stw.Write(index[i]);
                                    stw.Write(", ");
                                }
                                stw.Write(index.Last());
                                stw.Write("]");


                                throw new ArithmeticException("Inf found at " + stw.ToString() + "-th entry.");
                            }
                        } else {
                            Problems++;
                        }
                    }
            });

            return Problems;
        }

        /// <summary>
        /// (re-) allocates the array
        /// </summary>
        /// <param name="__Lengths">
        /// The number of entries to be allocated for each dimension.
        /// </param>
        public void Allocate(params int[] __Lengths) {
            if (this.m_ShallowCopy == true)
                throw new NotSupportedException("(Re-)Allocation is not supported for shallow copies.");

            if (__Lengths.Length != this.Dimension) {
                throw new ArgumentOutOfRangeException("Mismatch between array dimension and number of given lengths.");
            }

            if (m_LockedForever)
                throw new ApplicationException("illegal call - object is locked.");
            for (int i = 0; i < this.m_Dimension; i++) {
                this.SetLength(i, __Lengths[i]);
            }


            int absLen = 1;
            for (int i = 0; i < __Lengths.Length; i++) {
                absLen *= __Lengths[i];
            }

            if (absLen < 0) {
                throw new ArgumentOutOfRangeException("Array length smaller than 0 makes no sense.");
            }

            if (this.m_Storage != null && this.m_Storage.Length == absLen)
                Array.Clear(this.m_Storage, 0, absLen);
            else
                m_Storage = new double[absLen];
            m_ShallowCopy = false;
            CreateCycles();

            m_Offset = 0;
        }

        /// <summary>
        /// Uses <see cref="SetCycle"/> and <see cref="GetLength"/> to set up
        /// <see cref="m_StorageLayout"/>
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private void CreateCycles() {
            SetCycle(m_Dimension - 1, 1);
            for (int i = m_Dimension - 2; i >= 0; i--) {
                SetCycle(i, GetCycle(i + 1) * GetLength(i + 1));
            }
        }

        /// <summary>
        /// initializes this array as a copy of a 2d-array
        /// </summary>
        /// <param name="source"></param>
        public void InitializeFrom(double[,] source) {
            int[] NewLen = new int[] { source.GetLength(0), source.GetLength(1) };
            Allocate(NewLen);
            //CopyFrom(source);
            Set2DSubArray(source, -1, -1);
        }

        ///// <summary>
        ///// Copies values from a 2-dimensional .NET-array into this array; This
        ///// array must be 2-dimensional to support his operation (see
        ///// <see cref="Dimension"/>).
        ///// </summary>
        ///// <param name="source">The source array</param>
        //public void CopyFrom(double[,] source) {
        //    //CopyFrom(source, 0);
        //    Set(source, -1, -1);
        //}

        ///// <summary>
        ///// Copies values from a 2-dimensional .NET-array into this array; This
        ///// array must be 2-dimensional to support his operation (see
        ///// <see cref="Dimension"/>); An offset into the first index of this
        ///// array can be specified: <paramref name="source"/>[0,0] will be
        ///// copied to this[<paramref name="i0"/>,0];
        ///// </summary>
        ///// <param name="source">The source array</param>
        ///// <param name="i0">offset into this array, 1st index;</param>
        //public void CopyFrom(double[,] source, int i0) {
        //    if (this.m_Lengths.Length != 2) {
        //        throw new ArgumentOutOfRangeException("dimension mismatch");
        //    }

        //    if (m_IsContinious) {
        //        Buffer.BlockCopy(
        //            source,
        //            0,
        //            m_Storage,
        //            Index(i0, 0) * sizeof(double),
        //            source.Length * sizeof(double));
        //    } else {
        //        for (int i = 0; i < source.GetLength(0); i++) {
        //            for (int j = 0; j < source.GetLength(1); j++) {
        //                this[i + i0, j] = source[i, j];
        //            }
        //        }
        //    }
        //}

        ///// <summary>
        ///// Copies values from a 2-dimensional .NET-array into this array; This
        ///// array must be 2-dimensional to support his operation (see
        ///// <see cref="Dimension"/>); An offset into this array can be
        ///// specified: <paramref name="source"/>[0,0] will be copied to
        ///// this[<paramref name="i0"/>,<paramref name="i1"/>];
        ///// </summary>
        ///// <param name="source">The source array</param>
        ///// <param name="i0">offset into this array, 1st index;</param>
        ///// <param name="i1">offset into this array, 2nd index;</param>
        ///// <remarks>
        ///// In comparison to <see cref="CopyFrom(double[,])"/> and
        ///// <see cref="CopyFrom(double[,], int)"/>, this operation is rather
        ///// expensive for large sources.
        ///// </remarks>
        //public void CopyFrom(double[,] source, int i0, int i1) {
        //    if (this.m_Lengths.Length != 2)
        //        throw new ArgumentOutOfRangeException("dimension mismatch");

        //    int I = source.GetLength(0), J = source.GetLength(1);
        //    for (int i = 0; i < I; i++) {
        //        for (int j = 0; j < J; j++) {
        //            this[i + i0, j + i1] = source[i, j];
        //        }
        //    }
        //}

        /// <summary>
        /// Copies values from this array into a 2-dimensional .NET-array;
        /// This array must be 2-dimensional to support his operation (<see cref="Dimension"/>);
        /// </summary>
        /// <param name="dest">
        /// target for the copy operation;
        /// the [0,0]-entry of this array will end up in the [<paramref name="i0"/>,<paramref name="i1"/>]-entry
        /// of <paramref name="dest"/>;
        /// </param>
        /// <param name="i0">offset into this array, 1st index;</param>
        /// <param name="i1">offset into this array, 2nd index;</param>
        public void CopyTo(double[,] dest, int i0, int i1) {
            if (this.Dimension != 2)
                throw new ArgumentOutOfRangeException("dimension mismatch");

            int I = this.GetLength(0), J = this.GetLength(1);
            for (int i = 0; i < I; i++) {
                for (int j = 0; j < J; j++) {
                    dest[i + i0, j + i1] = this[i, j];
                }
            }
        }

        /// <summary>
        /// This array must be 2-dimensional to support his operation (<see cref="Dimension"/>);
        /// </summary>
        /// <returns>
        /// a 2d .NET-array which holds the values of this array;
        /// </returns>
        public double[,] To2DArray() {
            if (this.Dimension != 2)
                throw new ApplicationException("this array must be 2d");
            double[,] ret = new double[GetLength(0), GetLength(1)];
            CopyTo(ret, 0, 0);
            return ret;
        }

        /// <summary>
        /// This array must be 1-dimensional to support his operation (<see cref="Dimension"/>);
        /// </summary>
        /// <returns>
        /// a 1d .NET-array which holds the values of this array;
        /// </returns>
        public double[] To1DArray() {
            if (this.Dimension != 1)
                throw new ApplicationException("this array must be 1d");
            double[] ret = new double[GetLength(0)];

            int L = GetLength(0);
            for (int i = 0; i < L; i++)
                ret[i] = this[i];
            return ret;
        }

        /// <summary>
        /// empty constructor for De-Serialization.
        /// </summary>
        private MultidimensionalArray() {
        }


        /// <summary>
        /// creates a new multidimensional array with a fixed number of
        /// dimensions, but does not allocate Memory
        /// </summary>
        /// <param name="NumberOfdimensions">
        /// The number of dimensions of the array.
        /// </param>
        /// <returns></returns>
        public MultidimensionalArray(int NumberOfdimensions) {
            if (NumberOfdimensions <= 0)
                throw new ArgumentOutOfRangeException();
            if (NumberOfdimensions > 8)
                throw new NotSupportedException();

            this.m_Dimension = NumberOfdimensions;
        }

        /// <summary>
        /// creates a new multidimensional array and allocates memory, nD - version
        /// </summary>
        /// <param name="Lengths"></param>
        /// <returns></returns>
        public static MultidimensionalArray Create(params int[] Lengths) {
            MultidimensionalArray ret = new MultidimensionalArray(Lengths.Length);
            ret.Allocate(Lengths);
            return ret;
        }

        /// <summary>
        /// Creates a wrapper for the given one-dimensional array. That is,
        /// the result will be a shallow 'copy' of
        /// <paramref name="wrappedArray"/>.
        /// </summary>
        /// <param name="wrappedArray">
        /// The array to be wrapped by the newly created object
        /// </param>
        /// <param name="lengths">
        /// The dimensions of the wrapper.
        /// </param>
        /// <returns>
        /// A <see cref="MultidimensionalArray"/> with
        /// <see cref="MultidimensionalArray.Storage"/> =
        /// <paramref name="wrappedArray"/>.
        /// </returns>
        public static MultidimensionalArray CreateWrapper(double[] wrappedArray, params int[] lengths) {
            Debug.Assert(lengths.All((i) => i > 0), "All lengths must be positive");
            //Debug.Assert(
            //    lengths.Aggregate((i, j) => i * j) == wrappedArray.Length,
            //    "Length of wrapped array must be equal to length of wrapper");

            MultidimensionalArray wrapper = new MultidimensionalArray(lengths.Length);
            wrapper.m_Storage = wrappedArray;
            for (int i = 0; i < lengths.Length; i++)
                wrapper.SetLength(i, lengths[i]);
            wrapper.m_ShallowCopy = false;
            wrapper.m_Offset = 0;
            wrapper.CreateCycles();

            if (wrapper.Length > wrappedArray.Length)
                throw new ArgumentException();

            return wrapper;
        }

        ///// <summary>
        ///// creates a new multidim array and allocates memory, 2D - version
        ///// </summary>
        ///// <param name="I"></param>
        ///// <param name="J"></param>
        ///// <returns></returns>
        //public static MultidimensionalArray Create(int I, int J) {
        //    return Create(new int[] { I, J });
        //}

        ///// <summary>
        ///// creates a new multidim array and allocates memory, 1D - version
        ///// </summary>
        ///// <param name="I"></param>
        ///// <returns></returns>
        //public static MultidimensionalArray Create(int I) {
        //    return Create(new int[] { I });
        //}

        ///// <summary>
        /////  creates a new multidim array and allocates memory, 3D - version
        ///// </summary>
        ///// <param name="I"></param>
        ///// <param name="J"></param>
        ///// <param name="K"></param>
        ///// <returns></returns>
        //public static MultidimensionalArray Create(int I, int J, int K) {
        //    return Create(new int[] { I, J, K});
        //}

        #region IMatrix Members

        /// <summary>
        /// the 1st length, i.e. equal to <see cref="GetLength"/>(1);
        /// </summary>
        public int NoOfCols {
            get {
                if (Dimension != 2)
                    throw new NotSupportedException("array is not 2D and cannot be represented like a matrix");
                return GetLength(1);
            }
        }

        /// <summary>
        /// the 0-th length, i.e. equal to <see cref="GetLength"/>(0)
        /// </summary>
        public int NoOfRows {
            get {
                if (Dimension != 2)
                    throw new NotSupportedException("array is not 2D and cannot be represented like a matrix");
                return GetLength(0);
            }
        }

        /// <summary>
        /// copies a part of this array to <paramref name="x"/>;
        /// </summary>
        /// <param name="x">destination</param>
        /// <param name="SubArrayIdx"></param>
        public void CopyTo<T>(T x, params int[] SubArrayIdx) where T : IList<double> {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);

            if (sub.Dimension != 1)
                throw new ArgumentException("illegal subarray selection - not a row");
            if ((x.Count) < sub.Length)
                throw new ArgumentException("destination vector is to short");

            int L = sub.Length;
            for (int i = 0; i < L; i++)
                x[i] = sub[i];
        }

        /// <summary>
        /// copies a part of this array to <paramref name="x"/>;
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        /// <param name="xOffset">offset into <paramref name="x"/>, where the first element will be inserted</param>
        public void CopyTo<T>(T x, int[] Istart, int[] Iend, int xOffset = 0) where T : IList<double> {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);

            if (sub.Dimension != 1)
                throw new ArgumentException("illegal subarray selection - not a row");
            if ((x.Count - xOffset) < sub.Length)
                throw new ArgumentException("destination vector is to short");

            int L = sub.Length;
            for (int i = 0; i < L; i++)
                x[i + xOffset] = sub[i];
        }

        /// <summary>
        /// Copies to a one-dimensional array/list
        /// </summary>
        /// <param name="array"></param>
        /// <param name="RowWise">
        /// if true, elements are taken row-by-row (column index rotating fastest) and copied to <paramref name="array"/>;
        /// if false, elements are taken column-by-column;
        /// </param>
        /// <param name="arrayoffset"></param>
        /// <remarks>
        /// This array must be 2-dimensional to support his operation (<see cref="Dimension"/>);
        /// </remarks>
        public void CopyTo<T>(T array, bool RowWise, int arrayoffset)
            where T : IList<double> {
            if (Dimension != 2)
                throw new NotSupportedException("array is not 2D and cannot be represented like a matrix");

            int I = NoOfRows, J = NoOfCols;
            int targetind = arrayoffset;

            if (RowWise) {
                // concatenate rows
                // ++++++++++++++++

                if (IsContinious && typeof(T).Equals(typeof(double[]))) {
                    // optimized version

                    double[] _array = array as double[];
                    Array.Copy(this.m_Storage, this.Index(0, 0), _array, arrayoffset, this.Length);
                } else {

                    if (IsContinious) {
                        // still somewhat optimized

                        int srcInd = this.Index(0, 0);
                        int IJ = I * J;
                        for (int ij = 0; ij < IJ; ij++)
                            array[ij] = m_Storage[ij + srcInd];

                    } else {
                        // loop over rows...
                        for (int i = 0; i < I; i++) {
                            // loop over columns...
                            for (int j = 0; j < J; j++) {
                                array[targetind] = this[i, j];
                                targetind++;
                            }
                        }
                    }
                }
            } else {
                // concatenate columns
                // +++++++++++++++++++

                // loop over columns...
                for (int j = 0; j < J; j++) {
                    // loop over rows...
                    for (int i = 0; i < I; i++) {
                        array[targetind] = this[i, j];
                        targetind++;
                    }
                }

            }

        }

        #endregion

        #region ICloneable Members

        /// <summary>
        /// creates a non-shallow copy of this object
        /// </summary>
        virtual public object Clone() {
            return CloneAs();
        }

        #endregion

        /// <summary>
        /// Creates a non-shallow copy of this object; there is one difference, allthough:
        /// if this object is locked (<see cref="IsLocked"/>), the clone is not.
        /// </summary>
        public MultidimensionalArray CloneAs() {
            MultidimensionalArray ret = MultidimensionalArray.Create(this.Lengths);
            if (this.IsContinious) {
                int iSrc0 = this.m_Offset;
                int Len = this.Length;
                Array.Copy(this.m_Storage, iSrc0, ret.m_Storage, 0, Len);
            } else {
                CloneRecursive(ret, new int[this.Dimension], 0);
            }
            return ret;
        }

        /// <summary>
        /// If this object is a shallow copy, recursively clones all parents
        /// </summary>
        /// <param name="ret"></param>
        /// <param name="idx"></param>
        /// <param name="d"></param>
        private void CloneRecursive(MultidimensionalArray ret, int[] idx, int d) {
            if (d + 1 >= this.Dimension) {
                for (idx[d] = this.GetLength(d) - 1; idx[d] >= 0; idx[d]--) {
                    ret[idx] = this[idx];
                }
            } else {
                for (idx[d] = this.GetLength(d) - 1; idx[d] >= 0; idx[d]--) {
                    CloneRecursive(ret, idx, d + 1);
                }
            }
        }

        /// <summary>
        /// calling this method locks the array values for the whole lifetime
        /// of the object, i.e.  any further assess to <see cref="this[int]"/>
        /// will cause an exception. Of course, locking cannot be 100 percent
        /// water-tight, e.g. direct manipulation of <see cref="Storage"/> is
        /// still possible. This method is intended mainly for debugging
        /// purpose and may disappear in future.
        /// </summary>
        public void LockForever() {
            m_LockedForever = true;
        }

        /// <summary>
        /// Indicates whether <see cref="LockForever"/> has been called in the livetime of this
        /// object or not.
        /// </summary>
        public bool IsLocked {
            get {
                return m_LockedForever;
            }
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        public void SetSubVector<T>(T x, int[] Istart, int[] Iend) where T : IList<double> {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.Clear();
            sub.AccVector(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        public void SetSubMatrix(IMatrix x, int[] Istart, int[] Iend) {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.Clear();
            sub.AccMatrix(1.0, x);
        }

        /// <summary>
        /// Sets a part of this array to <paramref name="x"/>
        /// ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        public void Set2DSubArray(double[,] x, int[] Istart, int[] Iend) {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.Clear();
            sub.Acc2DArray(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="Istart">start indices (including) into this array</param>
        /// <param name="Iend">end indices (including) into this array</param>
        public void SetSubArray(MultidimensionalArray x, int[] Istart, int[] Iend) {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.Clear();
            sub.Acc(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        public void SetVector<T>(T x) where T : IList<double> {
            this.Clear();
            this.AccVector(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="SubArrayIdx"></param>
        public void SetSubVector<T>(T x, params int[] SubArrayIdx) where T : IList<double> {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.Clear();
            sub.AccVector(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="SubArrayIdx"></param>
        public void SetSubMatrix(IMatrix x, params int[] SubArrayIdx) {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.Clear();
            sub.AccMatrix(1.0, x);
        }

        /// <summary>
        /// sets a this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        public void SetMatrix(IMatrix x) {
            this.Clear();
            this.AccMatrix(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="SubArrayIdx"></param>
        public void Set2DSubArray(double[,] x, params int[] SubArrayIdx) {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.Clear();
            sub.Acc2DArray(1.0, x);
        }

        /// <summary>
        /// sets this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        public void Set2DArray(double[,] x) {
            this.Clear();
            this.Acc2DArray(1.0, x);
        }

        /// <summary>
        /// sets a part of this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="SubArrayIdx"></param>
        public void SetSubArray(MultidimensionalArray x, params int[] SubArrayIdx) {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.Clear();
            sub.Acc(1.0, x);
        }

        /// <summary>
        /// sets this array to <paramref name="x"/> ('Copy-From-Operation');
        /// </summary>
        /// <param name="x">new values</param>
        public void Set(MultidimensionalArray x) {
            if (x.m_Dimension != this.m_Dimension)
                throw new ArgumentException("Number of dimensions must be equal.");
            for (int i = this.m_Dimension - 1; i >= 0; i--) {
                if (x.GetLength(i) != this.GetLength(i))
                    throw new ArgumentException("Mismatch in dimension " + i);
            }


            if (this.IsContinious && x.IsContinious) {
                Array.Copy(x.m_Storage, x.m_Offset, this.Storage, this.m_Offset, this.Length);
            } else {
                int[] insert = new int[this.Dimension];
                insert.SetAll(-1);
                SetSubArray(x, insert);
            }
        }

        /// <summary>
        /// sets all entries in a part of this array to <paramref name="x"/>.
        /// </summary>
        /// <param name="x">new value</param>
        /// <param name="Istart">start indices (including)</param>
        /// <param name="Iend">end indices (including)</param>
        public void SetAll(double x, int[] Istart, int[] Iend) {
            var sub = this.ExtractSubArrayShallow(Istart, Iend);
            sub.SetAll(x);
        }

        /// <summary>
        /// sets all entries in a part of this array to <paramref name="x"/>.
        /// </summary>
        /// <param name="x">new values</param>
        /// <param name="SubArrayIdx"></param>
        public void SetAll(double x, params int[] SubArrayIdx) {
            var sub = this.ExtractSubArrayShallow(SubArrayIdx);
            sub.SetAll(x);
        }

        /// <summary>
        /// throws an exception, if the lengths of this array do not match <paramref name="l"/>.
        /// </summary>
        /// <param name="l"></param>
        [Conditional("DEBUG")]
        public void CheckLengths(params int[] l) {
            if (this.Dimension != l.Length)
                throw new ApplicationException("mismatch in array dimension: expecting a " + l.Length + " - dimensional array, this array is " + this.Dimension + " - dimensional.");
            for (int i = 0; i < this.Dimension; i++) {
                if (this.GetLength(i) != l[i])
                    throw new ApplicationException("mismatch in array length, in " + i + "-th dimension.");
            }
        }

        /// <summary>
        /// <see cref="Equals(MultidimensionalArray)"/>
        /// </summary>
        public override bool Equals(object obj) {
            if (obj is MultidimensionalArray) {
                return Equals((MultidimensionalArray)obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Taken from the first two entries and the length of this array.
        /// </summary>
        public override int GetHashCode() {
            if (this.Length <= 0)
                return 0;

            double entry0 = this.m_Storage[m_Offset];

            unsafe {
                int* pentry = (int*)(&entry0);
                int p1, p2;
                p1 = *pentry;
                pentry++;
                p2 = *pentry;

                return (p1 ^ p2 + this.Length);
            }
        }

        #region IEquatable<MultidimensionalArray> Members

        /// <summary>
        /// Checks all entries for equality
        /// </summary>
        /// <param name="other">Array to compare to</param>
        /// <returns>
        /// True, if the entries of <paramref name="other"/> are pairwise equal
        /// to the entries of this array.
        /// </returns>
        public bool Equals(MultidimensionalArray other) {
            if (object.ReferenceEquals(this, other))
                return true;

            if (this.Dimension != other.Dimension)
                return false;

            if (this.Length != other.Length)
                return false;

            if (this.Length == 0)
                return true;

            if (other.m_Dimension != this.m_Dimension)
                return false;

            for (int i = this.m_Dimension - 1; i >= 0; i--) {
                if (other.GetLength(i) != this.GetLength(i)) {
                    return false;
                }
            }

            if (this.IsContinious && other.IsContinious) {
                // accelerated version
                unsafe {
                    fixed (double* pThis = &this.m_Storage[0], pOthr = &other.m_Storage[0]) {
                        int L = this.Length;

                        // Interpret doubles as longs for faster comparison
                        long* pX = (long*)(pThis + this.m_Offset);
                        long* pY = (long*)(pOthr + this.m_Offset);

                        for (int l = L - 1; l >= 0; l--) {
                            if (*pX != *pY) {
                                return false;
                            }

                            pX++;
                            pY++;
                        }

                        return true;
                    }
                }
            } else {
                // Reference implementation
                bool eq = true;
                this.ApplyAll(delegate(int[] index, ref double entry) {
                    if (other[index] != entry)
                        eq = false;
                });
                return eq;
            }
        }

        #endregion
    }
}
