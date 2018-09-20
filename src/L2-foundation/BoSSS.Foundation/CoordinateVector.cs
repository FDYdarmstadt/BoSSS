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
using System.Collections;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;

namespace BoSSS.Foundation {

    /// <summary>
    /// By using a <see cref="CoordinateMapping"/>, this class
    /// presents the DG coordinates of a list of <see cref="DGField"/>-objects
    /// as one one-dimensional vector, without allocating additional memory.
    /// </summary>
    public class CoordinateVector : IList<double> {

        bool m_PresentExternal = false;

        /// <summary>
        /// if true, also the DG coordinates of external cells (aka. ghost cells) are
        /// presented in this vector;
        /// This property can only be set at construction time.
        /// </summary>
        public bool PresentExternal {
            get {
                return m_PresentExternal;
            }
        }

        /// <summary>
        /// scales all fields in this mapping (see <see cref="CoordinateMapping.Fields"/>) by a factor <paramref name="a"/>;
        /// </summary>
        /// <param name="a"></param>
        public void Scale(double a) {
            foreach (DGField f in m_mapping.m_Fields)
                f.Scale(a);
        }

        /// <summary>
        /// accumulates a vector to the coordinates of all fields in this mapping (see <see cref="CoordinateMapping.Fields"/>);
        /// </summary>
        /// <param name="x">
        /// Length of <paramref name="x"/> vector must be equal the update length (<see cref="Count"/> or <see cref="UnsetteledCoordinateMapping.LocalLength"/>);
        /// </param>
        /// <param name="alpha">
        /// accumulation factor for <paramref name="x"/>.
        /// </param>
        public void Acc<T>(double alpha, T x) where T : IList<double> {
            bool onlyLocal = true;
            if (x.Count != m_mapping.LocalLength)
                throw new ArgumentException("length of x must be either NUpdate", "x");

            if (!m_mapping.CellDepLength) {
                // constant number of DOF per cell
                for (int f = 0; f < m_mapping.m_Fields.Length; f++) {
                    int Skip = m_mapping.m_j0CoordinateIndex[f];
                    int Stride = m_mapping.MaxTotalNoOfCoordinatesPerCell;
                    m_mapping.m_Fields[f]._Acc(alpha, x, Skip, Stride, onlyLocal);
                }
            } else {
                GenericBlas.daxpy(x.Count, alpha, x, 1, this, 1);
            }
        }

        /// <summary>
        /// <see cref="Mapping"/>
        /// </summary>
        CoordinateMapping m_mapping;

        /// <summary>
        /// the mapping of this vector to the list of DG fields;
        /// </summary>
        public CoordinateMapping Mapping {
            get {
                return m_mapping;
            }
        }

        /// <summary>
        /// Constructs a new vector, based on the list of DG fields,
        /// <paramref name="fields"/>;
        /// </summary>
        /// <param name="fields"></param>
        /// <param name="presentExternal">
        /// value for the <see cref="PresentExternal"/>-property;
        /// </param>
        public CoordinateVector(bool presentExternal, params DGField[] fields) :
            this(new CoordinateMapping(fields), presentExternal) {
        }

        /// <summary>
        /// Constructs a new vector, based on the list of DG fields,
        /// <paramref name="fields"/>;
        /// </summary>
        /// <remarks>
        /// <see cref="PresentExternal"/> is initialized to false.
        /// </remarks>
        public CoordinateVector(params DGField[] fields) :
            this(false, fields) {
        }

        /// <summary>
        /// Constructs a new vector, based on the mapping <paramref name="mapping"/>;
        /// </summary>
        /// <param name="mapping">
        /// The mapping to transform indices; This parameter is
        /// used to initialize <see cref="Mapping"/>.
        /// </param>
        /// <param name="presentExternal">
        /// value for the <see cref="PresentExternal"/>-property;
        /// </param>
        public CoordinateVector(CoordinateMapping mapping, bool presentExternal) {
            m_mapping = mapping;
            m_PresentExternal = presentExternal;
            if (m_PresentExternal)
                m_Count = m_mapping.Ntotal;
            else
                m_Count = m_mapping.LocalLength;

        }

        /// <summary>
        /// Constructs a new vector, based on the mapping <paramref name="mapping"/>;
        /// </summary>
        /// <param name="mapping">
        /// The mapping to transform indices; This parameter is
        /// used to initialize <see cref="Mapping"/>.
        /// </param>
        /// <remarks>
        /// <see cref="PresentExternal"/> is initialized to false.
        /// </remarks>
        public CoordinateVector(CoordinateMapping mapping) :
            this(mapping, false) {
        }

        /// <summary>
        /// Copies the content of the vector <paramref name="array"/> into the 
        /// <see cref="DGField"/>'s in this mapping (<see cref="CoordinateMapping.Fields"/>);
        /// External cells are also copied, but no network update occurs.
        /// </summary>
        /// <param name="array"></param>
        /// <param name="arrayIndex"></param>
        public void CopyFrom<T>(T array, int arrayIndex) where T : IList<double> {
            Copy(array, arrayIndex, false);
        }

        /// <summary>
        /// this = this + <paramref name="x"/>*<paramref name="xScale"/>;
        /// using this function changes the coordinates of the fields in this mapping (<see cref="CoordinateMapping.Fields"/>);
        /// Also external cells are affected;
        /// </summary>
        /// <param name="x"></param>
        /// <param name="xScale"></param>
        /// <typeparam name="VectorType"></typeparam>
        public void axpy<VectorType>(VectorType x, double xScale) where VectorType : IList<double> {
            int nc;
            if (x.Count == m_mapping.LocalLength)
                nc = m_mapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            else if (x.Count == this.Count)
                nc = m_mapping.GridDat.iLogicalCells.Count;
            else
                throw new ArgumentException("length of x must be either NUpdate or Count.");

            int nf = m_mapping.m_Fields.Length;

            int[] len = new int[nf];
            bool slowVersion = false;
            for (int n = 0; n < len.Length; n++) {
                int l_max = m_mapping.m_Fields[n].Basis.MaximalLength;
                int l_min = m_mapping.m_Fields[n].Basis.MinimalLength;
                if (l_max == l_min)
                    len[n] = l_max;
                else {
                    len[n] = -1;
                    slowVersion = true;
                }
            }

            IMatrix[] Coordinates = new IMatrix[m_mapping.m_Fields.Length];
            for (int f = 0; f < Coordinates.Length; f++)
                Coordinates[f] = m_mapping.m_Fields[f].Coordinates;

            if (slowVersion) {
                // for at least one field, the number of DG coordinates 
                // may vary from cell to cell
                for (int c = 0; c < nc; c++) {
                    for (int f = 0; f < nf; f++) {
                        int L = len[f];
                        if (L < 0) {
                            L = m_mapping.m_Fields[f].Basis.GetLength(c);
                        }

                        //int Lmax = m_mapping.m_Fields[f].Basis.MaximalLength;
                        for (int n = 0; n < L; n++) {
                            Coordinates[f][c, n] += x[m_mapping.LocalUniqueCoordinateIndex(f, c, n)] * xScale;
                        }

                    }
                }
            } else {
                // optimized version: the same number of DG coordinates in each cell
                //                    for all fields
                int i = 0;
                for (int c = 0; c < nc; c++) {
                    for (int f = 0; f < nf; f++) {
                        for (int n = 0; n < len[f]; n++) {
                            Coordinates[f][c, n] += x[i] * xScale;
                            i++;
                        }

                    }
                }
            }
        }

        ///// <summary>
        ///// copys from vector <paramref name="source"/>, which indices are defined by <paramref name="SourceMapping"/>
        ///// to another vector <paramref name="dest"/>, which indices are defined by <paramref name="DestMapping"/> .
        ///// the source mapping must be "a subset" of the the destination mapping, i.e. all <see cref="Field"/>'s of the 
        ///// <paramref name="SourceMapping"/> (see <see cref="CoordinateMapping.Fields"/>) must be part of fields-collection in the
        ///// paramref name="DestMapping"/>-mapping.
        ///// External cells are also copied, but no network update occurs.
        ///// </summary>
        ///// <param name="SourceMapping"></param>
        ///// <param name="source"></param>
        ///// <param name="DestMapping"></param>
        ///// <param name="dest"></param>
        //public static void Copy(CoordinateMapping SourceMapping, double[] source, CoordinateMapping DestMapping, double[] dest) {
        //    // check input
        //    if(SourceMapping.m_Context != DestMapping.m_Context)
        //        throw new ArgumentException("source and destination mapping must refer top the same context");
        //    if (SourceMapping.Count != source.Length)
        //        throw new ArgumentException("size mismatch between source mapping and source array");
        //    if (DestMapping.Count != dest.Length)
        //        throw new ArgumentException("size mismatch between destination mapping and destination array");

        //    int[] Src2DstFieldIndices = new int[SourceMapping.m_Fields.Length];
        //    for (int f = 0; f < SourceMapping.m_Fields.Length; f++) {

        //        Src2DstFieldIndices[f] = Array.IndexOf<Field>(DestMapping.m_Fields, SourceMapping.m_Fields[f]);
        //        if(Src2DstFieldIndices[f] < 0)
        //            throw new ArgumentException("source mapping must be a subset of destination mapping.");
        //    }

        //    // cache some values ...
        //    Context ctx = SourceMapping.m_Context;
        //    int NoOfCells = ctx.GridDat.NoOfLocalUpdatedCells;
        //    int NoOfFields = SourceMapping.m_Fields.Length;

        //    int[] Len = new int[NoOfFields];
        //    int[] j0Dst = new int[NoOfFields];
        //    for (int f = 0; f < NoOfFields; f++) {
        //        Len[f] = SourceMapping.m_Fields[f].NoOfCoordinatesPerCell;
        //        j0Dst[f] = DestMapping.m_j0CoordinateIndex[Src2DstFieldIndices[f]];
        //    }

        //    int TotalNoCoordsDst = DestMapping.MaxTotalNoOfCoordinatesPerCell;

        //    // copy
        //    int iSrc, iDst;
        //    iSrc = 0;
        //    for( int j = 0; j < NoOfCells; j++) {
        //        for( int f = 0; f < NoOfFields; f++) {
        //            iDst = TotalNoCoordsDst*j + j0Dst[f];
        //            Array.Copy(source,iSrc,dest,iDst,Len[f]);
        //            iSrc += Len[f];
        //        }
        //    }


        //}



        /// <summary>
        /// implementation of the <see cref="CopyTo"/> and <see cref="CopyFrom"/>
        /// function
        /// </summary>
        /// <param name="array"></param>
        /// <param name="arrayIndex"></param>
        /// <param name="direction">
        /// true: <see cref="CopyTo"/>-mode; false: <see cref="CopyFrom"/>-mode;
        /// </param>
        private void Copy<T>(T array, int arrayIndex, bool direction)
            where T : IList<double> {
            // check array size
            if (array.Count < this.Count + arrayIndex)
                throw new ArgumentException("array to short");

            // cache some values for faster processing ...
            // --------------------------------------------
            int NoOfCells = m_mapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells;
            int NoOfFields = m_mapping.m_Fields.Length;
            IMatrix[] Storages = new IMatrix[NoOfFields];
            int[] Len = new int[NoOfFields];
            int[] n0 = new int[NoOfFields];
            for (int f = 0; f < NoOfFields; f++) {
                Storages[f] = m_mapping.m_Fields[f].Coordinates;

                int l_min = m_mapping.m_Fields[f].Basis.MinimalLength;
                int l_max = m_mapping.m_Fields[f].Basis.MaximalLength;
                if (l_max == l_min)
                    Len[f] = l_max;
                else {
                    Len[f] = -1;
                }

                n0[f] = this.m_mapping.LocalUniqueCoordinateIndex(f, 0, 0);
            }
            int Lmax = this.m_mapping.MaxTotalNoOfCoordinatesPerCell;

            // do copy
            // -------


            // copy
            int iScr, iDst;
            if (direction) {
                // CopyTo - Modus
                iDst = arrayIndex;
                for (int j = 0; j < NoOfCells; j++) {
                    
                    for (int f = 0; f < NoOfFields; f++) {
                        int L = Len[f];
                        if (L < 0)
                            L = m_mapping.m_Fields[f].Basis.GetLength(j);
                        int _iDst = iDst + n0[f];

                        for (int n = 0; n < L; n++) {
                            array[_iDst] = Storages[f][j, n];
                            _iDst++;
                        }
                        //iScr = m_mapping.m_Fields[f].m_Coordinates.Index(j, 0);
                        //Array.Copy(Storages[f], iScr, array, iDst, Len[f]);
                        //iDst += Len[f];
                    }

                    iDst += Lmax;
                }
            } else {
                // CopyFrom - Modus
                iScr = arrayIndex;
                for (int j = 0; j < NoOfCells; j++) {
                    for (int f = 0; f < NoOfFields; f++) {
                        int L = Len[f];
                        if (L < 0)
                            L = m_mapping.m_Fields[f].Basis.GetLength(j);
                        int _iScr = iScr + n0[f];
                        
                        for (int n = 0; n < L; n++) {
                            Storages[f][j, n] = array[_iScr];
                            _iScr++;
                        }
                        //iDst = m_mapping.m_Fields[f].m_Coordinates.Index(j, 0);
                        //Array.Copy(array, iScr, Storages[f], iDst, Len[f]);
                        //iScr += Len[f];
                    }
                    iScr += Lmax;
                }
            }
        }

        #region IList<double> Member

        /// <summary>
        /// implemented from interface, but not supported;
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public int IndexOf(double item) {
            throw new NotSupportedException();
        }

        /// <summary>
        /// not supported - array is static;
        /// </summary>
        /// <param name="index"></param>
        /// <param name="item"></param>
        public void Insert(int index, double item) {
            throw new NotSupportedException();
        }

        /// <summary>
        /// not supported - array is static;
        /// </summary>
        /// <param name="index"></param>
        public void RemoveAt(int index) {
            throw new NotSupportedException();
        }


        /// <summary>
        /// access all coordinates of all <see cref="DGField"/>'s as a linear array
        /// </summary>
        /// <param name="index">a local coordinate index</param>
        /// <returns></returns>
        public double this[int index] {
            get {
                if (index < 0 || index >= m_Count)
                    throw new IndexOutOfRangeException();
                DGField f;
                int j, n;
                m_mapping.LocalFieldCoordinateIndex(index, out f, out j, out n);
                return f.Coordinates[j, n];
            }
            set {
                if (index < 0 || index >= m_Count)
                    throw new IndexOutOfRangeException();
                DGField f;
                int j, n;
                m_mapping.LocalFieldCoordinateIndex(index, out f, out j, out n);
                f.Coordinates[j, n] = value;
            }
        }

        #endregion

        #region ICollection<double> Member

        /// <summary>
        /// not supported - array is static (i.e length cannot be changed);
        /// </summary>
        /// <param name="item"></param>
        public void Add(double item) {
            throw new NotSupportedException();
        }

        /// <summary>
        /// sets all entries to 0.0;
        /// </summary>
        public void Clear() {
            for (int f = 0; f < m_mapping.m_Fields.Length; f++) {
                m_mapping.m_Fields[f].Coordinates.Clear();
            }
        }

        /// <summary>
        /// not supported;
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Contains(double item) {
            throw new NotSupportedException();
        }

        /// <summary>
        /// Copies the content of the <see cref="DGField"/>'s in this mapping (<see cref="CoordinateMapping.Fields"/>)
        /// to a continuous region of memory;
        /// External cells are also copied, but no network update occurs.
        /// </summary>
        /// <param name="array">the destination array</param>
        /// <param name="arrayIndex">offset index into the destination array</param>
        public void CopyTo<T>(T array, int arrayIndex) where T : IList<double> {
            Copy(array, arrayIndex, true);
        }

        /// <summary>
        /// Copies the content of the <see cref="DGField"/>'s in this mapping (<see cref="CoordinateMapping.Fields"/>)
        /// to a continuous region of memory;
        /// External cells are also copied, but no network update occurs.
        /// </summary>
        /// <param name="array">the destination array</param>
        /// <param name="arrayIndex">offset index into the destination array</param>
        public void CopyTo(double[] array, int arrayIndex) {
            Copy<double[]>(array, arrayIndex, true);
        }

        /// <summary>
        /// allocates an array and copies all DG coordiantes into it;
        /// </summary>
        /// <returns></returns>
        public double[] ToArray() {
            double[] ret = new double[this.Count];
            CopyTo(ret, 0);
            return ret;
        }


        int m_Count;

        /// <summary>
        /// length of this vector, see <see cref="PresentExternal"/>
        /// </summary>
        public int Count {
            get {
                return m_Count;
            }
        }

        /// <summary>
        /// Alias for this.Count
        /// length of this vector, see <see cref="PresentExternal"/>
        /// </summary>
        public int Length
        {
            get
            {
                return this.Count;
            }
        }


        ///// <summary>
        ///// total number of coordinates in all processors
        ///// </summary>
        //public long GlobalCount {
        //    get {
        //        return m_Context.GridDat.GlobalNoOfCells * m_TotalNoOfCoordinatesPerCell;
        //    }
        //}

        /// <summary>
        /// false
        /// </summary>
        public bool IsReadOnly {
            get {
                return false;
            }
        }

        /// <summary>
        /// not supported - array is static;
        /// </summary>
        /// <param name="item"></param>
        /// <returns></returns>
        public bool Remove(double item) {
            throw new NotSupportedException();
        }

        #endregion



        #region IEnumerable<double> Member

        /// <summary>
        /// returns an generic enumerator;
        /// </summary>
        private class MyGenericEnumerator : MyEnumerator, IEnumerator<double> {

            public MyGenericEnumerator(CoordinateVector owner)
                : base(owner) {
            }


            #region IEnumerator<double> Member

            public new double Current {
                get {
                    return m_Owner[i];
                }
            }

            #endregion

            #region IDisposable Member

            public void Dispose() {
            }

            #endregion
        }


        /// <summary>
        /// returns an enumerator;
        /// </summary>
        /// <returns></returns>
        public IEnumerator<double> GetEnumerator() {
            return new MyGenericEnumerator(this);
        }

        #endregion

        #region IEnumerable Member

        /// <summary>
        /// Implementation of the IEnumerator-Interface
        /// </summary>
        private class MyEnumerator : IEnumerator {

            protected CoordinateVector m_Owner;

            protected int i = -1;

            public MyEnumerator(CoordinateVector owner) {
                m_Owner = owner;
            }

            #region IEnumerator Member

            public object Current {
                get {
                    return m_Owner[i];
                }
            }

            public bool MoveNext() {
                i++;
                if (i >= m_Owner.Count)
                    return false;
                else
                    return true;
            }

            public void Reset() {
                i = 0;
            }

            #endregion
        }


        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
            return new MyEnumerator(this);
        }

        #endregion
    }
}
