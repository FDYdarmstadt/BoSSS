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
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Foundation.XDG {
    
    /// <summary>
    /// a sparse data structure, used by e.g. <see cref="XDGField"/>
    /// </summary>
    public class FieldStorage : IMatrix, ICloneable {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="_NoOfCells"></param>
        /// <param name="_Nmin"></param>
        /// <param name="_Nmax"></param>
        internal FieldStorage(int _NoOfCells, int _Nmin, int _Nmax) {
            m_NMin = _Nmin;
            m_NoOfCells = _NoOfCells;
            m_NMax = _Nmax;

            m_BaseStorage = new double[m_NMin * m_NoOfCells];
            BaseStorageMda = MultidimensionalArray.CreateWrapper(m_BaseStorage, m_NoOfCells, m_NMin);
            m_ExtendedStorage = new double[NoOfCells][];
        }
              

        int m_NoOfCells;

        int m_NMin;


        /// <summary>
        /// 
        /// </summary>
        public int NMin { get { return m_NMin; } }

        int m_NMax;


        /// <summary>
        /// maximum number of supported DG coordinates
        /// </summary>
        public int NMax { get { return m_NMax; } }


        /// <summary>
        /// initiation of resize operation (<see cref="Resize"/> and <see cref="FinishResize"/>);
        /// </summary>
        /// <param name="NMax_New"></param>
        /// <returns></returns>
        public void BeginResize(int NMax_New) {
            m_jCnt = -1;
            //m_NMin_New = NMin_New;
            m_NMax = NMax_New;

            //if (m_NMin != m_NMin_New)
            //    m_BaseStorageNew = new double[m_NoOfCells * m_NMin_New];
        }

        int m_jCnt = int.MinValue;


        /// <summary>
        /// resizes the total number of XDG coordinates for cell <paramref name="j"/>;
        /// Must be called in between of <see cref="BeginResize"/> and <see cref="FinishResize"/>;
        /// </summary>
        /// <param name="j"></param>
        /// <param name="NewN"></param>
        public void Resize(int j, int NewN) {
            if (j < 0 || j > m_NoOfCells)
                throw new IndexOutOfRangeException("param j");
            if (j != (m_jCnt + 1))
                throw new ArgumentException("Resize must be called with j going sequentially from 0 to NumberOfRows-1", "j");
            m_jCnt++;

            if (NewN == m_NMin) {
                m_ExtendedStorage[j] = null;
            } else {
                int ln = NewN - m_NMin;
                double[] extOld = m_ExtendedStorage[j];
                if (extOld == null) {
                    if (ln > 0)
                        m_ExtendedStorage[j] = new double[ln];
                } else {
                    if (ln <= 0)
                        m_ExtendedStorage[j] = null;
                    else if (ln != extOld.Length)
                        Array.Resize<double>(ref m_ExtendedStorage[j], ln);
                }
            }
        }


        /// <summary>
        /// finalization of resize operation (<see cref="Resize"/> and <see cref="BeginResize"/>);
        /// </summary>
        public void FinishResize() {
            if (m_jCnt != (m_NoOfCells - 1))
                throw new ApplicationException("Resize must be called with j going sequentially from 0 to NumberOfRows-1");
            m_jCnt = int.MinValue;
        }


        /// <summary>
        /// initializes this to be a non-shallow copy of <paramref name="other"/>;
        /// </summary>
        internal void CopyFrom(FieldStorage other) {
            if (other.m_NMin != this.m_NMin || other.m_NoOfCells != this.m_NoOfCells)
                throw new ArgumentException();

            int J = m_NoOfCells;
            int NMin = m_NMin;

            Array.Copy(other.m_BaseStorage, this.m_BaseStorage, other.m_BaseStorage.Length);

            this.BeginResize(other.m_NMax);
            for (int j = 0; j < J; j++) {
                double[] exto = other.m_ExtendedStorage[j];

                this.Resize(j, (exto != null) ? (exto.Length + NMin) : 0);
            }
            this.FinishResize();
            this.m_NMax = other.m_NMax;


            for (int j = 0; j < J; j++) {
                double[] exto = other.m_ExtendedStorage[j];
                double[] extt = this.m_ExtendedStorage[j];


                if(exto != null)
                    Array.Copy(exto, extt, exto.Length);
            }
        }


        /// <summary>
        /// 
        /// </summary>
        public int NoOfCells { get { return m_NoOfCells; } }


        /// <summary>
        /// shallow copy of <see cref="m_BaseStorage"/>.
        /// </summary>
        public MultidimensionalArray BaseStorageMda {
            get;
            private set;
        }

        
        internal double[] m_BaseStorage;

        internal double[][] m_ExtendedStorage;

        /// <summary>
        /// returns the number of DG coordinates in cell <paramref name="j"/>
        /// </summary>
        /// <param name="j"></param>
        /// <returns></returns>
        public int GetNoOfCoordinates(int j) {
            int r = m_NMin;
            if (m_ExtendedStorage[j] != null)
                r += m_ExtendedStorage[j].Length;
            return r;
        }

        //public int AllocatedLength(int j) {
        //    int r = m_NMin;
        //    if (m_ExtendedStorage[j] != null)
        //        r += m_ExtendedStorage[j].Length;
        //    return r;
        //}


        #region IMatrix Members
      
        /// <summary>
        /// sets all entries to 0.0;
        /// </summary>
        public void Clear() {
            Array.Clear(m_BaseStorage, 0, m_BaseStorage.Length);
            foreach (double[] ext in m_ExtendedStorage) {
                if (ext != null) {
                    Array.Clear(ext, 0, ext.Length);
                }
            }
        }
        
        /// <summary>
        /// equal to <see cref="NMax"/>
        /// </summary>
        public int NoOfCols {
            get { return NMax;  }
        }

        /// <summary>
        /// equal to <see cref="NoOfCells"/>
        /// </summary>
        public int NoOfRows {
            get { return NoOfCells; }
        }


        /// <summary>
        /// computes the index into the <see cref="m_BaseStorage"/>-array,
        /// for some (<paramref name="j"/>,<paramref name="n"/>)-pair;
        /// </summary>
        /// <param name="j"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        internal int IndBase(int j, int n) {
            if (j < 0 || j >= m_NoOfCells)
                throw new IndexOutOfRangeException();
            if (n < 0 || n >= m_NMin)
                throw new IndexOutOfRangeException();

            return (j * m_NMin + n);
        }

        /// <summary>
        /// sets/get an entry;
        /// "Getting" is supported for all entries (in the matrix) (non-allocated entries are returned as 0.0),
        /// "Setting" is only supported for allocated entries;
        /// </summary>
        /// <param name="j"></param>
        /// <param name="n">DG coordinate index;
        /// For the "get" operation values between 0 (including) 
        /// and <see cref="NMax"/> (excluding) are valid.
        /// For the "set" operation, only values between o (including) 
        /// and <see cref="GetNoOfCoordinates"/>(<paramref name="j"/>) (excluding) are valid.
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// setting an entry for which no memory is allocated
        /// will throw an exception;
        /// </remarks>
        public double this[int j, int n] {
            get {
                if (j < 0 || j >= m_NoOfCells)
                    throw new IndexOutOfRangeException();
                if (n < 0 || n >= m_NMax)
                    throw new IndexOutOfRangeException();

                if (n < m_NMin) {
                    return m_BaseStorage[IndBase(j, n)];
                } else {
                    int nExt = n - m_NMin;

                    double[] ext = m_ExtendedStorage[j];
                    if (ext == null || nExt >= ext.Length)
                        return 0.0;
                    else
                        return ext[nExt];
                }
            }
            set {
                if (j < 0 || j >= m_NoOfCells)
                    throw new IndexOutOfRangeException();
                if (n < 0 || n >= m_NMax)
                    throw new IndexOutOfRangeException();

                if (n < m_NMin) {
                    m_BaseStorage[IndBase(j, n)] = value;
                } else {
                    int nExt = n - m_NMin;

                    double[] ext = m_ExtendedStorage[j];
                    if (ext == null || nExt >= ext.Length) {
                        return; // send value to nirvana!
                        //throw new IndexOutOfRangeException("trying to set non-allocated entry");
                    }
                    //} else {
                    //    //if (nExt >= ext.Length)
                    //    //    throw new IndexOutOfRangeException("trying to set non-allocated entry");
                    //}
                        
                    m_ExtendedStorage[j][nExt] = value;    
                }
            }
        }

        #endregion

        #region IMatrix Members

        /// <summary>
        /// accumulates <paramref name="a"/>*<paramref name="M"/> to this matrix
        /// </summary>
        public void Acc(double a, IMatrix M) {
            if (M.NoOfCols != this.NoOfCols || M.NoOfRows != this.NoOfRows)
                throw new ArgumentException("M must have same size and length as this matrix");
            int _NoOfRows = this.NoOfRows;
            int _NoOfCols = this.NoOfCols;

            for (int i = 0; i < _NoOfRows; i++)
                for (int j = 0; j < _NoOfCols; j++)
                    this[i, j] += a * M[i, j];
        }

        /// <summary>
        /// multiplies all entries by number <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        public void Scale(double a) {
            BLAS.dscal(m_BaseStorage.Length, a, m_BaseStorage, 1);

            for (int j = m_ExtendedStorage.Length - 1; j >= 0; j--) {
                double[] st = m_ExtendedStorage[j];
                if (st != null) {
                    BLAS.dscal(st.Length, a, st, 1);
                }
            }
        }

        #endregion

        #region ICloneable Members

        /// <summary>
        /// creates a non-shallow copy of this object;
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            FieldStorage r = new FieldStorage(this.NoOfCells, m_NMin, m_NMax);

            Array.Copy(this.m_BaseStorage, r.m_BaseStorage, this.m_BaseStorage.Length);
            for (int j = 0; j < m_NoOfCells; j++) {
                double[] st = m_ExtendedStorage[j];
                if (st != null) {
                    r.m_ExtendedStorage[j] = (double[])m_ExtendedStorage[j].Clone();
                }
            }

            return r;
        }

        #endregion

        #region IMatrix Members

        /// <summary>
        /// Copies values from this array into a 2-dimensional .NET-array;
        /// </summary>
        /// <param name="dest">
        /// target for the copy operation;
        /// the [0,0]-entry of this matrix will end up in the [<paramref name="i0"/>,<paramref name="i1"/>]-entry
        /// of <paramref name="dest"/>;
        /// </param>
        /// <param name="i0">offset into this array, 1st index;</param>
        /// <param name="i1">offset into this array, 2nd index;</param>
        public void CopyTo(double[,] dest, int i0, int i1) {

            int I = this.NoOfRows, J = this.NoOfCols;
            for( int i = 0; i < I; i++)
                for (int j = 0; j < J; j++) {
                    dest[i + i0, j + i0] = this[i, j];
                }

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
        public void CopyTo<T>(T array, bool RowWise, int arrayoffset) where T : IList<double> {

            int I = this.NoOfRows, J = this.NoOfCols;
            if (RowWise) {

                // loop over rows...
                for (int i = 0; i < I; i++)
                    // loop over columns..
                    for (int j = 0; j < J; j++) {
                        array[arrayoffset] = this[i, j];
                        arrayoffset++;
                    }
            } else {
                // loop over columns...
                for (int j = 0; j < J; j++)
                    // loop over rows...
                    for (int i = 0; i < I; i++) {
                        array[arrayoffset] = this[i, j];
                        arrayoffset++;
                    }
            }
        }

        #endregion
    }
}
