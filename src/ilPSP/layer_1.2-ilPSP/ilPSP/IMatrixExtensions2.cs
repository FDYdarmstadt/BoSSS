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
using ilPSP.LinSolvers;
using System.Collections.Generic;
using ilPSP.Utils;

namespace ilPSP {
	
    /// <summary>
    /// Container for extension methods for the <see cref="IMatrix"/> interface
    /// </summary>
    public static partial class IMatrixExtensions  {

        /// <summary>
        /// Converts the given a matrix to an <see cref="MsrMatrix"/>
        /// </summary>
        /// <param name="M">The matrix to be converted</param>
        /// <returns>
        /// A MsrMatrix with the same content as <paramref name="M"/>
        /// </returns>
        public static MsrMatrix ToMsrMatrix(this IMatrix M) {
            MsrMatrix result = new MsrMatrix(M.NoOfRows, M.NoOfCols, 1, 1);

            for (int i = 0; i < M.NoOfRows; i++) {
                for (int j = 0; j < M.NoOfCols; j++) {
                    result[i, j] = M[i, j];
                }
            }

            return result;
        }


        /// <summary>
        /// extracts the <paramref name="RowNo"/>-th row from
        /// <paramref name="inp"/> and copies it to <paramref name="outp"/>,
        /// version without memory allocation;
        /// </summary>
        /// <param name="inp">
        /// input matrix
        /// </param>
        /// <param name="RowNo">
        /// row which should be extracted
        /// </param>
        /// <param name="outp">
        /// an array with length greater or equal to 1st length of <paramref name="inp"/>;
        /// </param>
        /// <param name="i0">
        /// offset into <paramref name="outp"/>, coping starts from this index
        /// </param>
        /// <param name="inc">
        /// skip between to consecutive elements that are written in <paramref name="outp"/>
        /// </param>
        public static void GetRow<T>(this IMatrix inp, int RowNo, T outp, int i0 = 0, int inc = 1) where T : IList<double> {

            if (RowNo < 0 || RowNo >= inp.NoOfRows)
                throw new IndexOutOfRangeException("RowNo out of range");
            if ((outp.Count - i0) < inp.NoOfCols*inc)
                throw new ArgumentException("output array to short", "outp");

            int I1 = inp.NoOfCols;
            for (int i = 0; i < I1; i++)
                outp[i * inc + i0] = inp[RowNo, i];
        }


        /// <summary>
        /// sets the <paramref name="RowNo"/>-th row from
        /// <paramref name="inp"/> 
        /// to values provided by <paramref name="row"/>.
        /// </summary>
        /// <param name="inp">
        /// matrix that should be altered
        /// </param>
        /// <param name="RowNo">
        /// row index of the row to set
        /// </param>
        /// <param name="row">
        /// an array with length greater or equal to 1st length of <paramref name="inp"/>;
        /// </param>
        /// <param name="ReadOffset">
        /// offset into <paramref name="row"/>, coping starts from this index
        /// </param>
        /// <param name="ReadInc">
        /// skip between to consecutive elements that are taken from <paramref name="row"/>
        /// </param>
        public static void SetRow<T>(this IMatrix inp, int RowNo, T row, int ReadOffset = 0, int ReadInc = 1) where T : IEnumerable<double> {

            if (RowNo < 0 || RowNo >= inp.NoOfRows)
                throw new IndexOutOfRangeException("RowNo out of range");
            //if ((row.Count - i0) < inp.NoOfCols * inc)
            //    throw new ArgumentException("array to short", "row");

            //int I1 = Math.Max(inp.NoOfCols, (row.Count - i0) / inc);
            //for (int i = 0; i < I1; i++)
            //    inp[RowNo, i] = row[i * inc + i0];


            

            int i = 0;
            if (ReadInc == 1 && ReadOffset == 0) {
                foreach (double row_i in row) {
                    inp[RowNo, i] = row_i;
                    i++;
                }
            } else {
                int j = 0;
                foreach (double row_i in row) {
                    if (i >= ReadOffset && ((j - ReadOffset) % ReadInc == 0)) {
                        inp[RowNo, i] = row_i;
                        i++;
                    }
                    j++;
                }
            }
        } 
        
        /*
        /// <summary>
        /// sets the <paramref name="ColNo"/>-th column from
        /// <paramref name="inp"/> 
        /// to values provided by <paramref name="col"/>.
        /// </summary>
        /// <param name="inp">
        /// matrix that should be altered
        /// </param>
        /// <param name="ColNo">
        /// column index of the column to set
        /// </param>
        /// <param name="col">
        /// an array with length greater or equal to 1st length of <paramref name="inp"/>;
        /// </param>
        /// <param name="i0">
        /// offset into <paramref name="col"/>, coping starts from this index
        /// </param>
        /// <param name="inc">
        /// skip between to consecutive elements that are taken from <paramref name="col"/>
        /// </param>
        public static void SetColumn<T>(this IMatrix inp, int ColNo, T col, int i0 = 0, int inc = 1) where T : IList<double> {

            if (ColNo < 0 || ColNo >= inp.NoOfCols)
                throw new IndexOutOfRangeException("RowNo out of range");
            if ((col.Count - i0) < inp.NoOfRows * inc)
                throw new ArgumentException("array to short", "row");

            int I1 = inp.NoOfRows;
            for (int i = 0; i < I1; i++)
                inp[i, ColNo] = col[i * inc + i0];
        }
        */
        


        /// <summary>
        /// extracts the <paramref name="RowNo"/>-th row from
        /// <paramref name="inp"/>.
        /// </summary>
        /// <param name="inp">
        /// input matrix
        /// </param>
        /// <param name="RowNo">
        /// row which should be extracted
        /// </param>
        /// <returns>
        /// an array with length equal to 2nd length of <paramref name="inp"/>, containing the
        /// <paramref name="RowNo"/>-th row of <paramref name="inp"/>
        /// </returns>
        public static double[] GetRow(this IMatrix inp, int RowNo) {
            double[] ret = new double[inp.NoOfCols];
            GetRow(inp, RowNo, ret);
            return ret;
        }


        /// <summary>
        /// extracts <paramref name="N"/> items from the <paramref name="RowNo"/>-th row of matrix <paramref name="inp"/>, starting at column <paramref name="j0"/>
        /// </summary>
        public static double[] GetRowPart(this IMatrix inp, int RowNo, int j0, int N) {
            double[] ret = new double[N];
            for (int i = 0; i < N; i++)
                ret[i] = inp[RowNo, i + j0];
            return ret;
        }


        /// <summary>
        /// extracts the <paramref name="ColNo"/>-th column from
        /// <paramref name="inp"/> and copies it to <paramref name="outp"/>,
        /// version without memory allocation;
        /// </summary>
        /// <param name="inp">
        /// input matrix
        /// </param>
        /// <param name="ColNo">
        /// column index which should be extracted
        /// </param>
        /// <param name="outp">
        /// an array with length greater or equal to 1st length of <paramref name="inp"/>;
        /// </param>
        /// <param name="WriteOffset">
        /// offset into <paramref name="outp"/>, coping starts from this index
        /// </param>
        /// <param name="WriteInc">
        /// skip between to consecutive elements that are written in <paramref name="outp"/>
        /// </param>
        /// <param name="NoOfElm">
        /// Number of elements to take; if negative, ignored and the number of rows is extracted.
        /// </param>
        public static void GetColumn<T>(this IMatrix inp, int ColNo, T outp, int WriteOffset = 0, int WriteInc = 1, int NoOfElm = -1) where T : IList<double> {

            int I1 = NoOfElm >= 0 ? NoOfElm : inp.NoOfRows;

            if (ColNo < 0 || ColNo >= inp.NoOfCols)
                throw new IndexOutOfRangeException("ColNo out of range");
            if (outp.Count - WriteOffset < I1*WriteInc)
                throw new ArgumentException("output array to short", "outp");

            for (int i = 0; i < I1; i++)
                outp[i*WriteInc + WriteOffset] = inp[i, ColNo];
        }


        /// <summary>
        /// extracts the <paramref name="ColNo"/>-th column from
        /// <paramref name="inp"/>.
        /// </summary>
        /// <param name="inp">
        /// input matrix
        /// </param>
        /// <param name="ColNo">
        /// column index which should be extracted
        /// </param>
        /// <returns>
        /// an array with length equal to 1st length of <paramref name="inp"/>, containing the
        /// <paramref name="ColNo"/>-th column of <paramref name="inp"/>
        /// </returns>
        /// <param name="NoOfElm">
        /// Number of elements to take; if negative, ignored and the number of rows is extracted.
        /// </param>        
        public static double[] GetColumn(this IMatrix inp, int ColNo, int NoOfElm = -1) {
            double[] ret = new double[NoOfElm >= 0 ? NoOfElm : inp.NoOfRows];
            GetColumn(inp, ColNo, ret, NoOfElm:NoOfElm);
            return ret;
        }

		/// <summary>
        /// copies all entries of <paramref name="o"/> to <paramref name="inp"/>
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="o"></param>
        public static void CopyFrom(this IMatrix inp, IMatrix o) {
            if (inp.NoOfCols != o.NoOfCols)
                throw new ArgumentException("mismatch in number of columns.", "o");
            if (inp.NoOfRows != o.NoOfRows)
                throw new ArgumentException("mismatch in number of rows.", "o");

            int I = inp.NoOfRows, J = inp.NoOfCols;
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    inp[i, j] = o[i, j];
        }

        /// <summary>
        /// Eigenvalues for a symmetrical matrix
        /// </summary>
        /// <param name="Inp">
        /// On entry some symmetrical matrix, unchanged on exit. 
        /// Symmetry is not checked.
        /// Only the lower matrix triangle is considered, the upper part is ignored.
        /// </param>
        /// <returns></returns>
        public static double[] EigsSymm(this IMatrix Inp) {
            var ret = EigenspaceSymmInternal(Inp, false);
            return ret.EigenVals;
        }


        /// <summary>
        /// Eigenvalues and Eigenvectors of a symmetrical matrix
        /// </summary>
        /// <param name="Inp">
        /// On entry some symmetrical matrix, unchanged on exit. 
        /// Symmetry is not checked.
        /// Only the lower matrix triangle is considered, the upper part is ignored.
        /// </param>
        /// <returns>
        /// Eigenvalues and Eigenvectors
        /// </returns>
        public static (double[] EigenVals, MultidimensionalArray EigenVect) EigenspaceSymm(this IMatrix Inp) {
            var ret = EigenspaceSymmInternal(Inp, true);
            return ret;
        }

        static (double[] EigenVals, MultidimensionalArray EigenVect) EigenspaceSymmInternal(this IMatrix Inp, bool ComputeVectors) {
            if(Inp.NoOfCols != Inp.NoOfRows) {
                throw new ArgumentException("Not supported for non-symmetrical matrices.");
            }
            int N = Inp.NoOfRows;

            int JOBZ = ComputeVectors ? 'V' : 'N';
            // 'N':  Compute eigenvalues only;
            // 'V':  Compute eigenvalues and eigenvectors.

            int UPLO = 'L';
            // 'U':  Upper triangle of A is stored;
            // 'L':  Lower triangle of A is stored.

            unsafe {
                double[] InpBuffer = TempBuffer.GetTempBuffer(out int RefInp, N * N);
                double[] Eigis = new double[N];
                MultidimensionalArray EigiVect = ComputeVectors ? MultidimensionalArray.Create(N, N) : null;

                fixed(double* pInp = InpBuffer, pEigis = Eigis) {
                    CopyToUnsafeBuffer(Inp, pInp, true);

                    int LDA = N;
                    int info;

                    // phase 1: work size estimation
                    double WorkSize;
                    int LWORK = -1; // triggers estimation 
                    LAPACK.F77_LAPACK.DSYEV_(ref JOBZ, ref UPLO, ref N, pInp, ref LDA, pEigis, &WorkSize, ref LWORK, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(RefInp);
                        throw new ArithmeticException("LAPACK DSYEV (symmetrical matrix eigenvalues) returned info " + info);
                    }

                    LWORK = (int)WorkSize;

                    // phase 2: computation
                    double[] WorkBuffer = TempBuffer.GetTempBuffer(out int RefWork, LWORK * 1);
                    fixed(double* pWork = WorkBuffer) {

                        LAPACK.F77_LAPACK.DSYEV_(ref JOBZ, ref UPLO, ref N, pInp, ref LDA, pEigis, pWork, ref LWORK, out info);
                        TempBuffer.FreeTempBuffer(RefWork);
                        if(info != 0) {
                            TempBuffer.FreeTempBuffer(RefInp);
                            throw new ArithmeticException("LAPACK DSYEV (symmetrical matrix eigenvalues) returned info " + info);
                        }

                        if(EigiVect != null)
                            CopyFromUnsafeBuffer(EigiVect, pInp, true);

                    }
                }
                TempBuffer.FreeTempBuffer(RefInp);


                return (Eigis, EigiVect);
            }



        }


    }
}

