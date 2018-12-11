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
        /// <param name="i0">
        /// offset into <paramref name="row"/>, coping starts from this index
        /// </param>
        /// <param name="inc">
        /// skip between to consecutive elements that are taken from <paramref name="row"/>
        /// </param>
        public static void SetRow<T>(this IMatrix inp, int RowNo, T row, int i0 = 0, int inc = 1) where T : IList<double> {

            if (RowNo < 0 || RowNo >= inp.NoOfRows)
                throw new IndexOutOfRangeException("RowNo out of range");
            if ((row.Count - i0) < inp.NoOfCols * inc)
                throw new ArgumentException("array to short", "row");

            int I1 = inp.NoOfCols;
            for (int i = 0; i < I1; i++)
                inp[RowNo, i] = row[i * inc + i0];
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
        /// <param name="i0">
        /// offset into <paramref name="outp"/>, coping starts from this index
        /// </param>
        /// <param name="inc">
        /// skip between to consecutive elements that are written in <paramref name="outp"/>
        /// </param>
        public static void GetColumn<T>(this IMatrix inp, int ColNo, T outp, int i0 = 0, int inc = 1) where T : IList<double> {

            if (ColNo < 0 || ColNo >= inp.NoOfCols)
                throw new IndexOutOfRangeException("ColNo out of range");
            if (outp.Count - i0 < inp.NoOfRows*inc)
                throw new ArgumentException("output array to short", "outp");

            int I1 = inp.NoOfRows;
            for (int i = 0; i < I1; i++)
                outp[i*inc + i0] = inp[i, ColNo];
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
        public static double[] GetColumn(this IMatrix inp, int ColNo) {
            double[] ret = new double[inp.NoOfRows];
            GetColumn(inp, ColNo, ret);
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

    }
}

