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

namespace ilPSP.Utils {


    /// <summary>
    /// some static functions, which should help manipulating matrixes and vectors;
    /// not very high-performant, but helpful;
    /// </summary>
    static public class ArrayTools {

        /// <summary>
        /// Deletes the <paramref name="i"/>-th entry from <paramref name="A"/>.
        /// </summary>
        public static void RemoveAt<T>(ref T[] A, int i) {
            if (i < 0 || i >= A.Length) {
                throw new IndexOutOfRangeException();
            }
            T[] newA = new T[A.Length - 1];
            int j = 0;
            int k = 0;
            while (j < A.Length) {
                if (j != i) {
                    newA[k] = A[j];
                    k++;
                }
                j++;
            }
            A = newA;
        }

        


        /// <summary>
        /// returns an array of length <paramref name="l"/>, with each entry set to <paramref name="c"/>.
        /// </summary>
        public static double[] Const(int l, double c = 1.0) {
            double[] r = new double[l];
            r.SetAll(c);
            return r;
        }


        /// <summary>
        /// Resize vector (1-dimensional array) to a matrix (2-dimensional array),
        /// version without memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="outp">
        /// </param>
        /// <param name="ColumnByColumn">
        /// if false, values from <paramref name="inp"/> are copied row-by-row (
        /// i.e. output contains the first row of <paramref name="inp"/>,
        /// second row of <paramref name="inp"/>, and so on (2nd index rotationg fastest,
        /// C-order, ...));
        /// if true, values from <paramref name="inp"/> are copied column-by-column (
        /// i.e. output contains the first column of <paramref name="inp"/>,
        /// second column of <paramref name="inp"/>, and so on (1st index rotationg fastest,
        /// FORTRAN order, ...));
        /// </param>
        public static void Resize<T>(this T[] inp, T[,] outp, bool ColumnByColumn) {


            int cnt = 0;
            int I1 = outp.GetLength(0);
            int I2 = outp.GetLength(1);

            if (I1 * I2 < inp.Length)
                throw new ArgumentException("output array is too small;", "outp");

            if (!ColumnByColumn) {
                for (int i = 0; i < I1; i++) {
                    for (int j = 0; j < I2; j++) {
                        outp[i, j] = inp[cnt];
                        cnt++;
                    }
                }
            } else {
                for (int j = 0; j < I2; j++) {
                    for (int i = 0; i < I1; i++) {
                        outp[i, j] = inp[cnt];
                        cnt++;
                    }
                }
            }
        }


        /// <summary>
        /// Resize vector (1-dimensional array) to a matrix (2-dimensional array),
        /// version wit memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        /// <param name="RetNoOfColumns">
        /// length of first index of return value;
        /// </param>
        /// <param name="RetNoOfRows">
        /// length of second index of return value;
        /// </param>
        /// <param name="ColumnByColumn">
        /// if false, values from <paramref name="inp"/> are copied row-by-row (
        /// i.e. output contains the first row of <paramref name="inp"/>,
        /// second row of <paramref name="inp"/>, and so on (2nd index rotationg fastest,
        /// C-order, ...));
        /// if true, values from <paramref name="inp"/> are copied column-by-column (
        /// i.e. output contains the first column of <paramref name="inp"/>,
        /// second column of <paramref name="inp"/>, and so on (1st index rotationg fastest,
        /// FORTRAN order, ...));
        /// </param>
        public static T[,] Resize<T>(this T[] inp, bool ColumnByColumn, int RetNoOfRows, int RetNoOfColumns) {
            T[,] ret = new T[RetNoOfRows, RetNoOfColumns];
            Resize(inp, ret, ColumnByColumn);
            return ret;
        }




        /// <summary>
        /// Resize  matrix (2-dimensional array) to a vector (1-dimensional array),
        /// version without memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="output"></param>
        /// <param name="ColumnByColumn">
        /// if false, values from <paramref name="inp"/> are copied row-by-row (
        /// i.e. output contains the first row of <paramref name="inp"/>,
        /// second row of <paramref name="inp"/>, and so on (2nd index rotationg fastest,
        /// C-order, ...));
        /// if true, values from <paramref name="inp"/> are copied column-by-column (
        /// i.e. output contains the first column of <paramref name="inp"/>,
        /// second column of <paramref name="inp"/>, and so on (1st index rotationg fastest,
        /// FORTRAN order, ...));
        /// </param>
        public static void Resize<T>(this T[,] inp, T[] output, bool ColumnByColumn) {

            int I1 = inp.GetLength(0);
            int I2 = inp.GetLength(1);
            if (output.Length < I1 * I2)
                throw new ArgumentException("output array to short", "output");

            int cnt = 0;
            if (!ColumnByColumn) {
                for (int i = 0; i < I1; i++)
                    for (int j = 0; j < I2; j++) {
                        output[cnt] = inp[i, j];
                        cnt++;
                    }
            } else {
                for (int j = 0; j < I2; j++)
                    for (int i = 0; i < I1; i++) {
                        output[cnt] = inp[i, j];
                        cnt++;
                    }
            }

        }


        /// <summary>
        /// Resize  matrix (2-dimensional array) to a vector (1-dimensional array),
        /// version with memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="ColumnByColumn">
        /// if false, values from <paramref name="inp"/> are copied row-by-row (
        /// i.e. output contains the first row of <paramref name="inp"/>,
        /// second row of <paramref name="inp"/>, and so on (2nd index rotationg fastest,
        /// C-order, ...));
        /// if true, values from <paramref name="inp"/> are copied column-by-column (
        /// i.e. output contains the first column of <paramref name="inp"/>,
        /// second column of <paramref name="inp"/>, and so on (1st index rotationg fastest,
        /// FORTRAN order, ...));
        /// </param>
        /// <returns>
        /// an array of equal total number of elements like <paramref name="inp"/>
        /// </returns>
        public static T[] Resize<T>(this T[,] inp, bool ColumnByColumn) {
            T[] ret = new T[inp.GetLength(0) * inp.GetLength(1)];
            Resize(inp, ret, ColumnByColumn);
            return ret;
        }


        /// <summary>
        /// copies <paramref name="col"/> to the <paramref name="colind"/>-th 
        /// column of <paramref name="mtx"/>
        /// </summary>
        /// <param name="col">
        /// length must be equal to 1st length of <paramref name="mtx"/>
        /// </param>
        /// <param name="colind">
        /// column index (1nd index) int <paramref name="mtx"/>, where <paramref name="col"/> should be inserted;
        /// </param>
        /// <param name="mtx"></param>
        public static void SetColumn<T,V>(this T[,] mtx, V col, int colind)
            where V : IList<T>    
        {

            int I = mtx.GetLength(0);

            if (I != col.Count)
                throw new ArgumentException("length of col must be equal to 1st length of mtx", "col");

            for (int i = 0; i < I; i++)
                mtx[i, colind] = col[i];
        }

        /// <summary>
        /// copies <paramref name="row"/> to the <paramref name="rowind"/>-th 
        /// row of <paramref name="mtx"/>
        /// </summary>
        /// <param name="row">
        /// length must be equal to 2nd length of <paramref name="mtx"/>
        /// </param>
        /// <param name="rowind">
        /// row index (2nd index) into <paramref name="mtx"/>, where <paramref name="row"/> should be inserted;
        /// </param>
        /// <param name="mtx"></param>
        public static void SetRow<T, V>(this T[,] mtx, int rowind, V row) 
            where V : IList<T>
        {

            int J = mtx.GetLength(1);

            if (J != row.Count)
                throw new ArgumentException("length of row must be equal to 2nd length of mtx", "row");

            for (int j = 0; j < J; j++)
                mtx[rowind, j] = row[j];
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
        public static void GetColumn<T>(this T[,] inp, int ColNo, T[] outp) {

            if (ColNo < 0 || ColNo >= inp.GetLength(1))
                throw new IndexOutOfRangeException("ColNo out of range");
            if (outp.Length < inp.GetLength(0))
                throw new ArgumentException("output array to short", "outp");

            int I1 = inp.GetLength(0);
            for (int i = 0; i < I1; i++)
                outp[i] = inp[i, ColNo];
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
        public static T[] GetColumn<T>(this T[,] inp, int ColNo) {
            T[] ret = new T[inp.GetLength(0)];
            GetColumn(inp, ColNo, ret);
            return ret;
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
        public static void GetRow<T>(this T[,] inp, int RowNo, T[] outp) {

            if (RowNo < 0 || RowNo >= inp.GetLength(0))
                throw new IndexOutOfRangeException("RowNo out of range");
            if (outp.Length < inp.GetLength(1))
                throw new ArgumentException("output array to short", "outp");

            int I1 = inp.GetLength(1);
            for (int i = 0; i < I1; i++)
                outp[i] = inp[RowNo, i];
        }

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
        public static T[] GetRow<T>(this T[,] inp, int RowNo) {
            T[] ret = new T[inp.GetLength(1)];
            GetRow(inp, RowNo, ret);
            return ret;
        }



        /// <summary>
        /// transposes a matrix, input and output matrixes are different
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="outp"></param>
        /// <returns></returns>
        public static void Transpose<T>(T[,] inp, T[,] outp) {
            if (outp.GetLength(0) != inp.GetLength(1))
                throw new ArgumentException("1st length of inp must be equal to 2nd length of outp");
            if (outp.GetLength(1) != inp.GetLength(0))
                throw new ArgumentException("2nd length of inp must be equal to 1st length of outp");

            int I = outp.GetLength(0);
            int J = outp.GetLength(1);
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    outp[i, j] = inp[j, i];
        }

        /// <summary>
        /// transposes a quadratic matrix, input and output matrixes are the same array
        /// </summary>
        /// <param name="mtx"></param>
        public static void Transpose<T>(this T[,] mtx) {
            if (mtx.GetLength(0) != mtx.GetLength(1))
                throw new ArgumentException("matrix has to be quadratic", "mtx");

            int I = mtx.GetLength(0);
            for (int i = 0; i < I; i++) {
                for (int j = i + 1; j < I; j++) {
                    T buf = mtx[i, j];
                    mtx[i, j] = mtx[j, i];
                    mtx[j, i] = buf;
                }

            }

        }

        /// <summary>
        /// Extracts a sub-matrix from a matrix,
        /// version without memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="RowInd">row indices which sould be taken from <paramref name="inp"/></param>
        /// <param name="ColInd">column indices which should be taken from <paramref name="inp"/></param>
        /// <param name="outp"></param>
        public static void GetSubMatrix<T>(this T[,] inp, int[] RowInd, int[] ColInd, T[,] outp) {
            if (outp.GetLength(0) != RowInd.Length)
                throw new ArgumentException("mismatch between number of rows in outp and length of RowInd-array", "RowInd,outp");

            int I1 = RowInd.Length;
            int I2 = ColInd.Length;

            for (int i = 0; i < I1; i++) {
                for (int j = 0; j < I2; j++) {
                    outp[i, j] = inp[RowInd[i], ColInd[j]];
                }
            }
        }

        /// <summary>
        /// extracts a sub-matrix from a matrix,
        /// version with memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="RowInd">row indices which should be taken from <paramref name="inp"/></param>
        /// <param name="ColInd">column indices which should be taken from <paramref name="inp"/></param>
        /// <returns>
        /// </returns>
        public static T[,] GetSubMatrix<T>(this T[,] inp, int[] RowInd, int[] ColInd) {
            T[,] outp = new T[RowInd.Length, ColInd.Length];
            GetSubMatrix(inp, RowInd, ColInd, outp);
            return outp;
        }


        /// <summary>
        /// extracts a submatrix from a matrix,
        /// version without memory allocation;
        /// Number of rows an columns that should be copied is determined by the size of the
        /// <paramref name="outp"/>-Matrix
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="RowStart">row to start with coping;</param>
        /// <param name="ColStart">column to start with coping;</param>
        /// <param name="outp"></param>
        public static void GetSubMatrix<T>(this T[,] inp, int RowStart, int ColStart, T[,] outp) {
            int I1 = outp.GetLength(0);
            int I2 = outp.GetLength(1);

            for (int i = 0; i < I1; i++)
                for (int j = 0; j < I2; j++) {
                    outp[i, j] = inp[RowStart + i, ColStart + j];
                }
        }

        /// <summary>
        /// Extracts a part of a vector.
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="i0">start index</param>
        /// <param name="Length">number of elements to copy</param>
        /// <returns>
        /// an array containing the entries <paramref name="i0"/> (including) to <paramref name="i0"/>+<paramref name="Length"/> (excluding)
        /// of the input <paramref name="inp"/>
        /// </returns>
        public static T[] GetSubVector<T>(this IList<T> inp, int i0, int Length) {
            T[] ret = new T[Length];
            for (int i = 0; i < Length; i++) {
                ret[i] = inp[i0 + i];
            }
            return ret;
        }

        



        /// <summary>
        /// Extracts a sub-vector from <paramref name="inp"/>.
        /// The <paramref name="idxSrc"/>[i]-th entry of <paramref name="inp"/>
        /// is copied to
        /// the <paramref name="idxTarg"/>[i]-th entry of the return value
        /// </summary>
        /// <param name="inp">source.</param>
        /// <param name="idxSrc">list of indices, to pich values from the source array</param>
        /// <param name="idxTarg">
        /// Optional, can be null; if this is the case, it is assumed to be {0,1,2,3, ... }.
        /// </param>
        /// <param name="Lout">
        /// alternative length specifier for the output array; ignored if negative;
        /// </param>
        /// <returns></returns>
        public static O[] GetSubVector<I1, I2, O>(this IList<O> inp, I1 idxSrc, I2 idxTarg = default(I2), int Lout = -1) 
            where I1: IList<int>
            where I2: IList<int>
        {
            if(Lout < 0) {
                if(idxTarg != null)
                    Lout = idxTarg.Max();
                else
                    Lout = idxSrc.Count;
            }
            O[] ret = new O[Lout];
            GetSubVector(inp, ret, idxSrc, idxTarg);
            return ret;
        }
        

        /// <summary>
        /// Extracts a sub-vector from <paramref name="inp"/>.
        /// The <paramref name="idxSrc"/>[i]-th entry of <paramref name="inp"/>
        /// is copied to
        /// the <paramref name="idxTarg"/>[i]-th entry of <paramref name="output"/>.
        /// </summary>
        /// <param name="inp">source.</param>
        /// <param name="output">destination.</param>
        /// <param name="idxSrc">list of indices, to pich values from the source array</param>
        /// <param name="idxTarg">
        /// Optional, can be null; if this is the case, it is assumed to be {0,1,2,3, ... }.
        /// </param>
        public static void GetSubVector<I1, I2, O>(this IList<O> inp, IList<O> output, I1 idxSrc, I2 idxTarg = default(I2)) 
            where I1: IList<int>
            where I2: IList<int>
//            where O: struct
        {
            if (idxSrc.Count > output.Count)
                throw new ArgumentException("Cannot copy more elements than length of output array.");
                        
            if (idxTarg != null) {
                if (idxTarg.Count != idxSrc.Count)
                    throw new ArgumentException("If idxTarg is provided, its lenth must match the length of idxSrc", "idxTarg");
                
                if (idxTarg.Count > output.Count)
                    throw new ArgumentException();
            }

            if (idxTarg == null) {
                int L = output.Count;
                for (int l = 0; l < L; l++) {
                    int iTarg;
                    iTarg = l;
                    int iScr = idxSrc[l];

                    output[iTarg] = inp[iScr];
                }
            }
            else {
                int L = idxTarg.Count;
                for (int l = 0; l < L; l++) {
                    int iTarg;
                        iTarg = idxTarg[l];

                    int iScr = idxSrc[l];

                    output[iTarg] = inp[iScr];
                }
            }
        }

        /// <summary>
        /// Sets a specific portion of a 2D-array.
        /// </summary>
        /// <param name="target">copy destination</param>
        /// <param name="RowStart">offset into the 1st index of <paramref name="target"/></param>
        /// <param name="ColStart">offset into the 2nd index of <paramref name="target"/></param>
        /// <param name="src">copy source</param>
        public static void SetSubMatrix<T>(this T[,] target, int RowStart, int ColStart, T[,] src) {
            for (int i = 0; i < src.GetLength(0); i++) {
                for (int j = 0; j < src.GetLength(1); j++) {
                    target[i + RowStart, j + ColStart] = src[i, j];
                }
            }
        }



        /// <summary>
        /// extracts a submatrix from a matrix,
        /// version with memory allocation;
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="RowStart">row to start with coping;</param>
        /// <param name="ColStart">column to start with coping;</param>
        /// <param name="NoOfCols">
        /// no of columns to copy;
        /// </param>
        /// <param name="NoOfRows">
        /// no of rows to copy;
        /// </param>
        /// <returns>
        /// a matrix of size <paramref name="NoOfRows"/>x<paramref name="NoOfCols"/>, containing
        /// the specified submatrix
        /// </returns>
        public static T[,] GetSubMatrix<T>(this T[,] inp, int RowStart, int NoOfRows, int ColStart, int NoOfCols) {
            T[,] outp = new T[NoOfRows, NoOfCols];
            GetSubMatrix(inp, RowStart, ColStart, outp);
            return outp;
        }

        /// <summary>
        /// concatenates two vectors <paramref name="a"/> and <paramref name="b"/>
        /// </summary>
        public static T[] Cat<T>(this IList<T> a, IList<T> b) {
            T[] ret = new T[a.Count + b.Count];
            a.CopyTo(ret, 0);
            b.CopyTo(ret, a.Count);
            return ret;
        }

        /// <summary>
        /// concatenates two vectors <paramref name="a"/> and <paramref name="b"/>
        /// </summary>
        public static T[] Cat<T>(this T[] a, T[] b) {
            T[] ret = new T[a.Length + b.Length];
            Array.Copy(a, ret, a.Length);
            Array.Copy(b, 0, ret, a.Length, b.Length);
            return ret;
        }

        /*
        /// <summary>
        /// concatenates two or more arrays
        /// </summary>
        public static T[] Cat<T>(params IList<T>[] some) {
            int L = 0;
            for( int i = 0; i < some.Length; i++)
                L+=some.Length;

            T[] ret = new T[L];

            int i0 = 0;
            for (int i = 0; i < some.Length; i++) {
                some[i].CopyTo(ret, i0);
                i0 += some[i].Count;
            }

            return ret;
        }
        */

        /*
        /// <summary>
        /// concatenates list, arrays and single objects to an array of type <typeparamref name="T"/>
        /// </summary>
        public static T[] Cat<T>(T erstes, params object[] andere) {
            return Cat<T>(new T[] { erstes }, andere);
        }*/

        /// <summary>
        /// concatenates list, arrays and single objects to an array of type <typeparamref name="T"/>
        /// </summary>
        public static T[] Cat<T>(System.Collections.IEnumerable erstes, params object[] andere) {
            return Cat<T>(Cat<T>(new T[0], erstes), andere);
        }

        /// <summary>
        /// concatenates list, arrays and single objects to an array of type <typeparamref name="T"/>
        /// </summary>
        public static T[] Cat<T>(this IEnumerable<T> erstes, params object[] andere) {
            return Cat_<T>(erstes.ToArray(), andere);
        }

        static T[] Cat_<T>(IList<T> erstes, params object[] andere) {
            List<T> ret = new List<T>();
            int i0 = 0;
            for (; i0 < erstes.Count; i0++)
                ret.Add(erstes[i0]);

            for (int i = 0; i < andere.Length; i++) {

                if (andere[i] is IEnumerable<T>) {
                    IEnumerable<T> ol = andere[i] as IEnumerable<T>;

                    foreach (T o in ol) {
                        ret.Add(o);
                    }
                //} else if (andere[i] is System.Collections.IEnumerable) {
                //    var ol = andere[i] as System.Collections.IEnumerable;

                //    foreach (object o in ol) {
                //        ret.Add((T)o);
                //    }

                } else {
                    ret.Add((T)(andere[i]));
                }
            }

            return ret.ToArray();
        }

        /*

        /// <summary>
        /// vertical concatenation of two vectors, <paramref name="a"/> and <paramref name="b"/>,
        /// of length n into a 2-by-n - matrix;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatVert<T>(T[] a, T[] b) {
            if (a.Length != b.Length)
                throw new ArgumentException("length of vectors must be equal");

            int J = a.Length;
            T[,] ret = new T[2, J];
            for (int j = 0; j < J; j++) {
                ret[0, j] = a[j];
                ret[1, j] = b[j];
            }
            return ret;
        }

        /// <summary>
        /// vertical concatenation of matrix <paramref name="b"/> and
        /// vector <paramref name="a"/>;
        /// the 2nd length of <paramref name="a"/> ("number of columns") must be
        /// equal to the length of <paramref name="b"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatVert<T>(T[,] a, T[] b) {
            if (b.Length != a.GetLength(1))
                throw new ArgumentException("number of columns must be equal for both arrays");

            int I1 = a.GetLength(0);
            int J = a.GetLength(1);

            T[,] ret = new T[1 + I1, J];


            for (int i = 0; i < I1; i++)
                for (int j = 0; j < J; j++)
                    ret[i, j] = a[i, j];

            for (int j = 0; j < J; j++)
                ret[I1, j] = b[j];

            return ret;
        }

        /// <summary>
        /// vertical concatenation of vector <paramref name="a"/> and
        /// matrix <paramref name="b"/>;
        /// the 2nd length of <paramref name="b"/> ("number of columns") must be
        /// equal to the length of <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatVert<T>(T[] a, T[,] b) {
            if (a.Length != b.GetLength(1))
                throw new ArgumentException("number of columns must be equal for both arrays");

            int I2 = b.GetLength(0);
            int J = b.GetLength(1);

            T[,] ret = new T[1 + I2, J];

            for (int j = 0; j < J; j++)
                ret[0, j] = a[j];

            for (int i = 0; i < I2; i++)
                for (int j = 0; j < J; j++)
                    ret[i + 1, j] = b[i, j];

            return ret;
        }
        */


        /// <summary>
        /// concatenates two arrays vertically;
        /// the 2nd length of the two arguments ("number of columns") must be equal;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>
        /// Vertical concatenation of matrices <paramref name="a"/> and <paramref name="b"/>;
        /// </returns>
        public static T[,] CatVert<T>(T[,] a, T[,] b) {
            if (a.GetLength(1) != b.GetLength(1))
                throw new ArgumentException("number of columns must be equal for both arrays");

            int I1 = a.GetLength(0);
            int I2 = b.GetLength(0);
            int J = a.GetLength(1);

            T[,] ret = new T[I1 + I2, J];

            for (int i = 0; i < I1; i++)
                for (int j = 0; j < J; j++)
                    ret[i, j] = a[i, j];
            for (int i = 0; i < I2; i++)
                for (int j = 0; j < J; j++)
                    ret[i + I1, j] = b[i, j];

            return ret;
        }

        /*
        /// <summary>
        /// horizontal concatenation of two vectors, <paramref name="a"/> and <paramref name="b"/>,
        /// of length n int a n-by-2 - matrix;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatHoriz<T>(T[] a, T[] b) {
            if (a.Length != b.Length)
                throw new ArgumentException("length of vectors must be equal");

            int J = a.Length;
            T[,] ret = new T[J, 2];
            for (int j = 0; j < J; j++) {
                ret[j, 0] = a[j];
                ret[j, 1] = b[j];
            }
            return ret;
        }

        /// <summary>
        /// horizontal concatenation of matrix <paramref name="a"/> and
        /// vector <paramref name="b"/>;
        /// the 1st length of <paramref name="a"/> ("number of rows") must be
        /// equal to the length of <paramref name="b"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatHoriz<T>(T[,] a, T[] b) {
            if (b.Length != a.GetLength(0))
                throw new ArgumentException("number of rows must be equal for both arrays");

            int I = a.GetLength(0);
            int J1 = a.GetLength(1);

            T[,] ret = new T[I, J1 + 1];


            for (int i = 0; i < I; i++)
                for (int j = 0; j < J1; j++)
                    ret[i, j] = a[i, j];

            for (int i = 0; i < I; i++)
                ret[i, J1] = b[i];

            return ret;
        }

        /// <summary>
        /// horizontal concatenation of vector <paramref name="a"/> and
        /// matrix <paramref name="b"/>;
        /// the 1st length of <paramref name="b"/> ("number of rows") must be
        /// equal to the length of <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static T[,] CatHoriz<T>(T[] a, T[,] b) {
            if (a.Length != b.GetLength(0))
                throw new ArgumentException("number of rows must be equal for both arrays");

            int I = b.GetLength(0);
            int J2 = b.GetLength(1);

            T[,] ret = new T[I, 1 + J2];

            for (int i = 0; i < I; i++)
                ret[i, 0] = a[i];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J2; j++)
                    ret[i, j + 1] = b[i, j];

            return ret;
        }*/


        /// <summary>
        /// converts a vector into a matrix with one row
        /// </summary>
        public static T[,] ToRowVec<T>(this T[] vec) {
            var R = new T[1, vec.Length];
            for(int i = 0; i < vec.Length; i++) {
                R[0, i] = vec[i];
            }
            return R;
        }

        /// <summary>
        /// converts a vector into a matrix with one column
        /// </summary>
        public static T[,] ToColVec<T>(this T[] vec) {
            var R = new T[vec.Length, 1];
            for(int i = 0; i < vec.Length; i++) {
                R[i, 0] = vec[i];
            }
            return R;
        }

        /// <summary>
        /// concatenates two matrices horizontally;
        /// the 1st length of the two arguments ("number of rows") must be equal;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>
        /// Horizontal concatenation of matrices <paramref name="a"/> and <paramref name="b"/>;
        /// </returns>
        public static T[,] CatHoriz<T>(T[,] a, T[,] b) {
            if (a.GetLength(0) != b.GetLength(00))
                throw new ArgumentException("number of columns must be equal for both arrays");

            int I = a.GetLength(0);
            int J1 = a.GetLength(1);
            int J2 = b.GetLength(1);

            T[,] ret = new T[I, J1 + J2];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J1; j++)
                    ret[i, j] = a[i, j];

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J2; j++)
                    ret[i, j + J1] = b[i, j];

            return ret;
        }



        ///// <summary>
        ///// converts a generic list to an array
        ///// </summary>
        ///// <typeparam name="T">element type of the returned array</typeparam>
        ///// <param name="a"></param>
        ///// <returns></returns>
        //public static T[] List2Array<T>(IList<T> a) {
        //    T[] ret = new T[a.Count];
        //    a.CopyTo(ret, 0);
        //    return ret;
        //}

        ///// <summary>
        ///// converts a generic enumerable to a static array
        ///// </summary>
        ///// <typeparam name="T">element type of the returned array</typeparam>
        ///// <param name="a"></param>
        ///// <returns></returns>
        //public static T[] Enum2Array<T>(System.Collections.IEnumerable a) {
        //    int cnt = 0;
        //    foreach (var entry in a) {
        //        cnt++;
        //    }
        //    T[] ret = new T[cnt];
        //    cnt = 0;
        //    foreach (var entry in a) {
        //        ret[cnt] = (T)entry;
        //        cnt++;
        //    }
        //    return ret;
        //}
        
        ///// <summary>
        ///// converts a generic enumerable to a static array
        ///// </summary>
        ///// <typeparam name="T">element type of the returned array</typeparam>
        ///// <param name="a"></param>
        ///// <returns></returns>
        //public static T[] Enum2Array<T>(IEnumerable<T> a) {
        //    int cnt = 0;
        //    foreach (var entry in a) {
        //        cnt++;
        //    }
        //    T[] ret = new T[cnt];
        //    cnt = 0;
        //    foreach (var entry in a) {
        //        ret[cnt] = entry;
        //        cnt++;
        //    }
        //    return ret;
        //}

        /// <summary>
        /// converts a sublist of a generic list <paramref name="a"/>,
        /// starting at <paramref name="i0"/>, being <paramref name="len"/> elements long,
        /// to an array
        /// </summary>
        /// <typeparam name="T">element type of the returned array</typeparam>
        /// <param name="a"></param>
        /// <param name="i0">index into <paramref name="a"/>, defining the start of the sub-list</param>
        /// <param name="len">length of the sub-list</param>
        /// <returns></returns>
        public static T[] List2Array<T>(IList<T> a, int i0, int len) {
            T[] ret = new T[len];
            for (int i = 0; i < len; i++)
                ret[i] = a[i + i0];
            return ret;
        }

        /// <summary>
        /// initializes all entries of <paramref name="a"/> with value
        /// <paramref name="val"/>;
        /// </summary>
        public static void SetAll<Q>(this Q[,] a, Q val) {

            int Rank = a.Rank;
            int[] ind = new int[Rank];

            SetRecursive(a, val, ind, 0);
        }
        
        ///// <summary>
        ///// initializes all entries of <paramref name="a"/> with value
        ///// <paramref name="val"/>;
        ///// </summary>
        ///// <param name="a"></param>
        ///// <param name="val"></param>
        //public static void SetAll<T>(this T[] a, T val) {
        //    SetAllL<T[], T>(a, val);
        //}

        /// <summary>
        /// initializes all entries of <paramref name="a"/> with value
        /// <paramref name="val"/>;
        /// </summary>
        /// <param name="a"></param>
        /// <param name="val"></param>
        public static void SetAll<T,V>(this T a, V val) where T: IList<V> {

            for (int i = a.Count - 1; i >= 0; i--)
                a[i] = val;
        }

        /// <summary>
        /// Vector Addition
        /// </summary>
        public static double[] Plus(this double[] a, double[] b) {
            if (a.Length != b.Length)
                throw new ArgumentException();
            int D = a.Length;

            double[] R = a.CloneAs();
            R.AccV(1.0, b);
            return R;
        }


        /// <summary>
        /// Vector Difference
        /// </summary>
        public static double[] Minus(this double[] a, double[] b) {
            if (a.Length != b.Length)
                throw new ArgumentException();
            int D = a.Length;

            double[] R = a.CloneAs();
            R.AccV(-1.0, b);
            return R;
        }

        /// <summary>
        /// Multiplication with a scalar
        /// </summary>
        public static double[] Mul(this double[] a, double b) {
            double[] R = a.CloneAs();
            R.ScaleV(b);
            return R;
        }

        /// <summary>
        /// implementation of <see cref="SetAll"/>
        /// </summary>
        static void SetRecursive(Array a, object val, int[] indices, int dim) {
            int rank = a.Rank;

            int l0 = a.GetLowerBound(dim);
            int l1 = a.GetUpperBound(dim);

            if (dim == (rank - 1)) {
                for (indices[dim] = l0; indices[dim] <= l1; indices[dim]++)
                    a.SetValue(val, indices);
            } else {
                for (indices[dim] = l0; indices[dim] <= l1; indices[dim]++)
                    SetRecursive(a, val, indices, dim + 1);
            }
        }

        /// <summary>
        /// Checks two one-dimensional arrays for equality using the default
        /// equality comparer for the given type.
        /// </summary>
        /// <typeparam name="T">
        /// The type of the objects in both arrays
        /// </typeparam>
        /// <param name="first">First array</param>
        /// <param name="second">Second array</param>
        /// <returns>
        /// If any of the operands is null, false is returned. Otherwise, all
        /// elements are checked for equality using <paramref name="EqualityFunc"/> (if provided), otherwise
        /// <see cref="EqualityComparer{T}.Equals"/>
        /// </returns>
        /// <param name="EqualityFunc">
        /// optional comparison function
        /// </param>
        public static bool ListEquals<T>(this IEnumerable<T> first, IEnumerable<T> second, Func<T,T,bool> EqualityFunc = null) {
            if (ReferenceEquals(first, second)) {
                return true;
            }

            if (first == null && second == null) {
                return true;
            }

            if (first == null || second == null) {
                return false;
            }

            if (first.Count() != second.Count()) {
                return false;
            }

            if (EqualityFunc == null) {
                EqualityComparer<T> comp = EqualityComparer<T>.Default;
                EqualityFunc = comp.Equals;
            }

            var Ea = first.GetEnumerator();
            var Eb = second.GetEnumerator();
            while (Ea.MoveNext() && Eb.MoveNext()) {
                var first_i = Ea.Current;
                var second_i = Eb.Current;

                if (!EqualityFunc(first_i, second_i)) {
                    return false;
                }
            }

            if(Ea.MoveNext() || Eb.MoveNext())
                return false;

            return true;
        }

        /// <summary>
        /// Checks two  arrays for equality using  using 
        /// either the default equality comparer for the given type or a custom comparison (<paramref name="EqualityFunc"/>).
        /// </summary>
        /// <typeparam name="T">
        /// The type of the objects in both arrays
        /// </typeparam>
        /// <param name="first">First array</param>
        /// <param name="second">Second array</param>
        /// <returns>
        /// If any of the operands is null, false is returned. Otherwise, all
        /// elements are checked for equality using <paramref name="EqualityFunc"/> (if provided), otherwise
        /// <see cref="EqualityComparer{T}.Equals"/>
        /// </returns>
        /// <param name="EqualityFunc">
        /// optional comparison function
        /// </param>
        public static bool AreEqual<T>(T[] first, T[] second, Func<T, T, bool> EqualityFunc = null) {
            return ListEquals(first, second, EqualityFunc);
        }


        /// <summary>
        /// Checks two two-dimensional arrays for equality using 
        /// either the default equality comparer for the given type or a custom comparison (<paramref name="EqualityFunc"/>).
        /// </summary>
        /// <typeparam name="T">
        /// The type of the objects in both arrays
        /// </typeparam>
        /// <param name="first">First array</param>
        /// <param name="second">Second array</param>
        /// <returns>
        /// If any of the operands is null, false is returned. Otherwise, all
        /// elements are checked for equality using <paramref name="EqualityFunc"/> (if provided), otherwise
        /// <see cref="EqualityComparer{T}.Equals"/>
        /// </returns>
        /// <param name="EqualityFunc">
        /// optional comparison function
        /// </param>
        public static bool AreEqual<T>(T[,] first, T[,] second, Func<T, T, bool> EqualityFunc = null) {
            if (ReferenceEquals(first, second)) {
                return true;
            }

            if (first == null || second == null) {
                return false;
            }

            if (first.GetLength(0) != second.GetLength(0)) {
                return false;
            }

            if (first.GetLength(1) != second.GetLength(1)) {
                return false;
            }

            if (EqualityFunc == null) {
                EqualityComparer<T> comp = EqualityComparer<T>.Default;
                EqualityFunc = comp.Equals;
            }

            for (int i = 0; i < first.GetLength(0); i++) {
                for (int j = 0; j < first.GetLength(1); j++) {
                    if (!EqualityFunc(first[i, j], second[i, j])) {
                        return false;
                    }
                }

            }

            return true;
        }

        /// <summary>
        /// Checks two three-dimensional arrays for equality using  using 
        /// either the default equality comparer for the given type or a custom comparison (<paramref name="EqualityFunc"/>).
        /// </summary>
        /// <typeparam name="T">
        /// The type of the objects in both arrays
        /// </typeparam>
        /// <param name="first">First array</param>
        /// <param name="second">Second array</param>
        /// <returns>
        /// If any of the operands is null, false is returned. Otherwise, all
        /// elements are checked for equality using <paramref name="EqualityFunc"/> (if provided), otherwise
        /// <see cref="EqualityComparer{T}.Equals"/>
        /// </returns>
        /// <param name="EqualityFunc">
        /// optional comparison function
        /// </param>
        public static bool AreEqual<T>(T[, ,] first, T[, ,] second, Func<T, T, bool> EqualityFunc = null) {
            if (ReferenceEquals(first, second)) {
                return true;
            }

            if (first == null || second == null) {
                return false;
            }

            if (first.GetLength(0) != second.GetLength(0)) {
                return false;
            }

            if (first.GetLength(1) != second.GetLength(1)) {
                return false;
            }

            if (first.GetLength(2) != second.GetLength(2)) {
                return false;
            }

            if (EqualityFunc == null) {
                EqualityComparer<T> comp = EqualityComparer<T>.Default;
                EqualityFunc = comp.Equals;
            }


            for (int i = 0; i < first.GetLength(0); i++) {
                for (int j = 0; j < first.GetLength(1); j++) {
                    for (int k = 0; k < first.GetLength(2); k++) {
                        if (!EqualityFunc(first[i, j, k], second[i, j, k])) {
                            return false;
                        }
                    }
                }

            }

            return true;
        }

        /// <summary>
        /// adds some object at the end of an array.
        /// </summary>
        public static void AddToArray<T>(this T e, ref T[] A) {
            if (A == null)
                A = new T[1];
            else
                Array.Resize(ref A, A.Length + 1);

            A[A.Length - 1] = e;
        }
               

        /// <summary>
        /// clones an enumeration and all entries in it.
        /// </summary>
        public static IEnumerable<T> CloneNonshallow<T>(this IEnumerable<T> L) where T : ICloneable {
            return L.Select(e => (T)e.Clone());
        }

        /// <summary>
        /// Clears all entries in <paramref name="A"/>.
        /// </summary>
        public static void Clear<T>(this T[] A) {
            Array.Clear(A, 0, A.Length);
        }
    }
}
