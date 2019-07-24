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

using ilPSP.LinSolvers;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

namespace ilPSP {

    /// <summary>
    /// Augmenting 2d <see cref="MultidimensionalArray"/>'s with some matrix operations.
    /// </summary>
    public static partial class IMatrixExtensions {

        /// <summary>
        /// Copies from a multidimensional array to a matrix.<br/>
        /// 1st index of multidimensional array is interpreted as row;<br/>
        /// 2nd index of multidimensional array is interpreted as column
        /// </summary>
        /// <param name="m">a multidimensional array with exactly two dimensions</param>
        /// <param name="M">matrix</param>
        static public void Set(this IMatrix M, MultidimensionalArray m) {
            if(m.Dimension != 2)
                throw new ArgumentException("m.Dimension != 2", "m");
            if(m.GetLength(0) != M.NoOfRows)
                throw new ArgumentException("m.GetLength(0) != _NoOfRows", "m");
            if(m.GetLength(1) != M.NoOfCols)
                throw new ArgumentException("m.GetLength(1) != _NoOfCols", "m");


            for(int i = 0; i < M.NoOfRows; i++)
                for(int j = 0; j < M.NoOfCols; j++)
                    M[i, j] = m[i, j];

        }

        /// <summary>
        /// writes this matrix into a text file;
        /// </summary>
        /// <param name="name">path to the text file</param>
        /// <param name="M">matrix</param>
        /// <param name="fm">file creation mode</param>
        static public void SaveToTextFile(this IMatrix M, string name, FileMode fm = FileMode.Create) {
            NumberFormatInfo nfi = NumberFormatInfo.InvariantInfo;

            using(var txt = new StreamWriter(new FileStream(name, fm))) {
                SaveToStream(M, txt);
                txt.Flush();
            }
        }

        /// <summary>
        /// Writes this matrix to a stream, in text format.
        /// </summary>
        /// <param name="txt">path to the text file</param>
        /// <param name="M">matrix</param>
        static public void SaveToStream(this IMatrix M, TextWriter txt) {
            NumberFormatInfo nfi = NumberFormatInfo.InvariantInfo;

            
            int _NoOfRows = M.NoOfRows;
            int _NoOfCols = M.NoOfCols;

            for(int i = 0; i < _NoOfRows; i++) {
                for(int j = 0; j < _NoOfCols; j++) {

                    txt.Write(M[i, j].ToString(nfi));

                    if(j < _NoOfCols - 1)
                        txt.Write(" ");

                }
                txt.WriteLine();
            }

            txt.Flush();

        }


        /// <summary>
        /// reads a matrix from a text file;
        /// </summary>
        /// <param name="name">path to the text file</param>
        /// <param name="M">pre-allocated matrix to store the loaded data</param>
        static public void LoadFromTextFile(this IMatrix M, string name) {
            using(StreamReader txt = new StreamReader(name)) {
                try {
                    LoadFromStream(M, txt);
                } finally {
                    txt.Close();
                }
            }
        }

        /// <summary>
        /// reads a matrix from a text file;
        /// </summary>
        /// <param name="txt">Text reade on the stream to load from.</param>
        /// <param name="M">pre-allocated matrix to store the loaded data</param>
        static public void LoadFromStream(this IMatrix M, TextReader txt) {
            NumberFormatInfo nfi = NumberFormatInfo.InvariantInfo;

            char[] splits = new char[] { ' ', '\t', '\n', '\r', '\t', ',', ';' };
            
            for(int i = 0; i < M.NoOfRows; i++) {
                var line = txt.ReadLine();
                if(line == null)
                    throw new NotSupportedException("mismatch in number of rows");
                var items = line.Split(splits, StringSplitOptions.RemoveEmptyEntries);

                if(items.Length != M.NoOfCols)
                    throw new NotSupportedException(string.Format("trying to load matrix with {0} columns, but this matrix has {1} columns", items.Length, M.NoOfCols));

                for(int j = 0; j < M.NoOfCols; j++) {
                    try {
                        M[i, j] = double.Parse(items[j], nfi);
                    } catch(FormatException fe) {
                        var itm = items[j].ToLowerInvariant();

                        if(itm.Contains("nan"))
                            M[i, j] = double.NaN;
                        else if(itm.Contains("+inf"))
                            M[i, j] = double.PositiveInfinity;
                        else if(itm.Contains("-inf"))
                            M[i, j] = double.NegativeInfinity;
                        else if(itm.Contains("inf"))
                            M[i, j] = double.PositiveInfinity;
                        else
                            throw fe;
                    }
                }
            }
            if(txt.ReadLine() != null)
                throw new NotSupportedException("mismatch in number of rows");
        }

        /// <summary>
        /// Reads a matrix from a text file.
        /// </summary>
        /// <param name="name">path to the text file.</param>
        static public MultidimensionalArray LoadFromTextFile(string name) {
            using(StreamReader txt = new StreamReader(name)) {
                try {
                    return LoadFromStream(txt);
                } finally {
                    txt.Close();
                }
            }
        }
        
        /// <summary>
        /// Reads a matrix from a text stream.
        /// </summary>
        /// <param name="txt"></param>
        static public MultidimensionalArray LoadFromStream(TextReader txt) {
            NumberFormatInfo nfi = NumberFormatInfo.InvariantInfo;

            char[] splits = new char[] { ' ', '\t', '\n', '\r', '\t', ',', ';' };

            List<double[]> tempMatrix = new List<double[]>();

            int NoOfCols = -1;
            int iLine = -1;
            for(string line = txt.ReadLine(); line != null; line = txt.ReadLine()) {
                iLine++;
                if(line == null)
                    throw new NotSupportedException("mismatch in number of rows");
                var items = line.Split(splits, StringSplitOptions.RemoveEmptyEntries);

                if(NoOfCols < 0) {
                    NoOfCols = items.Length;
                } else {
                    if(NoOfCols != items.Length)
                        throw new NotSupportedException(string.Format("Found {0} columns in line {1}, but first line has {2} columns.", items.Length, iLine, NoOfCols));
                }

                double[] M = new double[items.Length];
                for(int j = 0; j < items.Length; j++) {
                    try {
                        M[j] = double.Parse(items[j], nfi);
                    } catch(FormatException fe) {
                        var itm = items[j].ToLowerInvariant();

                        if(itm.Contains("nan"))
                            M[j] = double.NaN;
                        else if(itm.Contains("+inf"))
                            M[j] = double.PositiveInfinity;
                        else if(itm.Contains("-inf"))
                            M[j] = double.NegativeInfinity;
                        else if(itm.Contains("inf"))
                            M[j] = double.PositiveInfinity;
                        else
                            throw fe;
                    }
                }

                tempMatrix.Add(M);
            }

            // convert to multidimensional array
            // =================================

            MultidimensionalArray R = MultidimensionalArray.Create(tempMatrix.Count, NoOfCols);
            for(int i = 0; i < R.NoOfRows; i++) {
                for(int j  =0; j < R.NoOfCols; j++) {
                    R[i, j] = tempMatrix[i][j];
                }
            }
            return R;
        }


        /// <summary>
        /// Finds those row of <paramref name="mda"/>, where the L2-distance between the row and <paramref name="Row"/> is minimal. 
        /// </summary>
        /// <returns>
        /// the row index of the minimum-distance row.
        /// </returns>
        static public int MindistRow(this IMatrix mda, double[] Row) {
            double Dmax;
            int iDmax;
            MindistRow(mda, Row, out Dmax, out iDmax);
            return iDmax;
        }

        /// <summary>
        /// Finds those row <paramref name="iDmax"/> of <paramref name="mda"/>, where the L2-distance between the row and <paramref name="Row"/> is minimal. 
        /// </summary>
        /// <param name="mda">some matrix</param>
        /// <param name="Row">some row</param>
        /// <param name="Dmax">minimum L2 distance over all rows</param>
        /// <param name="iDmax">index of minimum L2-distance row</param>
        /// <returns></returns>
        static public void MindistRow(this IMatrix mda, double[] Row, out double Dmax, out int iDmax) {
            if(Row.Length != mda.NoOfCols)
                throw new ArgumentException();


            int N = mda.NoOfCols;
            int M = mda.NoOfRows;

            Dmax = double.MaxValue;
            iDmax = int.MinValue;

            for(int i = 0; i < M; i++) {
                double dist = 0.0;
                for(int j = 0; j < N; j++) {
                    double a = Row[j] - mda[i, j];
                    dist += a * a;
                }

                if(Dmax > dist) {
                    Dmax = dist;
                    iDmax = i;
                }
            }

            Dmax = Math.Sqrt(Dmax);
        }


        /// <summary>
        /// Finds those row of <paramref name="mda"/>, where the L2-distance between the row and <paramref name="Row"/> is 
        /// below <paramref name="tol"/>; returns a negative number otherwise.
        /// </summary>
        static public int FindRow(this IMatrix mda, double[] Row, double tol) {
            double Dmax;
            int iDmax;
            MindistRow(mda, Row, out Dmax, out iDmax);
            if(Dmax < tol)
                return iDmax;
            else
                return -1;
        }

        /// <summary>
        /// Computes the L2-distance between any two different rows of <paramref name="mda"/>, and returns the minimum.
        /// </summary>
        static public double MindistBetweenRows(this IMatrix mda) {

            int N = mda.NoOfCols;
            int M = mda.NoOfRows;

            double Dmax = double.MaxValue;

            for(int i = 0; i < M; i++) {

                for(int k = i + 1; k < M; k++) {

                    double dist = 0.0;
                    for(int j = 0; j < N; j++) {
                        double a = mda[i, j] - mda[k, j];
                        dist += a * a;
                    }

                    if(Dmax > dist) {
                        Dmax = dist;
                    }
                }
            }

            Dmax = Math.Sqrt(Dmax);
            return Dmax;
        }

        /// <summary>
        /// Computes the L2-distance between any two different rows of <paramref name="mda"/>, and returns the maximum.
        /// </summary>
        static public double MaxdistBetweenRows(this IMatrix mda) {

            int N = mda.NoOfCols;
            int M = mda.NoOfRows;

            double Dmin = 0;

            for (int i = 0; i < M; i++) {

                for (int k = i + 1; k < M; k++) {

                    double dist = 0.0;
                    for (int j = 0; j < N; j++) {
                        double a = mda[i, j] - mda[k, j];
                        dist += a * a;
                    }

                    if (Dmin < dist) {
                        Dmin = dist;
                    }
                }
            }

            Dmin = Math.Sqrt(Dmin);
            return Dmin;
        }

        static unsafe void CopyToUnsafeBuffer<T>(T M, double* buffer, bool BufferInFortranOrder) where T : IMatrix {
            int I = M.NoOfRows, J = M.NoOfCols;

            if (BufferInFortranOrder) {
                // buffer in FORTRAN order
                for (int i = 0; i < I; i++) {
                    for (int j = 0; j < J; j++) {
                        buffer[i + j * I] = M[i, j];
                    }
                }
            } else {
                // buffer in C-order
                for (int i = 0; i < I; i++) {
                    for (int j = 0; j < J; j++) {
                        buffer[i * J + j] = M[i, j];
                    }
                }
            }
        }

        static unsafe void CopyFromUnsafeBuffer<T>(T M, double* buffer, bool BufferInFortranOrder) where T : IMatrix {
            int I = M.NoOfRows, J = M.NoOfCols;

            if (BufferInFortranOrder) {
                // buffer in FORTRAN order
                for (int i = 0; i < I; i++) {
                    for (int j = 0; j < J; j++) {
                        M[i, j] = buffer[i + j * I];
                    }
                }
            } else {
                // buffer in C-order
                for (int i = 0; i < I; i++) {
                    for (int j = 0; j < J; j++) {
                        M[i, j] = buffer[i * J + j];
                    }
                }
            }
        }


        /// <summary>
        /// Copies values from a array, either in FORTRAN or C-order
        /// </summary>
        /// <param name="src">
        /// source to copy from
        /// </param>
        /// <param name="iOffset">
        /// offset into <paramref name="src"/>;
        /// </param>
        /// <param name="FortranOrder">
        /// if true, values are taken in FORTRAN order, otherwise they are taken in c-order
        /// </param>
        /// <param name="b">
        /// some other matrix
        /// </param>
        static public void SetFromBuffer(this IMatrix b, double[] src, int iOffset, bool FortranOrder) {
            int cnt = iOffset;

            int m_NoOfCols = b.NoOfCols;
            int m_NoOfRows = b.NoOfRows;

            for (int j = 0; j < m_NoOfCols; j++) {
                for (int i = 0; i < m_NoOfRows; i++) {
                    double val = src[cnt];
                    cnt++;

                    if (FortranOrder)
                        b[i, j] = val;
                    else
                        b[j, i] = val;
                }
            }
        }

        /// <summary>
        /// The 1-Norm (maximum absolute column sum norm) of this matrix;
        /// </summary>
        /// <returns></returns>
        /// <param name="b">the matrix</param>
        static public double OneNorm(this IMatrix b) {
            double norm = 0;

            int m_NoOfCols = b.NoOfCols;
            int m_NoOfRows = b.NoOfRows;

            for (int j = 0; j < m_NoOfCols; j++) {
                double colnrm = 0;

                for (int i = 0; i < m_NoOfRows; i++) {
                    colnrm += Math.Abs(b[i, j]);
                }

                if (colnrm > norm)
                    norm = colnrm;
            }

            return norm;
        }

        /// <summary>
        /// The Infinity-Norm (maximum absolute row sum norm) of a matrix;
        /// </summary>
        static public double InfNorm(this IMatrix b) {
            double norm = 0;

            int m_NoOfCols = b.NoOfCols;
            int m_NoOfRows = b.NoOfRows;

            for (int i = 0; i < m_NoOfRows; i++) {
                double rownrm = 0;

                for (int j = 0; j < m_NoOfCols; j++) {
                    rownrm += Math.Abs(b[i, j]);
                }

                if (rownrm > norm)
                    norm = rownrm;
            }

            return norm;
        }

        /// <summary>
        /// detects NAN's and INF's 
        /// </summary>
        /// <param name="TestNaN">turn NAN detection on/off</param>
        /// <param name="TestInf">turn INF detection on/off</param>
        /// <returns>true, if a suspicious value was found</returns>
        /// <param name="b">the matrix</param>
        static public bool ContainsNanOrInf(this IMatrix b, bool TestNaN = true, bool TestInf = true) {
            int m_NoOfCols = b.NoOfCols;
            int m_NoOfRows = b.NoOfRows;

            for (int i = 0; i < m_NoOfRows; i++) {

                for (int j = 0; j < m_NoOfCols; j++) {
                    if (TestNaN && double.IsNaN(b[i, j])) {
                        return true;
                    }

                    if (TestInf && double.IsInfinity(b[i, j])) {
                        return true;
                    }
                }
            }

            return false;
        }

        /// <summary>
        /// General matrix vector product
        /// <paramref name="y"/> = <paramref name="y"/>*<paramref name="yScaling"/>
        /// + <paramref name="M"/>*<paramref name="x"/>*<paramref name="xScaling"/>
        /// </summary>
        /// <param name="xScaling">scaling of <paramref name="x"/></param>
        /// <param name="x">the vector to multiply this matrix with</param>
        /// <param name="yScaling">
        /// prescaling (before accumulation) of <paramref name="y"/>
        /// </param>
        /// <param name="y">
        /// vector to accumulate the result of the operation
        /// </param>
        /// <param name="transpose">true for transpose multiplication</param>
        /// <param name="M">the matrix</param>
        static public void gemv<MatrixType, VectorType1, VectorType2>(this MatrixType M, double xScaling, VectorType1 x, double yScaling, VectorType2 y, bool transpose = false)
            where MatrixType : IMatrix
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {
            //

            // Make sure $x is "cloned" if it is also the target of the
            // operation
            IList<double> xCopy = x;
            if (ReferenceEquals(x, y)) {
                xCopy = x.ToList();
            }

            int m_NoOfCols = M.NoOfCols, m_NoOfRows = M.NoOfRows;

            if (!transpose) {
                if (m_NoOfCols != x.Count)
                    throw new ArgumentException("Length of vector x must be equal to number of columns.");
                if (m_NoOfRows != y.Count)
                    throw new ArgumentException("Length of vector y must be equal to number of rows.");


                for (int i = 0; i < m_NoOfRows; i++) {
                    double yi = 0;

                    for (int j = 0; j < m_NoOfCols; j++)
                        yi += M[i, j] * xCopy[j];

                    y[i] = y[i] * yScaling + yi * xScaling;
                }

            } else {
                if (m_NoOfRows != x.Count)
                    throw new ArgumentException("Length of vector x must be equal to number of rows (transpose matrix multiply was selected!).");
                if (m_NoOfCols != y.Count)
                    throw new ArgumentException("Length of vector y must be equal to number of columns (transpose matrix multiply was selected!).");


                for (int i = 0; i < m_NoOfCols; i++) {
                    double yi = 0;

                    for (int j = 0; j < m_NoOfRows; j++)
                        yi += M[j, i] * xCopy[j];

                    y[i] = y[i] * yScaling + yi * xScaling;
                }
            }

        }

        /// <summary>
        /// General matrix/matrix multiplication:
        /// </summary>
        /// <returns>
        /// <paramref name="A"/>*<paramref name="B"/>
        /// </returns>
        static public MultidimensionalArray GEMM<Matrix2, Matrix3>(this Matrix2 A, Matrix3 B)
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix //
        {
            MultidimensionalArray R = MultidimensionalArray.Create(A.NoOfRows, B.NoOfCols);
            R.GEMM(1.0, A, B, 0.0);
            return R;
        }
        
        /// <summary>
        /// General matrix/matrix/matrix multiplication:
        /// </summary>
        /// <returns>
        /// <paramref name="A"/>*<paramref name="B"/>*<paramref name="C"/>
        /// </returns>
        static public MultidimensionalArray GEMM<Matrix2, Matrix3, Matrix4>(this Matrix2 A, Matrix3 B, Matrix4 C)
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix 
            where Matrix4 : IMatrix //
        {
            return GEMM(GEMM(A, B), C);
        }

        static MultidimensionalArray.MultiplyProgram GEMM_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ik", "kj");

        /// <summary>
        /// General matrix/matrix multiplication:
        /// <paramref name="M"/> = <paramref name="alpha"/>*<paramref name="A"/>*<paramref name="B"/> + <paramref name="beta"/>*<paramref name="M"/>;
        /// </summary>
        static public void GEMM<Matrix1, Matrix2, Matrix3>(this Matrix1 M, double alpha, Matrix2 A, Matrix3 B, double beta)
            where Matrix1 : IMatrix
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix //
        {
            if (A.NoOfCols != B.NoOfRows)
                throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
            if (A.NoOfRows != M.NoOfRows)
                throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
            if (B.NoOfCols != M.NoOfCols)
                throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");

            if (A is MultidimensionalArray && B is MultidimensionalArray && M is MultidimensionalArray) {
                MultidimensionalArray _A = A as MultidimensionalArray;
                MultidimensionalArray _B = B as MultidimensionalArray;
                MultidimensionalArray _M = M as MultidimensionalArray;

                _M.Multiply(alpha, _A, _B, beta, ref GEMM_Prog);
            } else {

                int K = A.NoOfCols;

                for (int i = M.NoOfRows - 1; i >= 0; i--) {
                    for (int j = M.NoOfCols - 1; j >= 0; j--) {
                        double r = 0;


                        for (int k = K - 1; k >= 0; k--)
                            r += A[i, k] * B[k, j];

                        r *= alpha;
                        M[i, j] = M[i, j] * beta + r;
                    }

                }
            }
        }

        ///// <summary>
        ///// General matrix/matrix multiplication:
        ///// <paramref name="C"/> = <paramref name="alpha"/>*<paramref name="A"/>*<paramref name="B"/> + <paramref name="beta"/>*<paramref name="C"/>;
        ///// </summary>
        //public static void gemm(double alpha, FullMatrix A, FullMatrix B, double beta, FullMatrix C) {
        //    if (A.NoOfCols != B.NoOfRows)
        //        throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
        //    if (A.NoOfRows != C.NoOfRows)
        //        throw new ArgumentException("A.NoOfRows != C.NoOfRows", "A,C");
        //    if (B.NoOfCols != C.NoOfCols)
        //        throw new ArgumentException("B.NoOfCols != C.NoOfCols", "B,C");

        //    int K = A.NoOfCols;

        //    for (int i = C.NoOfRows - 1; i >= 0; i--) {
        //        for (int j = C.NoOfCols - 1; j >= 0; j--) {
        //            double r = 0;


        //            for (int k = K - 1; k >= 0; k--)
        //                r += A[i, k] * B[k, j];

        //            r *= alpha;
        //            C[i, j] = C[i, j] * beta + r;
        //        }

        //    }
        //}

        /// <summary>
        /// total absolute sum of all entries
        /// </summary>
        /// <returns></returns>
        public static double TotalAbsoluteSum(this IMatrix M) {
            double r = 0.0;
            for (int i = 0; i < M.NoOfRows; i++) {
                for (int j = 0; j < M.NoOfCols; j++) {
                    r += Math.Abs(M[i, j]);
                }
            }
            return r;
        }

        /// <summary>
        /// Transposes this matrix without allocating additional memory. Thus,
        /// this matrix will be changed by the operation!
        /// </summary>
        /// <param name="M">the matrix</param>
        /// <returns>This object after transposition/reference-equal to <paramref name="M"/></returns>
        static public TM TransposeInPlace<TM>(this TM M) where TM : IMatrix {
            if (M.NoOfRows != M.NoOfCols)
                throw new NotSupportedException("in-place transposition is only supported for quadratic matrices.");


            M.TransposeTo(M);

            return M;
        }

        /// <summary>
        /// Returns a new matrix containing the transpose of this matrix.
        /// </summary>
        /// <param name="M">the matrix</param>
        /// <returns>The transpose of this matrix</returns>
        static public MultidimensionalArray Transpose<TM>(this TM M) where TM : IMatrix {
            MultidimensionalArray result = MultidimensionalArray.Create(M.NoOfCols, M.NoOfRows);

            M.TransposeTo(result);

            return result;
        }

        /// <summary>
        /// writes the transpose of matrix <paramref name="source"/> to <paramref name="target"/>;
        /// matrix <paramref name="source"/> is unchanged, if not reference-equal to <paramref name="target"/>
        /// </summary>
        /// <param name="source">original matrix</param>
        /// <param name="target">result of the operation</param>
        static public void TransposeTo<M1, M2>(this M1 source, M2 target)
            where M1 : IMatrix
            where M2 : IMatrix {


            if (target.NoOfCols != source.NoOfRows || target.NoOfRows != source.NoOfCols) {
                throw new ArgumentException("wrong size of target.");
            }

            int NoOfRows = source.NoOfRows, NoOfCols = source.NoOfCols;

            if (object.ReferenceEquals(source, target)) {
                if (NoOfRows != NoOfCols)
                    throw new NotSupportedException("in-place transposition is only supported for quadratic matrices.");


                for (int i = 0; i < NoOfRows; i++) {
                    for (int j = 0; j < i; j++) {
                        double temp = source[i, j];
                        target[i, j] = source[j, i];
                        target[j, i] = temp;
                    }
                }
            } else {
                for (int i = 0; i < source.NoOfRows; i++) {
                    for (int j = 0; j < source.NoOfCols; j++) {
                        target[j, i] = source[i, j];
                    }
                }
            }
        }

        /*
        /// <summary>
        /// standard gemm: 
        /// <paramref name="C"/> = <paramref name="beta"/>*<paramref name="C"/>
        /// + <paramref name="alpha"/>*<paramref name="A"/>*<paramref name="B"/>;
        /// </summary>
        /// <param name="A"></param>
        /// <param name="alpha"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        /// <param name="beta"></param>
        static public void Gemm(IMatrix A, double alpha, IMatrix B, IMatrix C, double beta) {
            if (A.NoOfCols != B.NoOfRows)
                throw new ArgumentException("Matrix size mismatch.");
            if (A.NoOfRows != C.NoOfRows)
                throw new ArgumentException("Matrix size mismatch.");
            if (B.NoOfCols != C.NoOfCols)
                throw new ArgumentException("Matrix size mismatch.");

            int M = A.NoOfRows, N = A.NoOfCols, K = B.NoOfCols;

            for (int i = 0; i < M; i++) {
                for (int j = 0; j < K; j++) {
                    double r = 0;

                    for (int k = 0; k < N; k++)
                        r += A[i, k] * B[k, j];

                    C[i, j] = C[i, j] * beta + r * alpha;
                }
            }
        }*/

        /// <summary>
        /// only supported for 1x1, 2x2, 3x3 and 4x4 - matrices;
        /// </summary>
        /// <returns></returns>
        static public double Determinant<T>(this T M) where T : IMatrix {
            int m_NoOfCols = M.NoOfCols, m_NoOfRows = M.NoOfRows;
            if (m_NoOfCols != m_NoOfRows)
                throw new NotSupportedException("Determinate only defined for quadratic matrices.");
            if (m_NoOfCols > 4) {
                throw new NotImplementedException("Only implemented for 1x1, 2x2, 3x3 and 4x4 matrices.");
            }

            switch (m_NoOfCols) {
                case 1:
                    return M[0, 0];

                case 2:
                    return (M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]);

                case 3:
                    return (M[0, 0] * M[1, 1] * M[2, 2] - M[0, 0] * M[1, 2] * M[2, 1]
                             - M[1, 0] * M[0, 1] * M[2, 2] + M[1, 0] * M[0, 2] * M[2, 1]
                             + M[2, 0] * M[0, 1] * M[1, 2] - M[2, 0] * M[0, 2] * M[1, 1]);

                case 4:
                    return (M[0, 0] * M[1, 1] * M[2, 2] * M[3, 3]
                            - M[0, 0] * M[1, 1] * M[2, 3] * M[3, 2]
                            - M[0, 0] * M[2, 1] * M[1, 2] * M[3, 3]
                            + M[0, 0] * M[2, 1] * M[1, 3] * M[3, 2]
                            + M[0, 0] * M[3, 1] * M[1, 2] * M[2, 3]
                            - M[0, 0] * M[3, 1] * M[1, 3] * M[2, 2]
                            - M[1, 0] * M[0, 1] * M[2, 2] * M[3, 3]
                            + M[1, 0] * M[0, 1] * M[2, 3] * M[3, 2]
                            + M[1, 0] * M[2, 1] * M[0, 2] * M[3, 3]
                            - M[1, 0] * M[2, 1] * M[0, 3] * M[3, 2]
                            - M[1, 0] * M[3, 1] * M[0, 2] * M[2, 3]
                            + M[1, 0] * M[3, 1] * M[0, 3] * M[2, 2]
                            + M[2, 0] * M[0, 1] * M[1, 2] * M[3, 3]
                            - M[2, 0] * M[0, 1] * M[1, 3] * M[3, 2]
                            - M[2, 0] * M[1, 1] * M[0, 2] * M[3, 3]
                            + M[2, 0] * M[1, 1] * M[0, 3] * M[3, 2]
                            + M[2, 0] * M[3, 1] * M[0, 2] * M[1, 3]
                            - M[2, 0] * M[3, 1] * M[0, 3] * M[1, 2]
                            - M[3, 0] * M[0, 1] * M[1, 2] * M[2, 3]
                            + M[3, 0] * M[0, 1] * M[1, 3] * M[2, 2]
                            + M[3, 0] * M[1, 1] * M[0, 2] * M[2, 3]
                            - M[3, 0] * M[1, 1] * M[0, 3] * M[2, 2]
                            - M[3, 0] * M[2, 1] * M[0, 2] * M[1, 3]
                            + M[3, 0] * M[2, 1] * M[0, 3] * M[1, 2]);

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// multiplies matrix <paramref name="M"/> with <paramref name="alpha"/>; entries of this
        /// matrix will be overwritten;
        /// </summary>
        static public void Scale<T>(this T M, double alpha) where T : IMatrix {
            if (M is MultidimensionalArray) {
                MultidimensionalArray _M = M as MultidimensionalArray;
                _M.Scale(alpha);
            } else {
                int NoRows = M.NoOfRows, Nocols = M.NoOfCols;
                for (int i = 0; i < NoRows; i++)
                    for (int j = 0; j < Nocols; j++)
                        M[i, j] *= alpha;
            }
        }

        /// <summary>
        /// Calculates the inverse of this matrix and stores the result in a
        /// newly allocated matrix. Use <see cref="Invert(FullMatrix)"/> to
        /// avoid memory allocation.
        /// </summary>
        /// <returns>
        /// A matrix containing the inverse of this matrix.
        /// </returns>
        static public MultidimensionalArray GetInverse<T>(this T Mtx) where T : IMatrix {
            var target = MultidimensionalArray.Create(Mtx.NoOfRows, Mtx.NoOfCols);
            InvertTo(Mtx, target);
            return target;
        }

        /// <summary>
        /// Cholesky factorization with LAPACK procedure DPOTRF.
        /// </summary>
        static public void Cholesky<T>(this T M) where T : IMatrix {
            if(M.NoOfCols != M.NoOfRows)
                throw new NotSupportedException("Expecting a square matrix.");
            int N = M.NoOfCols;

            // clear upper triangular part (not really necessary, but anyway)...
            for(int i = 0; i < N; i++) // loop over rows
                for(int j = i + 1; j < N; j++) // loop over upper-triangular columns
                    M[i, j] = 0.0;

            unsafe {
                int i0;
                double[] _this_entries =  TempBuffer.GetTempBuffer(out i0, M.NoOfCols * M.NoOfRows);
                fixed(double* this_entries = _this_entries) {
                    CopyToUnsafeBuffer(M, this_entries, true);

                    int info, UPLO = 'L';
                    LAPACK.F77_LAPACK.DPOTRF_(ref UPLO, ref N, this_entries, ref N, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK DPOTRF (Cholesky factorization) returned info " + info);
                    }
                    
                    CopyFromUnsafeBuffer(M, this_entries, true);
                }
                TempBuffer.FreeTempBuffer(i0);
            }

            

        }

        /// <summary>
        /// Inversion of matrix <paramref name="M"/> (in-place) but only if
        /// <paramref name="M"/> is symmetrical and positive definite, by using
        /// a Cholesky factorization.
        /// </summary>
        static public void InvertSymmetrical<T>(this T M) where T : IMatrix {
            if (M.NoOfCols != M.NoOfRows)
                throw new NotSupportedException("can't invert non-square matrix.");
            int N = M.NoOfCols;

            // clear upper triangular part (not really necessary, but anyway)...
            for (int i = 0; i < N; i++) // loop over rows
                for (int j = i + 1; j < N; j++) // loop over upper-triangular columns
                    M[i, j] = 0.0;

            unsafe {
                int i0;
                double[] _this_entries =  TempBuffer.GetTempBuffer(out i0, M.NoOfCols * M.NoOfRows);
                fixed(double* this_entries = _this_entries) {
                //fixed (double* this_entries = TempBuffer.GetTempBuffer(out i0, M.NoOfCols * M.NoOfRows)) {
                    CopyToUnsafeBuffer(M, this_entries, true);

                    int info, UPLO = 'L';
                    LAPACK.F77_LAPACK.DPOTRF_(ref UPLO, ref N, this_entries, ref N, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK DPOTRF (Cholesky factorization) returned info " + info);
                    }
                    LAPACK.F77_LAPACK.DPOTRI_(ref UPLO, ref N, this_entries, ref N, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK DPOTRI (Cholesky inversion) returned info " + info);
                    }

                    CopyFromUnsafeBuffer(M, this_entries, true);
                }
                TempBuffer.FreeTempBuffer(i0);
            }

            // set upper triangular part 
            for (int i = 0; i < N; i++) // loop over rows
                for (int j = i + 1; j < N; j++) // loop over upper-triangular columns
                    M[i, j] = M[j, i];

        }

        /// <summary>
        /// on exit, <paramref name="Left"/>*<paramref name="Q"/>*<paramref name="Right"/> = EYE.
        /// </summary>
        static public void SymmetricLDLInversion<M1, M2, M3>(this M1 Q, M2 Left, M3 Right)
            where M1 : IMatrix
            where M2 : IMatrix
            where M3 : IMatrix {
            double[] diag = new double[Q.NoOfRows];
            Q.SymmetricLDLInversion(Right, diag);

            Right.TransposeTo(Left);
            for (int i = 0; i < diag.Length; i++) {
                Left.RowScale(i, diag[i]);
            }
        }

        /// <summary>
        /// L1 distance between the upper and the lower diagonal entries
        /// </summary>
        static public double SymmetryError<T>(this T M) where T : IMatrix {
            if (M.NoOfCols != M.NoOfRows)
                throw new NotSupportedException();
            double R = 0;
            int N = M.NoOfRows;

            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    R += Math.Abs(M[i, j] - M[j, i]);
                }
            }

            return R;
        }


        /// <summary>
        /// Inversion of an upper or lower triangular matrix.
        /// </summary>
        /// <param name="Mtx">input</param>
        /// <param name="B">output; will be overwritten</param>
        /// <param name="Upper">If true, upper triangular invert; otherwise, lower.</param>
        static public void TriangularInvert<T1, T2>(this T1 Mtx, T2 B, bool Upper = true)
            where T1 : IMatrix
            where T2 : IMatrix //
        {

            if (Mtx.NoOfCols != Mtx.NoOfRows)
                throw new NotSupportedException("must be symmetrical");
            if (B.NoOfRows != Mtx.NoOfRows || B.NoOfCols != Mtx.NoOfCols)
                throw new ArgumentOutOfRangeException("output matrix must have the same size as this matrix.");
            int N = Mtx.NoOfCols;

#if DEBUG
            MultidimensionalArray MtxClone = MultidimensionalArray.Create(Mtx.NoOfRows, Mtx.NoOfCols);
            MtxClone.SetMatrix(Mtx);
#endif

            unsafe
            {
                int i0;
                double[] _B_entries = TempBuffer.GetTempBuffer(out i0, B.NoOfCols * B.NoOfRows);
                fixed (double* B_entries = _B_entries) {
                    CopyToUnsafeBuffer(Mtx, B_entries, true);


                    // clear lower triangular part 
                    for (int i = 0; i < N; i++) // loop over rows
                        for (int j = 0; j < i; j++) // loop over lower-triangular columns
                            B_entries[i + j * N] = 0.0;

                    int UPLO = Upper ? 'U' : 'L';
                    int DIAG = 'N'; // not unit-triangular
                    int info;
                    LAPACK.F77_LAPACK.DTRTRI_(ref UPLO, ref DIAG, ref N, B_entries, ref N, out info);

                    CopyFromUnsafeBuffer(B, B_entries, true);

                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK DTRTRI (triangular inversion) returned info " + info);
                    }


                }
                TempBuffer.FreeTempBuffer(i0);
            }
        }



        /// <summary>
        /// experimental LDL inversion/factorisation
        /// </summary>
        /// <remarks>
        /// on exit, Diag(<paramref name="Diag"/>)*<paramref name="B"/>^T*<paramref name="Mtx"/>*<paramref name="B"/> = EYE;
        /// </remarks>
        /// <param name="Mtx">input; remains unchanged</param>
        /// <param name="B">output; will be overwritten</param>
        /// <param name="Diag">
        /// output; can be null, if just a Choleski factorization is required -- in that case,
        /// an exception will occur if this matrix is not positive definite.
        /// </param>
        static public void SymmetricLDLInversion<T1, T2>(this T1 Mtx, T2 B, double[] Diag)
            where T1 : IMatrix
            where T2 : IMatrix //
        {

            if (Mtx.NoOfCols != Mtx.NoOfRows)
                throw new NotSupportedException("must be symmetrical");
            if (B.NoOfRows != Mtx.NoOfRows || B.NoOfCols != Mtx.NoOfCols)
                throw new ArgumentOutOfRangeException("output matrix must have the same size as this matrix.");
            if (Diag != null && Diag.Length != Mtx.NoOfRows)
                throw new ArgumentOutOfRangeException();

            int N = Mtx.NoOfRows;

            if (Diag != null)
                Diag.SetAll(1.0);

#if DEBUG
            MultidimensionalArray MtxClone = MultidimensionalArray.Create(Mtx.NoOfRows, Mtx.NoOfCols);
            MtxClone.SetMatrix(Mtx);
#endif

            unsafe {
                int i0;
                double[] _B_entries =  TempBuffer.GetTempBuffer(out i0, B.NoOfCols * B.NoOfRows);
                fixed(double* B_entries = _B_entries) {
                    CopyToUnsafeBuffer(Mtx, B_entries, true);

                    // at first, try with LAPACK Cholesky
                    int info;
                    {

                        int UPLO = 'U';
                        LAPACK.F77_LAPACK.DPOTRF_(ref UPLO, ref N, B_entries, ref N, out info);
                        if(info < 0) {
                            TempBuffer.FreeTempBuffer(i0);
                            throw new ArithmeticException("LAPACK DPOTRF (Cholesky factorisation) returned info " + info);
                        }
                    }
                    //info = 1;

                    if (info == 0) {
                        // cholesky successful: go on with it

                        // clear lower triangular part 
                        for (int i = 0; i < N; i++) // loop over rows
                            for (int j = 0; j < i; j++) // loop over lower-triangular columns
                                B_entries[i + j * N] = 0.0;

                        int UPLO = 'U', DIAG = 'N';
                        LAPACK.F77_LAPACK.DTRTRI_(ref UPLO, ref DIAG, ref N, B_entries, ref N, out info);

                        CopyFromUnsafeBuffer(B, B_entries, true);

                        if(info != 0) {
                            TempBuffer.FreeTempBuffer(i0);
                            throw new ArithmeticException("LAPACK DTRTRI (triangular inversion) returned info " + info);
                        }
                    } else {

                        if (Diag != null) {
                            // fall back to our homebrew LDL factorization (if desired by user)
                            B.Clear();
                            B.Acc(1.0, Mtx);
                            Ha(Mtx, B, Diag);
                            //GramSchmidt(this, B, Diag, N);
                        } else {
                            //Mtx.SaveToTextFile("C:\\temp\\Mtx.txt");
                            TempBuffer.FreeTempBuffer(i0);
                            throw new ArithmeticException("Not positive definite: LAPACK DPOTRF (Cholesky factorisation) returned info " + info);
                        }
                    }

                    // finally, do a test 
                    {
#if DEBUG
                        double Nrm = Math.Max(MtxClone.InfNorm(), B.InfNorm());
                        var Test = GEMM(GEMM(B.Transpose(), MtxClone), B);
                        if(Diag != null) {
                            for(int n = 0; n < N; n++) {
                                Test.RowScale(n, Diag[n]);
                            }
                        }

                        Test.AccEye(-1.0);

                        double err = Test.InfNorm();
                        if(err / Nrm > 1.0e-5) {
                            TempBuffer.FreeTempBuffer(i0);
                            throw new ArithmeticException();
                        }
#endif
                    }
                }
                TempBuffer.FreeTempBuffer(i0);
            }
        }

        static private void GramSchmidt<FullMatrix1, FullMatrix2>(FullMatrix1 Mtx, FullMatrix2 B, double[] Diag)
            where FullMatrix1 : IMatrix
            where FullMatrix2 : IMatrix {
            B.Clear();
            int N = B.NoOfRows;

            for (int m = 0; m < N; m++) {
                B[m, m] = 1.0;

                for (int n = 0; n <= (m - 1); n++) {
                    double zahler = 0;
                    double nenner = 0;
                    for (int o = 0; o <= n; o++) {
                        zahler = zahler + B[o, n] * Mtx[o, m];
                    }
                    for (int o = 0; o <= n; o++) {
                        for (int p = 0; p <= n; p++) {
                            nenner = nenner + B[o, n] * B[p, n] * Mtx[o, p];
                        }
                    }

                    double alpha_mn = zahler / nenner;

                    for (int o = 0; o <= n; o++) {
                        B[o, m] = B[o, m] - B[o, n] * alpha_mn;
                    }
                }

                // normalisation
                double _nenner = 0;
                for (int o = 0; o <= m; o++) {
                    for (int p = 0; p <= m; p++) {
                        _nenner = _nenner + B[o, m] * B[p, m] * Mtx[o, p];
                    }
                }
                double oneOverNorm = 1.0;
                if (_nenner > 0.0) {
                    oneOverNorm = Math.Sqrt(1.0 / _nenner);
                } else {
                    oneOverNorm = Math.Sqrt(1.0 / Math.Abs(_nenner));
                    Diag[m] = -1.0;
                    //fprintf(1,'negative norm for vector %d \n',m);
                }
                //oneOverNorm = 1.0/sqrt(nenner);
                for (int n = 0; n <= m; n++) {
                    B[n, m] = B[n, m] * oneOverNorm;
                }

            }
        }

        static private void Ha<FullMatrix1, FullMatrix2>(FullMatrix1 _Mtx, FullMatrix2 B, double[] Diag)
            where FullMatrix1 : IMatrix
            where FullMatrix2 : IMatrix {
            var Mtx = MultidimensionalArray.Create(_Mtx.NoOfRows, _Mtx.NoOfCols);
            Mtx.SetMatrix(_Mtx);
            FullMatrix2 LOW = B;
            LOW.Clear();
            LOW.AccEye(1.0);
            int N = B.NoOfRows;

            for (int n = 0; n < N; n++) { // loop over columns

                // pivotisierung:
                double m_nn = Mtx[n, n];
                int imax = n;
                for (int i = n + 1; i < N; i++) {
                    double m_ii = Mtx[i, i];
                    if (Math.Abs(m_ii) > Math.Abs(m_nn)) {
                        m_nn = m_ii;
                        imax = i;
                    }
                }
                if (imax != n) {
                    Mtx.SwapRows(imax, n);
                    Mtx.SwapColumns(imax, n);
                    Debug.Assert(Mtx[n, n] == m_nn);
                    LOW.SwapRows(imax, n);
                }


                m_nn = Mtx[n, n];
                if (m_nn < 0.0) {
                    // indefinite entry:
                    //Mtx.RowMul(n, -1.0);
                    Diag[n] *= -1.0;
                    m_nn = Mtx[n, n];
                }
                double sign_mnn = Math.Sign(m_nn);
                if (sign_mnn == 0.0)
                    throw new ArithmeticException();
                double oo_mnn = Math.Sqrt(Math.Abs(1.0 / m_nn));


                Mtx.RowScale(n, oo_mnn);
                LOW.RowScale(n, oo_mnn);

                for (int i = n + 1; i < N; i++) { // elimination loop: we eliminate row n+1 to N in column n
                    double oldVal = Mtx[i, n];
                    double fac = -oo_mnn * oldVal * sign_mnn;

                    Mtx.RowAdd(n, i, fac);
                    LOW.RowAdd(n, i, fac);

                    Debug.Assert(Math.Abs(Mtx[i, n]) <= Math.Max(Math.Abs(oldVal) * 1.0e-14, 1.0e-10));
                    Mtx[i, n] = 0.0;
                }

                Mtx[n, n] = 1.0;
                for (int j = n + 1; j < N; j++) { // elimination loop: we eliminate col n+1 to N in row n 
                    Mtx[n, j] = 0.0;
                }

            }

            // we produced a lower triangular matrix, but we want upper triangular
            LOW.TransposeInPlace();
        }

        /// <summary>
        /// Calculates the inverse of this matrix and stores the result in <paramref name="target"/>.
        /// </summary>
        /// <param name="source">the original matrix</param>
        /// <param name="target">the return value</param>
        static public void InvertTo<M1, M2>(this M1 source, M2 target)
            where M1 : IMatrix
            where M2 : IMatrix {
            int m_NoOfCols = source.NoOfCols, m_NoOfRows = source.NoOfRows;

            if (m_NoOfCols != m_NoOfRows)
                throw new NotSupportedException("can't invert non-square matrix.");
            if (target.NoOfRows != source.NoOfRows || target.NoOfCols != source.NoOfCols)
                throw new ArgumentException("target matrix must have the same size as source matrix.");

            switch (m_NoOfCols) {
                case 0: {
                    return;
                }

                case 1: {
                        double det = source.Determinant();
                        target[0, 0] = 1.0 / det;
                        return;
                    }

                case 2: {
                        double det = source.Determinant();
                        target[0, 0] = source[1, 1];
                        target[0, 1] = -source[0, 1];
                        target[1, 0] = -source[1, 0];
                        target[1, 1] = source[0, 0];
                        target.Scale(1.0 / det);
                        return;
                    }

                case 3: {
                        double det = source.Determinant();
                        target[0, 0] = source[1, 1] * source[2, 2] - source[1, 2] * source[2, 1];
                        target[0, 1] = -source[0, 1] * source[2, 2] + source[0, 2] * source[2, 1];
                        target[0, 2] = source[0, 1] * source[1, 2] - source[0, 2] * source[1, 1];

                        target[1, 0] = -source[1, 0] * source[2, 2] + source[1, 2] * source[2, 0];
                        target[1, 1] = source[0, 0] * source[2, 2] - source[0, 2] * source[2, 0];
                        target[1, 2] = -source[0, 0] * source[1, 2] + source[0, 2] * source[1, 0];

                        target[2, 0] = source[1, 0] * source[2, 1] - source[1, 1] * source[2, 0];
                        target[2, 1] = -source[0, 0] * source[2, 1] + source[0, 1] * source[2, 0];
                        target[2, 2] = source[0, 0] * source[1, 1] - source[0, 1] * source[1, 0];

                        target.Scale(1.0 / det);
                        return;
                    }

                default:
                    unsafe {

                        // use LAPACK
                        //
                        int N = source.NoOfCols;
                        int i0;
                        double[] _target_entries = TempBuffer.GetTempBuffer(out i0, N * N);
                        fixed (double* target_entries = _target_entries) {
                            int* IPIV = stackalloc int[N];
                            CopyToUnsafeBuffer(source, target_entries, true);

                            int info;
                            int lwork = -1;
                            if (N <= 64)
                                lwork = N * N;
                            else {
                                double twork = 0;
                                LAPACK.F77_LAPACK.DGETRI_(ref N, target_entries, ref N, IPIV, &twork, ref lwork, out info);
                                lwork = (int)twork;
                                if(info != 0) {
                                    TempBuffer.FreeTempBuffer(i0);
                                    throw new ArithmeticException("error in LAPACK routing DGETRI: info = " + info);
                                }
                            }

                            int i1;
                            double[] _work = TempBuffer.GetTempBuffer(out i1, lwork);
                            fixed(double* work = _work) {

                                int LDA = N;
                                LAPACK.F77_LAPACK.DGETRF(ref N, ref N, target_entries, ref LDA, IPIV, out info);
                                if(info != 0) {
                                    TempBuffer.FreeTempBuffer(i0);
                                    TempBuffer.FreeTempBuffer(i1);
                                    throw new ArithmeticException("error in LAPACK routing DGETRF: info = " + info);
                                }
                                LAPACK.F77_LAPACK.DGETRI_(ref N, target_entries, ref LDA, IPIV, work, ref lwork, out info);
                                if(info != 0) {
                                    TempBuffer.FreeTempBuffer(i0);
                                    TempBuffer.FreeTempBuffer(i1);
                                    throw new ArithmeticException("error in LAPACK routing DGETRI: info = " + info);
                                }
                                //IntPtr mem = Marshal.AllocHGlobal(N * N * 8 * 2);
                                //Marshal.Copy(target.Entries, 0, mem, target.entries.Length);
                                //unsafe {
                                //    LAPACK.DGETRI(ref N, (double*) mem, ref LDA, IPIV, work, ref lwork, out info);
                                //}
                                //Marshal.Copy(mem, target.entries, 0, target.entries.Length);

                                CopyFromUnsafeBuffer(target, target_entries, true);

                                
                            }
                            TempBuffer.FreeTempBuffer(i1);
                        }
                        TempBuffer.FreeTempBuffer(i0);

                        return;
                    }
            }
        }

        /// <summary>
        /// Calculates the inverse of this matrix and stores the result in <paramref name="target"/>.
        /// </summary>
        /// <param name="source">the original matrix</param>
        static public MultidimensionalArray InvertTo<M1>(this M1 source)
            where M1 : IMatrix // 
        {
            if (source.NoOfCols != source.NoOfRows)
                throw new NotSupportedException("can't invert non-square matrix.");
            MultidimensionalArray R = MultidimensionalArray.Create(source.NoOfRows, source.NoOfCols);

            InvertTo(source, R);

            return R;
        }
        


        /// <summary>
        /// Calculates the inverse of this matrix and stores it in-place.
        /// </summary>
        /// <param name="Mtx"></param>
        /// <param name="target">the return value</param>
        static public void Invert<M1>(this M1 Mtx)
            where M1 : IMatrix {
            int m_NoOfCols = Mtx.NoOfCols, m_NoOfRows = Mtx.NoOfRows;

            if (m_NoOfCols != m_NoOfRows)
                throw new NotSupportedException("can't invert non-square matrix.");

            switch (m_NoOfCols) {
                case 0:
                {
                    return;
                }

                case 1:
                {
                    Mtx[0, 0] = 1.0 / Mtx[0, 0];
                    return;
                }

                case 2:
                {
                    double Mtx00 = Mtx[0, 0];
                    double Mtx01 = Mtx[0, 1];
                    double Mtx10 = Mtx[1, 0];
                    double Mtx11 = Mtx[1, 1];

                    double det = Mtx.Determinant();
                    Mtx[0, 0] = Mtx11;
                    Mtx[0, 1] = -Mtx01;
                    Mtx[1, 0] = -Mtx10;
                    Mtx[1, 1] = Mtx00;
                    Mtx.Scale(1.0 / det);
                    return;
                }

                case 3:
                {
                    double Mtx00 = Mtx[0, 0];
                    double Mtx01 = Mtx[0, 1];
                    double Mtx02 = Mtx[0, 2];
                    double Mtx10 = Mtx[1, 0];
                    double Mtx11 = Mtx[1, 1];
                    double Mtx12 = Mtx[1, 2];
                    double Mtx20 = Mtx[2, 0];
                    double Mtx21 = Mtx[2, 1];
                    double Mtx22 = Mtx[2, 2];

                    double det = Mtx.Determinant();
                    Mtx[0, 0] = Mtx11 * Mtx22 - Mtx12 * Mtx21;
                    Mtx[0, 1] = -Mtx01 * Mtx22 + Mtx02 * Mtx21;
                    Mtx[0, 2] = Mtx01 * Mtx12 - Mtx02 * Mtx11;

                    Mtx[1, 0] = -Mtx10 * Mtx22 + Mtx12 * Mtx20;
                    Mtx[1, 1] = Mtx00 * Mtx22 - Mtx02 * Mtx20;
                    Mtx[1, 2] = -Mtx00 * Mtx12 + Mtx02 * Mtx10;

                    Mtx[2, 0] = Mtx10 * Mtx21 - Mtx11 * Mtx20;
                    Mtx[2, 1] = -Mtx00 * Mtx21 + Mtx01 * Mtx20;
                    Mtx[2, 2] = Mtx00 * Mtx11 - Mtx01 * Mtx10;

                    Mtx.Scale(1.0 / det);
                    return;
                }

                default:
                unsafe
                {

                    // use LAPACK
                    //
                    int N = Mtx.NoOfCols;
                    int i0;
                    double[] _target_entries = TempBuffer.GetTempBuffer(out i0, N * N);
                    fixed (double* target_entries = _target_entries)
                    {
                        int* IPIV = stackalloc int[N];
                        CopyToUnsafeBuffer(Mtx, target_entries, true);

                        int info;
                        int lwork = -1;
                        if (N <= 64)
                            lwork = N * N;
                        else {
                            double twork = 0;
                            LAPACK.F77_LAPACK.DGETRI_(ref N, target_entries, ref N, IPIV, &twork, ref lwork, out info);
                            lwork = (int)twork;
                            if (info != 0) {
                                TempBuffer.FreeTempBuffer(i0);
                                throw new ArithmeticException("error in LAPACK routing DGETRI: info = " + info);
                            }
                        }

                        int i1;
                        double[] _work = TempBuffer.GetTempBuffer(out i1, lwork);
                        fixed (double* work = _work)
                        {

                            int LDA = N;
                            LAPACK.F77_LAPACK.DGETRF(ref N, ref N, target_entries, ref LDA, IPIV, out info);
                            if (info != 0) {
                                TempBuffer.FreeTempBuffer(i0);
                                TempBuffer.FreeTempBuffer(i1);
                                throw new ArithmeticException("error in LAPACK routing DGETRF: info = " + info);
                            }
                            LAPACK.F77_LAPACK.DGETRI_(ref N, target_entries, ref LDA, IPIV, work, ref lwork, out info);
                            if (info != 0) {
                                TempBuffer.FreeTempBuffer(i0);
                                TempBuffer.FreeTempBuffer(i1);
                                throw new ArithmeticException("error in LAPACK routing DGETRI: info = " + info);
                            }
                            //IntPtr mem = Marshal.AllocHGlobal(N * N * 8 * 2);
                            //Marshal.Copy(target.Entries, 0, mem, target.entries.Length);
                            //unsafe {
                            //    LAPACK.DGETRI(ref N, (double*) mem, ref LDA, IPIV, work, ref lwork, out info);
                            //}
                            //Marshal.Copy(mem, target.entries, 0, target.entries.Length);

                            CopyFromUnsafeBuffer(Mtx, target_entries, true);


                        }
                        TempBuffer.FreeTempBuffer(i1);
                    }
                    TempBuffer.FreeTempBuffer(i0);

                    return;
                }
            }
        }



        /// <summary>
        /// <paramref name="M"/> = <paramref name="M"/> + <paramref name="o"/>*<paramref name="sclo"/>
        /// </summary>
        static public void Acc<T>(this T M, double sclo, IMatrix o) where T : IMatrix {
            if (o.NoOfCols != M.NoOfCols || o.NoOfRows != M.NoOfRows)
                throw new ArgumentException("o must have same size and length as M-matrix");
            int m_NoOfRows = M.NoOfRows, m_NoOfCols = M.NoOfCols;

            for (int i = 0; i < m_NoOfRows; i++)
                for (int j = 0; j < m_NoOfCols; j++)
                    M[i, j] += sclo * o[i, j];
        }

        /// <summary>
        /// <paramref name="M"/> = <paramref name="M"/> + ID*<paramref name="sclo"/>, where ID is the identity matrix;
        /// </summary>
        static public void AccEye<T>(this T M, double sclo) where T : IMatrix {
            if (M.NoOfCols != M.NoOfRows)
                throw new NotSupportedException("must be quadratic.");

            int m_NoOfRows = M.NoOfRows;
            for (int i = 0; i < m_NoOfRows; i++)
                M[i, i] += sclo;
        }

        /// <summary>
        /// Least-squares-solve (LAPACK function DGELSY) with multiple right-hand-side vectors.
        /// </summary>
        /// <param name="Mtx">
        /// The matrix, size \f$ M \times N \f$, i.e. \f$ M \f$ is the number of equations anf \f$ N \f$ the number of unknowns. </param>
        /// <param name="B">
        /// Input for multiple RHS vectors, output for the corresponding solutions, size is \f$ \max\{ M, N \} \times L \f$,
        /// i.e. \f$ L \f$ is the number of right-hand-sides.
        /// </param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        static public void LeastSquareSolve<T, R>(this T Mtx, R B, double RCOND = 1.0e-14)
            where T : IMatrix
            where R : IMatrix {

            int NRHS = B.NoOfCols;

            if (B.NoOfRows != Math.Max(Mtx.NoOfCols, Mtx.NoOfRows))
                throw new ArgumentException("B must have max(I,J) number of rows, where I,J are the dimensions of Mtx.");

            unsafe {
                int M = Mtx.NoOfRows;
                int N = Mtx.NoOfCols;

                int i0, i1;
                double[] __RHS = TempBuffer.GetTempBuffer(out i1, B.NoOfRows * B.NoOfCols);
                double[] _this_entries = TempBuffer.GetTempBuffer(out i0, N * M);
                fixed(double* _RHS = __RHS, this_entries = _this_entries) {
                    CopyToUnsafeBuffer(B, _RHS, true);

                    int LDB = B.NoOfRows;
                    CopyToUnsafeBuffer(Mtx, this_entries, true);

                    LAPACK.F77_LAPACK.DGELSY(M, N, this_entries, _RHS, LDB, NRHS, RCOND);

                    CopyFromUnsafeBuffer(B, _RHS, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
            }
        }

        /// <summary>
        /// Least-squares-solve (LAPACK function DGELSY) with multiple right-hand-side vectors.
        /// </summary>
        /// <param name="Mtx">
        /// The matrix, size \f$ M \times N \f$, i.e. \f$ M \f$ is the number of equations anf \f$ N \f$ the number of unknowns. </param>
        /// <param name="B">
        /// Input for the right-hand-side, size is \f$  M  \times L \f$,
        /// i.e. \f$ L \f$ is the number of right-hand-sides.
        /// </param>
        /// <param name="X">
        /// Output for the solutions, size is \f$  N  \times L \f$,
        /// i.e. \f$ L \f$ is the number of right-hand-sides.
        /// </param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        static public void LeastSquareSolve<T, R, S>(this T Mtx, R X, S B, double RCOND = 1.0e-14)
            where T : IMatrix
            where R : IMatrix
            where S : IMatrix {

            int NRHS = B.NoOfCols;

            if (B.NoOfRows != Mtx.NoOfRows)
                throw new ArgumentException("Illegal number of rows for B.");
            if (X.NoOfRows != Mtx.NoOfCols)
                throw new ArgumentException("Illegal number of rows for X.");
            if (X.NoOfCols != NRHS)
                throw new ArgumentException("Illegal number of columns for X.");
            
            unsafe
            {
                int M = Mtx.NoOfRows;
                int N = Mtx.NoOfCols;

                int i0, i1;
                double[] __RHS = TempBuffer.GetTempBuffer(out i1, NRHS*Math.Max(M, N));
                double[] _this_entries = TempBuffer.GetTempBuffer(out i0, N * M);
                fixed (double* _RHS = __RHS, this_entries = _this_entries)
                {
                    CopyToUnsafeBuffer(B, _RHS, true);

                    int LDB = B.NoOfRows;
                    CopyToUnsafeBuffer(Mtx, this_entries, true);

                    LAPACK.F77_LAPACK.DGELSY(M, N, this_entries, _RHS, LDB, NRHS, RCOND);

                    CopyFromUnsafeBuffer(X, _RHS, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
            }
        }



        /// <summary>
        /// least-squares-solve (LAPACK function DGELSY).
        /// </summary>
        /// <param name="Mtx">Input, matrix of the least-squares system.</param>
        /// <param name="x">on exit, the least-squares solution to the system</param>
        /// <param name="b">right hand side</param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        static public void LeastSquareSolve<T>(this T Mtx, double[] x, double[] b, double RCOND = 1.0e-14) where T : IMatrix {
            if (x.Length != Mtx.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Length != Mtx.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows");

            int i0, i1;
            double[] _RHS = TempBuffer.GetTempBuffer(out i1, Math.Max(x.Length, b.Length));
            Array.Copy(b, _RHS, b.Length);

            int M = Mtx.NoOfRows;
            int N = Mtx.NoOfCols;
            int LDB = b.Length;
            unsafe {
                double[] _this_entries = TempBuffer.GetTempBuffer(out i0, N * M);
                fixed (double* p_RHS = _RHS, this_entries = _this_entries) {
                    CopyToUnsafeBuffer(Mtx, this_entries, true);

                    LAPACK.F77_LAPACK.DGELSY(M, N, this_entries, p_RHS, _RHS.Length, 1, RCOND);

                }
            }
            Array.Copy(_RHS, x, x.Length);
            TempBuffer.FreeTempBuffer(i0);
            TempBuffer.FreeTempBuffer(i1);
        }

        /// <summary>
        /// Solves the linear equation system:
        /// 
        /// <paramref name="M"/>*<paramref name="x"/> = <paramref name="b"/>.
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system.</param>
        /// <param name="b">Right-hand-side of the equation system.</param>
        /// <param name="M">General quadratic, non-singular matrix.</param>
        static public void Solve<T>(this T M, double[] x, double[] b) where T : IMatrix {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve nonquadratic matrix.");
            if (x.Length != M.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Length != M.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows");
            unsafe {

                int L = M.NoOfCols;

                Array.Copy(b, x, x.Length);

                int* ipiv = stackalloc int[L];
                int i0;
                double[] _this_Entries = TempBuffer.GetTempBuffer(out i0, L * L);
                fixed (double* this_Entries = _this_Entries) {

                    CopyToUnsafeBuffer(M, this_Entries, true);

                    int info;
                    LAPACK.F77_LAPACK.DGETRF(ref L, ref L, this_Entries, ref L, ipiv, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        string infostring;
                        if (info < 0) {
                            infostring = String.Format("the {0}-th argument had an illegal value", info);
                        }
                        else {
                            infostring = "U("+info+@""","""+info+
                                ") is exactly zero. The factorization \n has been completed, but the factor U is exactly \n singular, and division by zero will occur if it is used \n to solve a system of equations.";
                        }

                        throw new ArithmeticException("LAPACK dgetrf info: " + infostring);
                    }
                    //         TRANS, N, NRHS, A,            LDA, IPIV, B, LDB
                    char transp = 'N';
                    int eins = 1;
                    LAPACK.F77_LAPACK.DGETRS(ref transp, ref L, ref eins, this_Entries, ref L, ipiv, x, ref L, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK dgetrs info: " + info);
                    }
                }
                TempBuffer.FreeTempBuffer(i0);
            }
        }


        /// <summary>
        /// Solves the symmetric, positive definite linear equation system:
        /// 
        /// <paramref name="M"/>*<paramref name="x"/> = <paramref name="b"/>.
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system.</param>
        /// <param name="b">Right-hand-side of the equation system.</param>
        /// <param name="M">General quadratic, symmetric, positive definite matrix.</param>
        static public void SolveSymmetric<T>(this T M, double[] x, double[] b) where T : IMatrix {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve nonquadratic matrix.");
            if (x.Length != M.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Length != M.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows");
            unsafe
            {

                int L = M.NoOfCols;

                Array.Copy(b, x, x.Length);

                int* ipiv = stackalloc int[L];
                int i0;
                double[] _this_Entries = TempBuffer.GetTempBuffer(out i0, L * L);
                fixed (double* this_Entries = _this_Entries, _x = x)
                {
                    CopyToUnsafeBuffer(M, this_Entries, true);

                    int uplo = 'U', nrhs = 1, info;
                    LAPACK.F77_LAPACK.DPOSV_(ref uplo, ref L, ref nrhs, this_Entries, ref L, _x, ref L, out info);
                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK dposv info: " + info);
                    }
                }
                TempBuffer.FreeTempBuffer(i0);
            }
        }

        /// <summary>
        /// Computes the QR factorization of <paramref name="M"/> , i.e.
        /// determines an orthogonal matrix <paramref name="Q"/> and an upper
        /// triangular matrix <paramref name="R"/> such that
        /// <paramref name="M"/> = <paramref name="Q"/>*<paramref name="R"/>
        /// </summary>
        /// <typeparam name="T">
        /// The type of the matrix <paramref name="M"/>
        /// </typeparam>
        /// <param name="M">Input matrix to be factorized</param>
        /// <param name="Q">An orthogonal matrix</param>
        /// <param name="R">An upper triangular matrix</param>
        static public void GetQRFactorization<T>(this T M, out IMatrix Q, out IMatrix R) where T : IMatrix {
            int N = M.NoOfRows;
            if (N != M.NoOfCols) {
                throw new NotImplementedException(
                    "Only implemented for symmetric matrices so far");
            }

            Q = MultidimensionalArray.Create(N, N);
            R = MultidimensionalArray.Create(N, N);
            double[] TAU = new double[N];
            unsafe {
                int i0;
                double[] _M = TempBuffer.GetTempBuffer(out i0, N * N);
                fixed (double* pM = _M, pTAU = &TAU[0]) {
                    CopyToUnsafeBuffer(M, pM, true);

                    LAPACK.F77_LAPACK.DGEQRF(ref N, ref N, pM, ref N, pTAU);
                    CopyFromUnsafeBuffer(R, pM, true);

                    LAPACK.F77_LAPACK.DORGQR(ref N, ref N, ref N, pM, ref N, pTAU);
                    CopyFromUnsafeBuffer(Q, pM, true);
                }
                TempBuffer.FreeTempBuffer(i0);
            }

            // Get rid of entries below diagonal (they're meaningless now)
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < i; j++) {
                    R[i, j] = 0.0;
                }
            }
        }

        /// <summary>
        /// Computes the QR factorization of <paramref name="M"/> using column
        /// pivoting, i.e. determines an orthogonal matrix
        /// <paramref name="Q"/>, an upper triangular matrix
        /// <paramref name="R"/> and a permutation vector <paramref name="P"/>
        /// that defines the reordering of the columns of <paramref name="M"/>
        /// such that <paramref name="M"/>(P,:) =
        /// <paramref name="Q"/>*<paramref name="R"/> in Matlab notation. 
        /// </summary>
        /// <typeparam name="T">
        /// The type of the matrix <paramref name="M"/>
        /// </typeparam>
        /// <param name="M">Input matrix to be factorized</param>
        /// <param name="P">A permutation vector</param>
        /// <param name="Q">An orthogonal matrix</param>
        /// <param name="R">An upper triangular matrix</param>
        /// <param name="permuteColumnToFront">
        /// If null or all entries are 'false', the permutation is chosen such
        /// that the diagonal elements of <paramref name="R"/> are decreasing
        /// in magnitude.<br />
        /// If the i-th entry is 'true', enforces a permutation that maps
        /// column i of <paramref name="M"/> onto one of the first columns
        /// (depending on how many entries are 'true') of
        /// <paramref name="Q"/>*<paramref name="R"/>. All unrestricted columns
        /// will still be ordered according the magnitude of the diagonal
        /// elements of <paramref name="R"/>.
        /// </param>
        static public void GetPivotedQRFactorization<T>(this T M, out int[] P, out IMatrix Q, out IMatrix R, bool[] permuteColumnToFront = null) where T : IMatrix {
            int N = M.NoOfRows;
            if (N != M.NoOfCols) {
                throw new NotImplementedException(
                    "Only implemented for symmetric matrices so far");
            }

            P = new int[N];
            if (permuteColumnToFront != null) {
                if (permuteColumnToFront.Length != N) {
                    throw new ArgumentException(
                        "Length of 'premuteColumnToFront' must be equal to number of rows/columns",
                        "premuteColumnToFront");
                }

                for (int i = 0; i < N; i++) {
                    P[i] = permuteColumnToFront[i] ? 1 : 0;
                }
            }

            Q = MultidimensionalArray.Create(N, N);
            R = MultidimensionalArray.Create(N, N);
            double[] TAU = new double[N];
            unsafe {
                int i0;
                double[] _M = TempBuffer.GetTempBuffer(out i0, N * N);
                fixed (double* pM = _M, pTAU = &TAU[0]) {
                    fixed (int* pP = &P[0]) {
                        CopyToUnsafeBuffer(M, pM, true);

                        LAPACK.F77_LAPACK.DGEQP3(ref N, ref N, pM, ref N, pP, pTAU);
                        CopyFromUnsafeBuffer(R, pM, true);

                        LAPACK.F77_LAPACK.DORGQR(ref N, ref N, ref N, pM, ref N, pTAU);
                        CopyFromUnsafeBuffer(Q, pM, true);
                    }
                }
                TempBuffer.FreeTempBuffer(i0);
            }

            // Get rid of entries below diagonal (they're meaningless now)
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < i; j++) {
                    R[i, j] = 0.0;
                }
            }

            // Convert from 1-based indexing to 0-based indexing
            for (int i = 0; i < N; i++) {
                P[i]--;
            }
        }

        /// <summary>
        /// Converts an implicit subspace representation (given as the solution of a singular matrix <paramref name="Mtx"/>)
        /// into an explicit representation.
        /// </summary>
        /// <param name="Mtx">
        /// A matrix implicitly defines a subspace.
        /// </param>
        /// <returns>
        /// A matrix whose columns span the solution space of the input matrix <paramref name="Mtx"/>.
        /// </returns>
        public static MultidimensionalArray GetSolutionSpace<T>(this T Mtx) where T : IMatrix {
            int I = Mtx.NoOfRows;

            (var RRE, var pivots, var cols, int rank) = ReducedRowEchelonForm(Mtx);

            Console.WriteLine("The rank of the matrix coming from the Reduced Row Echelon Form: rank={0}", rank);

            if (cols.Length != Mtx.NoOfCols - rank)
                throw new ArithmeticException("Error in reduced row echelon form.");

            var S = MultidimensionalArray.Create(Mtx.NoOfCols, cols.Length);
            int j = 0;
            foreach (int qj in cols) {
                /*for(int i = 0; i < rank; i++) {
                    S[i, jj] = -RRE[i, j];
                }
                S[rank + jj, jj] = 1.0;*/
                S[qj, j] = 1.0;

                for (int i = 0; i < rank; i++) {
                    Debug.Assert(S[pivots[i], j] == 0);
                    S[pivots[i], j] = -RRE[i, qj];
                }


                j++;
            }

#if DEBUG
            var Test = Mtx.GEMM(S);
            double norm = Test.InfNorm();
            double thresh = Mtx.InfNorm();
            if(norm > thresh * 1.0e-8)
                throw new ArithmeticException("Gauss elimination is fucked.");
#endif

            return S;
        }


        /// <summary>
        /// Computes a reduced row echelon form
        /// </summary>
        /// <param name="Mtx">some matrix</param>
        /// <returns>
        /// A tuple, containing
        /// - the reduced row echelon form of <paramref name="Mtx"/>
        /// - the indices of the pivots
        /// - the indices of non-identity rows
        /// - the rank of <paramref name="Mtx"/>
        /// </returns>
        public static (MultidimensionalArray,int[], int[],int) ReducedRowEchelonForm<T>(this T Mtx) where T : IMatrix {
            var M = MultidimensionalArray.Create(Mtx.NoOfRows, Mtx.NoOfCols);
            M.Acc(1.0, Mtx);
            Mtx = default(T);
            int I = M.NoOfRows;
            int J = M.NoOfCols;
            double tol = BLAS.MachineEps * Math.Max(I, J) * M.InfNorm(); 

            var cols = new List<int>();
            var pivots = new List<int>();
            int i = 0; // row counter
            int rank = 0;
            for(int j = 0; j < J; j++) {
                // find the pivot row
                int i_pivot;
                double val_pivot;
                {
                    i_pivot = int.MinValue;
                    val_pivot = -1.0;
                    for (int ii = i; ii < I; ii++) {
                        if (Math.Abs(M[ii, j]) > val_pivot) {
                            i_pivot = ii;
                            val_pivot = Math.Abs(M[ii, j]);
                        }
                    }
                }

                if(val_pivot <= tol) {
                    // Skip column j, making sure the approximately zero terms are actually zero.
                    for (int ii = i; ii < I; ii++)
                        M[ii, j] = 0.0;

                    cols.Add(j);
                } else {
                    // found a pivot 
                    pivots.Add(j);
                    rank++;

                    // Swap current row and pivot row
                    for (int jj = j; jj < J; jj++) {
                        double a = M[i, jj];
                        M[i, jj] = M[i_pivot, jj];
                        M[i_pivot, jj] = a;
                    }

                    // normalize the pivot row
                    double al = 1.0 / M[i, j];
                    for(int jj = j; jj < J; jj++) {
                        M[i, jj] *= al;
                    }

                    // eliminate the current column
                    for(int ii = 0; ii < I; ii++) {
                        if(ii != i) {
                            // Row[ii] = Row[ii] - M[ii,j]*Row[i]
                            double b = M[ii, j];

                            for (int jj = j; jj < J; jj++) {
                                M[ii, jj] = M[ii, jj] - b * M[i, jj];
                            }
                        }
                    }

                    i++;
                    
                    if (i >= I) {
                        // finished
                        j++;
                        for(; j < J; j++) {
                            cols.Add(j);
                        }

                        break;
                    }
                }
            }

            Debug.Assert(rank == pivots.Count);
            Debug.Assert(cols.Count == M.NoOfCols - rank);
            return (M, pivots.ToArray(), cols.ToArray(), rank);
        }



        /// <summary>
        /// Computes the QR factorization of <paramref name="M"/> using column
        /// pivoting, i.e. determines an orthogonal matrix
        /// <paramref name="Q"/>, an upper triangular matrix
        /// <paramref name="R"/> and a permutation matrix <paramref name="P"/>
        /// such that <paramref name="P"/>*<paramref name="M"/> =
        /// <paramref name="Q"/>*<paramref name="R"/>. 
        /// </summary>
        /// <typeparam name="T">
        /// The type of the matrix <paramref name="M"/>
        /// </typeparam>
        /// <param name="M">Input matrix to be factorized</param>
        /// <param name="P">A permutation matrix</param>
        /// <param name="Q">An orthogonal matrix</param>
        /// <param name="R">An upper triangular matrix</param>
        /// <param name="permuteColumnToFront">
        /// If null or all entries are 'false', the permutation is chosen such
        /// that the diagonal elements of <paramref name="R"/> are decreasing
        /// in magnitude.<br />
        /// If the i-th entry is 'true', enforces a permutation that maps
        /// column i of <paramref name="M"/> onto one of the first columns
        /// (depending on how many entries are 'true') of
        /// <paramref name="Q"/>*<paramref name="R"/>. All unrestricted columns
        /// will still be ordered according the magnitude of the diagonal
        /// elements of <paramref name="R"/>.
        /// </param>
        static public void GetPivotedQRFactorization<T>(this T M, out IMatrix P, out IMatrix Q, out IMatrix R, bool[] permuteColumnToFront = null) where T : IMatrix {
            int[] PEconomy;
            M.GetPivotedQRFactorization(out PEconomy, out Q, out R, permuteColumnToFront);

            P = MultidimensionalArray.Create(M.NoOfCols, M.NoOfRows);
            for (int i = 0; i < PEconomy.Length; i++) {
                P[PEconomy[i], i] = 1.0;
            }
        }

        
        /// <summary>
        /// swaps rows <paramref name="i"/> and <paramref name="j"/>
        /// </summary>
        static public void SwapRows<T>(this T M, int i, int j) where T : IMatrix {
            M.CheckRowAndCol(i, 0);
            M.CheckRowAndCol(j, 0);
            double a;
            for (int l = 0; l < M.NoOfCols; l++) {
                a = M[i, l];
                M[i, l] = M[j, l];
                M[j, l] = a;
            }
        }

        /// <summary>
        /// swaps columns <paramref name="i"/> and <paramref name="j"/>
        /// </summary>
        static public void SwapColumns<T>(this T M, int i, int j) where T : IMatrix {
            M.CheckRowAndCol(0, i);
            M.CheckRowAndCol(0, j);
            double a;
            for (int l = 0; l < M.NoOfRows; l++) {
                a = M[l, i];
                M[l, i] = M[l, j];
                M[l, j] = a;
            }
        }

        /// <summary>
        /// multiplies row <paramref name="i"/> by a factor <paramref name="alpha"/>;
        /// </summary>
        static public void RowScale<T>(this T M, int i, double alpha) where T : IMatrix {
            M.CheckRowAndCol(i, 0);
            for (int l = 0; l < M.NoOfCols; l++)
                M[i, l] *= alpha;
        }

        /// <summary>
        /// multiplies column <paramref name="i"/> by a factor <paramref name="alpha"/>;
        /// </summary>
        static public void ColScale<T>(this T M, int i, double alpha) where T : IMatrix {
            M.CheckRowAndCol(0, i);
            for (int l = 0; l < M.NoOfRows; l++)
                M[l, i] *= alpha;
        }

        /// <summary>
        /// accumulates row <paramref name="iSrc"/> times
        /// <paramref name="alpha"/> to row <paramref name="iDst"/>;
        /// </summary>
        static public void RowAdd<T>(this T M, int iSrc, int iDst, double alpha) where T : IMatrix {
            M.CheckRowAndCol(iDst, 0);
            M.CheckRowAndCol(iSrc, 0);
            for (int l = 0; l < M.NoOfCols; l++)
                M[iDst, l] += M[iSrc, l] * alpha;
        }

        /// <summary>
        /// accumulates column <paramref name="iSrc"/> times
        /// <paramref name="alpha"/> to column <paramref name="iDst"/>;
        /// </summary>
        static public void ColAdd<T>(this T M, int iSrc, int iDst, double alpha) where T : IMatrix {
            M.CheckRowAndCol(0, iDst);
            M.CheckRowAndCol(0, iSrc);
            for (int l = 0; l < M.NoOfRows; l++)
                M[l, iDst] += M[l, iSrc] * alpha;
        }

        /// <summary>
        /// throws an exception if either column or row index are our of range
        /// </summary>
        static private void CheckRowAndCol<T>(this T M, int i, int j) where T : IMatrix {
            if (i < 0 || i >= M.NoOfRows)
                throw new IndexOutOfRangeException("row index out of range");
            if (j < 0 || j >= M.NoOfCols)
                throw new IndexOutOfRangeException("column index out of range");
        }
        
               

        /// <summary>
        /// Accumulates  <paramref name="alpha"/>*<paramref name="row"/> to the <paramref name="RowNo"/>-th row of <paramref name="inp"/>.
        /// </summary>
        /// <param name="inp">
        /// Matrix that should be altered.
        /// </param>
        /// <param name="RowNo">
        /// Row index of the row to set.
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
        /// <param name="alpha">
        /// Scaling factor.
        /// </param>
        public static void AccRow<T>(this IMatrix inp, int RowNo, double alpha, T row, int i0 = 0, int inc = 1) where T : IList<double> {

            if(RowNo < 0 || RowNo >= inp.NoOfRows)
                throw new IndexOutOfRangeException("RowNo out of range");
            if((row.Count - i0) < inp.NoOfCols * inc)
                throw new ArgumentException("array to short", "row");

            int I1 = inp.NoOfCols;
            for(int i = 0; i < I1; i++)
                inp[RowNo, i] += row[i * inc + i0]*alpha;
        }

        /// <summary>
        /// clears the <paramref name="RowNo"/>-th row from
        /// <paramref name="inp"/> .
        /// </summary>
        /// <param name="inp">
        /// matrix that should be altered
        /// </param>
        /// <param name="RowNo">
        /// row index of the row to set
        /// </param>
        public static void ClearRow(this IMatrix inp, int RowNo) {

            if(RowNo < 0 || RowNo >= inp.NoOfRows)
                throw new IndexOutOfRangeException("RowNo out of range");

            int I1 = inp.NoOfCols;
            for(int i = 0; i < I1; i++)
                inp[RowNo, i] = 0.0;
        }

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
        public static void SetColumn<T>(this IMatrix inp, int ColNo, T col, int i0 = 0, int inc = 1) where T : IEnumerable<double> {

            if(ColNo < 0 || ColNo >= inp.NoOfCols)
                throw new IndexOutOfRangeException("RowNo out of range");
            if((col.Count() - i0) < inp.NoOfRows * inc)
                throw new ArgumentException("array to short", "row");

            int I1 = inp.NoOfRows;
            int i = 0;
            if(inc == 1 && i0 == 0) {
                foreach(double col_i in col) {
                    inp[i, ColNo] = col_i;
                    i++;
                }
            } else {
                int j = 0;
                foreach(double col_i in col) {
                    if(i >= i0 && ((j - i0) % inc == 0)) {
                        inp[i, ColNo] = col_i;
                        i++;
                    }
                    j++;
                }
            }
        }


        
   
        /*
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
        */

            /*
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
        */

    }
}
