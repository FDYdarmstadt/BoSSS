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
using System.Runtime.CompilerServices;

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
        /// <param name="txt">Text read on the stream to load from.</param>
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
        /// Finds those row of <paramref name="mda"/>, where the L2-distance between the row and <paramref name="Row"/> is minimal. 
        /// </summary>
        /// <returns>
        /// the row index of the minimum-distance row.
        /// </returns>
        static public int MindistRow(this IMatrix mda, Vector Row) {
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
        /// Finds those row <paramref name="iDmax"/> of <paramref name="mda"/>, where the L2-distance between the row and <paramref name="Row"/> is minimal. 
        /// </summary>
        /// <param name="mda">some matrix</param>
        /// <param name="Row">some row</param>
        /// <param name="Dmax">minimum L2 distance over all rows</param>
        /// <param name="iDmax">index of minimum L2-distance row</param>
        /// <returns></returns>
        static public void MindistRow(this IMatrix mda, Vector Row, out double Dmax, out int iDmax) {
            if(Row.Dim != mda.NoOfCols)
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
        static public int FindRow(this IMatrix mda, Vector Row, double tol) {
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

        /// <summary>
        /// From matrix <paramref name="M"/> to <paramref name="buffer"/>
        /// </summary>
        /// <param name="M">input matrix</param>
        /// <param name="buffer">
        /// output;
        /// </param>
        /// <param name="BufferInFortranOrder">
        /// - true:  <paramref name="buffer"/> will be written in FORTRAN order (column-wise)
        /// - false:  <paramref name="buffer"/> will be written in C order (row-wise)
        /// </param>
        internal static unsafe void CopyToUnsafeBuffer<T>(T M, double* buffer, bool BufferInFortranOrder) where T : IMatrix {
#if DEBUG
            if(M.GetType().IsValueType)
                throw new NotSupportedException("CopyTo value type -- probably not the expected result! (Using vector struct in CopyTo(...) - operation?)");
#endif

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

        /// <summary>
        /// From <paramref name="buffer"/> to matrix <paramref name="M"/>
        /// </summary>
        /// <param name="M">output matrix</param>
        /// <param name="buffer">
        /// output;
        /// </param>
        /// <param name="BufferInFortranOrder">
        /// - true: <paramref name="buffer"/> will be read in FORTRAN order (column-wise)
        /// - false: <paramref name="buffer"/> will be read in C order (row-wise)
        /// </param>
        internal static unsafe void CopyFromUnsafeBuffer<T>(T M, double* buffer, bool BufferInFortranOrder) where T : IMatrix {
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
        static public void GEMV<MatrixType, VectorType1, VectorType2>(this MatrixType M, double xScaling, VectorType1 x, double yScaling, VectorType2 y, bool transpose = false)
            where MatrixType : IMatrix
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> //
        {
            int m_NoOfCols = M.NoOfCols, m_NoOfRows = M.NoOfRows;
            if (!transpose) {
                if (m_NoOfCols != x.Count)
                    throw new ArgumentException("Length of vector x must be equal to number of columns.");
                if (m_NoOfRows != y.Count)
                    throw new ArgumentException("Length of vector y must be equal to number of rows.");
            } else {
                if (m_NoOfRows != x.Count)
                    throw new ArgumentException("Length of vector x must be equal to number of rows (transpose matrix multiply was selected!).");
                if (m_NoOfCols != y.Count)
                    throw new ArgumentException("Length of vector y must be equal to number of columns (transpose matrix multiply was selected!).");
            }
            if (m_NoOfCols == 0 || m_NoOfRows == 0)
                return;
            if (ReferenceEquals(x, y))
                throw new ArgumentException("in-place matrix-vector product is not supported");


            if (typeof(MatrixType) == typeof(MultidimensionalArray) && m_NoOfRows >= 2 && m_NoOfCols >= 2) {
                // +++++++++++++++++
                // optimized version
                // +++++++++++++++++


                MultidimensionalArray mdaM = M as MultidimensionalArray;
                if (mdaM.Dimension != 2)
                    throw new ArgumentException("Multidimensional Array must have 2 dimensions to be a matrix.");

                int off = mdaM.Index(0, 0);
                int LD = mdaM.Index(1, 0) - off; // leading dimension
                int SD = mdaM.Index(0, 1) - off; // small dimension

                double[] _x;
                if (typeof(VectorType1) == typeof(double[]))
                    _x = x as double[];
                else
                    _x = x.ToArray();

                //double[] yCeck = y.ToArray();

                if (SD == 1 && m_NoOfCols * m_NoOfRows >= 32) {
                    // ++++++++++++++++
                    // BLAS can be used
                    // ++++++++++++++++


                    double[] _y;
                    bool backCopy = false;
                    if (typeof(VectorType1) == typeof(double[])) {
                        _y = y as double[];
                    } else {
                        _y = y.ToArray();
                        backCopy = true;
                    }


                    unsafe {
                        fixed (double* _pmdaM = mdaM.Storage, px = _x, py = _y) {
                            double* pmdaM = _pmdaM + off;
                            // because of FORTRAN-vs-C arrays, the transpose passed to BLAS is inverted:
                            BLAS.dgemv(transpose ? 'N' : 'T', m_NoOfCols, m_NoOfRows, xScaling, pmdaM, LD, px, 1, yScaling, py, 1);

                        }
                    }

                    

                    if (backCopy) {
                        int I = _y.Length;
                        for (int i = 0; i < I; i++) {
                            y[i] = _y[i];
                        }
                    }

                    /*
                    {


                        if (!transpose) {

                            for (int i = 0; i < m_NoOfRows; i++) {
                                double yi = 0;

                                for (int j = 0; j < m_NoOfCols; j++)
                                    yi += M[i, j] * x[j];

                                yCeck[i] = yCeck[i] * yScaling + yi * xScaling;
                            }

                        } else {

                            for (int i = 0; i < m_NoOfCols; i++) {
                                double yi = 0;

                                for (int j = 0; j < m_NoOfRows; j++)
                                    yi += M[j, i] * x[j];

                                yCeck[i] = yCeck[i] * yScaling + yi * xScaling;
                            }
                        }

                        double Rel = Math.Max(y.L2Norm(), yCeck.L2Norm())*1e-7;

                        if (yCeck.L2Distance(y) / Rel > 1.0e-7)
                            throw new ArithmeticException("BLAS GEMV is fucked.");
                    }
                    */

                } else {
                    // ++++++++++++++++++++++++++++++++++++++++++++
                    // use internal implementation
                    // (either because very small or spread matrix)
                    // ++++++++++++++++++++++++++++++++++++++++++++

                    unsafe {
                        fixed (double* _pmdaM = mdaM.Storage, _px = _x) {

                            
                            if (!transpose) {
                                double* pmdaM = _pmdaM + off;
                                for (int i = 0; i < m_NoOfRows; i++) {
                                    double yi = 0;
                                    double* px = _px;

                                    for (int j = 0; j < m_NoOfCols; j++) {
                                        yi += *pmdaM * *px;
                                        px++;
                                        pmdaM += SD;
                                    }

                                    y[i] = y[i] * yScaling + yi * xScaling;
                                }

                            } else {

                                
                                for (int i = 0; i < m_NoOfCols; i++) {
                                    double* pmdaM = _pmdaM + off + i*SD;
                                    double* px = _px;
                                    double yi = 0;

                                    for (int j = 0; j < m_NoOfRows; j++) {
                                        yi += *pmdaM * *px;
                                        px++;
                                        pmdaM += LD;
                                    }
                                    y[i] = y[i] * yScaling + yi * xScaling;
                                }
                                
                            }
                        }
                    }

                    /*
                    {


                        if (!transpose) {

                            for (int i = 0; i < m_NoOfRows; i++) {
                                double yi = 0;

                                for (int j = 0; j < m_NoOfCols; j++)
                                    yi += M[i, j] * x[j];

                                yCeck[i] = yCeck[i] * yScaling + yi * xScaling;
                            }

                        } else {

                            for (int i = 0; i < m_NoOfCols; i++) {
                                double yi = 0;

                                for (int j = 0; j < m_NoOfRows; j++)
                                    yi += M[j, i] * x[j];

                                yCeck[i] = yCeck[i] * yScaling + yi * xScaling;
                            }
                        }

                        double Rel = Math.Max(y.L2Norm(), yCeck.L2Norm()) * 1e-7;

                        if (yCeck.L2Distance(y) / Rel > 1.0e-7)
                            throw new ArithmeticException("own GEMV is fucked.");
                    }
                    */
                }
            } else {
                // +++++++++++++++++
                // Reference version
                // +++++++++++++++++


                

                if (!transpose) {
                    
                    for (int i = 0; i < m_NoOfRows; i++) {
                        double yi = 0;

                        for (int j = 0; j < m_NoOfCols; j++)
                            yi += M[i, j] * x[j];

                        y[i] = y[i] * yScaling + yi * xScaling;
                    }

                } else {
                    
                    for (int i = 0; i < m_NoOfCols; i++) {
                        double yi = 0;

                        for (int j = 0; j < m_NoOfRows; j++)
                            yi += M[j, i] * x[j];

                        y[i] = y[i] * yScaling + yi * xScaling;
                    }
                }
            }
        }

        /// <summary>
        /// alias for <see cref="GEMM{Matrix2, Matrix3}(Matrix2, Matrix3)"/>
        /// </summary>
        static public MultidimensionalArray MatMatMul<Matrix2, Matrix3>(this Matrix2 A, Matrix3 B)
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix //    
        {
            return A.GEMM(B);
        }

        /// <summary>
        /// Alias for <see cref="GEMV"/>
        /// </summary>
        static public void MatVecMul<MatrixType, VectorType1, VectorType2>(this MatrixType M, double xScaling, VectorType1 x, double yScaling, VectorType2 y, bool transpose = false)
            where MatrixType : IMatrix
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> //
        {
            M.GEMV(xScaling, x, yScaling, y, transpose);
        }

        /// <summary>
        /// Alias for <see cref="GEMV"/>
        /// </summary>
        static public double[] MatVecMul<MatrixType, VectorType1>(this MatrixType M, double xScaling, VectorType1 x, bool transpose = false)
            where MatrixType : IMatrix
            where VectorType1 : IList<double>//
        {
            double[] y = new double[transpose ? M.NoOfCols : M.NoOfRows];
            M.GEMV(xScaling, x, 0.0, y, transpose);
            return y;
        }


        /// <summary>
        /// Alias for <see cref="GEMV"/>
        /// </summary>
        static public void MatVecMulInplace<MatrixType, VectorType1>(this MatrixType M, double xScaling, VectorType1 x, bool transpose = false)
            where MatrixType : IMatrix
            where VectorType1 : IList<double> //
        {
            double[] xc = x.ToArray();
            M.GEMV(xScaling, xc, 0.0, x, transpose);
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

        static MultidimensionalArray.MultiplyProgram GEMMnn_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ik", "kj"); // A*B
        static MultidimensionalArray.MultiplyProgram GEMMtn_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ki", "kj"); // A^T * B
        static MultidimensionalArray.MultiplyProgram GEMMnt_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ik", "jk"); // A   * B^T
        static MultidimensionalArray.MultiplyProgram GEMMtt_Prog = MultidimensionalArray.MultiplyProgram.Compile("ij", "ki", "jk"); // A^T * B^T

        /// <summary>
        /// General matrix/matrix multiplication:
        /// 
        /// <paramref name="C"/> = <paramref name="alpha"/>*<paramref name="A"/>^x*<paramref name="B"/>^y + <paramref name="beta"/>*<paramref name="C"/>,
        /// 
        /// where x and y are either T (transpose) or 1, depending on <paramref name="transA"/> and <paramref name="transB"/>, respectively.
        /// </summary>
        static public void GEMM<Matrix1, Matrix2, Matrix3>(this Matrix1 C, double alpha, Matrix2 A, Matrix3 B, double beta, bool transA = false, bool transB = false)
            where Matrix1 : IMatrix
            where Matrix2 : IMatrix
            where Matrix3 : IMatrix //
        {

            if (!transA && !transB) {
                if (A.NoOfCols != B.NoOfRows)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfRows != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != M.NoOfRows", "A,M");
                if (B.NoOfCols != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != M.NoOfCols", "B,M");
                if (A.NoOfCols == 0)
                    return;
            } else if (transA && !transB) {
                if (A.NoOfRows != B.NoOfRows)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfCols != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != M.NoOfRows", "A,M");
                if (B.NoOfCols != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != M.NoOfCols", "B,M");
                if (A.NoOfRows == 0)
                    return;
            } else if (!transA && transB) {
                if (A.NoOfCols != B.NoOfCols)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfRows != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != M.NoOfRows", "A,M");
                if (B.NoOfRows != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != M.NoOfCols", "B,M");
                if (A.NoOfCols == 0)
                    return;
            } else if (transA && transB) {
                if (A.NoOfRows != B.NoOfCols)
                    throw new ArgumentException("A.NoOfCols != B.NoOfRows", "A,B");
                if (A.NoOfCols != C.NoOfRows)
                    throw new ArgumentException("A.NoOfRows != M.NoOfRows", "A,M");
                if (B.NoOfRows != C.NoOfCols)
                    throw new ArgumentException("B.NoOfCols != M.NoOfCols", "B,M");
                if (A.NoOfRows == 0)
                    return;
            }

            bool inPlace = false;
            if (object.ReferenceEquals(C, A))
                inPlace = true;
            if (object.ReferenceEquals(C, B))
                inPlace = true;
            if (C.NoOfCols == 0 || C.NoOfRows == 0)
                return;
             

            if (!inPlace && A is MultidimensionalArray _A && B is MultidimensionalArray _B && C is MultidimensionalArray _C) {
                int a00 = _A.Index(0, 0);
                int b00 = _B.Index(0, 0);
                int c00 = _C.Index(0, 0);

                if (_A.NoOfCols > 1 && (_A.Index(0, 1) - a00 == 1) && _B.NoOfCols > 1 && (_B.Index(0, 1) - b00 == 1) && _C.NoOfCols > 1 && (_C.Index(0, 1) - c00 == 1)) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Data layout is suitable to use BLAS DGEMM
                    // directly on MultidimenasionalArray storage.
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    /*
                    var __C = _C.CloneAs();
                    if (!transA && !transB) {
                        __C.Multiply(alpha, _A, _B, beta, ref GEMMnn_Prog);
                    } else if (transA && !transB) {
                        __C.Multiply(alpha, _A, _B, beta, ref GEMMtn_Prog);
                    } else if (!transA && transB) {
                        __C.Multiply(alpha, _A, _B, beta, ref GEMMnt_Prog);
                    } else if (transA && transB) {
                        __C.Multiply(alpha, _A, _B, beta, ref GEMMtt_Prog);
                    }
                    */

                    unsafe {
                        fixed(double* _pA = _A.Storage, _pB = _B.Storage, _pC = _C.Storage) {
                            double* pA = _pA + a00, pB = _pB + b00, pC = _pC + c00;

                            // C-order vs. FORTRAN-order
                            // Note that we are using FORTRAN BLAS, while BoSSS stores in C-order;
                            // therefore, we have to trick with A, B and the transposition

                           
                            int TRANSA = transB ? 't' : 'n';
                            int TRANSB = transA ? 't' : 'n';

                            int M = !transB ? _B.NoOfCols : _B.NoOfRows;
                            int N = !transA ? _A.NoOfRows : _A.NoOfCols;
                            int K = !transB ? _B.NoOfRows : _B.NoOfCols;

   
                            int LDA = (TRANSA == 'n') ? Math.Max(1, M) : Math.Max(1, K);
                            int LDB = (TRANSB == 'n') ? Math.Max(1, K) : Math.Max(1, N);
                            int LDC = Math.Max(1, M);


                            BLAS.dgemm(TRANSA, TRANSB, M, N, K, alpha, pB, _B.GetCycle(0), pA, _A.GetCycle(0), beta, pC, _C.GetCycle(0));

                        }
                    }

                    /*
                    var ERR = _C.CloneAs();
                    ERR.Acc(-1.0, __C);
                    double gemmErr = ERR.L2Norm();
                    double denom = __C.L2Norm();
                    if (gemmErr > 1e-6) {
                        Console.WriteLine("gemm error " + gemmErr + ",  denom = " + denom + ", rel = " + (gemmErr / denom));

                        A.SaveToTextFile("A.txt");
                        B.SaveToTextFile("B.txt");
                        _C.SaveToTextFile("C1.txt");
                        __C.SaveToTextFile("C2.txt");

                        unsafe {
                            fixed (double* _pA = _A.Storage, _pB = _B.Storage, _pC = _C.Storage) {
                                double* pA = _pA + a00, pB = _pB + b00, pC = _pC + c00;

                                // C-order vs. FORTRAN-order
                                // Note that we are using FORTRAN BLAS, while BoSSS stores in C-order;
                                // therefore, we have to trick with A, B and the transposition


                                int TRANSA = transB ? 't' : 'n';
                                int TRANSB = transA ? 't' : 'n';
                                
                                int M = !transB ? _B.NoOfCols : _B.NoOfRows;
                                int N = !transA ? _A.NoOfRows : _A.NoOfCols;
                                int K = !transB ? _B.NoOfRows : _B.NoOfCols;

                                //int M = transA ? A.NoOfCols : A.NoOfRows;
                                //int N = transB ? B.NoOfRows : B.NoOfCols;
                                //int K = transA ? A.NoOfRows : A.NoOfCols;

                                int LDA = transB ? Math.Max(1, K) : Math.Max(1, M);
                                int LDB = transA ? Math.Max(1, N) : Math.Max(1, K);

                                int _LDA = (TRANSA == 'n') ? Math.Max(1, M) : Math.Max(1, K);
                                int _LDB = (TRANSB == 'n') ? Math.Max(1, K) : Math.Max(1, N);

                                int LDC = Math.Max(1, M);


                                BLAS.dgemm(TRANSA, TRANSB, M, N, K, alpha, pB, _B.GetCycle(0), pA, _A.GetCycle(0), beta, pC, _C.GetCycle(0));

                            }
                        }

                    }
                    _C.Set(__C);
                    */

                    return;
                } else if(_C.NoOfRows*_C.NoOfCols <= 512 && _A.NoOfRows * _A.NoOfCols <= 512 && _B.NoOfRows * _B.NoOfCols <= 512) {
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // Multidimensional array are not suitable for direst use of BLAS DGEMM
                    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    // Note: (fk, 04jul22)
                    // simple test shows, below 512 entries (for nearly quadratic matrix C), the overhead from copying data 
                    // out-weights the performance gain from BLAS DGEMM
                    // Obviously, this is processor-dependent and also might be different for very tall or shallow matrix shapes,
                    // but this is good enough for BoSSS;
                    //

                    // probably, the 'multiply'-program is faster than 'BLAS-DGEMM + copy overhead'?

                    if (!transA && !transB) {
                        _C.Multiply(alpha, _A, _B, beta, ref GEMMnn_Prog);
                    } else if (transA && !transB) {
                        _C.Multiply(alpha, _A, _B, beta, ref GEMMtn_Prog);
                    } else if (!transA && transB) {
                        _C.Multiply(alpha, _A, _B, beta, ref GEMMnt_Prog);
                    } else if (transA && transB) {
                        _C.Multiply(alpha, _A, _B, beta, ref GEMMtt_Prog);
                    }

                    return;
                }
            }


            {
                // +++++++++++++++++++++++++++++++++++++++++++++++++
                // not returned yet;
                // data must be copied before we can use BLAS DGEMM
                // +++++++++++++++++++++++++++++++++++++++++++++++++

                
                unsafe {
                    int TRANSA = transA ? 't' : 'n';
                    int TRANSB = transB ? 't' : 'n';

                    int M = transA ? A.NoOfCols : A.NoOfRows;
                    int N = transB ? B.NoOfRows : B.NoOfCols;
                    int K = transA ? A.NoOfRows : A.NoOfCols;

                    int LDA = transA ? Math.Max(1, K) : Math.Max(1, M);
                    int LDB = transB ? Math.Max(1, N) : Math.Max(1, K);
                    int LDC = Math.Max(1, M);

                    int i0, i1, i2;
                    double[] __A = TempBuffer.GetTempBuffer(out i0, A.NoOfRows * A.NoOfCols);
                    double[] __B = TempBuffer.GetTempBuffer(out i1, B.NoOfRows * B.NoOfCols);
                    double[] __C = TempBuffer.GetTempBuffer(out i2, C.NoOfRows * C.NoOfCols);
                    fixed (double* pA = __A, pB = __B, pC = __C) {
                        CopyToUnsafeBuffer(A, pA, true);
                        CopyToUnsafeBuffer(B, pB, true);
                        CopyToUnsafeBuffer(C, pC, true);

                        BLAS.dgemm(TRANSA, TRANSB, M, N, K, alpha, pA, LDA, pB, LDB, beta, pC, LDC);

                        CopyFromUnsafeBuffer(C, pC, true);
                    }
                    TempBuffer.FreeTempBuffer(i0);
                    TempBuffer.FreeTempBuffer(i1);
                    TempBuffer.FreeTempBuffer(i2);
                }
            }
        }

       

        


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
        static public MultidimensionalArray TransposeTo<TM>(this TM M) where TM : IMatrix {
            MultidimensionalArray result = MultidimensionalArray.Create(M.NoOfCols, M.NoOfRows);

            M.TransposeTo(result);

            return result;
        }

        /*
        /// <summary>
        /// Horizontal matrix concatenation
        /// </summary>
        static public MultidimensionalArray CatHoriz<TM1, TM2>(this TM1 LeftBlock, TM2 RightBlock)
            where TM1 : IMatrix
            where TM2 : IMatrix // 
        {
            var C = new IMatrix[] { LeftBlock, RightBlock };
            return C.CatHoriz();
        }
        */

        /// <summary>
        /// Horizontal concatenation
        /// </summary>
        static public MultidimensionalArray CatHoriz<TM>(this IEnumerable<TM> Block) 
            where TM : IMatrix // 
        {
            var C = new IMatrix[1, Block.Count()];
            int i = 0;
            foreach(var M in Block) {
                C[0, i] = M;
                i++;
            }

            return C.Cat();
        }

        /// <summary>
        /// Horizontal concatenation
        /// </summary>
        static public MultidimensionalArray CatHoriz<TM1,TM2>(this TM1 LeftestBlock, params TM2[] otherBlocks) 
            where TM1 : IMatrix
            where TM2 : IMatrix // 
        {
            var C = new List<IMatrix>();
            C.Add(LeftestBlock);
            foreach(var M in otherBlocks)
                C.Add(M);

            return C.CatHoriz();
        }

        /// <summary>
        /// Vertical concatenation
        /// </summary>
        static public MultidimensionalArray CatVert<TM>(this IEnumerable<TM> Block) 
            where TM : IMatrix // 
        {
            var C = new IMatrix[Block.Count(), 1];
            int i = 0;
            foreach(var M in Block) {
                C[i, 0] = M;
                i++;
            }

            return C.Cat();
        }

        /// <summary>
        /// Horizontal concatenation
        /// </summary>
        static public MultidimensionalArray CatVert<TM1,TM2>(this TM1 TopBlock, params TM2[] otherBlocks) 
            where TM1 : IMatrix
            where TM2 : IMatrix // 
        {
            var C = new List<IMatrix>();
            C.Add(TopBlock);
            foreach(var M in otherBlocks)
                C.Add(M);

            return C.CatVert();
        }

        /// <summary>
        /// Concatenates a bloc of matrices
        /// </summary>
        static public MultidimensionalArray Cat<TM>(this TM[,] block) where TM : IMatrix {
            int I = block.GetLength(0);
            int J = block.GetLength(1);

            // determine lengths of row and columns
            // ====================================
            int[] RowLengths = new int[I];
            int[] ColLengths = new int[J];

            RowLengths.SetAll(-1);
            ColLengths.SetAll(-1);
            for(int i = 0; i < block.GetLength(0); i++) { // loop over block rows...
                for(int j = 0; j < block.GetLength(1); j++) { // loop over block columns...
                    
                    if(block[i, j] == null)
                        continue;

                    if(RowLengths[i] < 0) {
                        int IB = block[i, j].NoOfRows;
                        RowLengths[i] = IB;
                    }

                    if(ColLengths[j] < 0) {
                        int JB = block[i, j].NoOfCols;
                        ColLengths[j] = JB;
                    } else {
                        break;
                    }
                }
            }
            for(int i = 0; i < I; i++) {
                int IB = RowLengths[i];
                if(IB < 0) {
                    throw new ArgumentException($"Unable to determine width of block-row {i}; all entries in row are null.");
                }
            }
            for(int j = 0; j < J; j++) {
                int JB = ColLengths[j];
                if(JB < 0) {
                    throw new ArgumentException($"Unable to determine width of block-column {j}; all entries in column are null.");
                }
            }

            int[] RowOffsets = new int[I];
            int[] ColOffsets = new int[J];
            for(int i = 1; i < I; i++) {
                RowOffsets[i] = RowOffsets[i - 1] + RowLengths[i - 1];
            }
            for(int j = 1; j < J; j++) {
                ColOffsets[j] = ColOffsets[j - 1] + ColLengths[j - 1];
            }

            // check
            // =====
            for(int i = 0; i < block.GetLength(0); i++) { // loop over block rows...
                for(int j = 0; j < block.GetLength(1); j++) { // loop over block columns...
                    int _IB, _JB;
                    if(block[i, j] == null) {
                        //_IB = 0;
                        //_JB = 0;
                        continue;
                    } else {
                        _IB = block[i, j].NoOfRows;
                        _JB = block[i, j].NoOfCols;
                    }

                    int IB = RowLengths[i];
                    int JB = ColLengths[j];

                    if(_IB != IB) {
                        throw new ArgumentException($"Mismatch in number of rows for block {i},{j}: this block has {_IB} rows, while other blocks in block-row {i} have {IB} rows.");
                    }
                    if(_JB != JB) {
                        throw new ArgumentException($"Mismatch in number of columns for block {i},{j}: this block has {_JB} columns, while other blocks in block-column {j} have {JB} columns.");
                    }
                }
            }

            // copy data
            // =========
            int TotRows = RowOffsets[I - 1] + RowLengths[I - 1];
            int TotCols = ColOffsets[J - 1] + ColLengths[J - 1];
            var Ret = MultidimensionalArray.Create(TotRows, TotCols);
            for(int i = 0; i < block.GetLength(0); i++) { // loop over block rows...
                for(int j = 0; j < block.GetLength(1); j++) { // loop over block columns...
                    if(block[i, j] == null)
                        continue;

                    int i0 = RowOffsets[i];
                    int iE = i0 + RowLengths[i] - 1;
                    int j0 = ColOffsets[j];
                    int jE = j0 + ColLengths[j] - 1;

                    Ret.ExtractSubArrayShallow(new[] { i0, j0 }, new[] { iE, jE }).AccMatrix(1.0, block[i, j]);
                }
            }

            // return
            // ======
            return Ret;
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
        /// newly allocated matrix. Use <see cref="InvertInPlace{M1}(M1)"/> to
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
        /// <param name="M">input/output</param>
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
                        var Test = GEMM(GEMM(B.TransposeTo(), MtxClone), B);
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

        /// <summary>
        /// Performs a Gram-Schmidt orthonormalization on a matrix <paramref name="Mtx"/>
        /// </summary>
        /// <param name="Mtx"></param>
        /// <param name="B">
        /// on exit, an upper triangular matrix so that 
        /// <paramref name="B"/>^T * <paramref name="Mtx"/> * <paramref name="B"/> is an identity matrix.
        /// </param>
        /// <param name="Diag">
        /// Not used, if Gram-Schmid performs as normal, i.e. if <paramref name="Mtx"/> is positive definite,
        /// If the algorithm runs into some trouble, i.e. if the matrix is positive definite or close to indefinitenes, 
        /// a negative norm may occur; 
        /// </param>
        static public void GramSchmidt<FullMatrix1, FullMatrix2>(this FullMatrix1 Mtx, FullMatrix2 B, double[] Diag)
            where FullMatrix1 : IMatrix
            where FullMatrix2 : IMatrix //
        {
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
        /// Calculates the inverse of this matrix and returns it without modifying the <paramref name="source"/>.
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
        static public void InvertInPlace<M1>(this M1 Mtx)
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
        /// Accumulates some vector to the diagonal of the matrix.
        /// </summary>
        static public void AccDiag<T,V>(this T M, V diag) 
            where T : IMatrix 
            where V: IEnumerable<double> //
        {
            if (M.NoOfCols != M.NoOfRows)
                throw new NotSupportedException("must be quadratic.");

            int i = 0;
            foreach(double d in diag) {
                M[i, i] += d;
                i++;
            }
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
        /// The matrix, size \f$ M \times N \f$, i.e. \f$ M \f$ is the number of equations and \f$ N \f$ the number of unknowns. </param>
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
        /// <param name="b">right hand side</param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        static public double[] LeastSquareSolve<T>(this T Mtx, double[] b, double RCOND = 1.0e-14) where T : IMatrix {
            double[] x = new double[Mtx.NoOfCols];
            Mtx.LeastSquareSolve(x, b, RCOND);
            return x;
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

        //static public Stopwatch DGETRF_stopwatch;// = new Stopwatch();

        /// <summary>
        /// Solves the linear equation system:
        /// 
        /// <paramref name="M"/>*<paramref name="x"/> = <paramref name="b"/>.
        /// </summary>
        static public double[] Solve<T,W>(this T M, W b)
            where T : IMatrix
            where W : IList<double> //
        {
            double[] x = new double[b.Count];
            M.Solve(x, b);
            return x;
        }

        /// <summary>
        /// Solves the linear equation system:
        /// 
        /// <paramref name="M"/>*<paramref name="x"/> = <paramref name="b"/>.
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system.</param>
        /// <param name="b">Right-hand-side of the equation system.</param>
        /// <param name="M">General quadratic, non-singular matrix.</param>
        static public void Solve<T,V,W>(this T M, V x, W b) 
            where T : IMatrix
            where V : IList<double>
            where W : IList<double> //
        {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve nonquadratic matrix.");
            if (x.Count != M.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Count != M.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows");
            unsafe {

                int L = M.NoOfCols;

                double[] _x = x as double[];
                if(_x == null) {
                    _x = new double[x.Count];
                }
                _x.SetV(b);

                int* ipiv = stackalloc int[L];
                int i0;
                double[] _this_Entries = TempBuffer.GetTempBuffer(out i0, L * L);
                fixed (double* this_Entries = _this_Entries, p_x = _x) {

                    CopyToUnsafeBuffer(M, this_Entries, true);

                    int info;
                    LAPACK.F77_LAPACK.DGETRF(ref L, ref L, this_Entries, ref L, ipiv, out info);
                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        string infostring;
                        if (info < 0) {
                            infostring = String.Format("the {0}-th argument had an illegal value", info);
                        }
                        else {
                            infostring = "U(" + info + @""",""" + info + ") is exactly zero. The factorization \n has been completed, but the factor U is exactly \n singular, and division by zero will occur if it is used \n to solve a system of equations.";
                        }

                        throw new ArithmeticException("LAPACK dgetrf info: " + infostring);
                    }
                    //         TRANS, N, NRHS, A,            LDA, IPIV, B, LDB
                    char transp = 'N';
                    int eins = 1;
                    LAPACK.F77_LAPACK.DGETRS(ref transp, ref L, ref eins, this_Entries, ref L, ipiv, p_x, ref L, out info);
                    if(info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        throw new ArithmeticException("LAPACK dgetrs info: " + info);
                    }
                }
                TempBuffer.FreeTempBuffer(i0);

                if(!object.ReferenceEquals(x, _x))
                    x.SetV(_x);
            }
        }
        

        /// <summary>
        /// Solves the linear equation system with multiple right-hand-sides:
        /// 
        /// <paramref name="M"/>*<paramref name="X"/> = <paramref name="B"/>.
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system; each column is the solution for one right-hand-side</param>
        /// <param name="B">Matrix of right-hand-sides; each column is a independent right-hand-side.</param>
        /// <param name="M">General quadratic, non-singular matrix.</param>
        static public void SolveEx<T, V, W>(this T M, V x, W b)
            where T : IMatrix
            where V : IMatrix
            where W : IMatrix //
        {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve non-quadratic matrix.");
            if (x.NoOfCols != b.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (x.NoOfRows != M.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (b.NoOfRows != M.NoOfRows)
                throw new ArgumentException("number of rows in b must be equal to number of rows in M");
            unsafe {

                int L = M.NoOfCols;

                int NoRhs = b.NoOfCols;

                int* ipiv = stackalloc int[L];
                double[] _this_Entries = TempBuffer.GetTempBuffer(out int i0, L * L);
                double[] _xRhs_Entries = TempBuffer.GetTempBuffer(out int i1, L * NoRhs);
                fixed (double* this_Entries = _this_Entries, xrhs_Entries = _xRhs_Entries) {

                    CopyToUnsafeBuffer(M, this_Entries, true);
                    CopyToUnsafeBuffer(b, xrhs_Entries, true);

                    int info;
                    LAPACK.F77_LAPACK.DGETRF(ref L, ref L, this_Entries, ref L, ipiv, out info);
                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        TempBuffer.FreeTempBuffer(i1);
                        string infostring;
                        if (info < 0) {
                            infostring = String.Format("the {0}-th argument had an illegal value", info);
                        } else {
                            infostring = "U(" + info + @""",""" + info + ") is exactly zero. The factorization \n has been completed, but the factor U is exactly \n singular, and division by zero will occur if it is used \n to solve a system of equations.";
                        }

                        throw new ArithmeticException("LAPACK dgetrf info: " + infostring);
                    }
                    //         TRANS, N, NRHS, A,            LDA, IPIV, B, LDB
                    char transp = 'N';
                    
                    LAPACK.F77_LAPACK.DGETRS(ref transp, ref L, ref NoRhs, this_Entries, ref L, ipiv, xrhs_Entries, ref L, out info);
                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        TempBuffer.FreeTempBuffer(i1);
                        throw new ArithmeticException("LAPACK dgetrs info: " + info);
                    }

                    CopyFromUnsafeBuffer(x, xrhs_Entries, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
            }
        }
        

        /// <summary>
        /// Solves the linear equation system:
        /// 
        /// <paramref name="M"/>*<paramref name="x"/> = <paramref name="b"/>
        /// under the condition
        /// <paramref name="N"/>*<paramref name="x"/> = <paramref name="c"/>
        /// if no such exact solution exists the first set is solved exactly,
        /// in such a fashion to minimize the residual of the second set
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system.</param>
        /// <param name="b">Right-hand-side of the main equation system.</param>
        /// <param name="M">General matrix.</param>
        /// <param name="c">Right-hand-side of the side equation system.</param>
        /// <param name="N">General matrix.</param>
        /// The second system is used to solve
        static public void SolveWithCondition<T>(this T M, double[] x, double[] b, T N, double[] c) 
            where T : IMatrix //
        {
            if (M.NoOfCols != N.NoOfCols)
                throw new ApplicationException("Solutionspace of both systems has to be of equal dimension");
            if (x.Length != M.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Length != M.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows of M");
            if (c.Length != N.NoOfRows)
                throw new ArgumentException("length of c must be equal to number of rows of N");

            // compute the nullspace of M
            //MultidimensionalArray SRREF = M.GetSolutionSpace();
            MultidimensionalArray S = M.GetSolutionSpaceSVD();
            if (S == null) {
                throw new ApplicationException("Something went wrong");
            }

            // retrieve the particular solution for the first system
            // I believe this could be done while retrieving the nullspace
            M.LeastSquareSolve(x, b);

            // Solve the augmented second system N * S * xi = (c - N * v)
            // In such a manner to minimize ||N * S * xi - (c - N * v)||
            // x = S * xi + v

            // non-shallow copy of c
            double[] rhs = new double[c.Length];
            c.CopyTo(rhs, 0);
            double[] xi = new double[S.NoOfCols];

            // create the new rhs
            N.MatVecMul(-1.0, x, 1.0, rhs);
            GEMM(N, S).LeastSquareSolve(xi, rhs);

            // compute the final solution
            S.MatVecMul(1.0, xi, 1.0, x);
        }

        /// <summary>
        /// computes a LU-Factorization of <paramref name="M"/> and stores it in-place, using LAPACK function DGETRF
        /// </summary>
        /// <param name="M">(input, output) General quadratic, non-singular matrix; on exit, the LU-factorization</param>
        /// <param name="_ipiv">(output, allocated by caller) pivot indices, as computed by LAPACK</param>
        static public void FactorizeLU<T>(this T M, int[] _ipiv) where T : IMatrix {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve nonquadratic matrix.");
            if (_ipiv.Length != M.NoOfCols)
                throw new ArgumentException("length of ipiv must be equal to number of columns");
            unsafe {

                int L = M.NoOfCols;
                double[] _this_Entries = TempBuffer.GetTempBuffer(out int i0, L * L);


                fixed (int* ipiv = _ipiv) {
                    fixed (double* this_Entries = _this_Entries) {

                        CopyToUnsafeBuffer(M, this_Entries, true);

                        int info;
                        LAPACK.F77_LAPACK.DGETRF(ref L, ref L, this_Entries, ref L, ipiv, out info);


                        if (info != 0) {
                            TempBuffer.FreeTempBuffer(i0);
                            string infostring;
                            if (info < 0) {
                                infostring = String.Format("the {0}-th argument had an illegal value", info);
                            } else {
                                infostring = "U(" + info + @""",""" + info +
                                    ") is exactly zero. The factorization \n has been completed, but the factor U is exactly \n singular, and division by zero will occur if it is used \n to solve a system of equations.";
                            }

                            throw new ArithmeticException("LAPACK dgetrf info: " + infostring);
                        }


                        CopyFromUnsafeBuffer(M, this_Entries, true);
                    }
                }

                TempBuffer.FreeTempBuffer(i0);
            }
        }


        /// <summary>
        /// Performs the backward substitution which has been obtained through <see cref="FactorizeLU{T}(T, int[])"/>
        /// </summary>
        static public void BacksubsLU<T>(this T M, int[] _ipiv, double[] x, double[] b) 
            where T : IMatrix //
        {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve non-quadratic matrix.");
            if (_ipiv.Length != M.NoOfCols)
                throw new ArgumentException("length of ipiv must be equal to number of columns");
            if (x.Length != M.NoOfCols)
                throw new ArgumentException("length of x must be equal to number of columns");
            if (b.Length != M.NoOfRows)
                throw new ArgumentException("length of b must be equal to number of rows");

            int L = M.NoOfCols;
            unsafe {
                int iBuf;

                //double[] _this_Entries;
                //int BufOffset;
                //if (M is MultidimensionalArray Mda && Mda.IsContinious) {
                //    _this_Entries = Mda.Storage; // fk, 06apr22: does not work, because MultidimensionalArray is in C-order
                //    iBuf = -1;
                //    BufOffset = Mda.Index(0, 0);
                //} else {
                //    _this_Entries = TempBuffer.GetTempBuffer(out iBuf, L * L);
                //    BufOffset = 0;
                //}

                double[] _this_Entries = TempBuffer.GetTempBuffer(out iBuf, L * L);
                int BufOffset = 0;

               
                fixed (int* ipiv = _ipiv) {

                    double* xx = stackalloc double[L + 2];
                    for (int i = 0; i < L; i++) {
                        xx[i + 1] = b[i];
                    }
                    xx[0] = 123.456;
                    xx[L + 1] = -987.76;

                    fixed (double* __this_Entries = _this_Entries) {

                        double* this_Entries = __this_Entries + BufOffset;
                        if(iBuf >= 0)
                            CopyToUnsafeBuffer(M, this_Entries, true);

                        char transp = 'N';
                        int eins = 1;
                        LAPACK.F77_LAPACK.DGETRS(ref transp, ref L, ref eins, this_Entries, ref L, ipiv, xx + 1, ref L, out int info);
                        if (info != 0) {
                            if(iBuf >= 0)
                                TempBuffer.FreeTempBuffer(iBuf);
                            throw new ArithmeticException("LAPACK dgetrs info: " + info);
                        }

                    }


                    if (xx[0] != 123.456 || xx[L + 1] != -987.76)
                        throw new MemberAccessException("LAPACK.DGETRS accessed memory that it should not touch");
                    for (int i = 0; i < L; i++) {
                        x[i] = xx[i + 1];
                    }

                }


                if (iBuf >= 0)
                    TempBuffer.FreeTempBuffer(iBuf);
            }
        }


        /// <summary>
        /// Performs the backward substitution which has been obtained through <see cref="FactorizeLU{T}(T, int[])"/>
        /// for multiple right-hand-sides
        /// </summary>
        static public void BacksubsLU<T, V, W>(this T M, int[] _ipiv, V x, W b)
            where T : IMatrix
            where V : IMatrix
            where W : IMatrix //
        {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve non-quadratic matrix.");
            if (x.NoOfCols != b.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (x.NoOfRows != M.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (b.NoOfRows != M.NoOfRows)
                throw new ArgumentException("number of rows in b must be equal to number of rows in M");
            if (_ipiv.Length != M.NoOfCols)
                throw new ArgumentException("length of ipiv must be equal to number of columns");
            unsafe {

                int L = M.NoOfCols;

                int NoRhs = b.NoOfCols;


                double[] _this_Entries = TempBuffer.GetTempBuffer(out int i0, L * L);
                double[] _xRhs_Entries = TempBuffer.GetTempBuffer(out int i1, L * NoRhs);
                fixed (double* this_Entries = _this_Entries, xrhs_Entries = _xRhs_Entries) {
                    fixed (int* ipiv = _ipiv) {
                        CopyToUnsafeBuffer(M, this_Entries, true);
                        CopyToUnsafeBuffer(b, xrhs_Entries, true);

                        int info;
                        /*
                        LAPACK.F77_LAPACK.DGETRF(ref L, ref L, this_Entries, ref L, ipiv, out info);
                        if (info != 0) {
                            TempBuffer.FreeTempBuffer(i0);
                            TempBuffer.FreeTempBuffer(i1);
                            string infostring;
                            if (info < 0) {
                                infostring = String.Format("the {0}-th argument had an illegal value", info);
                            } else {
                                infostring = "U(" + info + @""",""" + info + ") is exactly zero. The factorization \n has been completed, but the factor U is exactly \n singular, and division by zero will occur if it is used \n to solve a system of equations.";
                            }

                            throw new ArithmeticException("LAPACK dgetrf info: " + infostring);
                        }*/
                        //         TRANS, N, NRHS, A,            LDA, IPIV, B, LDB
                        char transp = 'N';
                        LAPACK.F77_LAPACK.DGETRS(ref transp, ref L, ref NoRhs, this_Entries, ref L, ipiv, xrhs_Entries, ref L, out info);
                        if (info != 0) {
                            TempBuffer.FreeTempBuffer(i0);
                            TempBuffer.FreeTempBuffer(i1);
                            throw new ArithmeticException("LAPACK dgetrs info: " + info);
                        }

                        CopyFromUnsafeBuffer(x, xrhs_Entries, true);
                    }
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
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
                throw new ArgumentException("Cannot solve nonquadratic matrix.");
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
        /// Solves the symmetric, positive definite linear equation system with multiple right-hand-sides:
        /// 
        /// <paramref name="M"/>*<paramref name="X"/> = <paramref name="B"/>.
        /// </summary>
        /// <param name="x">On exit, the solution of the equation system; each column is the solution for one right-hand-side</param>
        /// <param name="B">Matrix of right-hand-sides; each column is a independent right-hand-side.</param>
        /// <param name="M">General quadratic, non-singular matrix.</param>
        static public void SolveSymmetricEx<T,TX,TB>(this T M, TX x, TB b)
            where T : IMatrix 
            where TX : IMatrix 
            where TB : IMatrix //
        {
            if (M.NoOfRows != M.NoOfCols)
                throw new ApplicationException("Cannot solve nonquadratic matrix.");
            if (x.NoOfCols != b.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (x.NoOfRows != M.NoOfCols)
                throw new ArgumentException("mismatch between number of columns in x and b");
            if (b.NoOfRows != M.NoOfRows)
                throw new ArgumentException("number of rows in b must be equal to number of rows in M");
            unsafe {

                int L = M.NoOfCols;
                int NoRhs = b.NoOfCols;
               
                int* ipiv = stackalloc int[L];
                double[] _this_Entries = TempBuffer.GetTempBuffer(out int i0, L * L);
                double[] _xrhs_Entries = TempBuffer.GetTempBuffer(out int i1, L * L);
                fixed (double* this_Entries = _this_Entries, xrhs_Entries = _xrhs_Entries) {
                    CopyToUnsafeBuffer(M, this_Entries, true);
                    CopyToUnsafeBuffer(b, xrhs_Entries, true);

                    int uplo = 'U', info;
                    LAPACK.F77_LAPACK.DPOSV_(ref uplo, ref L, ref NoRhs, this_Entries, ref L, xrhs_Entries, ref L, out info);
                    if (info != 0) {
                        TempBuffer.FreeTempBuffer(i0);
                        TempBuffer.FreeTempBuffer(i1);
                        throw new ArithmeticException("LAPACK dposv info: " + info);
                    }

                    CopyFromUnsafeBuffer(x, xrhs_Entries, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
            }
        }

        /// <summary>
        /// Condition number of a quadratic matrix 
        /// </summary>
        /// <param name="M"></param>
        /// <param name="NORM">
        /// Switch between one and infinity norm
        /// - char '1': the one norm
        /// - char 'I': the infinity norm
        /// </param>
        static public double Cond<T>(this T M, int NORM = '1') where T : IMatrix {
            if (M.NoOfRows != M.NoOfCols)
                throw new ArgumentException("Only supported for quadratic matrix.");
            int N = M.NoOfCols;

            double RCOND = double.NaN;
            int INFO = 0;

            if(NORM != '1' && NORM != 'I')
                throw new ArgumentException("Unknown norm specifier; expecting either '1' for one-norm or 'I' for infinity norm.");

            unsafe {
                double[] _this_Entries = TempBuffer.GetTempBuffer(out int i0, N*N);
                double[] _work = TempBuffer.GetTempBuffer(out int i1, 4*N);
                int* Iwork = stackalloc int[N];
                int* IPIV = stackalloc int[N];

                fixed(double* A = _this_Entries, work = _work) {
                    CopyToUnsafeBuffer(M, A, true);

                    int LDA = N;
                    double ANORM;

                    if(NORM == '1' || NORM == 'O')
                        ANORM = M.OneNorm();
                    else if(NORM == 'I')
                        ANORM = M.InfNorm();
                    else
                        throw new ArgumentException("Illegal Matrix Norm Specifier.");

                    LAPACK.F77_LAPACK.DGETRF(ref N, ref N, A, ref LDA, IPIV, out INFO);
                    if(INFO > 0)
                        return double.PositiveInfinity;

                    if(INFO != 0)
                        throw new ArithmeticException("LAPACK DGETRF info is " + INFO);
                    
                    LAPACK.F77_LAPACK.DGECON_(ref NORM, ref N, A, ref LDA, ref ANORM, ref RCOND, work, Iwork, ref INFO);
                    if(INFO != 0)
                        throw new ArithmeticException("LAPACK DGECON info is " + INFO);
                }

                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
            }

            

            return 1.0 / RCOND;
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
                fixed (double* pM = _M, pTAU = TAU) {
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
                fixed (double* pM = _M, pTAU = TAU) {
                    fixed (int* pP = P) {
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

            //Console.WriteLine("The rank of the matrix coming from the Reduced Row Echelon Form: rank={0}", rank);

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
        /// Alternative for <see cref="GetSolutionSpace{T}(T)"/> using Lapack routine dgesvd
        /// Converts an implicit subspace representation (given as the solution of a singular matrix <paramref name="Mtx"/>)
        /// into an explicit representation.
        /// </summary>
        /// <param name="Mtx">
        /// A matrix implicitly defines a subspace.
        /// </param>
        /// <returns>
        /// A matrix whose columns span the nullspace of the input matrix <paramref name="Mtx"/>.
        /// </returns>
        public static MultidimensionalArray GetSolutionSpaceSVD<T>(this T Mtx) where T : IMatrix
        {
            int M = Mtx.NoOfRows;
            int N = Mtx.NoOfCols;

            int JOBU = 'N';
            int JOBVT = 'A';

            int LDA = M >= 1 ? M : 1;
            int LDU = 1;
            switch (JOBU)
            {
                case 'A':
                case 'S':
                    LDU = M;
                    break;
                default:
                    LDU = 1;
                    break;
            }

            int LDVT = 1;
            switch (JOBVT)
            {
                case 'A':
                case 'S':
                    LDVT = N;
                    break;
                default:
                    LDVT = 1;
                    break;
            }
            int LDS = Math.Min(M, N);

            MultidimensionalArray S = MultidimensionalArray.Create(LDS,1);
            MultidimensionalArray VT = MultidimensionalArray.Create(LDVT, LDVT);

            unsafe
            {

                int i0, i1, i2, i3;
                double[] __A = TempBuffer.GetTempBuffer(out i0, M * N);
                double[] __S = TempBuffer.GetTempBuffer(out i1, LDS);
                double[] __U = TempBuffer.GetTempBuffer(out i2, LDU * LDU);
                double[] __VT = TempBuffer.GetTempBuffer(out i3, LDVT * LDVT);
                fixed (double* _A = __A, _S = __S, _U = __U, _VT = __VT)
                {
                    CopyToUnsafeBuffer(Mtx, _A, true);

                    LAPACK.F77_LAPACK.DGESVD(JOBU, JOBVT, M, N, _A, LDA, _S, _U, LDU, _VT, LDVT);

                    CopyFromUnsafeBuffer(S, _S, true);
                    CopyFromUnsafeBuffer(VT, _VT, true);
                }
                TempBuffer.FreeTempBuffer(i0);
                TempBuffer.FreeTempBuffer(i1);
                TempBuffer.FreeTempBuffer(i2);
                TempBuffer.FreeTempBuffer(i3);
            }

            // we need the original V not V^T
            VT.TransposeInPlace();

            int nullspaceDim = M < N ? N - M : 0;

            double thresh = 1.0E-14;
            for (int i = 0; i < LDS; i++)
            {
                if (S.Storage[i] < thresh)
                {
                    nullspaceDim++;
                }
            }

            MultidimensionalArray Nullspace = null;
            if (nullspaceDim > 0)
            {
                Nullspace = MultidimensionalArray.Create(N, nullspaceDim);
                Nullspace.Acc(1.0, VT.ExtractSubArrayShallow(new int[] { 0, LDVT - nullspaceDim }, new int[] { LDVT - 1, LDVT - 1 }));
            }
            else
            {
                Console.WriteLine("Warning: Nullspace is empty");
            }

            return Nullspace;
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
        /// <param name="ReadOffset">
        /// offset into <paramref name="col"/>, coping starts from this index
        /// </param>
        /// <param name="ReadInc">
        /// skip between to consecutive elements that are taken from <paramref name="col"/>
        /// </param>
        public static void SetColumn<T>(this IMatrix inp, int ColNo, T col, int ReadOffset = 0, int ReadInc = 1) where T : IEnumerable<double> {

            if(ColNo < 0 || ColNo >= inp.NoOfCols)
                throw new IndexOutOfRangeException("ColNo out of range");
            //if((col.Count() - ReadOffset) < inp.NoOfRows * ReadInc)
            //    throw new ArgumentException("array to short", "row");
            // there will be an index exception anyway...

            int i = 0;
            if(ReadInc == 1 && ReadOffset == 0) {
                foreach(double col_i in col) {
                    inp[i, ColNo] = col_i;
                    i++;
                }
            } else {
                int j = 0;
                foreach(double col_i in col) {
                    if(i >= ReadOffset && ((j - ReadOffset) % ReadInc == 0)) {
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
