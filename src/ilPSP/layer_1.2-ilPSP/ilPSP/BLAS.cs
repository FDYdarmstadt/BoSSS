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
using MPI.Wrappers.Utils;
using MPI.Wrappers;
using ilPSP;
using ilPSP.Utils;
using System.IO;
using System.Globalization;
using System.Diagnostics;
using System.Linq;

namespace ilPSP.Utils {

    /// <summary>
    /// MPI-parallel variants of some BLAS functions
    /// </summary>
    public static class ParallelBlas {

        ///// <summary>
        ///// for the summation of one int, double, ..., over all MPI processors
        ///// </summary>
        //static public T MPI_Sum<T>(T x, MPI_Comm comm) 
        //    where T : struct {

        //    MPI_Datatype type;
        //    Type t = typeof(T);
        //    if(t==typeof(int))
        //        type = csMPI.Raw.MPI_Datatype.MPI_INT;
        //    else if(t==typeof(int))
        //        type = csMPI.Raw.MPI_Datatype.MPI_DOUBLE;
        //    else if(t==typeof(int))
        //        type = csMPI.Raw.MPI_Datatype.MPI_FLOAT;
        //    else if(t==typeof(int))
        //        type = csMPI.Raw.MPI_Datatype.MPI_LONG_LONG;
        //    else 
        //        throw new NotSupportedException("this function is not supported for type " + t.Name);
            
        //    unsafe {
        //        T outp = default(T);
        //        csMPI.Raw.Allreduce((IntPtr)(&x), (IntPtr)(&outp), 1, csMPI.Raw.MPI_Datatype.MPI_DOUBLE, csMPI.Raw.MPI_OP.SUM, comm);
        //        return outp;
        //    }
            
        //}

        /// <summary>
        /// MPI-parallel scalar product
        /// </summary>
        static public double MPI_ddot<V,W>(this V vec1, W vec2, MPI_Comm comm) 
            where V : IList<double> 
            where W: IList<double> 
        {
            if (vec1.Count != vec2.Count)
                throw new ArgumentException("vectors mus have the same lenght");
            return ddot(vec1.Count, vec1, 1, vec2, 1, comm);
        }

        /// <summary>
        /// MPI-parallel scalar product (on world communicator)
        /// </summary>
        static public double MPI_ddot<V,W>(this V vec, W vec2)
            where V : IList<double>
            where W : IList<double> {
            return vec.MPI_ddot(vec2, csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// MPI - parallel scalar product
        /// </summary>
        static public double ddot<TX, TY>(int N, TX x, int incx, TY y, int incy, MPI.Wrappers.MPI_Comm comm)
            where TX : IList<double>
            where TY : IList<double> {

            if( incx*x.Count < N)
                throw new ArgumentException("vector too short.","x");
            if( incy*y.Count < N)
                throw new ArgumentException("vector too short.","y");
            
            double locRes = 0;

            double[] dx = x as double[];
            double[] dy = y as double[];
            if (dx != null && dy != null) {
                // use dnese BLAS
                locRes = BLAS.ddot(N, dx, incx, dy, incy);
            } else {

                ISparseVector<double> spx = x as ISparseVector<double>;
                ISparseVector<double> spy = y as ISparseVector<double>;
                IList<double> _y = y;
                if (spy != null) {
                    if (spx == null || spy.Sparsity < spx.Sparsity) {
                        // y is sparser than x, use y !
                        spx = spy;
                        spy = null;
                        _y = x;
                        x = default(TX);
                        y = default(TY);
                        int buffy = incx;
                        incx = incy;
                        incy = buffy;
                    }
                }

                if (spx != null) {
                    // sparse implementation

                    foreach (var entry in spx.SparseStruct) {
                        int m = entry.Key % incx;
                        if (m != 0)
                            // entry is skipped by x-increment
                            continue;

                        int n = entry.Key / incx;

                        locRes += entry.Value * _y[n * incy];
                    }
                } else {
                    // default IList<double> - implementation
                    for (int n = 0; n < N; n++) {
                        locRes += x[n * incx] * y[n * incy];
                    }

                }
            }

            double globRes = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&locRes), (IntPtr)(&globRes), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, comm);
            }
            return globRes;
        }

        /// <summary>
        /// MPI-parallel 2-norm
        /// </summary>
        static public double MPI_L2Norm<V>(this V vec, MPI_Comm comm) where V : IList<double> {
            return drnm2(vec.Count, vec, 1, comm);
        }

        /// <summary>
        /// MPI-parallel 2-norm (on world communicator)
        /// </summary>
        static public double MPI_L2Norm<V>(this V vec) where V : IList<double> {
            return drnm2(vec.Count, vec, 1, csMPI.Raw._COMM.WORLD);
        }


        /// <summary>
        /// MPI-parallel 2-norm
        /// </summary>
        static public double drnm2<TX>(int N, TX x, int incx, MPI_Comm comm) where TX : IList<double> {
            double locRes = 0;

            double[] dx = x as double[];
            if (dx != null) {
                // double[] - implementation
                locRes = BLAS.dnrm2(N, dx, incx);
                locRes = locRes * locRes;
            } else {
                ISparseVector<double> spx = x as ISparseVector<double>;

                if (spx != null) {
                    // sparse implementation

                    foreach (var entry in spx.SparseStruct) {
                        int m = entry.Key % incx;
                        if (m != 0)
                            // entry is skipped by x-increment
                            continue;

                        double xi = entry.Value;
                        locRes += xi * xi;
                    }
                } else {
                    // default implementation
                    for (int n = 0; n < N; n++) {
                        double xi = x[n * incx];
                        locRes += xi * xi;
                    }
                }
            }

            

            double globRes = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&locRes), (IntPtr)(&globRes), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, comm);
            }
            return Math.Sqrt(globRes);
        }
    }

    /// <summary>
    /// some utility BLAS-style functions
    /// </summary>
    public static class GenericBlas {

        
        /// <summary>
        /// Creates linear interpolated nodes between two points.
        /// </summary>
        /// <param name="a">minimum</param>
        /// <param name="b">maximum</param>
        /// <param name="n">number of nodes desired</param>
        /// <returns>an array of length <paramref name="n"/>,
        /// with first entry equal to <paramref name="a"/>, 
        /// last entry equal to <paramref name="b"/>, and 
        /// all other points linear interpolated in between.</returns>
        public static double[] Linspace(double a, double b, int n) {
            if (a >= b)
                throw new ArgumentException("minimum >= maximum");
            if (n <= 1)
                throw new ArgumentException("Number of nodes must be at least 2.");

            double[] r = new double[n];
            double dx = (b - a) / ((double)(n - 1));
            for (int i = 0; i < n; i++) {
                r[i] = a + dx * ((double)i);
            }
            return r;
        }

        /// <summary>
        /// Nodes between two points, compressed at both ends
        /// </summary>
        /// <param name="l">lower limit.</param>
        /// <param name="r">upper limit.</param>
        /// <param name="a">scaling between linear ans sinus-mapping</param>
        /// <param name="n">number of nodes</param>
        /// <returns></returns>
        public static double[] SinLinSpacing(double l, double r, double a, int n) {
            double[] linnodes = GenericBlas.Linspace(-Math.PI * 0.5, Math.PI * 0.5, n);
            double[] linnodes2 = GenericBlas.Linspace(-1, 1, n);
            double[] nodes = new double[n];

            for (int i = 0; i < n; i++)
                //nodes[i] = linnodes2[i] * (1 - a) + (1.0 - Math.Sin(linnodes[i])) * a;
                nodes[i] = linnodes2[i] * (1 - a) + Math.Sin(linnodes[i]) * a;

            for (int i = 0; i < n; i++)
                nodes[i] = nodes[i] * (r - l) * 0.5 + l;
            return nodes;
        }

        /// <summary>
        /// generic dswap
        /// </summary>
        static public void dswap<U, V>(int N, U DX, int INCX, V DY, int INCY) 
            where U : IList<double>
            where V : IList<double>
        {
            if(object.ReferenceEquals(DX, DY))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug.
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
           

            if((DX is double[]) && (DY is double[])) {
                BLAS.dswap(N, DX as double[], INCX, DY as double[], INCY);
                return;
            } else {

                for(int i = 0; i < N; i++) {
                    double a;
                    a = DX[i * INCX];
                    DX[i * INCX] = DY[i * INCY];
                    DY[i * INCY] = a;
                }
            }
        }


        /// <summary>
        /// blas dscal: <paramref name="X"/> = <paramref name="X"/>*<paramref name="alpha"/>;
        /// </summary>
        static public void dscal<T>(int N, double alpha, T X, int incx) where T : IList<double> {

            ISparseVector<double> spx = X as ISparseVector<double>;
            double[] arx = X as double[];


            if (arx != null) {
                // double - array -> use BLAS
                // ++++++++++++++++++++++++++
                BLAS.dscal(N, alpha, arx, incx);
                return;
            }

            if (spx != null) {
                // sparce vector -> compute only nonzero entries
                // +++++++++++++++++++++++++++++++++++++++++++++

                int[] idx = new int[spx.NonZeros];
                spx.SparseStruct.Keys.CopyTo(idx, 0);

                foreach (int i in idx) {
                    int m = i % incx;
                    if (m != 0)
                        // entry is skipped by x-increment
                        continue;

                    spx[i] *= alpha;
                }

                return;

            }


            {
                // reference version
                // +++++++++++++++++

                for (int i = 0; i < N; i += incx)
                    X[i] *= alpha;
            }
        }


        /// <summary>
        /// blas daxpy: <paramref name="Y"/> = <paramref name="Y"/> + <paramref name="alpha"/>*<paramref name="X"/>;
        /// </summary>
        static public void daxpy<TX, TY>(int N,
                         double alpha, TX X, int INCX,
                         TY Y, int INCY)
            where TX : IList<double>
            where TY : IList<double> 
        {
            if(object.ReferenceEquals(X, Y))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug.
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
           
            {
                // sparse vector branch
                // ++++++++++++++++++++
                ISparseVector<double> spx = X as ISparseVector<double>;
                if (spx != null) {


                    foreach (var entry in spx.SparseStruct) {
                        int m = entry.Key % INCX;
                        if (m != 0)
                            // entry is skipped by x-increment
                            continue;

                        int n = entry.Key / INCX;

                        double xi = entry.Value;
                        Y[n * INCY] += xi * alpha;
                    }

                    return;
                }
            }


            {
                // both arrays-branch -> use BLAS
                // ++++++++++++++++++++++++++++++


                double[] _XasDouble = X as double[];
                double[] _YasDouble = Y as double[];
                if (_XasDouble != null && _YasDouble != null) {
                    BLAS.daxpy(N, alpha, _XasDouble, INCX, _YasDouble, INCY);
                    return;
                }

            }
            {
                // default branch
                // ++++++++++++++

                for (int n = 0; n < N; n++)
                    Y[n * INCY] += X[n * INCX] * alpha;

                return;
            }
        }

        public static object Cartesian2DGrid(double[] v1, double[] v2) {
            throw new NotImplementedException();
        }


        /// <summary>
        /// L2-Norm
        /// </summary>
        static public double L2Norm<T>(this T a) where T : IList<double> {
            return Math.Sqrt(L2NormPow2(a));
        }


        /// <summary>
        /// scales some vector <paramref name="a"/> by scalar <paramref name="alpha"/>.
        /// </summary>
        static public void ScaleV<T>(this T a, double alpha) where T : IList<double> {
            GenericBlas.dscal(a.Count, alpha, a, 1);
        }

        /// <summary>
        /// scales some vector <paramref name="a"/> by scalar <paramref name="alpha"/>.
        /// </summary>
        static public void ScaleV<T,R>(this T a, double alpha, R index)
            where T : IList<double>
            where R : IList<int>
        {
            int L = index.Count;
            for (int i = 0; i < L; i++) {
                a[index[i]] *= alpha;
            }
        }

        /// <summary>
        /// for all indices <em>i</em>, this sets
        /// <paramref name="a"/>[<em>i</em>] = <paramref name="b"/>[<em>i</em>] 
        /// </summary>
        static public void CopyEntries<T, R>(this T a, R b)
            where T : IList<double>
            where R : IList<double>
        {
            if(object.ReferenceEquals(a, b))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug.
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
           
            a.ClearEntries();
            a.AccV(1.0, b);
        }

        /// <summary>
        /// clear all entries.
        /// </summary>
        static public void ClearEntries<T>(this T a) where T : IList<double> {
            int L = a.Count;
            if(a is Array) {
                // optimized for arrays
                Array.Clear(a as double[], 0, L);
                return;
            } 

            // default:
            for (int i = 0; i < L; i++) {
                a[i] = 0.0;
            }
        }

        /// <summary>
        /// clear all entries.
        /// </summary>
        static public void ClearEntries<T,R>(this T a, R index) 
            where T : IList<double>
            where R : IList<int>
        {
            int L = index.Count;
            for (int i = 0; i < L; i++) {
                a[index[i]] = 0.0;
            }
        }

        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/> + <paramref name="alpha"/>*<paramref name="B"/>
        /// </summary>
        static public void AccV<T, V>(this T a, double alpha, V B)
            where T : IList<double>
            where V : IList<double> //
        {
            int N = a.Count;
            if(B.Count != N)
                throw new ArgumentException("element count must be equal");
            if(object.ReferenceEquals(a, B))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug, i.e. the user called it un-intentionally. 
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
            
            GenericBlas.daxpy(N, alpha, B, 1, a, 1);
        }
        
        /// <summary>
        /// <paramref name="a"/> = <paramref name="alpha"/>*<paramref name="B"/>
        /// </summary>
        static public void SetV<T, V>(this T a, V B, double alpha = 1.0)
            where T : IList<double>
            where V : IList<double> //
        {
            int N = a.Count;
            if(B.Count != N)
                throw new ArgumentException("element count must be equal");
            if(object.ReferenceEquals(a, B))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug, i.e. the user called it un-intentionally. 
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
            a.ClearEntries();
            GenericBlas.daxpy(N, alpha, B, 1, a, 1);
        }

        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/>*<paramref name="beta"/> + <paramref name="alpha1"/>*<paramref name="B1"/>
        /// </summary>
        static public void SumV<T, V1>(this T a, double beta, double alpha1, V1 B1)
            where T : IList<double>
            where V1 : IList<double>
        {
            if (beta != 1.0)
                a.ScaleV(beta);
            a.AccV(alpha1, B1);
        }

        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/>*<paramref name="beta"/> + <paramref name="alpha1"/>*<paramref name="B1"/> + <paramref name="alpha2"/>*<paramref name="B2"/>
        /// </summary>
        static public void SumV<T, V1, V2>(this T a, double beta, double alpha1, V1 B1, double alpha2, V2 B2)
            where T : IList<double>
            where V1 : IList<double>
            where V2 : IList<double> 
        {
            if(beta != 1.0)
                a.ScaleV(beta);
            a.AccV(alpha1, B1);
            a.AccV(alpha2, B2);
        }

        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/>*<paramref name="beta"/> + <paramref name="alpha1"/>*<paramref name="B1"/> + <paramref name="alpha2"/>*<paramref name="B2"/> + <paramref name="alpha3"/>*<paramref name="B3"/>
        /// </summary>
        static public void SumV<T, V1, V2, V3>(this T a, double beta, double alpha1, V1 B1, double alpha2, V2 B2, double alpha3, V3 B3)
            where T : IList<double>
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double> 
        {
            if (beta != 1.0)
                a.ScaleV(beta);
            a.AccV(alpha1, B1);
            a.AccV(alpha2, B2);
            a.AccV(alpha3, B3);
        }

        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/>*<paramref name="beta"/> + <paramref name="alpha1"/>*<paramref name="B1"/> + <paramref name="alpha2"/>*<paramref name="B2"/> + <paramref name="alpha3"/>*<paramref name="B3"/> + <paramref name="alpha4"/>*<paramref name="B4"/>
        /// </summary>
        static public void SumV<T, V1, V2, V3, V4>(this T a, double beta, double alpha1, V1 B1, double alpha2, V2 B2, double alpha3, V3 B3, double alpha4, V4 B4)
            where T : IList<double>
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double>
            where V4 : IList<double>
        {
            if (beta != 1.0)
                a.ScaleV(beta);
            a.AccV(alpha1, B1);
            a.AccV(alpha2, B2);
            a.AccV(alpha3, B3);
            a.AccV(alpha4, B4);
        }


        /// <summary>
        /// <paramref name="a"/> = <paramref name="a"/>*<paramref name="beta"/> + <paramref name="alpha1"/>*<paramref name="B1"/> + <paramref name="alpha2"/>*<paramref name="B2"/> + <paramref name="alpha3"/>*<paramref name="B3"/> + <paramref name="alpha4"/>*<paramref name="B4"/> + <paramref name="alpha5"/>*<paramref name="B5"/>
        /// </summary>
        static public void SumV<T, V1, V2, V3, V4, V5>(this T a, double beta, double alpha1, V1 B1, double alpha2, V2 B2, double alpha3, V3 B3, double alpha4, V4 B4, double alpha5, V5 B5)
            where T : IList<double>
            where V1 : IList<double>
            where V2 : IList<double>
            where V3 : IList<double>
            where V4 : IList<double>
            where V5 : IList<double>
        {
            if (beta != 1.0)
                a.ScaleV(beta);
            a.AccV(alpha1, B1);
            a.AccV(alpha2, B2);
            a.AccV(alpha3, B3);
            a.AccV(alpha4, B4);
            a.AccV(alpha5, B5);
        }

        /// <summary>
        /// checks all entries for infinity or NAN - values, and
        /// throws an <see cref="ArithmeticException"/> if found;
        /// </summary>
        static public void CheckForNanOrInfV<T>(this T v, bool CheckForInf = true, bool CheckForNan = true, bool ExceptionIfFound = true)
            where T: IEnumerable<double>
        {
            int cnt = 0;
            foreach (double a in v) {
                
                if (CheckForNan)
                    if (double.IsNaN(a)) {
                        if (ExceptionIfFound)
                            throw new ArithmeticException("NaN found at " + cnt + "-th entry.");
                    }

                if (CheckForInf)
                    if (double.IsInfinity(a)) {
                        if (ExceptionIfFound)
                            throw new ArithmeticException("Inf found at " + cnt + "-th entry.");
                    }

                cnt++;
            }
        }

        /// <summary>
        /// accumulation of subvectors
        /// this[<paramref name="acc_index"/>[i]] = this[<paramref name="acc_index"/>[i] + <paramref name="acc_index_shift"/>] + <paramref name="alpha"/>*<paramref name="b"/>[<paramref name="b_index"/>[i] + <paramref name="acc_index_shift"/>]
        /// </summary>
        static public void AccV<T, V, R, S>(this T acc, double alpha, V b, R acc_index, S b_index, int acc_index_shift = 0, int b_index_shift = 0)
            where T : IList<double>
            where V : IList<double> 
            where R : IList<int>
            where S : IList<int>
        {
            if(object.ReferenceEquals(acc, b))
                // while formally there is no problem that this should work, it seems much more likely
                // that it's a bug.
                throw new ArgumentException("Illegal use: reference-equality of input vectors -- this might be a mis-use instead of intention.");
           

            if(acc_index != null && b_index != null) {
                if (acc_index.Count != b_index.Count)
                    throw new ArgumentOutOfRangeException("length of 'acc_index' and 'b_index' must match.");

                int N = acc_index.Count;
                for (int i = 0; i < N; i++) {
                    acc[acc_index[i] + acc_index_shift] += alpha*b[b_index[i] + b_index_shift];
                }

            } else if( acc_index != null && b_index == null) {

                int N = acc_index.Count;
                for (int i = 0; i < N; i++) {
                    acc[acc_index[i] + acc_index_shift] += alpha*b[i + b_index_shift];
                }
            } else if (acc_index == null && b_index != null) {

                int N = b_index.Count;
                for (int i = 0; i < N; i++) {
                    acc[i + acc_index_shift] += alpha*b[b_index[i] + b_index_shift];
                }
            } else if (acc_index == null && b_index == null) {
                int N = acc.Count;
                if (acc.Count != b.Count)
                    throw new ArgumentOutOfRangeException("length of 'acc' and 'b' must match.");

                for (int i = 0; i < N; i++) {
                    acc[i + acc_index_shift] += alpha * b[i + b_index_shift];
                }
            } else {
                throw new Exception("should never be reached");
            }


        }

        /// <summary>
        /// L2-distance
        /// </summary>
        static public double L2Dist<TX, TY>(this TX a, TY b)
            where TX : IList<double>
            where TY : IList<double> {
            return Math.Sqrt(L2DistPow2(a, b));
        }

        /// <summary>
        /// L2-Norm to the power of 2
        /// </summary>
        static public double L2NormPow2<T>(this T a) where T : IList<double> {
            double ret = 0;
            for (int i = a.Count - 1; i >= 0; i--) {
                double ai = a[i];
                ret += ai * ai;
            }
            return ret;
        }

        /// <summary>
        /// L2-distance to the power of 2
        /// </summary>
        static public double L2DistPow2<TX, TY>(this TX a, TY b)
            where TX : IList<double>
            where TY : IList<double> {
            Debug.Assert(a.Count == b.Count, "mismatch in vector length");
            double ret = 0;
            int I = a.Count;
            for (int i = 0; i < I; i++) {
                double ai = a[i] - b[i];
                ret += ai * ai;
            }
            return ret;
        }

        /// <summary>
        /// l2 inner product
        /// </summary>
        static public double InnerProd<TX, TY>(this TX a, TY b)
            where TX : IList<double>
            where TY : IList<double> {
            Debug.Assert(a.Count == b.Count, "mismatch in vector length");
            double ret = 0;
            for (int i = a.Count - 1; i >= 0; i--) {
                ret += a[i] * b[i];
            }
            return ret;
        }

        /// <summary>
        /// Computes the angle, in radians, between the vectors 
        /// (<paramref name="a"/>-<paramref name="c"/>) and (<paramref name="b"/>-<paramref name="c"/>).
        /// </summary>
        static public double Angle<TX, TC, TY>(TX a, TC c, TY b)
            where TX : IList<double>
            where TC : IList<double>
            where TY : IList<double>
        {
            Debug.Assert(a.Count == b.Count, "mismatch in vector length");
            Debug.Assert(a.Count == c.Count, "mismatch in vector length");
            
            double cos_alpha = 0, l1 = 0, l2 = 0;
            for (int i = a.Count - 1; i >= 0; i--) {
                double v1i = a[i] - c[i];
                double v2i = b[i] - c[i];

                l1 += v1i*v1i;
                l2 += v2i*v2i;

                cos_alpha += v1i*v2i;
            }
            cos_alpha *= 1.0/Math.Sqrt(l1*l2);

            return Math.Acos(cos_alpha);
        }

        /// <summary>
        /// Computes the angle, in radians, between the vectors 
        /// <paramref name="a"/> and <paramref name="b"/>.
        /// </summary>
        static public double Angle<TX, TY>(TX a, TY b)
            where TX : IList<double>
            where TY : IList<double> //
        {
            Debug.Assert(a.Count == b.Count, "mismatch in vector length");
            
            double cos_alpha = 0, l1 = 0, l2 = 0;
            for (int i = a.Count - 1; i >= 0; i--) {
                double v1i = a[i];
                double v2i = b[i];

                l1 += v1i*v1i;
                l2 += v2i*v2i;

                cos_alpha += v1i*v2i;
            }
            cos_alpha *= 1.0/Math.Sqrt(l1*l2);

            return Math.Acos(cos_alpha);
        }
    }
    
    /// <summary>
    /// subset of BLAS
    /// </summary>
    public sealed class UnsafeDBLAS : DynLibLoader {
		
		// workaround for .NET bug:
		// https://connect.microsoft.com/VisualStudio/feedback/details/635365/runtimehelpers-initializearray-fails-on-64b-framework
		static PlatformID[] Helper() {
			PlatformID[] p = new PlatformID[5];
			p[0] = PlatformID.Win32NT;
			p[1] = PlatformID.Unix;
			p[2] = PlatformID.Unix;
			p[3] = PlatformID.Unix;
			p[4] = PlatformID.Unix;
			return p;
		}

        /// <summary>
        /// ctor
        /// </summary>
        public UnsafeDBLAS() :
            base(new string[] { "libacml_dll.dll", "libacml.so", "libatlas.so", "libblas.so", "libopenblas.so" },
                  new string[5][][], 
                  new GetNameMangling[] { DynLibLoader.CAPITAL_LETTERS, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore },
                  Helper(), //new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix, PlatformID.Unix },
                  new int[] { -1, -1, -1, -1, -1 }) { }

        
        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate double _DDOT(ref int N, double* DX, ref int INCX, double* DY, ref int INCY);
        

#pragma warning disable        649
        _DDOT ddot;
        _DNRM2 dnrm2;
        _DSWAP dswap;
        _DGEMM dgemm;        
        _DAXPY daxpy;
        _DSCAL dscal;
#pragma warning restore       649

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DDOT DDOT {
            get { return ddot; }
        }

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate double _DSWAP(ref int N, double* DX, ref int INCX, double* DY, ref int INCY);


        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DSWAP DSWAP {
            get { return dswap; }
        }

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate void _DGEMM(ref int TRANSA, ref int TRANSB,
                                           ref int M, ref int N, ref int K,
                                           ref double ALPHA,
                                           double* A, ref int LDA,
                                           double* B, ref int LDB,
                                           ref double BETA,
                                           double* C, ref int LDC);


        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DGEMM DGEMM {
            get { return dgemm; }
        }

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate void _DAXPY(ref int N,
                                           ref double DA, double* DX, ref int INCX,
                                           double* DY, ref int INCY);



        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DAXPY DAXPY {
            get { return daxpy; }
        }

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate void _DSCAL(ref int n, ref double a, double* x, ref int incx);

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DSCAL DSCAL {
            get { return dscal; }
        }

        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe delegate double _DNRM2(ref int n, double* x, ref int incx);


        /// <summary> FORTRAN BLAS routine </summary>
        public unsafe _DNRM2 DNRM2 {
            get { return dnrm2; }
        }
    }



    /// <summary>
    /// some parts of the BLAS interface, which are used by BoSSS;
    /// </summary>
    static public class BLAS {

        /// <summary>
        /// the machine double accuracy;
        /// </summary>
        public static double MachineEps {
            get {
                double machEps = 1.0d;

                do {
                    machEps /= 2.0d;
                }
                while((double)(1.0 + machEps) != 1.0);

                return 2 * machEps;
            }
        }

        static UnsafeDBLAS m_BLAS;

        /// <summary>Cos
        /// most native BLAS interface available
        /// </summary>
        public static UnsafeDBLAS F77_BLAS { get { return m_BLAS; }}

        /// <summary>
        /// static ctor
        /// </summary>
        static BLAS() {
            m_BLAS = new UnsafeDBLAS();
        }

        /// <summary> FORTRAN-Style BLAS routine </summary>
        static public double DDOT(ref int N, double[] DX, ref int INCX, double[] DY, ref int INCY) {
            unsafe {
                fixed (double* pDX = &DX[0], pDY = &DY[0]) {
                    return m_BLAS.DDOT(ref N, pDX, ref INCX, pDY, ref INCY);
                }
            }
        }

        /// <summary> C-Style BLAS routine </summary>
        static public unsafe double ddot(int N, double* DX, int INCX, double* DY, int INCY) {
            return m_BLAS.DDOT(ref N, DX, ref INCX, DY, ref INCX);
        }

        /// <summary> C-Style BLAS routine </summary>
        static public double ddot(int N, double[] DX, int INCX, double[] DY, int INCY) {
            return DDOT(ref N, DX, ref INCX, DY, ref INCX);
        }

        /// <summary> FORTRAN-Style BLAS routine </summary>
        static public void DSWAP(ref int N, double[] DX, ref int INCX, double[] DY, ref int INCY) {
            unsafe {
                fixed (double* pDX = &DX[0], pDY = &DY[0]) {
                    m_BLAS.DSWAP(ref N, pDX, ref INCX, pDY, ref INCY);
                }
            }
        }

        /// <summary> C-Style BLAS routine </summary>
        static public void dswap(int N, double[] DX, int INCX, double[] DY, int INCY) {
            DSWAP(ref N, DX, ref INCX, DY, ref INCY);
        }

        /// <summary> C-Style BLAS routine </summary>
        unsafe static public void dswap(int N, double* DX, int INCX, double* DY, int INCY) {
            m_BLAS.DSWAP(ref N, DX, ref INCX, DY, ref INCY);
        }


        /// <summary> FORTRAN-Style BLAS routine </summary>
        static public void DGEMM(ref int TRANSA, ref int TRANSB,
                                 ref int M, ref int N, ref int K,
                                 ref double ALPHA,
                                 double[] A, ref int LDA,
                                 double[] B, ref int LDB,
                                 ref double BETA,
                                 double[] C, ref int LDC) {
            unsafe {
                fixed (double* pA = &A[0], pB = &B[0], pC = &C[0]) {
                    m_BLAS.DGEMM(ref TRANSA, ref TRANSB,
                                 ref M, ref N, ref K,
                                 ref ALPHA,
                                 pA, ref LDA,
                                 pB, ref LDB,
                                 ref BETA,
                                 pC, ref LDC);
                }
            }
        }

        /// <summary>
        /// native blas in C-stype
        /// (matrices are still in FORTRAN order)
        /// </summary>
        unsafe static public void dgemm(int TRANSA, int TRANSB,
                                        int M, int N, int K,
                                        double ALPHA,
                                        double* A, int LDA,
                                        double* B, int LDB,
                                        double BETA,
                                        double* C, int LDC) {
            unsafe {
                m_BLAS.DGEMM(ref TRANSA, ref TRANSB,
                             ref M, ref N, ref K,
                             ref ALPHA,
                             A, ref LDA,
                             B, ref LDB,
                             ref BETA,
                             C, ref LDC);
            }
        }


        /// <summary>
        /// C-style BLAS
        /// </summary>
        static public void daxpy(int N,
                                 double DA, double[] DX, int INCX,
                                 double[] DY, int INCY) {
            DAXPY(ref N, ref DA, DX, ref INCX, DY, ref INCY);
        }

        /// <summary>
        /// C-style BLAS
        /// </summary>
        unsafe static public void daxpy(int N,
                                 double DA, double* DX, int INCX,
                                 double* DY, int INCY) {
            m_BLAS.DAXPY(ref N, ref DA, DX, ref INCX, DY, ref INCY);
        }


        /// <summary>
        /// C-stype BLAS
        /// </summary>
        static public void dscal(int n, double a, double[] x, int incx) {
            DSCAL(ref n, ref a, x, ref incx);
        }

        /// <summary>
        /// C-stype BLAS
        /// </summary>
        unsafe static public void dscal(int n, double a, double* x, int incx) {
            m_BLAS.DSCAL(ref n, ref a, x, ref incx);
        }



        /// <summary>
        /// C-stype BLAS
        /// </summary>
        static public double dnrm2(int n, double[] x, int incx) {
            return DNRM2(ref n, x, ref incx);
        }
        
        /// <summary>
        /// FORTRAN-style BLAS
        /// </summary>
        static public void DAXPY(ref int N,
                                 ref double DA, double[] DX, ref int INCX,
                                 double[] DY, ref int INCY) {
            unsafe {
                fixed (double* pDX = &DX[0], pDY = &DY[0]) {
                    m_BLAS.DAXPY(ref N, ref DA, pDX, ref INCX, pDY, ref INCY);
                }
            }
        }


        /// <summary>
        /// FORTRAN-stype BLAS
        /// </summary>
        static public void DSCAL(ref int n, ref double a, double[] x, ref int incx) {
            unsafe {
                fixed (double* px = &x[0]) {
                    m_BLAS.DSCAL(ref n, ref a, px, ref incx);
                }
            }
        }


        /// <summary>
        /// FORTRAN-stype BLAS
        /// </summary>
        static public double DNRM2(ref int n, double[] x, ref int incx) {
            unsafe {
                fixed (double* px = &x[0]) {
                    return m_BLAS.DNRM2(ref n, px, ref incx);
                }
            }
        }
    }
}
