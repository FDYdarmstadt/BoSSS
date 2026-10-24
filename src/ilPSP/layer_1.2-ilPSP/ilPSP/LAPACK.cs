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
using System.Diagnostics;
using System.Runtime.InteropServices;
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {

    /// <summary>LAPACK convenience methods; provider and threading follow BLAS.</summary>
    public static class LAPACK {
        // Load imports only when LAPACK is first used, after native search-path bootstrap.
        static LAPACK() { }

        /// <summary>Sequential MKL import, for explicitly sequential calls.</summary>
        public static readonly UnsafeLAPACK seq_F77_LAPACK = new UnsafeLAPACK(Parallelism.SEQ);
        /// <summary>Threaded MKL import.</summary>
        public static readonly UnsafeLAPACK omp_F77_LAPACK = new UnsafeLAPACK(Parallelism.OMP);
        /// <summary>Current import, selected together with BLAS.</summary>
        public static UnsafeLAPACK F77_LAPACK => BLAS.ActiveLAPACK;

        /// <summary>Selected native DGECON_ entry point.</summary>
        public static UnsafeLAPACK._DGECON DGECON_ => F77_LAPACK.DGECON_;

        /// <summary>Selected native DPOSV_ entry point.</summary>
        public static UnsafeLAPACK._DPOSV DPOSV_ => F77_LAPACK.DPOSV_;

        /// <summary>Selected native DPOTRF_ entry point.</summary>
        public static UnsafeLAPACK._DPOTRF DPOTRF_ => F77_LAPACK.DPOTRF_;

        /// <summary>Selected native DPOTRI_ entry point.</summary>
        public static UnsafeLAPACK._DPOTRI DPOTRI_ => F77_LAPACK.DPOTRI_;

        /// <summary>Selected native DTRTRI_ entry point.</summary>
        public static UnsafeLAPACK._DTRTRI DTRTRI_ => F77_LAPACK.DTRTRI_;

        /// <summary>Selected native DSYTRF_ entry point.</summary>
        public static UnsafeLAPACK._DSYTRF DSYTRF_ => F77_LAPACK.DSYTRF_;

        /// <summary>Selected native DGETRF entry point.</summary>
        public static UnsafeLAPACK._DGETRF DGETRF => F77_LAPACK.DGETRF;

        /// <summary>Selected native DGETRS entry point.</summary>
        public static UnsafeLAPACK._DGETRS DGETRS => F77_LAPACK.DGETRS;

        /// <summary>Selected native DGETRI_ entry point.</summary>
        public static UnsafeLAPACK._DGETRI DGETRI_ => F77_LAPACK.DGETRI_;

        /// <summary>Selected native DGELSY_ entry point.</summary>
        public static UnsafeLAPACK._DGELSY DGELSY_ => F77_LAPACK.DGELSY_;

        /// <summary>Selected native DGELSS_ entry point.</summary>
        public static UnsafeLAPACK._DGELSS DGELSS_ => F77_LAPACK.DGELSS_;

        /// <summary>Selected native DGEQP3_ entry point.</summary>
        public static UnsafeLAPACK._DGEQP3 DGEQP3_ => F77_LAPACK.DGEQP3_;

        /// <summary>Selected native DGEQRF_ entry point.</summary>
        public static UnsafeLAPACK._DGEQRF DGEQRF_ => F77_LAPACK.DGEQRF_;

        /// <summary>Selected native DORGQR_ entry point.</summary>
        public static UnsafeLAPACK._DORGQR DORGQR_ => F77_LAPACK.DORGQR_;

        /// <summary>Selected native DGESVD_ entry point.</summary>
        public static UnsafeLAPACK._DGESVD DGESVD_ => F77_LAPACK.DGESVD_;

        /// <summary>Selected native DGEEV_ entry point.</summary>
        public static UnsafeLAPACK._DGEEV DGEEV_ => F77_LAPACK.DGEEV_;

        /// <summary>Selected native DSYEV_ entry point.</summary>
        public static UnsafeLAPACK._DSYEV DSYEV_ => F77_LAPACK.DSYEV_;

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public static void DPOTRF(char UPLO, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            int _UPLO = UPLO;

            unsafe {
                fixed (double* pa = A) {
                    DPOTRF_(ref _UPLO, ref N, pa, ref LDA, out INFO);
                }
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public static void DPOTRI(char UPLO, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            int _UPLO = UPLO;

            unsafe {
                fixed (double* pa = A) {
                    DPOTRI_(ref _UPLO, ref N, pa, ref LDA, out INFO);
                }
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public static void DTRTRI(int UPLO, int DIAG, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            if (DIAG != 'N' && DIAG != 'U')
                throw new ArgumentOutOfRangeException();

            int _UPLO = UPLO, _DIAG = DIAG;

            unsafe {
                fixed (double* pa = A) {
                    DTRTRI_(ref _UPLO, ref _DIAG, ref N, pa, ref LDA, out INFO);
                }
            }
        }

        /// <summary>
        /// factorization of a real symmetric matrix A
        /// </summary>
        unsafe public static void DSYTRF(int UPLO, int N, double* A, int LDA, int* IPIV, double* WORK, int* LWORK, out int INFO)
        {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();

            INFO = 0;
            DSYTRF_(ref UPLO, ref N, A, ref LDA, IPIV, WORK, LWORK, ref INFO);
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public static void DGETRFF(ref int M, ref int N, double[] A, ref int LDA, int[] IPIV, out int INFO) {
            unsafe {
                fixed (double* pa = A) {
                    fixed (int* pIPIV = IPIV) {
                        F77_LAPACK.DGETRF(ref M, ref N, pa, ref LDA, pIPIV, out INFO);
                    }
                }
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public static void DGETRI(ref int N, double[] A, ref int LDA, int[] IPIV, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pa = A, pwork = WORK) {
                    fixed (int* pIPIV = IPIV) {
                        F77_LAPACK.DGETRI_(ref N, pa, ref LDA, pIPIV, pwork, ref LWORK, out INFO);
                    }
                }
            }
        }

        /// <summary>
        /// Solves the linear least squares system
        /// <paramref name="A"/> x = <paramref name="B"/> using complete
        /// orthogonal factorization. Variant of
        /// <see cref="DGELSY(ref int, ref int, ref int, double[], ref int, double[], ref int, int[], ref double, out int, double[], ref int, out int)"/>
        /// with simplified interface/calling sequence.
        /// </summary>
        /// <param name="M">
        /// Number of rows of matrix <paramref name="A"/>
        /// </param>
        /// <param name="N">
        /// Number of columns of matrix <paramref name="A"/>
        /// </param>
        /// <param name="A">
        /// <paramref name="M"/>x<paramref name="N"/> system-matrix of the
        /// linear least squares system to be solved. Note that the entries
        /// have to be in FORTRAN-order
        /// </param>
        /// <param name="B">
        /// Right-hand side of the linear least squares system to be solved.
        /// On exit: Contains the result vector x
        /// </param>
        /// <param name="NRHS">
        /// Number of right hand sides.
        /// </param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        public static void DGELSY(int M, int N, double[] A, double[] B, int NRHS, double RCOND) {
            {
                Debug.Assert(
                    A.Length == M * N,
                    "A must be a MxN matrix");
                Debug.Assert(
                    B.Length == Math.Max(M, N) * NRHS,
                    "B must be a vector of length max(M,N) * NRHS");
                unsafe {
                    fixed (double* pA = A, pB = B) {
                        // In order to make results deterministic, the start address of
                        // A and B must not be different on different calls. For now,
                        // we ensure this by allocating memory by hand (thus, we avoid
                        // the non-deterministic influence of the garbage collector)
                        IntPtr _pA = Marshal.AllocHGlobal(A.Length * sizeof(double));
                        Marshal.Copy(A, 0, _pA, A.Length);

                        IntPtr _pB = Marshal.AllocHGlobal(B.Length * sizeof(double));
                        Marshal.Copy(B, 0, _pB, B.Length);

                        int LDB = B.Length / NRHS;
                        DGELSY(M, N, (double*)_pA, (double*)_pB, LDB, NRHS, RCOND);

                        Marshal.Copy(_pB, B, 0, B.Length);
                        Marshal.FreeHGlobal(_pA);
                        Marshal.FreeHGlobal(_pB);
                    }
                }
            }
        }

        /// <summary>
        /// Solves the linear least squares system
        /// A x = B (multiple right-hand sides) using complete orthogonal
        /// factorization, where matrix A starts at <paramref name="pA"/> and B
        /// starts at <paramref name="pB"/>. Variant of
        /// <see cref="DGELSY(ref int, ref int, ref int, double[], ref int, double[], ref int, int[], ref double, out int, double[], ref int, out int)"/>
        /// with simplified interface/calling sequence.
        /// </summary>
        /// <param name="M">
        /// Number of rows of matrix <paramref name="pA"/>
        /// </param>
        /// <param name="N">
        /// Number of columns of matrix <paramref name="pB"/>
        /// </param>
        /// <param name="pA">
        /// <paramref name="M"/>x<paramref name="N"/> system-matrix of the
        /// linear least squares system to be solved. Note that the entries
        /// have to be in FORTRAN-order
        /// </param>
        /// <param name="pB">
        /// Right-hand side of the linear least squares system to be solved.
        /// On exit: Contains the result vector x
        /// </param>
        /// <param name="LDB">Leading dimension of B</param>
        /// <param name="NRHS">
        /// Number of right hand sides.
        /// </param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        internal static unsafe void DGELSY(int M, int N, double* pA, double* pB, int LDB, int NRHS, double RCOND) {
            {
                int RANK;
                int LWORK = -1;
                int INFO;

                int* JPVT = stackalloc int[N];
                double* pLENGTH = stackalloc double[1];

                F77_LAPACK.DGELSY_(ref M, ref N, ref NRHS, (double*)pA, ref M, (double*)pB, ref LDB, JPVT, ref RCOND, out RANK, pLENGTH, ref LWORK, out INFO);

                if (INFO != 0) {
                    throw new ApplicationException("Failed to determine optimal work size");
                }

                // Don't use stackalloc for $work since its size may be so
                // large that the stack overflows!
                LWORK = (int)*pLENGTH;
                double[] work = new double[LWORK];
                fixed (double* pWORK = work) {
                    F77_LAPACK.DGELSY_(ref M, ref N, ref NRHS, (double*)pA, ref M, (double*)pB, ref LDB, JPVT, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                }

                if (INFO != 0) {
                    throw new ApplicationException("Failed to solve least squares system");
                }
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgelsy.f
        /// </summary>
        public static void DGELSY(ref int M, ref int N, ref int NRHS, double[] A, ref int LDA, double[] B, ref int LDB, int[] JPVT, ref double RCOND, out int RANK, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pA = A, pB = B, pWORK = WORK) {
                    fixed (int* pJPVT = JPVT) {
                        F77_LAPACK.DGELSY_(ref M, ref N, ref NRHS, pA, ref LDA, pB, ref LDB, pJPVT, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                    }
                }
            }
        }

        /// <summary>
        /// Solves the linear least squares system
        /// <paramref name="A"/> x = <paramref name="B"/> using singular value
        /// decomposition. Variant of
        /// <see cref="DGELSS(ref int, ref int, ref int, double[], ref int, double[], ref int, double[], ref double, out int, double[], ref int, out int)"/>
        /// with simplified interface/calling sequence.
        /// </summary>
        /// <param name="M">
        /// Number of rows of matrix <paramref name="A"/>
        /// </param>
        /// <param name="N">
        /// Number of columns of matrix <paramref name="A"/>
        /// </param>
        /// <param name="A">
        /// <paramref name="M"/>x<paramref name="N"/> system-matrix of the
        /// linear least squares system to be solved. Note that the entries
        /// have to be in FORTRAN-order
        /// </param>
        /// <param name="B">
        /// Right-hand side of the linear least squares system to be solved.
        /// On exit: Contains the result vector x
        /// </param>
        /// <param name="NRHS">
        /// Number of right hand sides.
        /// </param>
        /// <param name="RCOND">
        /// Condition that defines the size of the eigenvalues that are
        /// considered zero (e.g., to cope with numerical round-off). A
        /// negative number implies that all non-negative eigenvalues are taken
        /// into account
        /// </param>
        /// <returns>
        /// An estimate for the condition number of <paramref name="A"/>.
        /// </returns>
        public static unsafe double DGELSS(int M, int N, double[] A, double[] B, int NRHS, double RCOND) {
            {
                Debug.Assert(
                    A.Length == M * N,
                    "A must be a MxN matrix");
                Debug.Assert(
                    B.Length == Math.Max(M, N) * NRHS,
                    "B must be a vector of length max(M,N) * NRHS");

                fixed (double* pA = A, pB = B) {
                    // In order to make results deterministic, the start address of
                    // A and B must not be different on different calls. For now,
                    // we ensure this by allocating memory by hand (thus, we avoid
                    // the non-deterministic influence of the garbage collector)
                    IntPtr _pA = Marshal.AllocHGlobal(A.Length * sizeof(double));
                    Marshal.Copy(A, 0, _pA, A.Length);

                    IntPtr _pB = Marshal.AllocHGlobal(B.Length * sizeof(double));
                    Marshal.Copy(B, 0, _pB, B.Length);

                    int LDB = B.Length / NRHS;
                    int RANK;
                    int LWORK = -1;
                    int INFO;

                    double* pS = stackalloc double[Math.Min(M, N)];
                    double* pLENGTH = stackalloc double[1];

                    F77_LAPACK.DGELSS_(ref M, ref N, ref NRHS, (double*)_pA, ref M, (double*)_pB, ref LDB, pS, ref RCOND, out RANK, pLENGTH, ref LWORK, out INFO);

                    if (INFO != 0) {
                        throw new ApplicationException("Failed to determine optimal work size");
                    }

                    // Don't use stackalloc for $work since its size may be so
                    // large that the stack overflows!
                    LWORK = (int)*pLENGTH;
                    double[] work = new double[LWORK];
                    fixed (double* pWORK = work) {
                        F77_LAPACK.DGELSS_(ref M, ref N, ref NRHS, (double*)_pA, ref M, (double*)_pB, ref LDB, pS, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                    }

                    if (INFO != 0) {
                        throw new ApplicationException("Failed to solve least squares system");
                    }

                    Marshal.Copy(_pB, B, 0, B.Length);
                    Marshal.FreeHGlobal(_pA);
                    Marshal.FreeHGlobal(_pB);

                    return *pS / *(pS + Math.Min(M, N) - 1);
                }
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgelss.f
        /// </summary>
        public static void DGELSS(ref int M, ref int N, ref int NRHS, double[] A, ref int LDA, double[] B, ref int LDB, double[] S, ref double RCOND, out int RANK, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pA = A, pB = B, pS = S, pWORK = WORK) {
                    F77_LAPACK.DGELSS_(ref M, ref N, ref NRHS, pA, ref LDA, pB, ref LDB, pS, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                }
            }
        }

        /// <summary>
        /// DGEQP3 computes a QR factorization of a matrix <paramref name="A"/>
        /// with column pivoting (A*P = Q*R) using Level 3 BLAS
        /// </summary>
        /// <param name="M">
        /// The number of rows of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="N">
        /// The number of columns of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="A">
        /// <paramref name="M"/>x<paramref name="N"/> matrix to be factorized.
        /// Note that the entries have to be in FORTRAN-order
        /// </param>
        /// <param name="LDA">
        /// Leading dimension of <paramref name="A"/>
        /// </param>
        /// <param name="JPVT">
        /// On entry: If JPVT(J) &lt;&gt; 0, the J-th column of A is permuted
        /// to the front of A*P (a leading column); if JPVT(J)=0, the J-th
        /// column of A is a free column.<br/>
        /// On exit: If JPVT(J)=K, then the J-th column of A*P was the
        /// K-th column of A
        /// </param>
        /// <param name="TAU"></param>
        /// <remarks>
        /// See <see cref="DORGQR"/> for information about the contents of
        /// <paramref name="A"/> and <paramref name="TAU"/> on exit.
        /// </remarks>
        public static unsafe void DGEQP3(ref int M, ref int N, double* A, ref int LDA, int* JPVT, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double LENGTH;
            unsafe {
                F77_LAPACK.DGEQP3_(ref M, ref N, A, ref M, JPVT, TAU, &LENGTH, ref LWORK, out INFO);
            }
            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)LENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = work) {
                F77_LAPACK.DGEQP3_(ref M, ref N, A, ref M, JPVT, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine pivoted QR factorization");
            }
        }

        /// <summary>
        /// Computes a QR factorization of a real
        /// <paramref name="M"/>-by-<paramref name="N"/> matrix
        /// <paramref name="A"/>: A = Q * R
        /// </summary>
        /// <param name="M">
        /// The number of rows of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="N">
        /// The number of columns of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="A">
        /// On entry: The <paramref name="M"/>-by-<paramref name="N"/> matrix
        /// <paramref name="A"/>.<br />
        /// On exit: The elements on and above the diagonal of the array 
        /// contain the min(M,N)-by-N upper trapezoidal matrix R (R is upper
        /// triangular if M >= N); the elements below the diagonal, with the
        /// array <paramref name="TAU"/>, represent the orthogonal matrix Q as
        /// a product of min(m,n) elementary reflectors
        /// </param>
        /// <param name="LDA">
        /// The leading dimension of the array <paramref name="A"/>
        /// </param>
        /// <param name="TAU">
        /// The scalar factors of the elementary reflectors.
        /// </param>
        /// <remarks>
        /// See <see cref="DORGQR"/> for information about the contents of
        /// <paramref name="A"/> and <paramref name="TAU"/> on exit.
        /// </remarks>
        public static unsafe void DGEQRF(ref int M, ref int N, double* A, ref int LDA, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double* pLENGTH = stackalloc double[1];
            F77_LAPACK.DGEQRF_(ref M, ref N, A, ref M, TAU, pLENGTH, ref LWORK, out INFO);

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)*pLENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = work) {
                F77_LAPACK.DGEQRF_(ref M, ref N, A, ref M, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine QR factorization");
            }
        }

        /// <summary>
        /// Generates an <paramref name="M"/>-by-<paramref name="N"/> real
        /// matrix Q with orthonormal columns, which is defined as the first
        /// <paramref name="N"/> columns of a product of <paramref name="K"/>
        /// elementary reflectors of order <paramref name="M"/>
        /// </summary>
        /// <param name="M">
        /// The number of rows of the matrix Q
        /// </param>
        /// <param name="N">
        /// The number of columns of the matrix Q
        /// </param>
        /// <param name="K">
        /// The number of elementary reflectors whose product defines the
        /// matrix Q
        /// </param>
        /// <param name="A">
        /// On entry: The i-th column must contain the vector which defines the
        /// elementary reflector H(i), for i = 1,2,...,k, as returned by
        /// <see cref="DGEQP3"/> in the first k columns of its array argument
        /// <paramref name="A"/>.<br/>
        /// On exit: The <paramref name="M"/>-by-<paramref name="N"/> matrix Q.
        /// </param>
        /// <param name="LDA">
        /// Leading dimension of <paramref name="A"/>
        /// </param>
        /// <param name="TAU">
        /// The scalar factor of the elementary reflector H(i), as returned by
        /// <see cref="DGEQP3"/>
        /// </param>
        public static unsafe void DORGQR(ref int M, ref int N, ref int K, double* A, ref int LDA, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double* pLENGTH = stackalloc double[1];
            F77_LAPACK.DORGQR_(ref M, ref N, ref K, A, ref LDA, TAU, pLENGTH, ref LWORK, out INFO);

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)*pLENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = work) {
                F77_LAPACK.DORGQR_(ref M, ref N, ref K, A, ref LDA, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to construct matrix Q");
            }
        }

        /// <summary>
        /// DGESVD computes the singular value decomposition (SVD) of a real
        /// M-by-N matrix A, optionally computing the left and/or right singular
        /// vectors.
        /// </summary>
        /// <param name="JOBU">
        /// JOBU is CHARACTER*1
        /// Specifies options for computing all or part of the matrix U:
        ///   = 'A':  all M columns of U are returned in array U:
        ///   = 'S':  the first min(m, n) columns of U(the left singular
        ///            vectors) are returned in the array U;
        ///   = 'O':  the first min(m, n) columns of U(the left singular
        ///            vectors) are overwritten on the array A;
        ///   = 'N':  no columns of U(no left singular vectors) are
        ///          computed.
        /// </param>
        /// <param name="JOBVT">
        /// JOBVT is CHARACTER*1
        /// Specifies options for computing all or part of the matrix
        ///   V** T:
        ///   = 'A':  all N rows of V** T are returned in the array VT;
        ///   = 'S':  the first min(m, n) rows of V** T(the right singular
        ///            vectors) are returned in the array VT;
        ///   = 'O':  the first min(m, n) rows of V** T(the right singular
        ///            vectors) are overwritten on the array A;
        ///   = 'N':  no rows of V**T(no right singular vectors) are
        ///          computed.
        /// 
        ///  JOBVT and JOBU cannot both be 'O'.
        /// </param>
        /// <param name="M">
        /// The number of rows of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="N">
        /// The number of columns of the matrix <paramref name="A"/>
        /// </param>
        /// <param name="A">
        /// <paramref name="M"/>x<paramref name="N"/> matrix to be decomposited.
        /// Note that the entries have to be in FORTRAN-order
        /// </param>
        /// <param name="LDA">
        /// Leading dimension of <paramref name="A"/>
        /// </param>
        /// <param name="S">
        /// Singular values of <paramref name="A"/>, sorted in decending order
        /// </param>
        /// <param name="U">
        /// Array containing the left singular vectors of <paramref name="A"/>
        /// </param>
        /// <param name="LDU">
        /// Leading dimension of <paramref name="U"/>
        /// </param>
        /// <param name="VT">
        /// Array containing the right singular vectors of <paramref name="A"/>
        /// </param>
        /// <param name="LDVT">
        /// Leading dimension of <paramref name="VT"/>
        /// </param>
        public static unsafe void DGESVD(int JOBU, int JOBVT, int M, int N, double* A, int LDA, double* S, double* U, int LDU, double* VT, int LDVT)
        {
            throw new NotImplementedException("temporarily deactivated");
            /*
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double LENGTH;
            unsafe
            {
                dgesvd(ref JOBU, ref JOBVT, ref M, ref N, A, ref LDA, S, U, ref LDU, VT, ref LDVT, &LENGTH, ref LWORK, out INFO);
            }
            if (INFO != 0)
            {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)LENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = &work[0])
            {
                dgesvd(ref JOBU, ref JOBVT, ref M, ref N, A, ref LDA, S, U, ref LDU, VT, ref LDVT, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0)
            {
                throw new ApplicationException("Failed to compute Singular Value Decomposition");
         
            }
            */
        }

    }
}
