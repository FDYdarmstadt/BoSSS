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
using ilPSP.Tracing;

namespace ilPSP.Utils {

    /// <summary>
    /// some parts of the LAPACK interface (FORTRAN 77 - style -- interface), which are used by BoSSS;
    /// </summary>
    public class LAPACK : DynLibLoader {

        static LAPACK _F77_LAPACK;

        /// <summary>
        /// entry point to FORTRAN 77 - style LAPACK
        /// </summary>
        public static LAPACK F77_LAPACK {
            get {
                return LAPACK._F77_LAPACK;
            }
        }

        static LAPACK() {
            _F77_LAPACK = new LAPACK();
        }

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
        public LAPACK() :
            base(new string[] { "libacml_dll.dll", "libacml.so*", "libacml.so", "liblapack.so*", "liblapack.so" },
                 new string[5][][],
                 new GetNameMangling[] { DynLibLoader.CAPITAL_LETTERS, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.SmallLetters_TrailingUnderscore },
                 Helper(),
                 new int[] { -1, -1, -1, -1, -1 }) {
        }

#pragma warning disable        649
        _DGETRF dgetrf;
        _DGETRS dgetrs;
        _DGETRI dgetri;
        _DGELSY dgelsy;
        _DGELSS dgelss;
        _DPOTRF dpotrf;
        _DTRTRI dtrtri;
        _DPOTRI dpotri;
        _DGEQRF dgeqrf;
        _DGEQP3 dgeqp3;
        _DORGQR dorgqr;
        _DPOSV dposv;
#pragma warning restore 649

        /// <summary>
        /// Solution by Cholesky factorization
        /// </summary>
        public unsafe delegate void _DPOSV(ref int UPLO, ref int N, ref int NRHS, double* A, ref int LDA, double* B, ref int LDB, out int INFO);

        /// <summary>
        /// Solution by Cholesky factorization
        /// </summary>
        public unsafe _DPOSV DPOSV_ {
            get {
                return dposv;
            }
        }

        /// <summary>
        /// Cholesky factorization
        /// </summary>
        public unsafe delegate void _DPOTRF(ref int UPLO, ref int N, double* A, ref int LDA, out int INFO);

        /// <summary>
        /// Cholesky-factorization
        /// </summary>
        public unsafe _DPOTRF DPOTRF_ {
            get {
                return dpotrf;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public void DPOTRF(char UPLO, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            int _UPLO = UPLO;

            unsafe {
                fixed (double* pa = &A[0]) {
                    DPOTRF_(ref _UPLO, ref N, pa, ref LDA, out INFO);
                }
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public unsafe delegate void _DPOTRI(ref int UPLO, ref int N, double* A, ref int LDA, out int INFO);

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public unsafe _DPOTRI DPOTRI_ {
            get {
                return dpotri;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public void DPOTRI(char UPLO, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            int _UPLO = UPLO;

            unsafe {
                fixed (double* pa = &A[0]) {
                    DPOTRI_(ref _UPLO, ref N, pa, ref LDA, out INFO);
                }
            }
        }


        /// <summary>
        /// triangular matrix inversion
        /// </summary>
        public unsafe delegate void _DTRTRI(ref int UPLO, ref int DIAG, ref int N, double* A, ref int LDA, out int INFO);

        /// <summary>
        /// triangular matrix inversion
        /// </summary>
        public unsafe _DTRTRI DTRTRI_ {
            get {
                return dtrtri;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public void DTRTRI(char UPLO, char DIAG, int N, double[] A, int LDA, out int INFO) {
            if (UPLO != 'U' && UPLO != 'L')
                throw new ArgumentOutOfRangeException();
            if (DIAG != 'N' && DIAG != 'U')
                throw new ArgumentOutOfRangeException();

            int _UPLO = UPLO, _DIAG = DIAG;


            unsafe {
                fixed (double* pa = &A[0]) {
                    DTRTRI_(ref _UPLO, ref _DIAG, ref N, pa, ref LDA, out INFO);
                }
            }
        }


        /// <summary>
        /// FORTRAN-style LAPACK, matrices in FORTRAN order;
        /// </summary>
        public unsafe delegate void _DGETRF(ref int m, ref int n, double* a, ref int lda, int* ipiv, out int info);

        /// <summary>
        /// FORTRAN-style LAPACK, matrices in FORTRAN order;
        /// </summary>
        public unsafe _DGETRF DGETRF {
            get {
                return dgetrf;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public void DGETRFF(ref int M, ref int N, double[] A, ref int LDA, int[] IPIV, out int INFO) {
            unsafe {
                fixed (double* pa = A) {
                    fixed (int* pIPIV = IPIV) {
                        dgetrf(ref M, ref N, pa, ref LDA, pIPIV, out INFO);
                    }
                }
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        unsafe public delegate void _DGETRS(ref char transa, ref int n, ref int nrhs, double* a, ref int lda, int* ipiv, double[] b, ref int ldb, out int info);

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public unsafe _DGETRS DGETRS {
            get {
                return dgetrs;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public unsafe delegate void _DGETRI(ref int N, double* A, ref int LDA, int* IPIV, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public unsafe _DGETRI DGETRI_ {
            get {
                return dgetri;
            }
        }

        /// <summary>
        /// FORTRAN-style LAPACK
        /// </summary>
        public void DGETRI(ref int N, double[] A, ref int LDA, int[] IPIV, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pa = A, pwork = WORK) {
                    fixed (int* pIPIV = IPIV) {
                        dgetri(ref N, pa, ref LDA, pIPIV, pwork, ref LWORK, out INFO);
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
        public void DGELSY(int M, int N, double[] A, double[] B, int NRHS, double RCOND) {
            using (new FuncTrace()) {
                Debug.Assert(
                    A.Length == M * N,
                    "A must be a MxN matrix");
                Debug.Assert(
                    B.Length == Math.Max(M, N) * NRHS,
                    "B must be a vector of length max(M,N) * NRHS");
                unsafe {
                    fixed (double* pA = &A[0], pB = &B[0]) {
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
        unsafe internal void DGELSY(int M, int N, double* pA, double* pB, int LDB, int NRHS, double RCOND) {
            using (new FuncTrace()) {
                int RANK;
                int LWORK = -1;
                int INFO;

                int* JPVT = stackalloc int[N];
                double* pLENGTH = stackalloc double[1];

                dgelsy(ref M, ref N, ref NRHS, (double*)pA, ref M, (double*)pB, ref LDB, JPVT, ref RCOND, out RANK, pLENGTH, ref LWORK, out INFO);

                if (INFO != 0) {
                    throw new ApplicationException("Failed to determine optimal work size");
                }

                // Don't use stackalloc for $work since its size may be so
                // large that the stack overflows!
                LWORK = (int)*pLENGTH;
                double[] work = new double[LWORK];
                fixed (double* pWORK = &work[0]) {
                    dgelsy(ref M, ref N, ref NRHS, (double*)pA, ref M, (double*)pB, ref LDB, JPVT, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                }

                if (INFO != 0) {
                    throw new ApplicationException("Failed to solve least squares system");
                }
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgelsy.f
        /// </summary>
        /// <param name="M"></param>
        /// <param name="N"></param>
        /// <param name="NRHS"></param>
        /// <param name="A"></param>
        /// <param name="LDA"></param>
        /// <param name="B"></param>
        /// <param name="LDB"></param>
        /// <param name="JPVT"></param>
        /// <param name="RCOND"></param>
        /// <param name="RANK"></param>
        /// <param name="WORK"></param>
        /// <param name="LWORK"></param>
        /// <param name="INFO"></param>
        public void DGELSY(ref int M, ref int N, ref int NRHS, double[] A, ref int LDA, double[] B, ref int LDB, int[] JPVT, ref double RCOND, out int RANK, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pA = &A[0], pB = &B[0], pWORK = &WORK[0]) {
                    fixed (int* pJPVT = &JPVT[0]) {
                        dgelsy(ref M, ref N, ref NRHS, pA, ref LDA, pB, ref LDB, pJPVT, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                    }
                }
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgelsy.f
        /// </summary>
        /// <param name="M"></param>
        /// <param name="N"></param>
        /// <param name="NRHS"></param>
        /// <param name="A"></param>
        /// <param name="LDA"></param>
        /// <param name="B"></param>
        /// <param name="LDB"></param>
        /// <param name="JPVT"></param>
        /// <param name="RCOND"></param>
        /// <param name="RANK"></param>
        /// <param name="WORK"></param>
        /// <param name="LWORK"></param>
        /// <param name="INFO"></param>
        public unsafe delegate void _DGELSY(ref int M, ref int N, ref int NRHS, double* A, ref int LDA, double* B, ref int LDB, int* JPVT, ref double RCOND, out int RANK, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DGELSY"/>
        /// </summary>
        public unsafe _DGELSY DGELSY_ {
            get {
                return dgelsy;
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
        public unsafe double DGELSS(int M, int N, double[] A, double[] B, int NRHS, double RCOND) {
            using (new FuncTrace()) {
                Debug.Assert(
                    A.Length == M * N,
                    "A must be a MxN matrix");
                Debug.Assert(
                    B.Length == Math.Max(M, N) * NRHS,
                    "B must be a vector of length max(M,N) * NRHS");

                fixed (double* pA = &A[0], pB = &B[0]) {
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

                    dgelss(ref M, ref N, ref NRHS, (double*)_pA, ref M, (double*)_pB, ref LDB, pS, ref RCOND, out RANK, pLENGTH, ref LWORK, out INFO);

                    if (INFO != 0) {
                        throw new ApplicationException("Failed to determine optimal work size");
                    }

                    // Don't use stackalloc for $work since its size may be so
                    // large that the stack overflows!
                    LWORK = (int)*pLENGTH;
                    double[] work = new double[LWORK];
                    fixed (double* pWORK = &work[0]) {
                        dgelss(ref M, ref N, ref NRHS, (double*)_pA, ref M, (double*)_pB, ref LDB, pS, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
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
        public void DGELSS(ref int M, ref int N, ref int NRHS, double[] A, ref int LDA, double[] B, ref int LDB, double[] S, ref double RCOND, out int RANK, double[] WORK, ref int LWORK, out int INFO) {
            unsafe {
                fixed (double* pA = &A[0], pB = &B[0], pS = &S[0], pWORK = &WORK[0]) {
                    dgelss(ref M, ref N, ref NRHS, pA, ref LDA, pB, ref LDB, pS, ref RCOND, out RANK, pWORK, ref LWORK, out INFO);
                }
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgelss.f
        /// </summary>
        public unsafe delegate void _DGELSS(ref int M, ref int N, ref int NRHS, double* A, ref int LDA, double* B, ref int LDB, double* S, ref double RCOND, out int RANK, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DGELSS"/>
        /// </summary>
        public unsafe _DGELSS DGELSS_ {
            get {
                return dgelss;
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
        public unsafe void DGEQP3(ref int M, ref int N, double* A, ref int LDA, int* JPVT, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double* pLENGTH = stackalloc double[1];
            dgeqp3(ref M, ref N, A, ref M, JPVT, TAU, pLENGTH, ref LWORK, out INFO);

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)*pLENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = &work[0]) {
                dgeqp3(ref M, ref N, A, ref M, JPVT, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine pivoted QR factorization");
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgeqp3.f
        /// </summary>
        public unsafe delegate void _DGEQP3(ref int M, ref int N, double* A, ref int LDA, int* JPVT, double* TAU, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DGEQP3"/>
        /// </summary>
        public unsafe _DGEQP3 DGEQP3_ {
            get {
                return dgeqp3;
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
        public unsafe void DGEQRF(ref int M, ref int N, double* A, ref int LDA, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double* pLENGTH = stackalloc double[1];
            dgeqrf(ref M, ref N, A, ref M, TAU, pLENGTH, ref LWORK, out INFO);

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)*pLENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = &work[0]) {
                dgeqrf(ref M, ref N, A, ref M, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine QR factorization");
            }
        }

        /// <summary>
        /// See http://www.netlib.org/lapack/double/dgeqrf.f
        /// </summary>
        public unsafe delegate void _DGEQRF(ref int M, ref int N, double* A, ref int LDA, double* TAU, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DGEQRF"/>
        /// </summary>
        public unsafe _DGEQRF DGEQRF_ {
            get {
                return dgeqrf;
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
        public unsafe void DORGQR(ref int M, ref int N, ref int K, double* A, ref int LDA, double* TAU) {
            // Determine optimal workspace size
            int INFO;
            int LWORK = -1;
            double* pLENGTH = stackalloc double[1];
            dorgqr(ref M, ref N, ref K, A, ref LDA, TAU, pLENGTH, ref LWORK, out INFO);

            if (INFO != 0) {
                throw new ApplicationException("Failed to determine optimal work size");
            }

            // Don't use stackalloc for $work since its size may be so
            // large that the stack overflows!
            LWORK = (int)*pLENGTH;
            double[] work = new double[LWORK];

            // Actual computation
            fixed (double* pWORK = &work[0]) {
                dorgqr(ref M, ref N, ref K, A, ref LDA, TAU, pWORK, ref LWORK, out INFO);
            }

            if (INFO != 0) {
                throw new ApplicationException("Failed to construct matrix Q");
            }
        }

        /// <summary>
        /// See <see cref="http://www.netlib.org/lapack/double/dorgqr.f"/>
        /// </summary>
        public unsafe delegate void _DORGQR(ref int M, ref int N, ref int K, double* A, ref int LDA, double* TAU, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DORGQR"/>
        /// </summary>
        public unsafe _DORGQR DORGQR_ {
            get {
                return dorgqr;
            }
        }
    }
}
