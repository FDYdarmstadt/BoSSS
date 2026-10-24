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
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {

    /// <summary>Native FORTRAN-style LAPACK imports.</summary>
    public sealed class UnsafeLAPACK : DynLibLoader {
        /// <summary>Loads all required LAPACK exports for the chosen provider.</summary>
        public UnsafeLAPACK(Parallelism parallelism, BLASProvider provider = BLASProvider.MKL)
            : base(BLAS_LAPACK_Libraries.Get(provider, parallelism)) {
            Provider = provider;
        }

        internal UnsafeLAPACK(UnsafeDBLAS owner) : base(owner, BoSSS_Prefix) {
            Provider = owner.Provider;
        }

        /// <summary>The implementation behind this import.</summary>
        public BLASProvider Provider { get; }

#pragma warning disable        649
        _DGETRF dgetrf;
        public unsafe delegate void _DGETRS(ref char transa, ref int n, ref int nrhs, double* a, ref int lda, int* ipiv, double* b, ref int ldb, out int info);

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
        _DPOSV  dposv;
        //_DGESVD dgesvd;
        //_DSYTRF dsytrf;
        _DGECON dgecon;
        _DGEEV dgeev;
        _DSYEV dsyev;
#pragma warning restore 649

        /// <summary>
        /// Reciprocal condition number
        /// </summary>
        public unsafe delegate void _DGECON(ref int NORM, ref int N, double* A, ref int LDA, ref double ANORM, ref double RCOND, double* WORK, int* IWORK, ref int INFO);

        /// <summary>
        /// Reciprocal condition number
        /// </summary>
        public _DGECON DGECON_ {
            get {
                return dgecon;
            }
        }

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
        /// factorization of a real symmetric matrix A
        /// </summary>
        public unsafe delegate void _DSYTRF(ref int UPLO, ref int N, double* A, ref int LDA, int* IPIV, double* WORK, int* LWORK, ref int INFO);

        /// <summary>
        /// factorization of a real symmetric matrix A
        /// </summary>
        public unsafe _DSYTRF DSYTRF_ {
            get {
                return null;
                //return dsytrf;
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
        /// See http://www.netlib.org/lapack/double/dgelsy.f
        /// </summary>
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
        /// See: http://www.netlib.org/lapack/double/dorgqr.f"
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

        /// <summary>
        /// Singular Value Decomposition
        /// </summary>
        public unsafe delegate void _DGESVD(ref int JOBU, ref int JOBVT, ref int M, ref int N, double* A, ref int LDA, double* S, double* U, ref int LDU, double* VT, ref int LDVT, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// Singular Value Decomposition
        /// </summary>
        public unsafe _DGESVD DGESVD_
        {
            get
            {
                throw new NotImplementedException("temporarily deactivated");
                //return dgesvd;
            }
        }

        /// <summary>
        /// General matrix eigenvalues.
        /// </summary>
        public unsafe delegate void _DGEEV(ref int JOBVL, ref int JOBVR, ref int N, double* A, ref int LDA, double* WR, double* WI, double* VL, ref int LDVL, double* VR, ref int LDVR, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DGEEV"/>
        /// </summary>
        public unsafe _DGEEV DGEEV_ {
            get {
                return dgeev;
            }
        }

        /// <summary>
        /// Symmetrical matrix eigenvalues.
        /// </summary>
        public unsafe delegate void _DSYEV(ref int JOBZ, ref int UPLO, ref int N, double* A, ref int LDA, double* W, double* WORK, ref int LWORK, out int INFO);

        /// <summary>
        /// <see cref="_DSYEV"/>
        /// </summary>
        public unsafe _DSYEV DSYEV_ {
            get {
                return dsyev;
            }
        }
    }
}
