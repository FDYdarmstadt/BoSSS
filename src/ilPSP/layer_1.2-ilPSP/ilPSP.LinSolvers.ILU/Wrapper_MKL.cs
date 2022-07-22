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
using System.Text;
using MPI.Wrappers.Utils;
using ilPSP;
using System.Security;

namespace ilPSP.LinSolvers.ILU {

#pragma warning disable 1591

    /// <summary>
    /// raw function wrapper for loading the ILU pre-conditioner from the Intel MKL libraries,
    /// see https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/precondition-based-on-incomplete-lu-factorization.html
    /// </summary>
    /// <remarks>
    /// Licensing: despite being closed software, Intel MKL can be redistributed.
    /// </remarks>
    [SuppressUnmanagedCodeSecurity]
    public class Wrapper_MKL : DynLibLoader {

        /// <summary>
        /// Read from Environment which type of parallel library should be used.
        /// Returns a list of libraries in specific order to search for.
        /// </summary>
        static string[] SelectLibrary(Parallelism par) {
            string[] liborder;
            switch(par) {

                case Parallelism.SEQ:
                liborder = new string[] { "ILU.dll", "libBoSSSnative_seq.so", "libBoSSSnative_omp.so" };
                break;

                default:
                throw new NotSupportedException($"Unsupported level of parallelism {par} for Intel MKL ILU library.");
            }
            return liborder;
        }


        /// <summary>
        /// ctor
        /// </summary>
        public Wrapper_MKL() : base(
            SelectLibrary(Parallelism.SEQ),
            new string[3][][],
            new GetNameMangling[] { DynLibLoader.SmallLetters_TrailingUnderscore, DynLibLoader.BoSSS_Prefix, DynLibLoader.BoSSS_Prefix },
            new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix, PlatformID.Unix },
            new int[] { -1, -1, -1 }) {
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="n">matrix size</param>
        /// <param name="a"></param>
        /// <param name="ia"></param>
        /// <param name="ja"></param>
        /// <param name="bilu0"></param>
        /// <param name="ipar"></param>
        /// <param name="dpar"></param>
        /// <param name="ierr"></param>
        /// <returns></returns>
        public unsafe delegate int _dcsrilu0(int* n, double* a, int* ia, int* ja, double* bilu0, int* ipar, double* dpar, int* ierr);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="a"></param>
        /// <param name="ia"></param>
        /// <param name="ja"></param>
        /// <param name="bilut"></param>
        /// <param name="ibilut"></param>
        /// <param name="jbilut"></param>
        /// <param name="tol"></param>
        /// <param name="maxfil"></param>
        /// <param name="ipar"></param>
        /// <param name="dpar"></param>
        /// <param name="ierr"></param>
        /// <returns></returns>
        public unsafe delegate int _dcsrilut(int* n, double* a, int* ia, int* ja,double* bilut,int* ibilut, int* jbilut,double* tol, int* maxfil, int* ipar, double* dpar, int* ierr);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="operation">
        /// SPARSE_OPERATION_NON_TRANSPOSE      = 10,
        /// SPARSE_OPERATION_TRANSPOSE          = 11,
        /// SPARSE_OPERATION_CONJUGATE_TRANSPOSE= 12
        /// </param>
        /// <param name="alpha"></param>
        /// <param name="A"></param>
        /// <param name="descr"></param>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public unsafe delegate int _d_trsv(sparse_operation_t operation, double alpha, IntPtr A, MATRIX_DESCR descr, double* y, double* x);

        /// <summary>
        /// 
        /// </summary>
        /// <param name="A">output, CSR (compressed sparse row) format handle</param>
        /// <param name="indexing">c-style 0, fortran-style 1</param>
        /// <param name="rows">number of rows</param>
        /// <param name="cols">number of columns</param>
        /// <param name="rows_start">pointerB</param>
        /// <param name="rows_end">pointerE</param>
        /// <param name="col_indx">column indices</param>
        /// <param name="values">matrix entries</param>
        /// <returns></returns>
        public unsafe delegate int _d_crcsr(IntPtr* A, sparse_index_base_t indexing, int rows, int cols, int* rows_start, int* rows_end, int* col_indx, double* values);

        /// <summary>
        /// For debugging
        /// </summary>
        /// <param name="source"></param>
        /// <param name="indexing"></param>
        /// <param name="rows"></param>
        /// <param name="cols"></param>
        /// <param name="rows_start"></param>
        /// <param name="rows_end"></param>
        /// <param name="col_indx"></param>
        /// <param name="values"></param>
        /// <returns></returns>
        public unsafe delegate int _export(IntPtr source, out sparse_index_base_t indexing, out int rows, out int cols, out int* rows_start, out int* rows_end, out int* col_indx, out double* values);

        /// <summary>
        /// destroys MKL internal CSR format matrix
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public unsafe delegate int _destroy(IntPtr A);

        /// <summary>
        /// optimizes MKL internal CSR format matrix
        /// </summary>
        /// <param name="A"></param>
        /// <returns></returns>
        public unsafe delegate int _optimize(IntPtr A);

#pragma warning disable 649
        _dcsrilu0 dcsrilu0;
        _dcsrilut dcsrilut;
        _d_trsv mkl_sparse_d_trsv;
        _d_crcsr mkl_sparse_d_create_csr;
        _destroy mkl_sparse_destroy;
        _optimize mkl_sparse_optimize;
        _export mkl_sparse_d_export_csr;
#pragma warning restore 649

        public struct MATRIX_DESCR {
            public sparse_matrix_type_t type;
            public sparse_fill_mode_t mode;
            public sparse_diag_type_t diag;
        }

        public struct SPARSE_MATRIX_T {}

        public enum sparse_matrix_type_t {
            SPARSE_MATRIX_TYPE_GENERAL = 20,   /*    General case                    */
            SPARSE_MATRIX_TYPE_SYMMETRIC = 21,   /*    Triangular part of              */
            SPARSE_MATRIX_TYPE_HERMITIAN = 22,   /*    the matrix is to be processed   */
            SPARSE_MATRIX_TYPE_TRIANGULAR = 23,
            SPARSE_MATRIX_TYPE_DIAGONAL = 24,    /* diagonal matrix; only diagonal elements will be processed */
            SPARSE_MATRIX_TYPE_BLOCK_TRIANGULAR = 25,
            SPARSE_MATRIX_TYPE_BLOCK_DIAGONAL = 26    /* block-diagonal matrix; only diagonal blocks will be processed */
        }

        public enum sparse_fill_mode_t {
            SPARSE_FILL_MODE_LOWER = 40,           /* lower triangular part of the matrix is stored */
            SPARSE_FILL_MODE_UPPER = 41,            /* upper triangular part of the matrix is stored */
            SPARSE_FILL_MODE_FULL = 42            /* upper triangular part of the matrix is stored */
        }

        public enum sparse_diag_type_t {
            SPARSE_DIAG_NON_UNIT = 50,           /* triangular matrix with non-unit diagonal */
            SPARSE_DIAG_UNIT = 51            /* triangular matrix with unit diagonal */
        }

        public enum sparse_operation_t {
            SPARSE_OPERATION_NON_TRANSPOSE = 10,
            SPARSE_OPERATION_TRANSPOSE = 11,
            SPARSE_OPERATION_CONJUGATE_TRANSPOSE = 12
        }

        public enum sparse_index_base_t {
            ARSE_INDEX_BASE_ZERO = 0,           /* C-style */
            SPARSE_INDEX_BASE_ONE = 1            /* Fortran-style */
        }

        /// <summary>
        /// PARDISO interface
        /// </summary>
        public unsafe _dcsrilu0 ILU0 { get { return dcsrilu0; } }

        /// <summary>
        /// PARDISO interface
        /// </summary>
        public unsafe _dcsrilut ILUT { get { return dcsrilut; } }

        /// <summary>
        /// Enables performant back and forward substitution.
        /// </summary>
        public unsafe _d_trsv Substitute { get { return mkl_sparse_d_trsv; } }

        /// <summary>
        /// Internal MKL matrix format. Needed for the substitution routine.
        /// </summary>
        public unsafe _d_crcsr CreateCSRMatrix { get { return mkl_sparse_d_create_csr;  } }

        public unsafe _destroy Destroy { get { return mkl_sparse_destroy; } }

        public unsafe _optimize Optimize { get { return mkl_sparse_optimize; } }

        public unsafe _export Export { get { return mkl_sparse_d_export_csr; } }

        static string WinMKL_lp64_mangling(string nmn) {
            return "mkl_pds_lp64_" + nmn;
        }

        static string WinMKLmangling(string nmn) {
            return "mkl_pds_" + nmn;
        }


        /// <summary>
        /// converts PARDISO error code to an hopefully more-explaning error string (taken from the manual)
        /// </summary>
        public static string ILUerror2string(int error) {
            switch (error) {
                case 0: return "no error";
                case -101: return "Indicates that the routine was interrupted because of an error: the number of elements in some matrix row specified in the sparse format is equal to or less than 0.";
                case -102: return "Indicates that the routine was interrupted because the value of the computed diagonal element is less than the product of the given tolerance and the current matrix row norm, and it cannot be replaced as ipar(31)=0.";
                case -103: return "Indicates that the routine was interrupted because the element ia(i + 1) is less than or equal to the element ia(i)";
                case -104: return "Indicates that the routine was interrupted because the memory is insufficient for the internal work arrays.";
                case -105: return "Indicates that the routine was interrupted because the input value of maxfil is less than 0.";
                case -106: return "Indicates that the routine was interrupted because the size n of the input matrix is less than 0.";
                case -107: return "Indicates that the routine was interrupted because an element of the array ja is less than 1, or greater than n";
                case 101: return "The value of maxfil is greater than or equal to n. The calculation is performed with the value of maxfil set to (n-1).";
                case 102: return "The value of tol is less than 0. The calculation is performed with the value of the parameter set to (-tol)";
                case 103: return "The absolute value of tol is greater than value of dpar(31); it can result in instability of the calculation.";
                case 104: return "The value of dpar(31) is equal to 0. It can cause calculations to fail.";
                default: return "error code not specified in manual";
            }
        }

        public static string MKLstatus2string(int status) {
            return Enum.GetName(typeof(Wrapper_MKL.MKL_Status), status);
        }

        enum MKL_Status {
            SPARSE_STATUS_SUCCESS = 0,
            SPARSE_STATUS_NOT_INITIALIZED = 1,// The routine encountered an empty handle or matrix array.
            SPARSE_STATUS_ALLOC_FAILED = 2,// Internal memory allocation failed.
            SPARSE_STATUS_INVALID_VALUE = 3,// The input parameters contain an invalid value.
            SPARSE_STATUS_EXECUTION_FAILED = 4,// Execution failed.
            SPARSE_STATUS_INTERNAL_ERROR = 5,// An error in algorithm implementation occurred.
            SPARSE_STATUS_NOT_SUPPORTED = 6,// The requested operation is not supported.
        }

    }
}
