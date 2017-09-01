using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.ALMOND {
    static class Wrappers {

        //
        // Basic initialization, step 1/to be called exactly once in the lifetime of an application.
        // Initializes the occa device used by ALMOND.
        //
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_Init([MarshalAs(UnmanagedType.AnsiBStr)] string OccaMode, ref int occaPlatId, ref int occaDevId, out int ierr);

        //
        // returns the handle of the occa device which is used for ALMOND
	    // Args:
	    // occaDevptr: output; the occa device handle; all ALMOND objects refer to that device.
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_GetOccaDev(out IntPtr occaDevptr, out int ierr);


        //
        // Basic initialization, step 2a/to be called exactly once in the lifetime of an application.
        // compiles the ALMOND kernels for single precision.
        // Args:
        // almond_dir: optional specification of the ALMOND directory;
        //             occa kernels are expected to be in the subdirectory $(almond_dir)/kernels
        //             If NULL, ignored.
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_buildAlmondKernels_f([MarshalAs(UnmanagedType.AnsiBStr)] string almond_dir, out int ierr);

        //
        // Basic initialization, step 2b/to be called exactly once in the lifetime of an application.
        // compiles the ALMOND kernels for double precision.
        // Args:
        // almond_dir: optional specification of the ALMOND directory;
        //             occa kernels are expected to be in the subdirectory $(almond_dir)/kernels
        //             If NULL, ignored.
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_buildAlmondKernels_d([MarshalAs(UnmanagedType.AnsiBStr)] string almond_dir, out int ierr);

        //
        // create a (double precision) matrix object
        // Args:
        // m: no of rows
        // n: no of columns
        // nnz: number of non-zeros
        // ia: array of length 'nnz', row indices
        // ja: array of length 'nnz', column indices
        // aa: array of length 'nnz', non-zero matrix entries
        // _A: output: handle to the almond matrix object
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_DefineMatrix_d(ref int m, ref int n, ref int nnz, int[] ia, int[] ja, double[] aa, out IntPtr _A, out int error);

        //
        // delete a (double precision) matrix object
        //
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_DeleteMatrix_d(ref IntPtr _A, out int error);

        //
        // create a double-precision solver object
        // Args:
        // _solver: output: handle for the created solver object.
        //
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_CreateSolver_d(ref IntPtr _A, out IntPtr _solver, out int ierr);

        //
        // delete a double-precision solver object
        //
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_DeleteSolver_d(ref IntPtr _solver, out int ierr);

        //
        // perform one run of the double-precision solver
        //
        [DllImport("ALMOND_C")]
        public static extern void ALMOND_UseSolver_d(ref int Norows, ref int _krylovtype, ref int maxIt, ref double tol, ref IntPtr _solver, double[] X, double[] RHS, out int iter, out int ierr);

    }
}
