using ilPSP.Utils;
using MPI.Wrappers.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Text;

namespace ilPSP {
    internal class MKLservice : DynLibLoader {

        public MKLservice() :
            base(BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPrequesiteLibraries(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP).Length.ForLoop<GetNameMangling>( i => DynLibLoader.Identity),
                 BLAS_LAPACK_Libstuff.GetPlatformID(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPointerSizeFilter(Parallelism.OMP)) //
        { }



        /// <summary>
        /// FORTRAN-style LAPACK, matrices in FORTRAN order;
        /// </summary>
        public unsafe delegate void _BoSSS_set_num_threads(int* nth);


#pragma warning disable        649
        _BoSSS_set_num_threads BoSSS_set_num_threads;
#pragma warning restore 649


        public readonly static MKLservice instance = new MKLservice();


        public static void SetNumThreads(int nth) {
            unsafe {
                instance.BoSSS_set_num_threads(&nth);
            }
        }
    }
}
