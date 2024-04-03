using ilPSP.Utils;
using MPI.Wrappers.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Text;

namespace ilPSP {
    public class MKLservice : DynLibLoader {

        public MKLservice() :
            base(BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPrequesiteLibraries(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP).Length.ForLoop<GetNameMangling>( i => DynLibLoader.Identity),
                 BLAS_LAPACK_Libstuff.GetPlatformID(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPointerSizeFilter(Parallelism.OMP)) //
        { }



        /// <summary>
        /// 
        /// </summary>
        public unsafe delegate void _BoSSS_set_num_threads(int* nth);

        public unsafe delegate int _BoSSS_bind_omp_threads(int NumThreads, int* CPUindices);


#pragma warning disable        649
        _BoSSS_set_num_threads BoSSS_set_num_threads;
        _BoSSS_bind_omp_threads BoSSS_bind_omp_threads;
#pragma warning restore 649


        public readonly static MKLservice instance = new MKLservice();


        public static void SetNumThreads(int nth) {
            if (nth > 1000000) {
                unsafe {
                    instance.BoSSS_set_num_threads(&nth);
                }
            }
        }

        public static void BindOMPthreads(int[] CPUindices) {
            unsafe {
                fixed (int* cores = CPUindices) {
                    int NumCpus = CPUindices.Length;
                    instance.BoSSS_bind_omp_threads(NumCpus, cores);
                }
            }
        }

    }
}
