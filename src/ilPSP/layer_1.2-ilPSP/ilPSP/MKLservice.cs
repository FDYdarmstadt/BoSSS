using ilPSP.Utils;
using MPI.Wrappers.Utils;
using System;
using System.Collections.Generic;
using System.Data;
using System.Runtime.InteropServices;
using System.Text;

namespace ilPSP {


    public enum OMPbindingStrategy {
        /// <summary>
        /// <see cref="MKLservice.BindOMPthreads_1To1"/>
        /// </summary>
        OneToOne = 1,

        /// <summary>
        /// <see cref="MKLservice.BindOMPthreads_AllToAll"/>
        /// </summary>
        AllToAll = 2,

        /// <summary>
        /// <see cref="MKLservice.BindOMPthreads_CloneFromMainThread"/>
        /// </summary>
        CloneFromMain = 3

    }


    public class MKLservice : DynLibLoader {

        public MKLservice() :
            base(BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPrequesiteLibraries(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetLibname(Parallelism.OMP).Length.ForLoop<GetNameMangling>( i => DynLibLoader.Identity),
                 BLAS_LAPACK_Libstuff.GetPlatformID(Parallelism.OMP),
                 BLAS_LAPACK_Libstuff.GetPointerSizeFilter(Parallelism.OMP)) //
        { }




        public unsafe delegate void _BoSSS_set_num_threads(int* nth);

        public unsafe delegate int _BoSSS_bind_omp_threads(int NumThreads, int* CPUindices);

        public unsafe delegate int _BoSSS_set_dynamic(int boolDynThreads);

        public unsafe delegate void _BoSSS_get_dynamic(int* boolDynThreads);

        public unsafe delegate int _BoSSS_get_max_threads();


#pragma warning disable        649
        _BoSSS_set_num_threads BoSSS_set_num_threads;
        _BoSSS_bind_omp_threads BoSSS_bind_omp_threads;
        _BoSSS_set_dynamic BoSSS_set_dynamic;
        _BoSSS_get_dynamic BoSSS_get_dynamic;
        _BoSSS_get_max_threads BoSSS_get_max_threads;
#pragma warning restore 649


        public readonly static MKLservice instance = new MKLservice();

        /// <summary>
        /// Setting number of threads for MKL/OpenMP
        /// </summary>
        public static void SetNumThreads(int nth) {
            unsafe {
                instance.BoSSS_set_num_threads(&nth);
            }
        }

        //[DllImport("libiomp5md.dll")]
        //public static extern void omp_set_num_threads(int nth);

        //[DllImport("libiomp5md.dll")]
        //public static extern int omp_get_num_threads();

        /// <summary>
        /// Not really sure?
        /// </summary>
        public static int GetMaxThreads() {
            unsafe {
                return instance.BoSSS_get_max_threads();
            }
        }

        /// <summary>
        /// Setting/Getting the state of OpenMP dynamic thread allocation 
        /// (internal control variable `omp_get_dynamic`, overrifed the bahavor of OMP_DYNAMIC environment variable)
        /// </summary>
        public static bool Dynamic {
            get {
                unsafe {
                    int ret = 0;
                    instance.BoSSS_get_dynamic(&ret);
                    return ret != 0;
                }
            }
            set {
                instance.BoSSS_set_dynamic(value ? 1 : 0);
            }

        }

        public static void BindOMPthreads(int[] CPUindices, OMPbindingStrategy strat) {
            switch(strat) {
                case OMPbindingStrategy.AllToAll: BindOMPthreads_AllToAll(CPUindices); break;
                case OMPbindingStrategy.CloneFromMain: BindOMPthreads_CloneFromMainThread(CPUindices.Length); break;
                case OMPbindingStrategy.OneToOne: BindOMPthreads_1To1(CPUindices); break;
                default: throw new NotImplementedException();
            }
        }


        public static void BindOMPthreads_1To1(int[] CPUindices) {
            int ret;
            int NumCpus = CPUindices.Length;
            unsafe {
                int* __CPUindices = stackalloc int[NumCpus];
                for (int i = 0; i < NumCpus; i++)
                    __CPUindices[i] = CPUindices[i]; // we need to make a copy to write a proper error message in the fail case.
                ret = instance.BoSSS_bind_omp_threads(NumCpus, __CPUindices); // `__CPUindices` is used for input and output;

                for (int i = 0; i < NumCpus; i++) {
                    if(__CPUindices[i] != 0)
                        Console.Error.WriteLine($"Error binding OMP thread #{i} to CPU #{CPUindices[i]}: thread #{i} return code {__CPUindices[i]}");
                }

                if(ret != 0) {
                    Console.Error.WriteLine($"BoSSS_bind_omp_threads returned {ret}.");
                }

            }
        }

        public static void BindOMPthreads_AllToAll(int[] CPUindices) {
            int ret;
            int NumThreads = CPUindices.Length;
            //int[] _CPUindices = new int[(NumThreads + 1)*(NumThreads)];
            unsafe {
                //fixed (int* __CPUindices = _CPUindices) {
                int* __CPUindices = stackalloc int[(NumThreads + 1)*(NumThreads)];

                for (int i = 0; i < NumThreads; i++) {
                    __CPUindices[i*(NumThreads + 1)] = -(NumThreads + 1);

                    for (int j = 0; j < NumThreads; j++) {
                        __CPUindices[i*(NumThreads + 1) + j + 1] = CPUindices[j];
                    }
                }


                ret = instance.BoSSS_bind_omp_threads(NumThreads, __CPUindices); // `__CPUindices` is used for input and output;

                for (int i = 0; i < NumThreads; i++) {
                    if (__CPUindices[i*(NumThreads + 1) + 1] != 0)
                        Console.Error.WriteLine($"Error binding OMP thread #{i} to CPU #{CPUindices[i]}: thread #{i} return code {__CPUindices[i*(NumThreads + 1) + 1]}");
                }

                if (ret != 0) {
                    Console.Error.WriteLine($"BoSSS_bind_omp_threads returned {ret}.");
                }


            }
        }

        public static void BindOMPthreads_CloneFromMainThread(int NumThreads) {
            int ret;
            unsafe {
                int* __CPUindices = stackalloc int[NumThreads];

                for (int i = 0; i < NumThreads; i++) {
                    __CPUindices[i] = -1222333444;
                }


                ret = instance.BoSSS_bind_omp_threads(NumThreads, __CPUindices); // `__CPUindices` is used for input and output;

                for (int i = 0; i < NumThreads; i++) {
                    if (__CPUindices[i] != 0)
                        Console.Error.WriteLine($"Error binding OMP thread #{i} equally to main thread: thread #{i} return code {__CPUindices[i]}");
                }

                if (ret != 0) {
                    Console.Error.WriteLine($"BoSSS_bind_omp_threads returned {ret}.");
                }

            }
        }

    }
}
