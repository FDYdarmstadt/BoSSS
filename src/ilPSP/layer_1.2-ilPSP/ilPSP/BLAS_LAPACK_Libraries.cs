using System;
using MPI.Wrappers.Utils;

namespace ilPSP.Utils {
    /// <summary>Library descriptors shared by the BLAS and LAPACK imports.</summary>
    internal static class BLAS_LAPACK_Libraries {
        internal static (string[], string[][][], DynLibLoader.GetNameMangling[], PlatformID[], int[]) Get(BLASProvider provider, Parallelism parallelism) {
            switch(provider) {
                case BLASProvider.MKL:
                    return (BLAS_LAPACK_IntelMKL_Libstuff.GetLibname(parallelism),
                        BLAS_LAPACK_IntelMKL_Libstuff.GetPrerequesiteLibraries(parallelism),
                        BLAS_LAPACK_IntelMKL_Libstuff.GetGetNameMangling(parallelism),
                        BLAS_LAPACK_IntelMKL_Libstuff.GetPlatformID(parallelism),
                        BLAS_LAPACK_IntelMKL_Libstuff.GetPointerSizeFilter(parallelism));
                case BLASProvider.OpenBLAS:
                    return BLAS_OpenBLAS_Libstuff.GetLibstuff();
                default:
                    throw new ArgumentOutOfRangeException(nameof(provider));
            }
        }
    }
}
