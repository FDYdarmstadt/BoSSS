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
using System.Security;
using MPI.Wrappers.Utils;

namespace ilPSP.LinSolvers.UMFPACK {

    /// <summary>
    /// Native wrapper for the double/int UMFPACK interface.
    /// </summary>
    [SuppressUnmanagedCodeSecurity]
    public class UnsafeUMFPACK : DynLibLoader {

        static Parallelism SelectPar(Parallelism par) {
            if(ilPSP.Environment.MaxNumOpenMPthreads <= 1 && par == Parallelism.OMP)
                par = Parallelism.SEQ;
            if(ilPSP.Environment.OpenMPenabled == false && par == Parallelism.OMP)
                par = Parallelism.SEQ;

            return par;
        }

        static (string[] _LibNames, string[][][] PrequesiteLibraries, GetNameMangling[] NameMangling, PlatformID[] OsFilter, int[] PointerSizeFilter) SelectLibrary(Parallelism par) {
            switch(SelectPar(par)) {
                case Parallelism.SEQ:
                return (new string[] { "UMFPACK_seq.dll", "libBoSSSnative_seq.so" },
                    new string[2][][],
                    new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                    new[] { PlatformID.Win32NT, PlatformID.Unix },
                    new int[] { -1, -1 });

                case Parallelism.OMP:
                return (new string[] { "UMFPACK_omp.dll", "libBoSSSnative_omp.so" },
                    new string[2][][],
                    new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                    new[] { PlatformID.Win32NT, PlatformID.Unix },
                    new int[] { -1, -1 });

                default:
                throw new NotSupportedException($"Unsupported level of parallelism {par} for UMFPACK.");
            }
        }

        /// <summary>
        /// ctor
        /// </summary>
        public UnsafeUMFPACK(Parallelism par) : base(SelectLibrary(par)) {
        }

        /// <summary>
        /// UMFPACK symbolic factorization.
        /// </summary>
        public unsafe delegate int _umfpack_di_symbolic(int n_row, int n_col, int* Ap, int* Ai, double* Ax, void** Symbolic, double* Control, double* Info);

        /// <summary>
        /// UMFPACK numeric factorization.
        /// </summary>
        public unsafe delegate int _umfpack_di_numeric(int* Ap, int* Ai, double* Ax, void* Symbolic, void** Numeric, double* Control, double* Info);

        /// <summary>
        /// UMFPACK solve.
        /// </summary>
        public unsafe delegate int _umfpack_di_solve(int sys, int* Ap, int* Ai, double* Ax, double* X, double* B, void* Numeric, double* Control, double* Info);

        /// <summary>
        /// Frees UMFPACK symbolic storage.
        /// </summary>
        public unsafe delegate void _umfpack_di_free_symbolic(void** Symbolic);

        /// <summary>
        /// Frees UMFPACK numeric storage.
        /// </summary>
        public unsafe delegate void _umfpack_di_free_numeric(void** Numeric);

#pragma warning disable 649
        _umfpack_di_symbolic umfpack_di_symbolic;
        _umfpack_di_numeric umfpack_di_numeric;
        _umfpack_di_solve umfpack_di_solve;
        _umfpack_di_free_symbolic umfpack_di_free_symbolic;
        _umfpack_di_free_numeric umfpack_di_free_numeric;
#pragma warning restore 649

        /// <summary>
        /// UMFPACK symbolic factorization.
        /// </summary>
        public unsafe _umfpack_di_symbolic Symbolic { get { return umfpack_di_symbolic; } }

        /// <summary>
        /// UMFPACK numeric factorization.
        /// </summary>
        public unsafe _umfpack_di_numeric Numeric { get { return umfpack_di_numeric; } }

        /// <summary>
        /// UMFPACK solve.
        /// </summary>
        public unsafe _umfpack_di_solve Solve { get { return umfpack_di_solve; } }

        /// <summary>
        /// Frees UMFPACK symbolic storage.
        /// </summary>
        public unsafe _umfpack_di_free_symbolic FreeSymbolic { get { return umfpack_di_free_symbolic; } }

        /// <summary>
        /// Frees UMFPACK numeric storage.
        /// </summary>
        public unsafe _umfpack_di_free_numeric FreeNumeric { get { return umfpack_di_free_numeric; } }

        /// <summary>
        /// Converts UMFPACK status codes to text.
        /// </summary>
        public string UMFPACKerror2string(int error) {
            switch(error) {
                case 0: return "no error";
                case 1: return "matrix is singular";
                case -1: return "out of memory";
                case -3: return "invalid matrix";
                case -4: return "argument missing";
                case -5: return "invalid Numeric object";
                case -6: return "invalid Symbolic object";
                case -8: return "invalid permutation";
                case -11: return "different pattern";
                case -13: return "invalid system";
                case -15: return "invalid scale";
                default: return "error code not specified in UMFPACK manual";
            }
        }
    }
}
