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

namespace ilPSP.LinSolvers.PARDISO {

    
//    /// <summary>
//    /// experimental stuff
//    /// </summary>
//    class MKL_threadCtrl : DynLibLoader {

//        public MKL_threadCtrl()
//            : base(
//                new string[] { "mkl_intel_thread.dll", "libmkl_intel_thread.so" },
//                new GetNameMangling[] { DynLibLoader.Identity, LinuxMapping },
//                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
//                new int[] { -1, -1 }) { }


//        public delegate void _set_num_threads(int n);

//#pragma warning disable        649
//        public _set_num_threads mkl_serv_mkl_set_num_threads;
//#pragma warning restore        649

//        public override void Dispose() {
//            //base.Dispose();
//        }


//        static string LinuxMapping(string s) {
//            return ("mkl_serv" + s);
//        }
//    }
    

    /// <summary>
    /// wrapper for loading the PARDISO solver from the Intel MKL libraries
    /// </summary>
    /// <remarks>
    /// <b>IMPORTANT: Licencing issues:</b><br/>
    /// PARDISO does not ship with free licence, neither source nor 
    /// binaries compiled from it can be shiped with this software;<br/>
    /// PARDISO is ether distributed with the INTEL MKL library, or it may be downloaded
    /// from http://www.pardiso-project.org/;
    /// </remarks>
    public class Wrapper_MKL : DynLibLoader {

        /// <summary>
        /// ctor
        /// </summary>
        public Wrapper_MKL()
            : base(
                new string[] { "PARDISO.dll" },
                new string[1][][],
                new GetNameMangling[] { SmallLetters_TrailingUnderscore},
                new PlatformID[] { PlatformID.Win32NT },
                new int[] { -1 }) { }

        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe delegate int _pardisoinit(void* pt, int* mtype, int* iparm);

        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe delegate int _pardiso(void* pt, int* maxfct, int* mnum, int* mtype,
                                            int* phase, int* n,
                                            void* a, int* ia, int* ja,
                                            int* perm, int* nrhs, int* iparm,
                                            int* msglvl, double* b, double* x,
                                            int* error);

#pragma warning disable        649
        _pardisoinit pardisoinit;
        _pardiso pardiso;
#pragma warning restore        649


        /// <summary>
        /// PARDISO interface
        /// </summary>
        public unsafe _pardisoinit PARDISOINIT { get { return pardisoinit; } }

        /// <summary>
        /// PARDISO interface
        /// </summary>
        public unsafe _pardiso PARDISO { get { return pardiso; } }


        static string WinMKL_lp64_mangling(string nmn) {
            return "mkl_pds_lp64_" + nmn;
        }

        static string WinMKLmangling(string nmn) {
            return "mkl_pds_" + nmn;
        }


        /// <summary>
        /// converts PARDISO error code to an hopefully more-explaning error string (taken from the manual)
        /// </summary>
        public string PARDISOerror2string(int error) {
            switch (error) {
                case 0: return "no error";
                case -1: return "input inconsistent";
                case -2: return "not enough memory";
                case -3: return "reordering problem";
                case -4: return "zero pivot, numerical factorization or iterative refinement problem";
                case -5: return "unclassified (internal) error";
                case -6: return "preordering failed (matrix types 11, 13 only)";
                case -7: return "diagonal matrix problem";
                case -8: return "32-bit integer overflow problem";
                case -9: return "not enough memory for OOC";
                case -10: return "problems with opening OOC temporary files";
                case -11: return "read/write problems with the OOC data file";
                default: return "error code not specified in manual";
            }
        }
    }
}
