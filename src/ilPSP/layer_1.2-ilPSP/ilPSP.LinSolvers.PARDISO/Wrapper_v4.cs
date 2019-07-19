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
   
    
    /// <summary>
    /// wrapper for loading the PARDISO Version 4
    /// </summary>
    /// <remarks>
    /// <b>IMPORTANT: Licensing issues:</b><br/>
    /// PARDISO does not ship with free license, neither source nor 
    /// binaries compiled from it can be shipped with this software;<br/>
    /// PARDISO is ether distributed with the INTEL MKL library, or it may be downloaded
    /// from http://www.pardiso-project.org/;
    /// </remarks>
    public class Wrapper_v4 : DynLibLoader {

        /// <summary>
        /// ctor
        /// </summary>
        public Wrapper_v4()
            : base(
                new string[] { "libpardiso412-WIN-X86-64.dll", "libpardiso412-WIN-X86.dll" },
                new string[2][][],
                new GetNameMangling[] { Identity, Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Win32NT },
                new int[] { 8, 4 }) 
        {
            System.Environment.SetEnvironmentVariable("PARDISOLICMESSAGE", "1"); // surpress license message
        }

        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe delegate void _pardisoinit(void* pt, 
            int* mtype, int* solver, int* iparm, double* dparm, int* error);


        /// <summary>
        /// PARDISO interface, see PARDISO documentation
        /// </summary>
        public unsafe delegate void _pardiso(void* pt, int* maxfct, int* mnum, int* mtype, 
                                             int* phase, int* n, 
                                             void* a, int* ia, int* ja,
                                             int* perm, int* nrhs, int* iparm,
                                             int* msglvl, double* b, double* x,
                                             int* error, double* dparm);


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
                case 0: return "No error.";
                case -1: return "Input inconsistent.";
                case -2: return "Not enough memory.";
                case -3: return "Reordering problem.";
                case -4: return "Zero pivot, numerical fact. or iterative refinement problem.";
                case -5: return "Unclassified (internal) error.";
                case -6: return "Preordering failed (matrix types 11, 13 only).";
                case -7: return "Diagonal matrix problem.";
                case -8: return "32-bit integer overflow problem.";
                case -10: return "No license file pardiso.lic found.";
                case -11: return "License is expired.";
                case -12: return "Wrong username or hostname.";
                case -100: return "Reached maximum number of Krylov-subspace iteration in iterative solver.";
                case -101: return "No sucient convergence in Krylov-subspace iteration within 25 iterations.";
                case -102: return "Error in Krylov-subspace iteration.";
                case -103: return "Break-Down in Krylov-subspace iteration";
                default: return "error code not specified in manual";
            }
        }
    }
}
