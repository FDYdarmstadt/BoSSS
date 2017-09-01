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
using System.Runtime.InteropServices;
using MPI.Wrappers;
using MPI.Wrappers.Utils;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {
    
    /// <summary>
    /// wrappers for HYPRE ParaSails preconditioner
    /// </summary>
    class ParaSails : DynLibLoader {

        /// <summary>
        /// ctor
        /// </summary>
        ParaSails() : base(
                new string[] { "HYPRE.dll", "libHYPRE.so" },
                new string[2][][],
                new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                new int[] { -1, -1 }) {
        }

        internal static ParaSails my = new ParaSails();

        internal delegate int _HYPRE_ParaSailsSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_ParaSailsSetup HYPRE_ParaSailsSetup;
#pragma warning restore        649

        internal delegate int _HYPRE_ParaSailsSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_ParaSailsSolve HYPRE_ParaSailsSolve;
#pragma warning restore        649

        static public int __HYPRE_ParaSailsSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_ParaSailsSolve(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }

        static public int __HYPRE_ParaSailsSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_ParaSailsSetup(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }
                
        /// <summary>
        /// Create a solver object, 8 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParaSailsCreate")]
        static extern int Create8(ulong MPI_comm, out   T_Solver solver);
        
        /// <summary>
        /// Create a solver object, 4 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParaSailsCreate")]
        static extern int Create4(uint MPI_comm, out   T_Solver solver);

        /// <summary>
        /// Create a solver object
        /// </summary>
        public static int HYPRE_ParaSailsCreate(MPI_Comm MPI_Comm, out T_Solver solver) {
            ulong _com8;
            uint _com4;
            // we need to convert the MPI comm in ilPSP (which is a FORTRAN MPI comm)
            // to a C-MPI comm: can be either 4 or 8 bytes!
            int sz = csMPI.Raw.MPI_Comm_f2c(MPI_Comm, out _com4, out _com8);
            switch (sz) {
                case 4: return Create4(_com4, out solver);
                case 8: return Create8(_com8, out solver);
                default: throw new NotImplementedException();
            }
        }

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsDestroy(T_Solver solver);
        
        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsSetThresh(T_Solver solver, double thresh);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsGetThresh(T_Solver solver, out double thresh);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsSetSym(T_Solver solver, int sym);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsGetSym(T_Solver solver, out int sym);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsSetLoadbal(T_Solver solver, double loadbal);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsGetLoadbal(T_Solver solver, double loadbal);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsGetFilter(T_Solver solver, out double filter);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsSetFilter(T_Solver solver, double filter);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsGetLogging(T_Solver solver, out int logging);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsSetLogging(T_Solver solver, int logging);

        [DllImport("HYPRE")]
        public extern static int HYPRE_ParaSailsBuildIJMatrix(T_Solver solver, out T_IJMatrix pij_A);

        [DllImport("HYPRE")]
        public static extern int HYPRE_ParaSailsSetNlevels(T_Solver solver, int nlevels);

        [DllImport("HYPRE")]
        public static extern int HYPRE_ParaSailsGetNlevels(T_Solver solver, out int nlevels);

        [DllImport("HYPRE")]
        public static extern int HYPRE_ParaSailsSetReuse(T_Solver solver, int reuse);

        [DllImport("HYPRE")]
        public static extern int HYPRE_ParaSailsGetReuse(T_Solver solver, out int reuse);

    }
}
