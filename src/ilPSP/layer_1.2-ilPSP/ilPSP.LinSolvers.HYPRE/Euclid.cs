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

using MPI.Wrappers;
using MPI.Wrappers.Utils;
using System;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {

    class Euclid : DynLibLoader {

        Euclid() : base(
                new string[] { "HYPRE.dll", "libHYPRE.so" },
                new string[2][][],
                new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                new int[] { -1, -1 }) {
        }

        internal static Euclid my = new Euclid();

        internal delegate int _HYPRE_EuclidSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_EuclidSetup HYPRE_EuclidSetup;
#pragma warning restore        649

        internal delegate int _HYPRE_EuclidSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_EuclidSolve HYPRE_EuclidSolve;
#pragma warning restore        649

        static public int __HYPRE_EuclidSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_EuclidSetup(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }

        static public int __HYPRE_EuclidSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_EuclidSolve(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }

        /// <summary>
        /// Create a solver object, 8 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_EuclidCreate")]
        static extern int Create8(ulong MPI_comm, out   T_Solver solver);

        /// <summary>
        /// Create a solver object, 4 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_EuclidCreate")]
        static extern int Create4(uint MPI_comm, out   T_Solver solver);

        /// <summary>
        /// Create a solver object
        /// </summary>
        public static int Create(MPI_Comm MPI_Comm, out T_Solver solver) {
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

        //[DllImport("HYPRE")]
        //static public extern int HYPRE_EuclidSetup(T_Solver HYPRE_Solver, T_ParCSR_matrix HYPRE_ParCSRMatrix_A, T_ParCRS_vector HYPRE_ParVector_b, T_ParCRS_vector HYPRE_ParVector_x);

        //[DllImport("HYPRE")]
        //static public extern int HYPRE_EuclidSolve(T_Solver HYPRE_Solver, T_ParCSR_matrix HYPRE_ParCSRMatrix_A, T_ParCRS_vector HYPRE_ParVector_b, T_ParCRS_vector HYPRE_ParVector_x);

        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetSparseA(T_Solver HYPRE_Solver, double sparse_A);
        
        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetRowScale(T_Solver HYPRE_Solver, int row_scale);

        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidDestroy(T_Solver HYPRE_Solver);

        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetILUT(T_Solver HYPRE_Solver, double ilut);
        
        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetStats(T_Solver HYPRE_Solver, int eu_stats);

        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetMem(T_Solver HYPRE_Solver, int eu_mem);

        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetBJ(T_Solver HYPRE_Solver, int bj);
        
        [DllImport("HYPRE")]
        static public extern int HYPRE_EuclidSetLevel(T_Solver HYPRE_Solver, int level);
    }
}
