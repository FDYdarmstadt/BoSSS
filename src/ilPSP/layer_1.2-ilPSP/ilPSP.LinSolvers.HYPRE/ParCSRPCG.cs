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
using System.Runtime.InteropServices;
using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {
   

    class ParCSRPCG {

        /// <summary>
        /// Create a solver object, 8 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGCreate")]
        static extern int Create8(ulong MPI_comm, out   T_Solver solver);
        
        /// <summary>
        /// Create a solver object, 4 byte MPI communicator handle
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGCreate")]
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

        /// <summary>
        /// Destroy a solver object
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGDestroy")]
        public static extern int Destroy(T_Solver solver);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetup")]
        public static extern int Setup(T_Solver solver, IntPtr ParCSRMatrix_A, IntPtr ParVector_b, IntPtr ParVector_x);

        /// <summary> 
        /// Solve the system
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSolve")]
        public static extern int Solve(T_Solver solver, IntPtr ParCSRMatrix_A, IntPtr ParVector_b, IntPtr ParVector_x);

        /// <summary>
        /// (Optional) Set the relative convergence tolerance
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetTol")]
        public static extern int SetTol(T_Solver solver, double tol);

        /// <summary>
        /// (Optional) Set the absolute convergence tolerance (default is 0)
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetAbsoluteTol")]
        public static extern int SetAbsoluteTol(T_Solver solver, double tol);

        /// <summary>
        /// (Optional) Set maximum number of iterations
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetMaxIter")]
        public static extern int SetMaxIter(T_Solver solver, int max_iter);

        /// <summary>
        /// (Optional) Use the two-norm in stopping criteria
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetTwoNorm")]
        public static extern int SetTwoNorm(T_Solver solver, int two_norm);

        /// <summary>
        /// (Optional) Additionally require that the relative difference in successive iterates be small
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetRelChange")]
        public static extern int SetRelChange(T_Solver solver, int rel_change);


        /// <summary>
        /// (Optional) Set the preconditioner to use
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGSetPrecond")]
        public static extern int SetPrecond(T_Solver solver, IntPtr PtrToParSolverFcn_precond, IntPtr PtrToParSolverFcn_precond_setup, IntPtr precond_solver);

        /// <summary>
        /// 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGGetPrecond")]
        public static extern int GetPrecond(T_Solver solver, out IntPtr precond_data);

        /// <summary>
        /// Return the number of iterations taken
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGGetNumIterations")]
        public static extern int GetNumIterations(T_Solver solver, out int num_iterations);

        /// <summary>
        /// Return the norm of the final relative residual
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRPCGGetFinalRelativeResidualNorm")]
        public static extern int GetFinalRelativeResidualNorm(T_Solver solver, out double norm);

        /// <summary>
        /// Setup routine for diagonal preconditioning
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRDiagScaleSetup")]
        public static extern int ParCSRDiagScaleSetup(T_Solver solver, IntPtr ParCSRMatrix_A, IntPtr ParVector_y, IntPtr ParVector_x);

        /// <summary>
        /// Solve routine for diagonal preconditioning
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_ParCSRDiagScale")]
        public static extern int ParCSRDiagScale(T_Solver solver, IntPtr ParCSRMatrix_HA, IntPtr ParVector_Hy, IntPtr ParVector_Hx);
    }
}
