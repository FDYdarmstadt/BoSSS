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

namespace ilPSP.LinSolvers.HYPRE.Wrappers {

    /// <summary>
    /// Interface to  HYPRE_PCG functions;
    /// </summary>
    class PCG {
        
        /// <summary>
        /// Prepare to solve the system 
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetup")]
        public static extern int Setup(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);

        /// <summary>
        /// Solve the system
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSolve")]
        public static extern int Solve(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);

        /// <summary>
        /// (Optional) Set the relative convergence tolerance
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetTol")]
        public static extern int SetTol(T_Solver solver, double tol);

        /// <summary>
        /// (Optional) Set the absolute convergence tolerance (default is 0)
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetAbsoluteTol")]
        public static extern int SetAbsoluteTol(T_Solver solver, double a_tol);

        /// <summary>
        /// (Optional) Get the absolute convergence tolerance (default is 0)
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetAbsoluteTol")]
        public static extern int GetAbsoluteTol(T_Solver solver, out double a_tol);

        /// <summary>
        /// (Optional) Set maximum number of iterations
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetMaxIter")]
        public static extern int SetMaxIter(T_Solver solver, int max_iter);

        /// <summary>
        /// (Optional) Use the two-norm in stopping criteria
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetTwoNorm")]
        public static extern int SetTwoNorm(T_Solver solver, int two_norm);

        /// <summary>
        /// (Optional) Additionally require that the relative difference in successive iterates be small
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetRelChange")]
        public static extern int SetRelChange(T_Solver solver, int rel_change);

        /// <summary>
        /// (Optional) Set the preconditioner to use
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetPrecond")]
        public static extern int SetPrecond(T_Solver solver, 
            IntPtr PtrToSolverFcn_precond_solve, IntPtr PtrToSolverFcn_precond_setup, T_Solver Solver_precond_solver);

        ///// <summary>
        ///// (Optional) Set the preconditioner to use
        ///// </summary>
        //[DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetPrecond_Mod")]
        //public static extern int SetPrecond(T_Solver solver, IntPtr name, T_Solver Solver_precond_solver);


        /// <summary>
        /// (Optional) Set the amount of logging to do
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetLogging")]
        public static extern int SetLogging(T_Solver solver, int logging);

        ///// <summary>
        ///// (Optional) Set the amount of printing to do to the screen
        ///// </summary>
        //[DllImport("Platform_Native", EntryPoint = "HYPRE_PCGSetPrintLevel")]
        //public static extern int SetPrintLevel(IntPtr solver, int level);

        /// <summary>
        /// Return the number of iterations taken
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetNumIterations")]
        public static extern int GetNumIterations(T_Solver solver, out int num_iterations);

        /// <summary>
        /// Return the norm of the final relative residual
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetFinalRelativeResidualNorm")]
        public static extern int GetFinalRelativeResidualNorm(T_Solver solver, out double norm);

        ///// <summary>
        ///// Return the residual
        ///// </summary>
        //[DllImport("Platform_Native", EntryPoint = "HYPRE_PCGGetResidual")]
        //public static extern int GetResidual(IntPtr solver, void** residual);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetTol")]
        public static extern int GetTol(T_Solver solver, out double tol);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetMaxIter")]
        public static extern int GetMaxIter(T_Solver solver, out int max_iter);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetTwoNorm")]
        public static extern int GetTwoNorm(T_Solver solver, out int two_norm);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetRelChange")]
        public static extern int GetRelChange(T_Solver solver, out int rel_change);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetPrecond")]
        public static extern int GetPrecond(T_Solver solver, out IntPtr precond_data_ptr);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetLogging")]
        public static extern int GetLogging(T_Solver solver, out int level);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetPrintLevel")]
        public static extern int GetPrintLevel(T_Solver solver, out int level);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGSetPrintLevel")]
        public static extern int SetPrintLevel(T_Solver solver, int level);

        /// <summary> </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_PCGGetConverged")]
        public static extern int GetConverged(T_Solver solver, out int converged);


    }
}
