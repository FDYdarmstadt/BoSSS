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

namespace ilPSP.LinSolvers.HYPRE.Wrappers
{

    class GMRES
    {
        /// <summary>
        /// Prepare to solve the system.  The coefficient data in {\tt b} and {\tt x} is
        /// ignored here, but information about the layout of the data may be used.
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetup")]
        public static extern int Setup(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);

        ///<summary>
        /// Solve the system 
        ///</summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSolve")]
        public static extern int Solve(T_Solver solver, T_ParCSR_matrix Matrix_A, T_ParCRS_vector Vector_b, T_ParCRS_vector Vector_x);

        /// <summary>
        /// (Optional) Set the preconditioner to use
        /// </summary>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetPrecond")]
        public static extern int SetPrecond(T_Solver solver,
            IntPtr PtrToSolverFcn_precond_solve, IntPtr PtrToSolverFcn_precond_setup, T_Solver Solver_precond_solver);


        /// <summary>
        /// (Optional) Set the absolute convergence tolerance (default is 0). 
        /// If one desires
        /// the convergence test to check the absolute convergence tolerance {\it only}, then
        /// set the relative convergence tolerance to 0.0.  (The convergence test is 
        /// $\|r\| \leq$ max(relative$\_$tolerance$\ast \|b\|$, absolute$\_$tolerance).)
        /// </summary>
        /// <param name="solver">Pointer to a HYPRE_GMRES Solver</param>
        /// <param name="p">absolute tolerance value (default is 0.0)</param>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetAbsoluteTol")]
        public static extern int SetAbsoluteTolerance(T_Solver solver, double p);

        /// <summary>
        /// Set relative convergence tolerance
        /// </summary>
        /// <param name="solver">Pointer to HYPRE_GMRES Solver</param>
        /// <param name="p">relative tolerance value</param>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetTol")]
        public static extern int SetRelativeTolerance(T_Solver solver, double p);

        /// <summary>
        ///  Get relative convergence tolerance
        /// </summary>
        /// <param name="solver">Pointer to HYPRE_GMRES Solver</param>
        /// <param name="tolerance">placeholder for return value</param>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetTol")]
        public static extern int GetRelativeTolerance(T_Solver solver, out  double tolerance);

        /// <summary>
        ///  Get absolute convergence tolerance
        /// </summary>
        /// <param name="solver">Pointer to HYPRE_GMRES Solver</param>
        /// <param name="tolerance">placeholder for return value</param>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetAbsoluteTol")]
        public static extern int GetAbsoluteTolerance(T_Solver solver, out  double tolerance);

        /// <summary>
        ///  Get maximal number of iteration 
        /// </summary>
        /// <param name="solver">Pointer to HYPRE_GMRES Solver</param>
        /// <param name="_MaxIterations">placeholder for return value</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetMaxIter")]
        public static extern int GetMaxIterations(T_Solver solver, out int _MaxIterations);

        /// <summary>
        /// (Optional)Set maximal number of iterations
        /// </summary>
        /// <param name="solver">Pointer to GMRESsolver</param>
        /// <param name="num_iterations">Number of iteration to set</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetMaxIter")]
        public static extern int SetMaxIterations(T_Solver solver, int num_iterations);

        /// <summary>
        ///  Get minimal number of iteration 
        /// </summary>
        /// <param name="solver">Pointer to HYPRE_GMRES Solver</param>
        /// <param name="num_iterations">placeholder for return value</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetMinIter")]
        public static extern int GetMinIterations(T_Solver solver, out int num_iterations);

        /// <summary>
        /// (Optional)Set minimal number of iterations
        /// </summary>
        /// <param name="solver">Pointer to GMRESsolver</param>
        /// <param name="num_iterations">Number of iteration to set</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetMinIter")]
        public static extern int SetMinIterations(T_Solver solver, int num_iterations);
        /// <summary>
        /// (Optional) Set maximum size of Krylov  space
        /// </summary>
        /// <param name="solver">Pointer to GMRESsolver</param>
        /// <param name="k_Dim">Dimensions to set</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetKDim")]
        public static extern int SetKrylovSpaceDim(T_Solver solver, int k_Dim);

        /// <summary>
        /// (Optional) Get maximum size of Krylov  space
        /// </summary>
        /// <param name="solver">Pointer to GMRESsolver</param>
        /// <param name="k_Dim">pointer to result spaceholder</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetKDim")]
        public static extern int GetKrylovSpaceDim(T_Solver solver, out int k_Dim);

        /// <summary>
        ///  Return num of iterations taken
        /// </summary>
        /// <param name="solver">Pointer to GMRESSolver</param>
       /// <param name="num_iterations">pointer to result spaceholder</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetNumIterations")]
        public static extern int GetNumOfIterations(T_Solver solver, out int num_iterations);

        /// <summary>
        /// (Optional) Set the amount of printing to do to the screen
        /// </summary>
        /// <param name="solver">Pointer to GMRESSolver</param>
        /// <param name="verbosity_lvl">Level of verbosity</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetPrintLevel")]
        public static extern int SetPrintLevel(T_Solver solver, int verbosity_lvl);

        /// <summary>
        ///  (Optional) Additionally require that the relative difference in succesive iterates be small
        /// </summary>
        /// <param name="solver">Pointer to GMRESSolver</param>
        /// <param name="rel_change">level of changes</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESSetRelChange")]
        public static extern int SetRelChange(T_Solver solver, int rel_change);

        /// <summary>
        ///  (Optional) Additionally require that the relative difference in succesive iterates be small
        /// </summary>
        /// <param name="solver">Pointer to GMRESSolver</param>
        /// <param name="rel_change">pointer to result spaceholder</param>
        /// <returns></returns>
        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetRelChange")]
        public static extern int GetRelChange(T_Solver solver, out int rel_change);

        [DllImport("HYPRE", EntryPoint = "HYPRE_GMRESGetConverged")]
        public static extern int GetConverged(T_Solver m_GMRESSolver, out  int _Converged);
    }
}
