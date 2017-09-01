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
using MPI.Wrappers.Utils;

namespace ilPSP.LinSolvers.HYPRE.Wrappers {
        
    /// <summary>
    /// Wrappers to the HYPER_BoomerAMG Interface
    /// </summary>
    class BoomerAMG : DynLibLoader {

        BoomerAMG()
            : base(
                new string[] { "HYPRE.dll", "libHYPRE.so" },
                new string[2][][],
                new GetNameMangling[] { DynLibLoader.Identity, DynLibLoader.Identity },
                new PlatformID[] { PlatformID.Win32NT, PlatformID.Unix },
                new int[] { -1, -1 }) {
        }

        internal static BoomerAMG my = new BoomerAMG();

        internal delegate int _HYPRE_BoomerAMGSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_BoomerAMGSetup HYPRE_BoomerAMGSetup;
#pragma warning restore        649

        internal delegate int _HYPRE_BoomerAMGSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x);
#pragma warning disable        649
        internal _HYPRE_BoomerAMGSolve HYPRE_BoomerAMGSolve;
#pragma warning restore        649

        static public int __HYPRE_BoomerAMGSetup(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_BoomerAMGSetup(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }

        static public int __HYPRE_BoomerAMGSolve(T_Solver solver, T_ParCSR_matrix ParCSRMatrix_A, T_ParCRS_vector ParVector_b, T_ParCRS_vector ParVector_x) {
            return my.HYPRE_BoomerAMGSolve(solver, ParCSRMatrix_A, ParVector_b, ParVector_x);
        }

        /// <summary>
        /// Create a solver object
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGCreate(out T_Solver solver);

        /// <summary>
        /// Destroy a solver object
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGDestroy(T_Solver solver);

        /// <summary>
        /// Solve the transpose system AT x = b or apply AMG as a preconditioner to
        /// the transpose system 
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSolveT(T_Solver solver, IntPtr ParCSRMatrix_A, IntPtr ParVector_b, IntPtr ParVector_x);

        /// <summary>
        /// (Optional) Set the convergence tolerance, if BoomerAMG is used as a solver
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetTol(T_Solver solver, double tol);

        /// <summary>
        /// Gets the convergence tolerance, if BoomerAMG is used as a solver
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetTol(T_Solver solver, out double tol);

        /// <summary>
        /// (Optional) Sets maximum number of iterations, if BoomerAMG is used as a solver
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMaxIter(T_Solver solver, int max_iter);


        /// <summary>
        /// (Optional) Gets maximum number of iterations, if BoomerAMG is used as a solver
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetMaxIter(T_Solver solver, out int max_iter);



        /// <summary>
        /// (Optional) Sets maximum number of multigrid levels
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMaxLevels(T_Solver solver, int max_levels);
        
        /// <summary>
        /// (Optional) Gets maximum number of multigrid levels
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetMaxLevels(T_Solver solver, out int max_levels);

        /// <summary>
        /// (Optional) Sets AMG strength threshold
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetStrongThreshold(T_Solver solver, double strong_threshold);
        
        /// <summary>
        /// (Optional) Sets AMG strength threshold
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetStrongThreshold(T_Solver solver, out double strong_threshold);

        /// <summary>
        /// (Optional) Sets a parameter to modify the definition of strength for diagonal
        /// dominant portions of the matrix
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMaxRowSum(T_Solver solver, double max_row_sum);

        /// <summary>
        /// (Optional) Sets a parameter to modify the definition of strength for diagonal
        /// dominant portions of the matrix
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetMaxRowSum(T_Solver solver, out double max_row_sum);

        /// <summary>
        /// (Optional) Defines which parallel coarsening algorithm is used.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetCoarsenType(T_Solver solver, int coarsen_type);

        /// <summary>
        /// (Optional) Defines which parallel coarsening algorithm is used.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetCoarsenType(T_Solver solver, out int coarsen_type);

        /// <summary>
        /// (Optional) Defines whether local or global measures are used
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMeasureType(T_Solver solver, int measure_type);

        /// <summary>
        /// (Optional) Defines whether local or global measures are used
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetMeasureType(T_Solver solver, out int measure_type);

        /// <summary>
        /// (Optional) Defines the type of cycle
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetCycleType(T_Solver solver, int cycle_type);

        /// <summary>
        /// (Optional) Defines the type of cycle
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetCycleType(T_Solver solver, out int cycle_type);


        /// <summary>
        /// (Optional) Defines the number of sweeps for the fine and coarse grid, the up and down cycle.
        /// Note: This routine will be phased out!!!! Use HYPRE BoomerAMGSetNumSweeps or
        /// HYPRE BoomerAMGSetCycleNumSweeps instead.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNumGridSweeps(T_Solver solver, int[] num_grid_sweeps);

        ///// <summary>
        ///// Returnes the number of sweeps for the fine and coarse grid, the up and down cycle.
        ///// Note: This routine will be phased out!!!! Use HYPRE BoomerAMGSetNumSweeps or
        ///// HYPRE BoomerAMGSetCycleNumSweeps instead.
        ///// </summary>
        //[DllImport("HYPRE")]
        //public static extern int HYPRE_BoomerAMGGetNumGridSweeps(Solver solver, out int[] num_grid_sweeps);


        /// <summary>
        /// (Optional) Sets the number of sweeps
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNumSweeps(T_Solver solver, int num_sweeps);

      
        /// <summary>
        /// Get the number of sweeps (!not only) at a specified cycle
        /// Note: This function is used to output any(!) *NumSweeps BoomerAMG setup parameter
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetCycleNumSweeps(T_Solver solver, out int num_sweeps, int k);

        /// <summary>
        /// (Optional) Sets the number of sweeps at a specified cycle
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetCycleNumSweeps(T_Solver solver, int num_sweeps, int k);



        /// <summary>
        /// (Optional) Defines the smoother to be used
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetRelaxType(T_Solver solver, int relax_type);

        /// <summary>
        /// (Optional) Defines the smoother at a given cycle
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetCycleRelaxType(T_Solver solver, int relax_type, int k);

        /// <summary>
        /// Defines the smoother at a given cycle
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetCycleRelaxType(T_Solver solver, out int relax_type, int k);


        /// <summary>
        /// (Optional) Defines in which order the points are relaxed
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetRelaxOrder(T_Solver solver, int relax_order);

        /// <summary>
        /// (Optional) Defines in which order the points are relaxed
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetRelaxOrder(T_Solver solver, out int relax_order);

        /// <summary>
        /// (Optional) Defines in which order the points are relaxed
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGSetGridRelaxPoints(T_Solver solver, int** grid_relax_points);

        /// <summary>
        /// (Optional) Defines in which order the points are relaxed
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGGetGridRelaxPoints(T_Solver solver, int** grid_relax_points);


        /// <summary>
        /// (Optional) Defines the relaxation weight for smoothed Jacobi and hybrid
        /// SOR on all levels.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetRelaxWt(T_Solver solver, double relax_weight);

        /// <summary>
        /// (Optional) Defines the relaxation weight for smoothed Jacobi and hybrid
        /// SOR on the user defined level.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetLevelRelaxWt(T_Solver solver, out double relax_weight, int level);

        /// <summary>
        /// number of AMG levels
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetNumLevels(T_Solver solver, out int num_levels);

        /// <summary>
        /// (Optional) Defines the relaxation weight for smoothed Jacobi and hybrid
        /// SOR on the user defined level.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetLevelRelaxWt(T_Solver solver, double relax_weight, int level);

        ///// <summary>
        ///// (Optional) Defines the outer relaxation weight for hybrid SOR.
        ///// Note: This routine will be phased out!!!!
        ///// Use HYPRE\_BoomerAMGSetOuterWt or HYPRE\_BoomerAMGSetLevelOuterWt instead.
        ///// </summary>
        //[DllImport("HYPRE")]
        //unsafe public static extern int HYPRE_BoomerAMGSetOmega(IntPtr solver, double* omega);

        ///// <summary>
        ///// (Optional) Returns the outer relaxation weight for hybrid SOR.
        ///// Note: This routine will be phased out!!!!
        ///// Use HYPRE_BoomerAMGSetOuterWt or HYPRE_BoomerAMGSetLevelOuterWt instead.
        ///// </summary>
        //[DllImport("HYPRE")]
        //unsafe public static extern int HYPRE_BoomerAMGGetOmega(IntPtr solver, double** omega);

        /// <summary>
        /// (Optional) Defines the outer relaxation weight for hybrid SOR and SSOR on all levels.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetOuterWt(T_Solver solver, double omega);

        /// <summary>
        /// (Optional) Defines the outer relaxation weight for hybrid SOR or SSOR on
        /// the user defined level
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetLevelOuterWt(T_Solver solver, double omega, int level);

        /// <summary>
        /// (Optional)
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetDebugFlag(T_Solver solver, int debug_flag);

        /// <summary>
        /// Returns the residual 
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetResidual(T_Solver solver, out IntPtr HYPRE_ParVector_residual);

        /// <summary>
        /// Returns the number of iterations taken
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGGetNumIterations(T_Solver solver, int* num_iterations);
 
        /// <summary>
        /// Returns the norm of the final relative residual
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGGetFinalRelativeResidualNorm(T_Solver solver, double* rel_resid_norm);

        /// <summary>
        /// (Optional) Defines a truncation factor for the interpolation.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetTruncFactor(T_Solver solver, double trunc_factor);

        /// <summary>
        /// (Optional) Defines the maximal number of elements per row for the interpolation.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetPMaxElmts(T_Solver solver, int P_max_elmts);

        /// <summary>
        /// (Optional) Defines the largest strength threshold for which the strength matrix
        /// S uses the communication package of the operator A.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSCommPkgSwitch(T_Solver solver, double S_commpkg_switch);

        /// <summary>
        /// (Optional) Defines which parallel interpolation operator is used
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetInterpType(T_Solver solver, int interp_type);

        /// <summary>
        /// (Optional)
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMinIter(T_Solver solver, int min_iter);

        /// <summary>
        /// (Optional) Enables the use of more complex smoothers.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSmoothType(T_Solver solver, int smooth_type);

        /// <summary>
        /// (Optional) Sets the number of levels for more complex smoothers.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSmoothNumLevels(T_Solver solver, int smooth_num_levels);

        /// <summary>
        /// (Optional) Sets the number of sweeps for more complex smoothers.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSmoothNumSweeps(T_Solver solver, int smooth_num_sweeps);

        /// <summary>
        /// (Optional) Requests automatic printing of setup and solve information.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetPrintLevel(T_Solver solver, int print_level);

        /// <summary>
        /// (Optional) Requests additional computations for diagnostic and similar data
        /// to be logged by the user.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetLogging(T_Solver solver, int logging);

        /// <summary>
        /// (Optional) Sets the size of the system of PDEs, if using the systems version.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNumFunctions(T_Solver solver, int num_functions);

        /// <summary>
        /// (Optional) Sets whether to use the nodal systems version.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNodal(T_Solver solver, int nodal);


        /// <summary>
        /// (Optional) Sets the mapping that assigns the function to each variable, if
        /// using the systems version. 
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGSetDofFunc(T_Solver solver, int* dof_func);

        /// <summary>
        /// (Optional) Defines the number of levels of aggressive coarsening
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetAggNumLevels(T_Solver solver, int agg_num_levels);

        /// <summary>
        /// (Optional) Defines the degree of aggressive coarsening.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNumPaths(T_Solver solver, int num_paths);

        /// <summary>
        /// (Optional) Defines which variant of the Schwarz method is used.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetVariant(T_Solver solver, int variant);

        /// <summary>
        /// (Optional) Defines the overlap for the Schwarz method. 
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetOverlap(T_Solver solver, int overlap);

        /// <summary>
        /// (Optional) Defines the type of domain used for the Schwarz method.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetDomainType(T_Solver solver, int domain_type);


        /// <summary>
        /// (Optional) Defines a smoothing parameter for the additive Schwarz method
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSchwarzRlxWeight(T_Solver solver, double schwarz_rlx_weight);

        /// <summary>
        /// (Optional) Defines symmetry for ParaSAILS.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetSym(T_Solver solver, int sym);

        /// <summary>
        /// (Optional) Defines number of levels for ParaSAILS.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetLevel(T_Solver solver, int level);

        /// <summary>
        /// (Optional) Defines threshold for ParaSAILS.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetThreshold(T_Solver solver, double threshold);

        /// <summary>
        /// (Optional) Defines threshold for ParaSAILS.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetThreshold(T_Solver solver, out double threshold);

        /// <summary>
        /// (Optional) Defines filter for ParaSAILS.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetFilter(T_Solver solver, double filter);

        /// <summary>
        /// (Optional) Defines drop tolerance for PILUT.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetDropTol(T_Solver solver, double drop_tol);

        /// <summary>
        /// (Optional) Defines maximal number of nonzeros for PILUT.
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetMaxNzPerRow(T_Solver solver, int max_nz_per_row);

        /// <summary>
        /// (Optional) Defines name of an input file for Euclid parameters.
        /// </summary>
        [DllImport("HYPRE")]
        unsafe public static extern int HYPRE_BoomerAMGSetEuclidFile(T_Solver solver, char* euclidfile);

        /// <summary>
        /// (Optional) Specifies the use of GSMG - geometrically smooth coarsening and
        /// interpolation. 
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetGSMG(T_Solver solver, int gsmg);

        /// <summary>
        /// (Optional) Defines the number of sample vectors used in GSMG or LS interpolation
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGSetNumSamples(T_Solver solver, int num_samples);
        
        /// <summary>
        /// 
        /// </summary>
        [DllImport("HYPRE")]
        public static extern int HYPRE_BoomerAMGGetLevelOuterWt(T_Solver solver, out double omega, int level);
    }
}
