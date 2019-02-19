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
using System.Threading.Tasks;
using System.Runtime.Serialization;

namespace BoSSS.Solution.Control {

    [Serializable]
    public class LinearSolverConfig {



        public enum Code {
            
            /// <summary>
            /// Automatic choose of linear solver depending on nonlinear solver, problem size, etc.
            /// </summary>
            automatic = 666,

            //direct solvers

            /// <summary>
            /// Direct solver (<see cref="ilPSP.LinSolvers.PARDISO.PARDISOSolver"/>) without any pre-processing of the matrix.
            /// </summary>
            classic_pardiso = 0,

            /// <summary>
            /// Direct solver (<see cref="ilPSP.LinSolvers.MUMPS"/>) without any pre-processing of the matrix.
            /// </summary>
            classic_mumps = 1,

            /// <summary>
            /// Classic Multigrid approach, especially useful for predoncitioning
            /// </summary>
            exp_multigrid = 3,

            /// <summary>
            /// ILU decomposition with modification for saddle-point (highly experimental)
            /// </summary>
            exp_ILU = 4,

            /// <summary>
            /// Direct solver from new solver framework, using a dense LU decomposition from Lapack.
            /// </summary>
            exp_direct_lapack = 7,

            //Schwarz: Domain decomposition (direct Solver)

            /// <summary>
            /// Schwarz method with METIS blocking and direct solver for coarse solve
            /// </summary>
            exp_schwarz_directcoarse = 10,

            /// <summary>
            /// Schwarz method with METIS blocking and direct solver for coarse solve with an overlap of 1
            /// </summary>
            exp_schwarz_directcoarse_overlap = 11,

            /// <summary>
            /// Schwarz method with Multigridblocking on the coarsest level and coarse solve
            /// </summary>
            exp_schwarz_Kcycle_directcoarse = 12,

            /// <summary>
            /// Schwarz method with Multigridblocking on the coarsest level and coarse solve with an overlap of 1
            /// </summary>
            exp_schwarz_Kcycle_directcoarse_overlap = 13,

            //GMRES (iterative Solver)

            /// <summary>
            /// GMRES without any preconditioning
            /// </summary>
            exp_softgmres = 20,

            /// <summary>
            /// GMRES with schwarz precoditioning using METIS blocking without overlap
            /// </summary>
            exp_softgmres_schwarz_directcoarse = 21,

            /// <summary>
            /// GMRES with schwarz precoditioning using METIS blocking with overlap
            /// </summary>
            exp_softgmres_schwarz_directcoarse_overlap = 22,

            /// <summary>
            /// GMRES with schwarz precoditioning using MG blocking with overlap
            /// </summary>
            exp_softgmres_schwarz_Kcycle_directcoarse_overlap = 23,

            //GMRES Solver with different experimental Preconditioner

            exp_Schur = 24,
            exp_localPrec = 25,
            exp_Simple = 26,
            exp_AS_1000 = 27,
            exp_AS_5000 = 28,
            exp_AS_10000 = 29,
            exp_AS_MG = 30,

            //CG versions

            /// <summary>
            /// Conjugate gradient (<see cref="ilPSP.LinSolvers.monkey.CG"/>) without any preconditioning.
            /// </summary>
            classic_cg = 40,

            /// <summary>
            /// Multiple levels of additive Schwarz, in a Krylov multi-grid cycle.
            /// </summary>
            exp_Kcycle_schwarz = 41,

            /// <summary>
            /// Conjugate gradient, with multi-grid preconditioner.
            /// </summary>
            exp_softpcg_mg = 42,

            /// <summary>
            /// Conjugate gradient, with additive Schwarz preconditioner.
            /// </summary>
            exp_softpcg_schwarz = 43,

            /// <summary>
            /// Conjugate gradient, with additive Schwarz preconditioner, including a coarse-grid solver.
            /// </summary>
            exp_softpcg_schwarz_directcoarse = 44,

        }

        /// <summary>
        /// If iterative saddle-point solvers like GMRES or Orthonormalization are used, the maximum number of basis vectors
        /// that are used to construct the accelerated solution.
        /// </summary>
        [DataMember]
        public int MaxKrylovDim=30;

        /// <summary>
        /// If iterative solvers are used, the maximum number of iterations.
        /// </summary>
        [DataMember]
        public int MaxSolverIterations = 2000;


        /// <summary>
        /// If iterative solvers are used, the minimum number of iterations.
        /// </summary>
        [DataMember]
        public int MinSolverIterations = 2;

        /// <summary>
        /// Convergence criterion for linear solver.
        /// </summary>
        [DataMember]
        public double ConvergenceCriterion = 1e-8;

        /// <summary>
        /// Sets the algorithm to use for linear solving, e.g. MUMPS or GMRES.
        /// </summary>
        [DataMember]
        public LinearSolverConfig.Code SolverCode = LinearSolverConfig.Code.classic_mumps;

        /// <summary>
        /// Sets the number of Multigrid levels. Multigrid approach is used to get a Preconditioner for Krylov solvers, e.g. GMRES.
        /// </summary>
        [DataMember]
        public int NoOfMultigridLevels = 1;

        //-------------------------
        // These parameters have to be set only, if exp_localPrec is used. They can be deleted, if exp_localPrec is removed.
        /// <summary>
        /// The physical viscosity has to be written to <see cref="exp_localPrec_muA"/>, if the experimental linear solver <see cref="LinearSolverConfig.Code.exp_localPrec"/> is used.
        /// </summary>
        [DataMember]
        public int exp_localPrec_muA = 1;

        /// <summary>
        /// The minimum time step has to be written to <see cref="exp_localPrec_Min_dt"/>, if the experimental linear solver <see cref="LinearSolverConfig.Code.exp_localPrec"/> is used.
        /// </summary>
        [DataMember]
        public int exp_localPrec_Min_dt = 0;

        /// <summary>
        /// If any blocking is used (Schwarz, block Jacobi), a target for the block size.
        /// Tests show that the ideal block size may be around 10000, but this may depend on computer, DG polynomial order, etc.
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.ExclusiveLowerBound(99.0)]
        public int TargetBlockSize = 10000;
    }
}
