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

    public enum LinearSolverCode {

        /// <summary>
        /// Automatic choose of linear solver depending on nonlinear solver, problem size, etc.
        /// </summary>
        automatic = 0,

        //direct solvers

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.MUMPS"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_mumps = 1,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.PARDISO.PARDISOSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_pardiso = 2,

        //Schwarz: Domain decomposition (direct Solver)

        /// <summary>
        /// Schwarz method with METIS blocking and direct solver for coarse solve with an overlap of 1
        /// </summary>
        exp_AS = 11,

        exp_AS_MG = 12,

        //GMRES (iterative Solver)

        /// <summary>
        /// GMRES without any preconditioning
        /// </summary>
        exp_softgmres = 20,


        //GMRES Solver with different experimental Preconditioner

        exp_gmres_Schur = 24,
        exp_gmres_localPrec = 25,
        exp_gmres_Simple = 26,
        exp_gmres_AS = 27,
        exp_gmres_AS_MG = 30,

        //CG versions

        /// <summary>
        /// Conjugate gradient (from monkey library) without any preconditioning.
        /// </summary>
        classic_cg = 40,

        /// <summary>
        /// Multiple levels of additive Schwarz, in a Krylov multi-grid cycle.
        /// </summary>
        exp_Kcycle_schwarz = 41,

        exp_softpcg_jacobi_mg =42,

        /// <summary>
        /// Conjugate gradient, with additive Schwarz preconditioner.
        /// </summary>
        exp_softpcg_schwarz = 43,

        /// <summary>
        /// Conjugate gradient, with additive Schwarz preconditioner, including a coarse-grid solver.
        /// </summary>
        exp_softpcg_schwarz_directcoarse = 44,

        /// <summary>
        /// GMRES with p-multigrid on the same mesh level; direct solver is used for 
        /// </summary>
        exp_gmres_levelpmg = 47,

        exp_gmres_schwarz_pmg = 48,

        /// <summary>
        /// highly experimental shit
        /// </summary>
        exp_decomposedMG_OrthoScheme = 50,

        /// <summary>
        /// Orthonormalization Scheme with p-multigrid preconditioner
        /// </summary>
        exp_OrthoS_pMG = 51,

        /// <summary>
        /// Work-in-progress: experimental stuff for rheology solver
        /// </summary>
        exp_Kcycle_schwarz_4Rheology = 52,



        selfmade = 999,
    }

    public enum LinearSolverMode
    {

        /// <summary>
        /// Standard Mode, perform the simulation (solve the linear system)
        /// </summary>
        Solve = 0,

        //direct solvers

        /// <summary>
        /// Set RHS to zero and examine the error spectrum before and after solving
        /// </summary>
        SpectralAnalysis = 1,

    }

    /// <summary>
    /// The linear solver config
    /// </summary>
    [Serializable]
    public class LinearSolverConfig : ICloneable, IEquatable<LinearSolverConfig> {

        /// <summary>
        /// This will print out more information about iterations.
        /// </summary>
        public bool verbose = false;

        /// <summary>
        /// If iterative saddle-point solvers like GMRES or Orthonormalization are used, the maximum number of basis vectors
        /// that are used to construct the accelerated solution.
        /// </summary>
        [DataMember]
        public int MaxKrylovDim = 30;

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
        public double ConvergenceCriterion = 1e-10;

        /// <summary>
        /// Sets the algorithm to use for linear solving, e.g. MUMPS or GMRES.
        /// </summary>
        [DataMember]
        public LinearSolverCode SolverCode = LinearSolverCode.classic_pardiso;

        /// <summary>
        /// Sets the number of Multigrid levels. Multigrid approach is used to get a Preconditioner for Krylov solvers, e.g. GMRES.
        /// </summary>
        [DataMember]
        public int NoOfMultigridLevels = 1;

        /// <summary>
        /// Sets the mode for the solver to run in
        /// </summary>
        [DataMember]
        public LinearSolverMode SolverMode = LinearSolverMode.Solve;

        //-------------------------
        // Probably legacy code: can be deleted, if exp_localPrec is removed.
        /// <summary>
        /// The physical viscosity has to be written to <see cref="exp_localPrec_muA"/>, if the experimental linear solver <see cref="LinearSolverCode.exp_localPrec"/> is used.
        /// </summary>
        [DataMember]
        public int exp_localPrec_muA = 1;


        /// <summary>
        /// The minimum time step has to be written to <see cref="exp_localPrec_Min_dt"/>, if the experimental linear solver <see cref="LinearSolverCode.exp_localPrec"/> is used.
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

        /// <summary>
        /// Clones the LinearConfig
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var clone = new LinearSolverConfig() {
                verbose = this.verbose,
                MaxKrylovDim = this.MaxKrylovDim,
                MaxSolverIterations = this.MaxSolverIterations,
                MinSolverIterations = this.MinSolverIterations,
                NoOfMultigridLevels = this.NoOfMultigridLevels,
                SolverCode = this.SolverCode,
                TargetBlockSize = this.TargetBlockSize
            };
            return clone;
        }

        /// <summary>
        /// Compares value not ref!
        /// </summary>
        /// <param name="compareto"></param>
        /// <returns></returns>
        public bool Equals(LinearSolverConfig compareto) {
            if(compareto == null)
                return false;

            return this.verbose == compareto.verbose &&
                this.MaxKrylovDim == compareto.MaxKrylovDim &&
                this.MaxSolverIterations == compareto.MaxSolverIterations &&
                this.MinSolverIterations == compareto.MinSolverIterations &&
                this.NoOfMultigridLevels == compareto.NoOfMultigridLevels &&
                this.SolverCode == compareto.SolverCode &&
                this.TargetBlockSize == compareto.TargetBlockSize;
        }
    }
}
