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
using ilPSP.LinSolvers.PETSc;

namespace BoSSS.Solution.Control {

    
    /// <summary>
    /// Abbreviations for various pre-defined solver configurations (objects implementing <see cref="AdvancedSolvers.ISolverFactory"/>)
    /// to be used with the 'convenience method' <see cref="SolverConfigFactory.GetConfig"/>.
    /// </summary>
    public enum LinearSolverCode {

        /// <summary>
        /// Automatic choose of linear solver depending on nonlinear solver, problem size, etc.
        /// </summary>
        automatic = 0,


        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.MUMPS"/>) without any pre-processing of the matrix.
        /// </summary>
        direct_mumps = 1,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.PARDISO.PARDISOSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        direct_pardiso = 2,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.UMFPACK.UMFPACKSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        direct_umfpack = 3,

        /// <summary>
        /// PETSc direct LU using PETSc's default factorization package. The
        /// portable baseline supports this option on a single MPI rank.
        /// </summary>
        direct_petsc = 10,

        /// <summary>
        /// PETSc direct LU using MUMPS. The selected native PETSc library must
        /// have been built with MUMPS.
        /// </summary>
        direct_petsc_mumps = 11,

        /// <summary>
        /// PETSc direct LU using SuperLU_DIST. The selected native PETSc library
        /// must have been built with SuperLU_DIST.
        /// </summary>
        direct_petsc_superlu_dist = 12,

        /// <summary>
        /// PETSc direct LU using SuiteSparse UMFPACK. This is a single-rank
        /// backend and must be present in the selected native PETSc library.
        /// </summary>
        direct_petsc_umfpack = 13,

        /// <summary>
        /// PETSc direct Cholesky using PETSc's default factorization package.
        /// The matrix must be symmetric positive-definite.
        /// </summary>
        direct_petsc_cholesky = 14,


        /// <summary>
        /// Direct solver, the default option
        /// </summary>
        direct = 2,


        ///// <summary>
        ///// Conjugate gradient (from monkey library) without any preconditioning (<see cref="AdvancedSolvers.MonkeySolver.Config"/>
        ///// </summary>
        //cg = 40,

        /// <summary>
        /// Multiple levels of additive Schwarz, in a Krylov multi-grid cycle; <see cref="AdvancedSolvers.OrthoMGSchwarzConfig"/>
        /// </summary>
        exp_Kcycle_schwarz = 41,

        /// <summary>
        /// Equal to <see cref="exp_Kcycle_schwarz"/>, but hardcoded to employ <see cref="AdvancedSolvers.SchwarzImplementation.CoarseMesh"/>
        /// </summary>
        exp_Kcycle_schwarz_CoarseMesh = 42,

        /// <summary>
        /// Equal to <see cref="exp_Kcycle_schwarz"/>, but hardcoded to employ <see cref="AdvancedSolvers.SchwarzImplementation.PerProcess"/>
        /// </summary>
        exp_Kcycle_schwarz_PerProcess = 43,


        /// <summary>
        /// GMRES with p-multigrid on the same mesh level; direct solver is used for lowest polynomial level, <see cref="AdvancedSolvers.PTGconfig"/>
        /// </summary>
        exp_gmres_levelpmg = 47,

        
        /// <summary>
        /// a k-cycle (i.e. a tree of <see cref="AdvancedSolvers.OrthonormalizationMultigrid"/> solvers), with ILU-preconditioner (<see cref="AdvancedSolvers.CellILU"/>) at each level
        /// </summary>
        exp_Kcycle_ILU = 54,

        /// <summary>
        /// p-Multigrid (i.e. multigrid over DG polynomial degree), <see cref="AdvancedSolvers.PmgConfig"/>
        /// </summary>
        pMultigrid = 61,

        /// <summary>
        /// A highly experimental variant of Schur complement designed for saddle point problems, <see cref="AdvancedSolvers.SchurPrecond"/>
        /// </summary>
        SchurPrecond = 71,

		/// <summary>
		/// Flexible GMRES, allowing variable preconditioners for GMRES <see cref="AdvancedSolvers.FlexGMRES"/>
		/// </summary>
		FlexGMRES = 81,


        /// <summary>
        /// A Schur complement solver for Stokes/Navier-Stokes systems, <see cref="AdvancedSolvers.SchurComplementForStokes"/>
        /// </summary>
        SchurComplementForStokes = 101,
    }
    

    /// <summary>
    /// Extension methods
    /// </summary>
    static public class SolverConfigFactory {

        /// <summary>
        /// Convenience method: turns a simple configuration code (<paramref name="config"/>)
        /// into a full-fledged solver configuration.
        /// </summary>
        static public BoSSS.Solution.AdvancedSolvers.ISolverFactory GetConfig(this LinearSolverCode config) {

            switch (config) {
                case LinearSolverCode.automatic:
                case LinearSolverCode.direct_pardiso:
                    return new AdvancedSolvers.DirectSolver.Config() { WhichSolver = AdvancedSolvers.DirectSolver._whichSolver.PARDISO };

                case LinearSolverCode.direct_umfpack:
                    return new AdvancedSolvers.DirectSolver.Config() { WhichSolver = AdvancedSolvers.DirectSolver._whichSolver.UMFPACK };

                case LinearSolverCode.direct_mumps:
                    return new AdvancedSolvers.DirectSolver.Config() { WhichSolver = AdvancedSolvers.DirectSolver._whichSolver.MUMPS };

                case LinearSolverCode.direct_petsc:
                    return new AdvancedSolvers.PETScDirectSolverConfig();

                case LinearSolverCode.direct_petsc_mumps:
                    return new AdvancedSolvers.PETScDirectSolverConfig() {
                        Backend = PETScDirectBackend.Mumps
                    };

                case LinearSolverCode.direct_petsc_superlu_dist:
                    return new AdvancedSolvers.PETScDirectSolverConfig() {
                        Backend = PETScDirectBackend.SuperLuDist,
                        // SuperLU_DIST uses PETSc's assembled matrix path. This
                        // copies the local BoSSS CSR data and lets PETSc perform
                        // its normal distributed matrix assembly.
                        MatrixInput = PETScMatrixInput.CopiedCsr
                    };

                case LinearSolverCode.direct_petsc_umfpack:
                    return new AdvancedSolvers.PETScDirectSolverConfig() {
                        Backend = PETScDirectBackend.Umfpack,
                        // UMFPACK is sequential and requires PETSc MATSEQAIJ.
                        // CopiedCsr lets PETSc select that matrix type for this
                        // single-rank solver instead of forcing MPIAIJ.
                        MatrixInput = PETScMatrixInput.CopiedCsr
                    };

                case LinearSolverCode.direct_petsc_cholesky:
                    return new AdvancedSolvers.PETScDirectSolverConfig() {
                        Factorization = PETScDirectFactorization.Cholesky
                    };

                case LinearSolverCode.exp_Kcycle_schwarz:
                    return new AdvancedSolvers.OrthoMGSchwarzConfig();

                case LinearSolverCode.exp_Kcycle_schwarz_CoarseMesh: {
                        var c = new AdvancedSolvers.OrthoMGSchwarzConfig();
                        c.SchwarzImplementation = AdvancedSolvers.SchwarzImplementation.CoarseMesh;
                        return c;
                    }
                case LinearSolverCode.exp_Kcycle_schwarz_PerProcess: {
                        var c = new AdvancedSolvers.OrthoMGSchwarzConfig();
                        c.SchwarzImplementation = AdvancedSolvers.SchwarzImplementation.PerProcess;
                        return c;
                    }


                case LinearSolverCode.exp_gmres_levelpmg:
                    return new AdvancedSolvers.PTGconfig();


                case LinearSolverCode.exp_Kcycle_ILU:
                    return new AdvancedSolvers.OrthoMGILUconfig();

                case LinearSolverCode.pMultigrid:
                    return new AdvancedSolvers.PmgConfig();

				//case LinearSolverCode.SchurPrecond:
				//	return new AdvancedSolvers.SchurPrecondConfig();

				case LinearSolverCode.SchurComplementForStokes:
					return new AdvancedSolvers.SchurComplementForStokes.Config();

                case LinearSolverCode.FlexGMRES:
					return new AdvancedSolvers.FGMRESConfig();

				default:
                    throw new NotImplementedException("todo: " + config);
            }


        }

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


    /*
    /// <summary>
    /// The linear solver config
    /// </summary>
    [Serializable]
    public class LinearSolverConfig : IEquatable<LinearSolverConfig> {

        private bool m_verbose = false;

        /// <summary>
        /// This will print out more information about iterations.
        /// </summary>
        public bool verbose {
            get { return m_verbose; }
            set {
                m_verbose = value;
            }
        }

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
        public int MinSolverIterations = 1;

        /// <summary>
        /// Convergence criterion for linear solver.
        /// </summary>
        [DataMember]
        public double ConvergenceCriterion = 1e-10;

        private LinearSolverCode m_SolverCode= LinearSolverCode.classic_pardiso;

        /// <summary>
        /// Sets the algorithm to use for linear solving, e.g. MUMPS or GMRES.
        /// </summary>
        [DataMember]
        public LinearSolverCode SolverCode {
            get { return m_SolverCode; }
            set {
                m_SolverCode = value;  
            }
        }

        /// <summary>
        /// Sets the **maximum** number of Multigrid levels to be used.
        /// The numbers of levels which are actually used is probably much less, and detemeined via <see cref="TargetBlockSize"/>.
        /// Multigrid approach is used to get a Preconditioner for Krylov solvers, e.g. GMRES.
        /// </summary>
        [DataMember]
        public int NoOfMultigridLevels = 1000000;

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
        /// 
        /// This also determines the actual number of multigrid levels used in dependence of the problem size;
        /// As soon as the number of DOF's on a certain multigrid level fall below this threshold, a direct 
        /// solver is used and no further multigrid levels are allocated.
        /// </summary>
        [DataMember]
        [BoSSS.Solution.Control.ExclusiveLowerBound(99.0)]
        public int TargetBlockSize = 10000;

        /// <summary>
        /// Determines maximal DG order within coarse system of a p-Multigrid. Only applicable for p-two-grid, e.g. Schwarz with p-MG or PTG <see cref="exp_gmres_levelpmg"/> preconditioner.
        /// </summary>
        [DataMember]
        public int pMaxOfCoarseSolver = 1;


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

    */
}
