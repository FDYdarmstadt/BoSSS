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

                case LinearSolverCode.direct_mumps:
                    return new AdvancedSolvers.DirectSolver.Config() { WhichSolver = AdvancedSolvers.DirectSolver._whichSolver.MUMPS };

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
