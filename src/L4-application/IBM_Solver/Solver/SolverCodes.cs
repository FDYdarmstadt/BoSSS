using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.IBM_Solver {

    /// <summary>
    /// All supported nonlinear solvers.
    /// </summary>
    public enum NonlinearSolverCodes {

        /// <summary>
        /// NewtonKrylov GMRES (<see cref="BoSSS.Solution.Multigrid.Newton"/>) with linear solver (<see cref="LinearSolverCodes"/>) used as preconditioner for matrix-free GMRES 
        /// </summary>
        NewtonGMRES = 0,

        /// <summary>
        /// Picard fixpoint solver (<see cref="BoSSS.Solution.Multigrid.FixpointIterator"/>) with linear solver (<see cref="LinearSolverCodes"/>) for the linearized equation system
        /// </summary>
        Picard = 1,

    }

    /// <summary>
    /// All supported linear solvers
    /// </summary>
    public enum LinearSolverCodes {

        /// <summary>
        /// Automatic choose of linear solver depending on nonlinear solver, problem size, etc.
        /// </summary>
        automatic = 0,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.MUMPS"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_mumps = 1,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.PARDISO.PARDISOSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_pardiso = 2,

        /// <summary>
        /// Schwarz method with METIS blocking and direct solver for coarse solve
        /// </summary>
        exp_schwarz_directcoarse = 3,

        /// <summary>
        /// Schwarz method with METIS blocking and direct solver for coarse solve with an overlap of 1
        /// </summary>
        exp_schwarz_directcoarse_overlap = 4,

        /// <summary>
        /// Schwarz method with Multigridblocking on the coarsest level and coarse solve
        /// </summary>
        exp_schwarz_Kcycle_directcoarse = 5,

        /// <summary>
        /// Schwarz method with Multigridblocking on the coarsest level and coarse solve with an overlap of 1
        /// </summary>
        exp_schwarz_Kcycle_directcoarse_overlap = 6,

        /// <summary>
        /// GMRES without any preconditioning
        /// </summary>
        exp_softgmres = 7,

        exp_softgmres_schwarz_directcoarse = 8,

        exp_softgmres_schwarz_Kcycle_directcoarse_overlap = 9,



    }

}
