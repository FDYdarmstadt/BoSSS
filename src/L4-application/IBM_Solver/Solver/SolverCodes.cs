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

        classic_mumps = 1,

        classic_pardiso = 2,

        exp_schwarz_directcoarse = 3,

        exp_schwarz_Kcycle_directcoarse = 4,    

        exp_softgmres = 5,

        exp_softgmres_schwarz_directcoarse = 6,

    }

}
