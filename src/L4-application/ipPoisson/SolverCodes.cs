using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson {
    
    
    /// <summary>
    /// All variants for supported solvers
    /// </summary>
    public enum SolverCodes {

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.PARDISO.PARDISOSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_pardiso = 0,

        /// <summary>
        /// Direct solver (<see cref="ilPSP.LinSolvers.MUMPS.MUMPSSolver"/>) without any pre-processing of the matrix.
        /// </summary>
        classic_mumps = 1,

        /// <summary>
        /// Conjugate gradient (<see cref="ilPSP.LinSolvers.monkey.CG"/>) without any preconditioning.
        /// </summary>
        classic_cg = 2,

        /// <summary>
        /// Direct solver from new solver framework (causes some overhead in comparison to <see cref="classic_pardiso"/>).
        /// </summary>
        exp_direct = 3,


        /// <summary>
        /// Direct solver from new solver framework, using a dense LU decomposition from Lapack.
        /// </summary>
        exp_direct_lapack = 7,

        /// <summary>
        /// Conjugate gradient, with multi-grid preconditioner.
        /// </summary>
        exp_softpcg_mg = 4,


        /// <summary>
        /// Conjugate gradient, with additive Schwarz preconditioner.
        /// </summary>
        exp_softpcg_schwarz = 5,

        /// <summary>
        /// Conjugate gradient, with additive Schwarz preconditioner, including a coarse-grid solver.
        /// </summary>
        exp_softpcg_schwarz_directcoarse = 6,

        /// <summary>
        /// Multiple levels of additive Schwarz, in a Krylov multi-grid cycle.
        /// </summary>
        exp_Kcycle_schwarz = 8,

      


    }
}
