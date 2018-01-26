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

        exp_direct,


        exp_softpcg_mg,


        exp_softpcg_schwarz,

        exp_softpcg_schwarz_directcoarse
    }
}
