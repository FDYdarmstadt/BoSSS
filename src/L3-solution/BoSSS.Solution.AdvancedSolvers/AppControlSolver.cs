using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Control {

    /// <summary>
    /// PDE-solver-control object which defines configuration options for nonlinear and linear equation solvers
    /// </summary>
    [Serializable]
    [DataContract]
    public class AppControlSolver : AppControl {

        /// <summary>
        /// Linked to <see cref="LinearSolverConfig.NoOfMultigridLevels"/>.
        /// </summary>
        [DataMember]
        public override int NoOfMultigridLevels {
            get {
                return LinearSolver.NoOfMultigridLevels;
            }
            set {
                LinearSolver.NoOfMultigridLevels = value;
            }
        }

        /// <summary>
        /// Configuration of 'primary' linear solver, respectively preconditioner used for <see cref="NonLinearSolver"/>.
        /// </summary>
        [DataMember]
        public LinearSolverConfig LinearSolver = new LinearSolverConfig();

        /// <summary>
        /// Configuration of 'primary' nonlinear solver, if used in application
        /// </summary>
        [DataMember]
        public NonLinearSolverConfig NonLinearSolver = new NonLinearSolverConfig();


    }
}
