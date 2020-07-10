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

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {
            if(!base.Equals(obj))
                return false;
            var other = obj as AppControlSolver;
            if(other == null)
                return false;

            if(this.LinearSolver != null) {
                if(!this.LinearSolver.Equals(other.LinearSolver))
                    return false;
            } else {
                if(other.LinearSolver != null)
                    return false;
            }

            if(this.NonLinearSolver != null) {
                if(!this.NonLinearSolver.Equals(other.NonLinearSolver))
                    return false;
            } else {
                if(other.NonLinearSolver != null)
                    return false;
            }


            return true;
        }

    }
}
