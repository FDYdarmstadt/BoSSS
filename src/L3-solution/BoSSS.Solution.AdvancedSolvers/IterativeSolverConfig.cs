using System;
using System.Collections.Generic;
using System.Runtime.Serialization;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {
    
    /// <summary>
    /// Base-class for iterative solver configurations
    /// </summary>
    [Serializable]
    abstract public class IterativeSolverConfig : ISolverFactory {

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


        /// <summary>
        /// Termination criterion based on the threshold <see cref="ConvergenceCriterion"/>.
        /// </summary>
        /// <returns>
        /// - 1st item: true: solver should continue; false: terminate;
        /// - 2nd item: if first item true, either success (true, i.e. solver converged successfully) or fail (false, e.g. reached the <see cref="MaxSolverIterations"/>);
        /// </returns>
        public (bool bNotTerminate, bool bSuccess) DefaultTermination(int iter, double R0_l2, double R_l2) {
            if(iter <= MinSolverIterations)
                return (true, false); // keep running

            if(R_l2 < R0_l2 * ConvergenceCriterion + ConvergenceCriterion)
                return (false, true); // success

            if(iter > MaxSolverIterations)
                return (false, false); // fail

            return (true, false); // keep running
        }

        /// <summary>
        /// Name/Description of the solver configuration
        /// </summary>
        public abstract string Name { get; }
        
        /// <summary>
        /// short version of <see cref="Name"/>, e.g. to be used in plots or tables
        /// </summary>
        public abstract string Shortname { get; }

        /// <summary>
        /// <see cref="ISolverFactory.CreateInstance"/>
        /// </summary>
        public abstract ISolverSmootherTemplate CreateInstance(MultigridOperator level);
    }
}
