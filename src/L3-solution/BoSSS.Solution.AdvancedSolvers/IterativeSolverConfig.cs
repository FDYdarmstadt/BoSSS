using System;
using System.Collections.Generic;
using System.Reflection;
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
        public int MaxSolverIterations = 300;


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

        /// <summary>
        /// 
        /// </summary>
        public bool Equals(ISolverFactory other) {
            return EqualsImpl(other);
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {
            return EqualsImpl(obj);
        }

        private bool EqualsImpl(object o) {
            if(object.ReferenceEquals(o, this))
                return true;

            var other = o as IterativeSolverConfig;
            if(other == null)
                return false;

            // equality through reflection 

            var t = this.GetType();
            if(!t.Equals(other.GetType()))
                return false;

            var flds = t.GetFields(BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic | BindingFlags.FlattenHierarchy | BindingFlags.Instance);

            foreach(var fld in flds) {
                var fld_this = fld.GetValue(this);
                var fld_othr = fld.GetValue(other);

                if(fld_this == fld_othr)
                    continue;
                if(fld_othr == null && fld_this != null)
                    return false;
                if(fld_othr != null && fld_this == null)
                    return false;

                if(!fld_othr.Equals(fld_this))
                    return false;
            }


            return true;
        }

        public override int GetHashCode() {
            return 0;
        }
    }
}
