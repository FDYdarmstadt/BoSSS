using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json;
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
    public class AppControlSolver : AppControl
    {

        /// <summary>
        /// ctor
        /// </summary>
        public AppControlSolver() {
            this.NoOfMultigridLevels = 1000000;
        }

        /*
        /// <summary>
        /// Linked to <see cref="LinearSolverConfig.NoOfMultigridLevels"/>.
        /// </summary>
        [JsonIgnore]
        public override int NoOfMultigridLevels {
            get {
                return LinearSolver.NoOfMultigridLevels;
            }
            set {
                LinearSolver.NoOfMultigridLevels = value;
            }
        }
        */

        /// <summary>
        /// Configuration of 'primary' linear solver, respectively preconditioner used for <see cref="NonLinearSolver"/>.
        /// </summary>
        [DataMember]
        public AdvancedSolvers.ISolverFactory LinearSolver = new AdvancedSolvers.DirectSolver.Config();

        /// <summary>
        /// Configuration of 'primary' nonlinear solver, if used in application
        /// </summary>
        [DataMember]
        public NonLinearSolverConfig NonLinearSolver = new NonLinearSolverConfig();

        /// <summary>
        /// Cut-cell volume fraction threshold for cell agglomeration
        /// </summary>
        [DataMember]
        public double AgglomerationThreshold = 0.1;


        /// <summary>
        /// In the case of multi-step methods (e.g. BDF2 and higher), multiple initial values, resp. 
        /// </summary>
        [DataMember]
        public bool MultiStepInit = true;

        /// <summary>
        /// Kind of timestepping to use
        /// </summary>
        [DataMember]
        public TimeSteppingScheme TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

        [JsonIgnore]
        public override _TimesteppingMode TimesteppingMode {
            get {
                return base.TimesteppingMode;
            }
            set {
                base.TimesteppingMode = value;
                if (value == _TimesteppingMode.Steady) {
                    TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

                }
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        
        /// <summary>
        /// List of active AMR level indicators 
        /// </summary>
        [DataMember]
        public List<AMRLevelIndicator> activeAMRlevelIndicators = new List<AMRLevelIndicator>();


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

            if(other.AgglomerationThreshold != this.AgglomerationThreshold)
                return false;

            if(other.TimeSteppingScheme != this.TimeSteppingScheme)
                return false;

            if(other.MultiStepInit != this.MultiStepInit)
                return false;

            if(other.LevelSet_ConvergenceCriterion != this.LevelSet_ConvergenceCriterion)
                return false;

            if (!IEnumerableExtensions.SetEquals(this.activeAMRlevelIndicators, other.activeAMRlevelIndicators))
                return false;


            return true;
        }

        /// <summary>
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double LevelSet_ConvergenceCriterion = 1.0e-6;
    }
}
