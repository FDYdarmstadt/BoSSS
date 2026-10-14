using System;
using System.Runtime.Serialization;

namespace BoSSS.Solution.LevelSetTools.StokesExtension {

    /// <summary>
    /// Options for the Stokes-extension velocity solve
    /// </summary>
    [Serializable]
    [DataContract]
    public class StokesExtensionControl {

        /// <summary>
        /// Overrides the solver dispatch in <see cref="StokesExtension.SolveExtension"/>; null (default) is a plain direct solve.
        /// </summary>
        [DataMember]
        public BoSSS.Solution.AdvancedSolvers.ISolverFactory LinearSolver = null;

        /// <summary>
        /// If true, warm-starts the solve from the previous time step's converged extension velocity.
        /// </summary>
        [DataMember]
        public bool WarmStart = false;

        /// <summary>
        /// If true, combines velocity and pressure into a single Schur-complement preconditioning block.
        /// </summary>
        [DataMember]
        public bool UseSchurBlockPrec = false;

        /// <summary>
        /// If true, solves via pseudo-transient continuation instead of a single solve.
        /// </summary>
        [DataMember]
        public bool UsePseudoTimestepping = false;

        /// <summary>
        /// Only relevant if <see cref="UsePseudoTimestepping"/>: initial pseudo-time step size. Null
        /// (default) resolves to the flow solve's fixed timestep size (<see cref="BoSSS.Solution.AppControl.dtFixed"/>)
        /// if one is set, or 1.0 otherwise.
        /// </summary>
        [DataMember]
        public double? PseudoDtau0 = null;

        /// <summary>
        /// Only relevant if <see cref="UsePseudoTimestepping"/>: upper bound for the pseudo-time step size.
        /// </summary>
        [DataMember]
        public double PseudoDtauMax = 1.0e8;

        /// <summary>
        /// Only relevant if <see cref="UsePseudoTimestepping"/>: maximum number of pseudo-time steps.
        /// </summary>
        [DataMember]
        public int MaxPseudoTimeSteps = 20;

        /// <summary>
        /// Only relevant if <see cref="UsePseudoTimestepping"/>: convergence criterion, relative to ‖RHS‖.
        /// </summary>
        [DataMember]
        public double PseudoTimeConvergenceCriterion = 1.0e-10;

        /// <summary>
        /// If true, skips local preconditioning entirely, handing the raw matrix to the solver.
        /// </summary>
        [DataMember]
        public bool DisableLocalPreconditioning = false;

        /// <summary>
        /// Artificial viscosity of the StokesExtension equation itself -- not the
        /// physical flow's viscosity. Affects both the conditioning of the extension's linear system
        /// and the relative strength of the interface-matching penalty against the viscous term.
        /// </summary>
        [DataMember]
        public double Viscosity = 1.0;

        /// <summary>
        /// Amount of debug information to write (via <see cref="BoSSS.Solution.LevelSetTools.StokesExtension.StokesExtension.SolveExtension"/>'s
        /// <c>FuncTrace</c>): 0 writes nothing extra; 2 or higher additionally computes and logs,
        /// per pseudo-time step, the interface-velocity-mismatch and divergence norms.
        /// </summary>
        [DataMember]
        public uint WriteDebugInfoLevel = 0;

        /// <summary>
        /// Shallow copy; used to resolve externally-derived defaults (e.g. <see cref="PseudoDtau0"/>'s
        /// flow-timestep fallback) onto a private instance instead of mutating a shared control object.
        /// </summary>
        public StokesExtensionControl Clone() {
            return (StokesExtensionControl)this.MemberwiseClone();
        }
    }
}
