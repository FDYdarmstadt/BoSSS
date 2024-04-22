using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Solution.LevelSetTools.PhasefieldLevelSet.PhasefieldControl;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet {

    public class PhasefieldSettings {

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// diff = (kappa*lambda)/epsilon
        /// diffusivity
        /// </summary>
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        [DataMember]
        public double diff = 1e-3;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// Cahn´s Number, rule of thumb ~(1.0 / 4.164) * h * 8.0 / (p + 1);
        /// </summary>
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        [DataMember]
        public double cahn = 1.0;

        /// <summary>
        /// ContactAngle for walls
        /// </summary>
        [DataMember]
        public Formula theta = new Formula("X => Math.PI/2.0");

        /// <summary>
        /// Timestepping scheme, use RK as the implementation does not remember "old" timesteps (necessary for BDF)
        /// </summary>
        [DataMember]
        public XdgTimestepping.TimeSteppingScheme TimeSteppingScheme = XdgTimestepping.TimeSteppingScheme.RK_ImplicitEuler;

        /// <summary>
        /// Correction to e.g. conserve area
        /// </summary>
        [DataMember]
        public Correction CorrectionType = Correction.None;
    }

    /// <summary>
    /// Model Type of Phasefield equation, see Halperin (1977)
    /// </summary>
    public enum ModelType {
        /// <summary>
        /// Order Parameter is nonconserved
        /// </summary>
        modelA,

        /// <summary>
        /// Order Parameter is conserved
        /// </summary>
        modelB,

        /// <summary>
        /// Mass is conserved
        /// </summary>
        modelC
    }

    /// <summary>
    /// Type of algebraic correction that is performed
    /// </summary>
    public enum Correction {
        /// <summary>
        /// Mass of a phase is conserved
        /// </summary>
        Mass,

        /// <summary>
        /// Total Concentration is conserved
        /// </summary>
        Concentration,

        /// <summary>
        /// No algebraic correction
        /// </summary>
        None
    }


    public class PhasefieldControl : AppControlSolver, ILevSetControl {

        /// <summary>
        /// ctor
        /// </summary>
        public PhasefieldControl() : base() {
            NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            NonLinearSolver.MaxSolverIterations = 50;

            TimesteppingMode = _TimesteppingMode.Transient;

            savetodb = false;
        }     

        /// <summary>
        /// Set the <see cref="ModelType"/>
        /// </summary>
        public ModelType ModTyp = ModelType.modelB;        

        public Correction CorrectionType = Correction.None;

        ///// <summary>
        ///// Type of algebraic correction that is performed
        ///// </summary>
        //public enum CurvatureCorrection {
        //    /// <summary>
        //    /// Curvature by divergence form
        //    /// </summary>
        //    FullyCoupled,

        //    /// <summary>
        //    /// Direct evaluation from concentration field, evaluated in each Newton iteration
        //    /// </summary>
        //    DirectCoupledIterative,

        //    /// <summary>
        //    /// Direct evaluation from concentration field, only evaluated once per timestep
        //    /// </summary>
        //    DirectCoupledOnce,

        //    /// <summary>
        //    /// No curvature correction
        //    /// </summary>
        //    None
        //}

        //public CurvatureCorrection CurvatureCorrectionType = CurvatureCorrection.None;

        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double penalty_poisson = 2.6;

        /// <summary>
        /// Set if Phasefield constants shall be calculated according to some metric or are fixed
        /// </summary>
        public bool FixedConstants = false;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation, to adjust surface vs bulk diffusion
        /// </summary>
        public double lambda = 0.0;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// diff = (kappa*lambda)/epsilon
        /// </summary>
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double diff = 1.0;

        /// <summary>
        /// Some parameter of Cahn Hilliard equation
        /// Cahn´s Number
        /// </summary>
        [BoSSS.Solution.Control.ExclusiveLowerBound(0.0)]
        public double cahn = 1.0;

    }
}
