using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Solution.Control.AppControl;

namespace FreeXNSE {
    public class FreeXNSE_Control : SolverWithLevelSetUpdaterControl {

        public FreeXNSE_Control() : this(false) {
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        public FreeXNSE_Control(bool equal) {

            //base.CutCellQuadratureType = CutCellQuadratureMethod.OneStepGaussAndStokes;
            //shift of Solver Information
            base.NoOfMultigridLevels = 1;
            base.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig(); //LinearSolver
            base.NonLinearSolver.MaxSolverIterations = 10; //Solver_MaxIterations
            base.NonLinearSolver.MinSolverIterations = 2; //Solver_MinIterations
            base.NonLinearSolver.ConvergenceCriterion = 0.0; //Solver_ConvergenceCriterion: solve as accurate as possible. Don't change this, Grüße von FK!
            base.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton; //NonLinearSolver
            base.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            EqualOrder = equal;

            this.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
            this.Timestepper_LevelSetHandling = BoSSS.Solution.XdgTimestepping.LevelSetHandling.LieSplitting;
        }

        /// <summary>
        /// Equal order discretization pressure/velocity
        /// </summary>
        [DataMember]
        public bool EqualOrder {
            get;
            private set;
        }

        ParameterFunctionPair[] m_ParameterLevelSet;
        /// <summary>
        /// Values for ParameterLevelSet
        /// </summary>
        [DataMember]
        virtual public ParameterFunctionPair[] ParameterLevelSet {
            get { return m_ParameterLevelSet; }
            set {
                if(m_ParameterLevelSet == null) {
                    this.Option_LevelSetEvolution = BoSSS.Solution.LevelSetTools.LevelSetEvolution.CustomLevelSet;
                }
                m_ParameterLevelSet = value;
            }
        }

        [DataMember]
        public Tuple<double[], double[]>[] DualSplinePhi0Initial;

        public Func<double, double[]> SemiCircleSplinePhi0Initial;

        [DataMember]
        public int NoOfNodes { get; set; } = 30;

        bool m_UseImmersedBoundary;
        /// <summary>
        /// Activation of second level-set (fluid/solid boundary)
        /// </summary>
        [DataMember]
        virtual public bool UseImmersedBoundary {
            get { return m_UseImmersedBoundary; }
            set {
                if(value == true) {
                    throw new NotImplementedException();
                    this.Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.StokesExtension;
                } else {
                    this.Option_LevelSetEvolution2 = BoSSS.Solution.LevelSetTools.LevelSetEvolution.None;
                }
                m_UseImmersedBoundary = value;
            }
        }

        /// <summary>
        /// Artificially fix the interface psoition, used for investigation of asymptotic solutions
        /// </summary>
        [DataMember]
        public bool FixedInterface {
            get;
            set;
        } = false;

        /// <summary>
        /// Dependency of the slip length on level-set (approximate interface distance)
        /// \beta(\phi) = \beta * \phi^SlipScaling
        /// <see cref="DimensionlessNumbers.beta"/>>
        /// </summary>
        [DataMember]
        public double SlipScaling {
            get;
            set;
        } = 0.0;

        /// <summary>
        /// Dependency of the viscosity on distance to contact point
        /// \mu = 1/Re * r^ViscosityScaling
        /// <see cref="DimensionlessNumbers.Re"/>>
        /// </summary>
        [DataMember]
        public double ViscosityScaling {
            get;
            set;
        } = 0.0;

        /// <summary>
        /// Dependency of the deviation in contact angle to that in contactline velocity
        /// (1/We (cos(\theta) - cos(\theta_0)) = \alpha * sign((u - uB) \cdot \normal_L) * |((u - uB) \cdot \normal_L)|^ContactAngleScaling)
        /// </summary>
        [DataMember]
        public double ContactAngleScaling {
            get;
            set;
        } = 1.0;

        /// <summary>
        /// Type of <see cref="XNSE"/>.
        /// </summary>
        public override Type GetSolverType() {
            return typeof(FreeXNSE);
        }

        [DataMember]
        public int Degree {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        public override void SetDGdegree(int p) {
            if(p < 1 && !EqualOrder)
                throw new ArgumentOutOfRangeException("Pressure degree must be 0 at minimum.");

            FieldOptions.Add("Velocity*", new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.Pressure, new FieldOpts() {
                Degree = EqualOrder ? p : p - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetDG, new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCG, new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetDGidx(1), new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            FieldOptions.Add(VariableNames.LevelSetCGidx(1), new FieldOpts() {
                Degree = p,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            Degree = p;
        }

        /// <summary>
        /// Viscosity, density and surface tension.
        /// </summary>
        [DataMember]
        public DimensionlessNumbers DimensionlessNumbers = new DimensionlessNumbers() {
            Oh = 1.0
        };

        /// <summary>
        /// Viscosity, density and surface tension.
        /// </summary>
        [DataMember]
        public ActiveTerms ActiveTerms = new ActiveTerms() {
            VelocityDivergence = VelocityDivergence.Central,
            Temporal = Temporal.On,
            Convective = Convective.LaxFriedrich,
            PressureGradient = PressureGradient.Central,
            Viscous = Viscous.SIP,
            SurfaceTension = SurfaceTension.LaplaceBeltrami
        };



        [DataMember]
        public SolverSettings SolverSettings = new SolverSettings();

        /// <summary>
        /// Time dependent (component-wise) gravitational acceleration (either A or B).
        /// </summary>
        [DataMember]
        public Forcing? VolumeForce {
            get;
            set;
        }

        /// <summary>
        /// Configuring <see cref="AppControl._TimesteppingMode.Steady"/> sets the <see cref="TimeSteppingScheme.ImplicitEuler"/>
        /// </summary>
        [JsonIgnore]
        public override _TimesteppingMode TimesteppingMode {
            get {
                return base.TimesteppingMode;
            }
            set {
                base.TimesteppingMode = value;
                if(value == _TimesteppingMode.Steady)
                    this.TimeSteppingScheme = BoSSS.Solution.XdgTimestepping.TimeSteppingScheme.ImplicitEuler;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            return base.GetHashCode();
        }

        /// <summary>
        /// 
        /// </summary>
        public override bool Equals(object obj) {
            //System.Diagnostics. dbg_launch();
            if(!base.Equals(obj))
                return false;

            var other = obj as FreeXNSE_Control;
            if(other == null)
                return false;

            if(!this.SolverSettings.Equals(other.SolverSettings))
                return false;
            if(!this.DimensionlessNumbers.Equals(other.DimensionlessNumbers))
                return false;
            if(!this.EqualOrder.Equals(other.EqualOrder))
                return false;
            if(!this.VolumeForce.Equals(other.VolumeForce))
                return false;

            return true;
        }

    }

    /// <summary>
    /// Forcing term on the right hand side, requires an actual implementation for serialization
    /// </summary>
    public abstract class Forcing {
        public abstract double[] Evaluate(double[] X, double t);

        public abstract bool equals(Forcing other);

        public override bool Equals(object obj) {

            var other = obj as Forcing;
            if(other == null)
                return false;

            if(!this.equals(other))
                return false;

            return true;
        }
    }

    public class SolverSettings {
        
        /// <summary>
        /// Only for debugging purpose:
        /// solver is turned of and residual of initial value/exact solution is evaluated, used to 
        /// test the consistency of the implementation.
        /// </summary>
        [DataMember]
        public bool SkipSolveAndEvaluateResidual = false;

        /// <summary>
        /// Terminates the simulation if the linear or nonlinear solver fails to converge
        /// </summary>
        [DataMember]
        public bool FailOnSolverFail = true;

        public override bool Equals(object obj) {
            var other = obj as SolverSettings;
            if(other == null)
                return false;

            if(!this.FailOnSolverFail.Equals(other.FailOnSolverFail))
                return false;
            if(!this.SkipSolveAndEvaluateResidual.Equals(other.SkipSolveAndEvaluateResidual))
                return false;

            return true;
        }
    }
    public class DimensionlessNumbers {

        /// <summary>
        /// Reynolds number
        /// </summary>
        [DataMember]
        public double Re {
            get;
            private set;
        }

        /// <summary>
        /// Capillary number
        /// </summary>
        [DataMember]
        public double Ca {
            get;
            private set;
        }

        /// <summary>
        /// Laplace number
        /// </summary>
        [DataMember]
        public double La {
            get;
            private set;
        }

        /// <summary>
        /// Weber number
        /// </summary>
        public double We {
            get;
            private set;
        }

        /// <summary>
        /// Ohnesorge number - THIS IS THE MAIN SETTING FOR THIS SOLVER!!!
        /// For Oh < 1 a capillary non-dimenzionalization is chosen
        /// For Oh > 1 a viscous non-dimenzionalization is chosen
        /// </summary>
        [DataMember]
        public double Oh {
            get { return Math.Sqrt(1.0/this.La); }
            set {
                if(value < 0.0) throw new ArgumentException();
                if(value > 1.0) {
                    this.Re = 1.0;
                    this.La = 1.0 / (value * value);
                    this.Ca = 1.0 / this.La;
                    this.We = this.Re * this.Ca;
                } else if(value <= 1.0) {
                    this.La = 1.0 / (value * value);
                    this.Ca = value;
                    this.Re = 1.0 / value;
                    this.We = 1.0;
                } else {
                    throw new ArgumentException();
                }
            }
        }

        /// <summary>
        /// Froude number
        /// </summary>
        public double Fr {
            get;
            set;
        }

        /// <summary>
        /// Eötvos number
        /// </summary>
        public double Eo {
            get { return We / Fr; }
        }

        /// <summary>
        /// Bond number
        /// </summary>
        public double Bo {
            get { return Eo; }
        }

        /// <summary>
        /// Contactangle (static) number
        /// </summary>
        [DataMember]
        public double Theta {
            get;
            set;
        } = Math.PI/ 2;

        /// <summary>
        /// Contactangle (static advancing) number
        /// </summary>
        [DataMember]
        public double ThetaAdv {
            get;
            set;
        } = 0;

        /// <summary>
        /// Contactangle (static receding) number
        /// </summary>
        [DataMember]
        public double ThetaRec {
            get;
            set;
        } = Math.PI;

        /// <summary>
        /// Frictionparameter (base) (1/Re \tensr{P} \nabla u \normalB = \beta * f(\phi) * \tensr{P} (u - uB))
        /// Set 
        /// 0.0 for freeslip,
        /// Negative or double.PositiveInfinity for No-Slip
        /// </summary>
        [DataMember]
        public double beta {
            get;
            set;
        }

        /// <summary>
        /// Frictionparameter for the contactline (base) (1/We (cos(\theta) - cos(\theta_0)) = \alpha * (u - uB) \cdot \normal_L)
        /// Set 
        /// 0.0 for quasi-static contact angle,
        /// Negative or double.PositiveInfinity for free contactangle
        /// </summary>
        [DataMember]
        public double alpha {
            get;
            set;
        }

        public bool Equals(object obj) {
            var other = obj as DimensionlessNumbers;
            if(other == null)
                return false;

            if(!this.Re.Equals(other.Re))
                return false;
            if(!this.Ca.Equals(other.Ca))
                return false;
            if(!this.La.Equals(other.La))
                return false;
            if(!this.Oh.Equals(other.Oh))
                return false;
            if(!this.Theta.Equals(other.Theta))
                return false;
            if(!this.beta.Equals(other.beta))
                return false;
            if(!this.alpha.Equals(other.alpha))
                return false;

            return true;
        }

        public static bool IsPhysical(double N) {
            return (N > 0 && N.IsFinite());
        }

    }


    public enum Viscous {
        Off,
        SIP,
        NIP,
        IIP,
        OBB
    }
    public enum Temporal {
        Off,
        On
    }
    public enum Convective {
        Off,
        LaxFriedrich,
        Temam,
        ConservativeTemam
    }
    public enum PressureGradient {
        Off,
        Central
    }
    public enum VelocityDivergence {
        Off,
        Central
    }
    public enum SurfaceTension {
        Off,
        LaplaceBeltrami,
        LaplaceBeltrami_BoussinesqScriven
    }
    public class ActiveTerms {

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Temporal Temporal {
            get;
            set;
        }
        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Convective Convective {
            get;
            set;
        }
        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public VelocityDivergence VelocityDivergence {
            get;
            set;
        }
        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public PressureGradient PressureGradient {
            get;
            set;
        }

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Viscous Viscous {
            get;
            set;
        }

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public SurfaceTension SurfaceTension {
            get;
            set;
        }

        public override bool Equals(object obj) {
            var other = obj as ActiveTerms;
            if(other == null)
                return false;

            if(!this.Temporal.Equals(other.Temporal))
                return false;
            if(!this.Convective.Equals(other.Convective))
                return false;
            if(!this.VelocityDivergence.Equals(other.VelocityDivergence))
                return false;
            if(!this.PressureGradient.Equals(other.PressureGradient))
                return false;
            if(!this.Viscous.Equals(other.Viscous))
                return false;
            if(!this.SurfaceTension.Equals(other.SurfaceTension))
                return false;

            return true;
        }

    }

}
