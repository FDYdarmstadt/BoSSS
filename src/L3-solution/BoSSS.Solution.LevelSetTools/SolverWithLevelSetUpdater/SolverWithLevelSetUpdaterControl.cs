using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Solution.LevelSetTools.EllipticReInit;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.XdgTimestepping;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {

    /// <summary>
    /// PDE-solver-control object which defines configuration options for nonlinear and linear equation solvers
    /// </summary>
    [Serializable]
    [DataContract]
    public class SolverWithLevelSetUpdaterControl : AppControlSolver {

        /// <summary>
        /// ctor
        /// </summary>
        public SolverWithLevelSetUpdaterControl() {
            Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
        }


        /// <summary>
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        virtual public LevelSetHandling Timestepper_LevelSetHandling {
            get;
            set;
        }

        /// <summary>
        /// underrelaxation of the level set movement in case of coupled iterative
        /// </summary>
        [DataMember]
        public double LSunderrelax = 1.0;


        /// <summary>
        /// Evolution option for the first level-set (name <see cref="BoSSS.Solution.NSECommon.VariableNames.LevelSetCG"/>)
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

        /// <summary>
        /// Evolution option for the second level-set (name <see cref="BoSSS.Solution.NSECommon.VariableNames.LevelSetCGidx"/>(1))
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution2 = LevelSetEvolution.None;

        /// <summary>
        /// getting [<see cref="Option_LevelSetEvolution"/>, <see cref="Option_LevelSetEvolution2"/>][<paramref name="iLs"/>]
        /// </summary>
        public LevelSetEvolution Get_Option_LevelSetEvolution(int iLs) {
            switch(iLs) {
                case 0: return Option_LevelSetEvolution;
                case 1: return Option_LevelSetEvolution2;
                default: throw new IndexOutOfRangeException();
            }    
        }
                

        /// <summary>
        /// setting [<see cref="Option_LevelSetEvolution"/>, <see cref="Option_LevelSetEvolution2"/>][<paramref name="iLs"/>]
        /// </summary>
        public void Set_Option_LevelSetEvolution(int iLs, LevelSetEvolution val) {
            switch(iLs) {
                case 0: Option_LevelSetEvolution = val; break;
                case 1: Option_LevelSetEvolution2 = val; break;
                default: throw new IndexOutOfRangeException();
            }    
        }

        /// <summary>
        /// If set to <see cref="AppControl._TimesteppingMode.Steady"/>, turns of level-set evolution
        /// </summary>
        [JsonIgnore]
        override public _TimesteppingMode TimesteppingMode {
            get {
                return base.TimesteppingMode;
            }
            set {
                if(value == _TimesteppingMode.Steady) {
                    Timestepper_LevelSetHandling = LevelSetHandling.None;
                    Option_LevelSetEvolution = LevelSetEvolution.None;
                    //Option_LevelSetEvolution2 = LevelSetEvolution.None;
                    Timestepper_LevelSetHandling = LevelSetHandling.None;
                }
                base.TimesteppingMode = value;
            }
        }

#region HIGKLY_SPECIFIC_CONFIGURATIONS

        ///// <summary>
        ///// options for additional penalization terms for fast marching
        ///// </summary>
        //[DataMember]
        //public Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

        /// <summary>
        /// Options for the initialization of the Fourier Level-set
        /// </summary>
        [DataMember]
        public FourierLevSetControl FourierLevSetControl;

        /// <summary>
        /// Options for the initialization of the Parameterized Level-set
        /// </summary>
        [DataMember]
        public LevelSetTools.ParameterizedLevelSet.ParameterizedLevelSetControl ParameterizedLevelSetControl;

        /// <summary>
        /// Options for the initialization of the Phasefield Level-set
        /// </summary>
        [DataMember]
        public PhasefieldSettings PhasefieldControl;

        /// <summary>
        /// Control Options for ExtVel
        /// </summary>
        public EllipticExtVelAlgoControl EllipticExtVelAlgoControl = new EllipticExtVelAlgoControl();
        
        /// <summary>
        /// An explicit expression (y = f(x)) of the initial 0 Level-set. Used for <see cref="SplineLevelSet"/>
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double, double> Phi0Initial;

        /// <summary>
        /// Control Options for ReInit
        /// </summary>
        [DataMember]
        public EllipticReInitAlgoControl ReInitControl = new EllipticReInitAlgoControl();

        /// <summary>
        /// if false, Neumann boundary conditions are applied everywhere. If true, the boundary conditions from the BoundaryConditionMap are applied.
        /// </summary>
        [DataMember]
        public StokesExtentionBoundaryOption StokesExtentionUseBCmap = StokesExtentionBoundaryOption.FreeSlipAtWall;

        /// <summary>
        /// If set, the Stokes-extension velocity solve (<see cref="LevelSetEvolution.StokesExtension"/>) is
        /// performed through PETSc using this options string, instead of the default direct-solver fallback
        /// chain. Null/empty (the default) preserves the original behavior.
        /// </summary>
        [DataMember]
        public string StokesExtensionPetscOptions = null;

        /// <summary>
        /// If true, the Stokes-extension velocity solve (<see cref="LevelSetEvolution.StokesExtension"/>)
        /// is warm-started from the previous time step's converged extension velocity. False (the default)
        /// always starts the solve from zero, matching the original behavior.
        /// </summary>
        [DataMember]
        public bool StokesExtensionWarmStart = false;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionPetscOptions"/> is set: forwarded to
        /// <see cref="BoSSS.Solution.AdvancedSolvers.PetscGMRESConfig.UseGeometricMultigrid"/> for
        /// the Stokes-extension velocity solve. False (the default) preserves the original
        /// flat-matrix behavior.
        /// </summary>
        [DataMember]
        public bool StokesExtensionUseGeometricMultigrid = false;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionUseGeometricMultigrid"/> is true: forwarded
        /// to <see cref="BoSSS.Solution.AdvancedSolvers.PetscGMRESConfig.CoarseKickIn"/>.
        /// </summary>
        [DataMember]
        public int StokesExtensionCoarseKickIn = 90000;

        /// <summary>
        /// If true, the Stokes-extension velocity solve's (<see cref="LevelSetEvolution.StokesExtension"/>)
        /// velocity and pressure <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.ChangeOfBasisConfig"/>
        /// are combined into a single <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.SchurComplement"/>
        /// block, mirroring the <c>UseSchurBlockPrec</c> option used for the main flow solve (e.g.
        /// <c>XNSE_Control.UseSchurBlockPrec</c>). False (the default) preserves the original per-field
        /// configuration.
        /// </summary>
        [DataMember]
        public bool StokesExtensionUseSchurBlockPrec = false;

        /// <summary>
        /// If true, the Stokes-extension velocity solve (<see cref="LevelSetEvolution.StokesExtension"/>) is
        /// solved via pseudo-transient continuation (repeated solves of a progressively less-regularized
        /// system, switched-evolution-relaxation) instead of a single solve. False (the default) preserves
        /// the original single-solve behavior. See <c>StokesExtension.UsePseudoTimestepping</c>.
        /// </summary>
        [DataMember]
        public bool StokesExtensionUsePseudoTimestepping = false;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionUsePseudoTimestepping"/> is set: initial pseudo-time
        /// step size.
        /// </summary>
        [DataMember]
        public double StokesExtensionPseudoDtau0 = 1.0;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionUsePseudoTimestepping"/> is set: upper bound for the
        /// pseudo-time step size's growth.
        /// </summary>
        [DataMember]
        public double StokesExtensionPseudoDtauMax = 1.0e8;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionUsePseudoTimestepping"/> is set: maximum number of
        /// pseudo-time steps.
        /// </summary>
        [DataMember]
        public int StokesExtensionMaxPseudoTimeSteps = 20;

        /// <summary>
        /// Only relevant if <see cref="StokesExtensionUsePseudoTimestepping"/> is set: convergence criterion
        /// on the true (un-regularized) residual norm, RELATIVE to ‖RHS‖.
        /// </summary>
        [DataMember]
        public double StokesExtensionPseudoTimeConvergenceCriterion = 1.0e-10;

        /// <summary>
        /// If true, the Stokes-extension velocity solve (<see cref="LevelSetEvolution.StokesExtension"/>) skips
        /// BoSSS's own local preconditioning transform entirely (neither <see cref="StokesExtensionUseSchurBlockPrec"/>'s
        /// combined block nor the standard per-field configs are used) -- every field gets
        /// <see cref="BoSSS.Solution.AdvancedSolvers.MultigridOperator.Mode.Eye"/> (identity), handing the raw
        /// system matrix directly to PETSc/the fallback solver. See <c>StokesExtension.DisableLocalPreconditioning</c>.
        /// Default (false) preserves the original behavior.
        /// </summary>
        [DataMember]
        public bool StokesExtensionDisableLocalPreconditioning = false;


        #endregion

        /// <summary>
        /// Width of the narrow band.
        /// </summary>
        [DataMember]
        public int LS_TrackerWidth = 1;



        /// <summary>
        /// Reinitilization period for the Level Set updater,
        /// if <see cref="Option_LevelSetEvolution"/> or <see cref="Option_LevelSetEvolution2"/>
        /// is set to 
        /// <see cref="LevelSetEvolution.StokesExtension"/> or <see cref="LevelSetEvolution.FastMarching"/>.
        /// </summary>
        [DataMember]
        public int ReInitPeriod = 0;

        /// <summary>
        /// in case of a restart set for internal timestep number within the ls evolver
        /// </summary>
        [DataMember]
        public int ReInitTimestepIndex = 0;


        /// <summary>
        /// Controls the behavior of the <see cref="ContinuityProjection"/>, i.e. the algorithm which enforces continuity of the level-set
        /// </summary>
        [DataMember]
        public ContinuityProjectionOption LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

        /// <summary>
        /// Controls the behavior of the <see cref="ContinuityProjection"/>, i.e. the algorithm which enforces continuity of the second level-set 
        /// </summary>
        [DataMember]
        public ContinuityProjectionOption LSContiProjectionMethod2 = ContinuityProjectionOption.ConstrainedDG;
    }


    /// <summary>
    /// Configuration options for the Stokes extension;
    /// This is (also) intended to be a transitional development setting which will be removed from the code when a final conclusion on boundary conditions for the Stokes extension is reached.
    /// </summary>
    public enum StokesExtentionBoundaryOption {
        
        /// <summary>
        /// Zero flux, old option, used mainly by the FSI-Codes (chen Miao, Lauritz Beck)
        /// </summary>
        ZeroFlux,


        /// <summary>
        /// 
        /// </summary>
        FreeSlipAtWall,


        /// <summary>
        /// apply the normal boundary conditions from the physical code
        /// </summary>
        useBcMap
    }
}

