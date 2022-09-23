using BoSSS.Solution.Control;
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
        public double LSunderrelax = 1.0;


        /// <summary>
        /// Evolution option for the first level-set (name <see cref="BoSSS.Solution.NSECommon.VariableNames.LevelSetCG"/>)
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

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
                    Option_LevelSetEvolution2 = LevelSetEvolution.None;
                    Timestepper_LevelSetHandling = LevelSetHandling.None;
                }
                base.TimesteppingMode = value;
            }
        }


        /// <summary>
        /// options for additional penalization terms for fast marching
        /// </summary>
        [DataMember]
        public Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms FastMarchingPenaltyTerms = Solution.LevelSetTools.Smoothing.JumpPenalization.jumpPenalizationTerms.Jump;

        /// <summary>
        /// Options for the initialization of the Fourier Level-set
        /// </summary>
        [DataMember]
        public FourierLevSetControl FourierLevSetControl;

        /// <summary>
        /// Options for the initialization of the Phasefield Level-set
        /// </summary>
        [DataMember]
        public PhasefieldControl PhasefieldControl;


        /// <summary>
        /// An explicit expression (y = f(x)) of the initial 0 Level-set. Used for <see cref="SplineLevelSet"/>
        /// </summary>
        [NonSerialized]
        [JsonIgnore]
        public Func<double, double> Phi0Initial;

        /// <summary>
        /// Width of the narrow band.
        /// </summary>
        [DataMember]
        public int LS_TrackerWidth = 1;

        /// <summary>
        /// Reinitilization period for Fastmarching
        /// </summary>
        [DataMember]
        public int FastMarchingReInitPeriod = 0;

        /// <summary>
        /// Reinitilization period for the LevelSetUpdater
        /// </summary>
        [DataMember]
        public int ReInitPeriod = 0;

        /// <summary>
        /// Controls the behavior of the <see cref="ContinuityProjection"/>, i.e. the algorithm which enforces continuity of the level-set
        /// </summary>
        [DataMember]
        public ContinuityProjectionOption LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;


    }


}

