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
        }


        /// <summary>
        /// See <see cref="LevelSetHandling"/>
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

        /// <summary>
        /// underrelaxation of the level set movement in case of coupled iterative
        /// </summary>
        public double LSunderrelax = 1.0;


        /// <summary>
        /// See <see cref="LevelSetEvolution"/>.
        /// </summary>
        [DataMember]
        public LevelSetEvolution Option_LevelSetEvolution = LevelSetEvolution.FastMarching;


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
    }


}

