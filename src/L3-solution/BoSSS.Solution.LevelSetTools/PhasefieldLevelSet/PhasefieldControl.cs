using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    public class PhasefieldControl : ILevSetControl
    {

        /// <summary>
        /// ctor
        /// </summary>
        public PhasefieldControl()
        {
            NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            NonLinearSolver.MaxSolverIterations = 50;
        }


        /// <summary>
        /// Model Type of Phasefield equation, see Halperin (1977)
        /// </summary>
        public enum ModelType
        {
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
        /// Set the <see cref="ModelType"/>
        /// </summary>
        public ModelType ModTyp = ModelType.modelA;

        /// <summary>
        /// Type of algebraic correction that is performed
        /// </summary>
        public enum Correction
        {
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

        public Correction CorrectionType = Correction.None;

        /// <summary>
        /// Type of algebraic correction that is performed
        /// </summary>
        public enum CurvatureCorrection
        {
            /// <summary>
            /// Curvature by divergence form
            /// </summary>
            FullyCoupled,

            /// <summary>
            /// Direct evaluation from concentration field, evaluated in each Newton iteration
            /// </summary>
            DirectCoupledIterative,

            /// <summary>
            /// Direct evaluation from concentration field, only evaluated once per timestep
            /// </summary>
            DirectCoupledOnce,

            /// <summary>
            /// No curvature correction
            /// </summary>
            None
        }

        public CurvatureCorrection CurvatureCorrectionType = CurvatureCorrection.None;

        public LinearSolverConfig LinearSolver = new LinearSolverConfig();

        public NonLinearSolverConfig NonLinearSolver = new NonLinearSolverConfig();

    }
}
