/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using System;

namespace NSE_SIMPLE {

    /// <summary>
    /// Solver configuration for common options
    /// </summary>
    public class SIMPLEControl : AppControl {

        /// <summary>
        /// Version of the SIMPLE algorithm.
        /// </summary>
        public SolutionAlgorithms Algorithm;

        /// <summary>
        /// Physical application modus.
        /// </summary>
        public PhysicsMode PhysicsMode;

        /// <summary>
        /// Reynolds number.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double Reynolds;

        /// <summary>
        /// Option for approximation of predictor.
        /// </summary>
        public PredictorApproximations PredictorApproximation = PredictorApproximations.Identity;//Optional

        /// <summary>
        /// True, if the approximation of the predictor is constant.
        /// </summary>
        public bool PredictorApproximationIsConstant {
            get {
                return (((PredictorApproximation == PredictorApproximations.Identity) || (PredictorApproximation == PredictorApproximations.Identity_IP1))
                    && (Algorithm == SolutionAlgorithms.Steady_SIMPLE));
            }
        }

        /// <summary>
        /// Factor for pressure stabilization in continuity equation
        /// - needed for the equal order formulation.
        /// A factor of 0.0 means no pressure stabilization.
        /// </summary>
        [InclusiveLowerBound(0.0)]
        public double PressureStabilizationScaling = 0.0;

        /// <summary>
        /// Maximum number of SIMPLE steps.
        /// </summary>
        [InclusiveLowerBound(1)]
        public int MaxNoSIMPLEsteps;

        /// <summary>
        /// Convergence criterion for pressure correction.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double L2NormPressureCorrection;

        /// <summary>
        /// Convergence criterion for velocity.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double L2NormVelocityResidual;

        /// <summary>
        /// Under-relaxation factor for velocity.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        [InclusiveUpperBound(1.0)]
        public double RelexationFactorVelocity;

        /// <summary>
        /// Under-relaxation factor for pressure.
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        [InclusiveUpperBound(1.0)]
        public double RelaxationFactorPressure;

        /// <summary>
        /// Reference point for pressure.
        /// Only needed if there is no Dirichlet boundary condition for the pressure.
        /// </summary>
        public double[] PressureReferencePoint = null;

        /// <summary>
        /// Mean value of pressure.
        /// Only needed if there is no Dirichlet boundary condition for the pressure.
        /// </summary>
        public double PressureMeanValue = double.NaN;

        /// <summary>
        /// Source for gradient of the pressure.
        /// E.g. needed with periodic boundary conditions.
        /// </summary>
        public double[] PressureGradientSource = null;

        /// <summary>
        /// Additional scaling of the penalty for the viscous SIP fluxes
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double ViscousPenaltyScaling = 1.0;

        /// <summary>
        /// Additional scaling of the penalty for the pressure SIP fluxes
        /// </summary>
        [ExclusiveLowerBound(0.0)]
        public double PressureCorrectionPenaltyScaling = 1.0;//Optional

        /// <summary>
        /// Linear solver configuration predictor.
        /// </summary>
        [NotNull]
        public Func<ISparseSolver> PredictorSolverFactory;

        /// <summary>
        /// Linear solver configuration corrector.
        /// </summary>
        [NotNull]
        public Func<ISparseSolver> CorrectorSolverFactory;

        /// <summary>
        /// Print results of linear solvers.
        /// </summary>
        public bool PrintLinerSolverResults = false;

        /// <summary>
        /// Analytic solution VelocityX.
        /// </summary>
        public Func<double[], double> AnalyticVelocityX = null;

        /// <summary>
        /// Analytic solution VelocityY.
        /// </summary>
        public Func<double[], double> AnalyticVelocityY = null;

        /// <summary>
        /// Analytic solution VelocityZ.
        /// </summary>
        public Func<double[], double> AnalyticVelocityZ = null;

        /// <summary>
        /// Analytic solution Pressure.
        /// </summary>
        public Func<double[], double> AnalyticPressure = null;

        /// <summary>
        /// Edge tags for calculating drag and lift.
        /// </summary>
        public string[] EdgeTagsDragAndLift = null;

        [InclusiveLowerBound(1)]
        public int SavePeriodSIMPLE;

        [InclusiveLowerBound(1)]
        [InclusiveUpperBound(4)]
        public int TimeOrder = 1;

        /// <summary>
        /// <see cref="PredictorApproximationUpdateCycle"/>
        /// </summary>
        [InclusiveLowerBound(1)]
        public int PredictorApproximationUpdateCycle = int.MaxValue;
    }

    /// <summary>
    /// Options SIMPLE algorithm.
    /// </summary>
    public enum SolutionAlgorithms {

        /// <summary>
        /// Standard steady SIMPLE
        /// </summary>
        Steady_SIMPLE,

        /// <summary>
        /// Standard unsteady SIMPLE using BDF for time discretization.
        /// </summary>
        Unsteady_SIMPLE
    }

    /// <summary>
    /// Approximation \f$ \mathcal A^\vartheta\f$ 
    /// of predictor matrix \f$ \mathbf{A}^\vartheta_{C} - \frac{1}{Re} \mathbf{A}_{D}\f$  in correction step of the SIMPLE-algorithm.
    /// </summary>        
    public enum PredictorApproximations {

        /// <summary>
        /// Have a guess...
        /// </summary>
        Identity,

        /// <summary>
        /// Interior penalty discretization of Laplace operator for pressure correction,
        /// with penalty for pressure correction.
        /// </summary>
        Identity_IP1,

        /// <summary>
        /// Diagonal matrix of predictor.
        /// </summary>
        Diagonal,

        /// <summary>
        /// Block diagonal matrix of predictor.
        /// </summary>
        BlockDiagonal
    }

    /// <summary>
    /// Option for under-relaxation of scalar variable
    /// </summary>
    public enum RelaxationTypes {

        Implicit,

        Explicit
    }
}
