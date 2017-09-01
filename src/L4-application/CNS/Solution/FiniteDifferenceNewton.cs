using System;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Utils;

namespace CNS.Solution {

    /// <summary>
    /// A solver for nonlinear systems of equation based on an approximation of
    /// the gradient of a nonlinear system by means of a forward difference
    /// approximation of the gradients of the nonlinear fluxes
    /// (see <see cref="LinearizedFlux" />)
    /// </summary>
    /// <remarks>
    /// For a generic nonlinear system of equations
    /// <latex>
    /// \vec{F}(\vec{U}) = 0,
    /// </latex>
    /// this class computed approximate solutions via
    /// <latex>
    /// \vec{U}_{n+1} = (1 - \alpha) \vec{U}_n + \alpha \vec{U}*
    /// </latex>
    /// where <latex mode="inline">\alpha</latex> is damping factor. The
    /// intermediate solution <latex mode="inline">\vec{U}*</latex> is the
    /// solution of the linear system
    /// <latex>
    /// M(\vec{U}_n) \vec{U*} = \vec{b}(\vec{U}_n)
    /// </latex>
    /// which is computed from the Taylor expansion of
    /// <latex mode="inline">\vec{F}(\vec{U})</latex> w.r.t.
    /// <latex mode="inline">\vec{U}</latex> around
    /// <latex mode="inline">\vec{U}_n</latex>. As a result, we have
    /// <latex>
    /// M(\vec{U}_n) = \nabla_\vec{U} \vec{F}(\vec{U}_n)
    /// </latex>
    /// and
    /// <latex>
    /// \vec{b}(\vec{U}_n) = \nabla_\vec{U} \vec{F}(\vec{U}_n) \vec{U}_n - \vec{F}(\vec{U}_n).
    /// </latex>
    /// </remarks>
    public class FiniteDifferenceNewton : IIterativeImplicitScheme {

        /// <summary>
        /// Configuration for the linear solver that ought to be used in every
        /// Newton iteration
        /// </summary>
        private Configuration solverConfiguration;

        /// <summary>
        /// Damping factor (or under-relaxation factor) of the Newton
        /// iteration. Should be between 0.0 (exclusive) and 1.0.
        /// </summary>
        private double dampingFactor;

        /// <summary>
        /// The current approximation of
        /// <latex mode="inline">M(\vec{U})</latex>.
        /// </summary>
        private MsrMatrix currentMatrix;

        /// <summary>
        /// The current approximation of
        /// <latex mode="inline">\vec{b}(\vec{U})</latex>
        /// </summary>
        private double[] currentAffineOffset;

        /// <summary>
        /// The current instance the solver for the linear equation system
        /// </summary>
        /// <remarks>
        /// Has to be recreated in every iteration.
        /// </remarks>
        private ISparseSolver currentSolver;

        /// <summary>
        /// Initializes a new instance of the
        /// <see cref="FiniteDifferenceNewton"/> class.
        /// </summary>
        /// <param name="originalOperator">
        /// The original operator <latex mode="inline">\vec{F}</latex>
        /// </param>
        /// <param name="solverConfiguration">
        /// The solver configuration
        /// </param>
        /// <param name="coordinateMapping">
        /// The coordinate mapping
        /// </param>
        /// <param name="dampingFactor">
        /// Damping factor (or under-relaxation factor) of the Newton
        /// iteration. Should be between 0.0 (exclusive) and 1.0.
        /// </param>
        public FiniteDifferenceNewton(SpatialOperator originalOperator, Configuration solverConfiguration, CoordinateMapping coordinateMapping, double dampingFactor = 1.0) {
            Mapping = coordinateMapping;
            DGCoordinates = new CoordinateVector(coordinateMapping);
            this.solverConfiguration = solverConfiguration;
            this.dampingFactor = dampingFactor;

            if (originalOperator.ParameterVar.Count > 0) {
                throw new NotImplementedException(
                    "Operators with parameters are currently not supported");
            }

            Operator = new SpatialOperator(
                originalOperator.DomainVar,
                originalOperator.DomainVar.Select(x => x + "0").ToList(),
                originalOperator.CodomainVar);
            foreach (var component in originalOperator.EquationComponents) {
                var linearizableFluxes = component.Value.OfType<INonlinearFlux>();

                foreach (var flux in linearizableFluxes) {
                    Operator.EquationComponents[component.Key].Add(
                        new LinearizedFlux(flux));
                }

                foreach (var flux in component.Value.Except(linearizableFluxes)) {
                    Operator.EquationComponents[component.Key].Add(flux);
                }
            }

            Operator.Commit();

            if (Operator.ContainsNonlinear) {
                throw new ArgumentException(
                    "The given operator contains nonlinear fluxes that could not"
                        + " be linearized automatically. This feature is currently"
                        + " only supported for fluxes of type INonlinearFlux.",
                    "originalOperator");
            }
        }

        #region INonlinearSystemSolver Members

        /// <summary>
        /// Gets the current iteration.
        /// </summary>
        public int CurrentIteration {
            get;
            private set;
        }

        /// <summary>
        /// Gets the operator which represents a linearization of the original
        /// operator supplied in the constructor
        /// </summary>
        public SpatialOperator Operator {
            get;
            private set;
        }

        /// <summary>
        /// Performs a single Newton iteration
        /// </summary>
        public void PerformIteration() {
            UpdateLinearization();

            double[] rhs = new double[Mapping.NUpdate];
            currentMatrix.SpMV(1.0 - dampingFactor, DGCoordinates, 0.0, rhs);
            BLAS.daxpy(Mapping.NUpdate, -dampingFactor, currentAffineOffset, 1, rhs, 1);

            currentSolver.Solve(DGCoordinates, rhs);

            CurrentIteration++;
        }

        #endregion

        #region ITimeStepper Members

        /// <summary>
        /// Always returns zero.
        /// </summary>
        public double Time {
            get {
                return 0.0;
            }
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="NewTime">Irrelevant</param>
        public void ResetTime(double NewTime) {
        }

        /// <summary>
        /// Currently, just performs an iteration of the scheme. Once this
        /// solver is extended to unsteady equations, this will change.
        /// </summary>
        /// <param name="dt">Time-step size</param>
        public void Perform(double dt) {
            PerformIteration();
        }

        /// <summary>
        /// The affected DG fields.
        /// </summary>
        public CoordinateMapping Mapping {
            get;
            private set;
        }

        /// <summary>
        /// The affected DG coordinates.
        /// </summary>
        public CoordinateVector DGCoordinates {
            get;
            private set;
        }

        #endregion

        #region IDisposable Members

        /// <summary>
        /// Disposes the current linear solver
        /// </summary>
        public void Dispose() {
            if (currentSolver != null) {
                currentSolver.Dispose();
            }
            GC.SuppressFinalize(this);
        }

        #endregion

        /// <summary>
        /// Updates all components of the linearization, i.e. the solver,
        /// the matrices and the affine offset.
        /// </summary>
        private void UpdateLinearization() {
            if (currentSolver != null) {
                currentSolver.Dispose();
            }

            currentMatrix = new MsrMatrix(Mapping.NUpdate);
            currentAffineOffset = new double[Mapping.NUpdate];
            Operator.ComputeMatrixEx(
                Mapping,
                Mapping,
                Mapping,
                currentMatrix,
                currentAffineOffset);

            currentSolver = SolverFactory.CreateSolver<ISparseSolver>(solverConfiguration);
            currentSolver.DefineMatrix(currentMatrix);
        }
    }
}
