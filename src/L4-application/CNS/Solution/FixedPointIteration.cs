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
    /// \f$ 
    /// \vec{F}(\vec{U}) = 0,
    /// \f$ 
    /// this class computed approximate solutions via
    /// \f$ 
    /// \vec{U}_{n+1} = (1 - \alpha) \vec{U}_n + \alpha \vec{U}*
    /// \f$ 
    /// where \f$ \alpha\f$  is damping factor. The
    /// intermediate solution \f$ \vec{U}*\f$  is the
    /// solution of the linear system
    /// \f$ 
    /// M(\vec{U}_n) \vec{U*} = \vec{b}(\vec{U}_n)
    /// \f$ 
    /// which is computed from the Taylor expansion of
    /// \f$ \vec{F}(\vec{U})\f$  w.r.t.
    /// \f$ \vec{U}\f$  around
    /// \f$ \vec{U}_n\f$ . As a result, we have
    /// \f$ 
    /// M(\vec{U}_n) = \nabla_\vec{U} \vec{F}(\vec{U}_n)
    /// \f$ 
    /// and
    /// \f$ 
    /// \vec{b}(\vec{U}_n) = \nabla_{\vec{U}} \vec{F}(\vec{U}_n) \vec{U}_n - \vec{F}(\vec{U}_n)
    /// \f$.
    /// </remarks>
    public class FixedPointIteration : IIterativeImplicitScheme {

        /// <summary>
        /// Configuration for the linear solver that ought to be used in every
        /// fixed-point iteration
        /// </summary>
        private ISparseSolver solver;

        /// <summary>
        /// Damping factor (or under-relaxation factor) of the fixed-point
        /// iteration. Should be between 0.0 (exclusive) and 1.0.
        /// </summary>
        private double dampingFactor;

        /// <summary>
        /// The current approximation of
        /// \f$ M(\vec{U})\f$ .
        /// </summary>
        private MsrMatrix currentMatrix;

        /// <summary>
        /// The current approximation of
        /// \f$ \vec{b}(\vec{U})\f$ 
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
        /// Parameter mapping for the linearized operator; contains the
        /// parameters of the original, non-linear operator plus the current
        /// values of each field in <see cref="Mapping"/>.
        /// </summary>
        private CoordinateMapping parameterMapping;

        /// <summary>
        /// Initializes a new instance of the
        /// <see cref="FixedPointIteration"/> class.
        /// </summary>
        /// <param name="originalOperator">
        /// The original operator \f$ \vec{F}\f$ 
        /// </param>
        /// <param name="solver">
        /// The solver configuration
        /// </param>
        /// <param name="coordinateMapping">
        /// The coordinate mapping
        /// </param>
        /// <param name="parameterMapping">
        /// The mapping of parameter fields.
        /// </param>
        /// <param name="dampingFactor">
        /// Damping factor (or under-relaxation factor) of the fixed-point
        /// iteration. Should be between 0.0 (exclusive) and 1.0.
        /// </param>
        public FixedPointIteration(SpatialOperator originalOperator, ISparseSolver solver, CoordinateMapping coordinateMapping, CoordinateMapping parameterMapping = null, double dampingFactor = 1.0) {
            Mapping = coordinateMapping;
            if (parameterMapping == null) {
                this.parameterMapping = Mapping;
            } else {
                this.parameterMapping = new CoordinateMapping(
                    coordinateMapping.Fields.Concat(parameterMapping.Fields).ToArray());
            }
            DGCoordinates = new CoordinateVector(coordinateMapping);
            this.solver = solver;
            this.dampingFactor = dampingFactor;

            Operator = new SpatialOperator(
                originalOperator.DomainVar,
                originalOperator.DomainVar.Select(x => x + "0").Concat(
                    originalOperator.ParameterVar).ToList(),
                originalOperator.CodomainVar, QuadOrderFunc.MaxDegTimesTwo());
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
        /// Performs a single fixed-point iteration
        /// </summary>
        public void PerformIteration() {
            UpdateLinearization();

            double[] rhs = new double[Mapping.LocalLength];
            if (dampingFactor != 1.0) {
                // Should be replaced by something real at some point (or at
                // least by a check if #CPUs = 1)
                currentMatrix.SpMVpara(1.0 - dampingFactor, DGCoordinates, 0.0, rhs);
            }
            // rhs = -factor * currentOffset + rhs 
            BLAS.daxpy(Mapping.LocalLength, -dampingFactor, currentAffineOffset, 1, rhs, 1);

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
        public double Perform(double dt) {
            PerformIteration();
            return dt;
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

            currentMatrix = new MsrMatrix(Mapping.LocalLength);
            currentAffineOffset = new double[Mapping.LocalLength];
            Operator.ComputeMatrixEx(
                Mapping,
                parameterMapping,
                Mapping,
                currentMatrix,
                currentAffineOffset);

            currentSolver = solver;
            currentSolver.DefineMatrix(currentMatrix);
        }
    }
}
