using System;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Utils;

namespace CNS.Solution {

    public class NewtonSolver : INonlinearSystemSolver {

        private Configuration solverConfiguration;

        private uint updateInterval;

        private double relaxationFactor;

        private MsrMatrix currentMatrix;

        private double[] currentAffineOffset;

        private ISparseSolver currentSolver;

        public NewtonSolver(SpatialOperator originalOperator, Configuration solverConfiguration, CoordinateMapping coordinateMapping, uint updateInterval = 1, double relaxationFactor = 1.0) {
            CoordinateMapping = coordinateMapping;
            DGCoordinates = new CoordinateVector(coordinateMapping);
            this.solverConfiguration = solverConfiguration;
            this.updateInterval = updateInterval;
            this.relaxationFactor = relaxationFactor;

            if (originalOperator.ParameterVar.Count > 0) {
                throw new NotImplementedException(
                    "Operators with parameters are currently not supported");
            }

            Operator = new SpatialOperator(
                originalOperator.DomainVar, originalOperator.DomainVar.Select(x => x + "0").ToList(), originalOperator.CodomainVar);
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

        public int CurrentIteration {
            get;
            private set;
        }

        public CoordinateVector DGCoordinates {
            get;
            private set;
        }

        public CoordinateMapping CoordinateMapping {
            get;
            private set;
        }

        public SpatialOperator Operator {
            get;
            private set;
        }

        public void PerformIteration() {
            if (CurrentIteration % updateInterval == 0) {
                UpdateLinearization();
            }

            double[] rhs = new double[CoordinateMapping.NUpdate];
            currentMatrix.SpMV(1.0 - relaxationFactor, DGCoordinates, 0.0, rhs);
            BLAS.daxpy(CoordinateMapping.NUpdate, -relaxationFactor, currentAffineOffset, 1, rhs, 1);

            currentSolver.Solve(DGCoordinates, rhs);

            CurrentIteration++;
        }

        public void ReInitialize() {
            UpdateLinearization();
            CurrentIteration = 0;
        }

        #endregion

        #region IDisposable Members

        public void Dispose() {
            if (currentSolver != null) {
                currentSolver.Dispose();
            }
        }

        #endregion

        private void UpdateLinearization() {
            if (currentSolver != null) {
                currentSolver.Dispose();
            }

            currentMatrix = new MsrMatrix(CoordinateMapping.NUpdate);
            currentAffineOffset = new double[CoordinateMapping.NUpdate];
            Operator.ComputeMatrixEx(
                CoordinateMapping,
                CoordinateMapping,
                CoordinateMapping,
                currentMatrix,
                currentAffineOffset);

            currentSolver = SolverFactory.CreateSolver<ISparseSolver>(solverConfiguration);
            currentSolver.DefineMatrix(currentMatrix);
        }
    }
}
