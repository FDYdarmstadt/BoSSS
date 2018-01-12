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
using System.Collections.Generic;

using BoSSS.Foundation;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {
    
    /// <summary>
    /// the well-known implicit euler scheme
    /// </summary>
    public class ImplicitEuler : ImplicitTimeStepper {

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, params DGField[] fields)
            : this(solver, spatialOpMtx, spatialOpAffine, new CoordinateMapping(fields)) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, params DGField[] fields)
            : this(solver, AllTrue(fields.Length), spatialOpMtx, spatialOpAffine, new CoordinateMapping(fields)) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields)
            : this(solver, AllTrue(fields.Fields.Count), spatialOpMtx, spatialOpAffine, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields)
            : base(solver, temporalOp, spatialOpMtx, spatialOpAffine, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, SpatialOperator spatialOp, CoordinateMapping fields)
            : this(solver, AllTrue(fields.Fields.Count), spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, bool[] temporalOp, SpatialOperator spatialOp, CoordinateMapping fields)
            : base(solver, temporalOp, spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, SpatialOperator spatialOp, params DGField[] fields)
            : this(solver, AllTrue(fields.Length), spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public ImplicitEuler(ISparseSolverExt solver, bool[] temporalOp, SpatialOperator spatialOp, params DGField[] fields)
            : this(solver, temporalOp, spatialOp, new CoordinateMapping(fields)) { }

        /// <summary>
        /// Create an empty scheme
        /// </summary>
        public ImplicitEuler() {}

        /// <summary>
        /// Solves the linear system (diag(1/<paramref name="dt"/>) +
        /// <em>M</em>) * x =
        /// <see cref="ImplicitTimeStepper.CurrentState"/> /
        /// <paramref name="dt"/> -
        /// <see cref="ImplicitTimeStepper.m_AffineOffset1"/> and writes the
        /// result to <see cref="ImplicitTimeStepper.CurrentState"/>.
        /// </summary>
        /// <param name="dt">The length of the timestep</param>
        protected override void PerformTimeStep(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {


                int n = Mapping.LocalLength;
                int np = m_diagVecOneSec.Length;

                double[] diag = new double[n];
                for (int i = 0; i < n; i++) {
                    diag[i] = m_diagVecOneSec[i % np];
                }
                BLAS.dscal(n, 1.0 / dt, diag, 1);

                double[] rhs = (double[])m_AffineOffset1.Clone();
                BLAS.dscal(n, -1.0, rhs, 1);

                for (int i = 0; i < n; i++) {
                    rhs[i] += diag[i] * CurrentState[i];
                }

                tr.Info("Calling solver");
                LastSolverResult = m_Solver.Solve<double[], CoordinateVector, double[]>(1.0, diag, CurrentState, rhs);

            }
        }

        /// <summary>
        /// contains the sparse solver statistics acquired in the last call to <see cref="ImplicitTimeStepper.Perform"/>
        /// </summary>
        public SolverResult LastSolverResult;
    }
}