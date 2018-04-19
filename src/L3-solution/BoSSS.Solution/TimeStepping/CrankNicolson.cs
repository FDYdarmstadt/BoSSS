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
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP.LinSolvers;

namespace BoSSS.Solution.Timestepping {
    
    /// <summary>
    /// The Crank-Nicolson scheme (see e.g.
    /// http://de.wikipedia.org/wiki/Crank-Nicolson for the basic idea)
    /// </summary>
    public class CrankNicolson : ImplicitTimeStepper {

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, params DGField[] fields)
            : this(solver, spatialOpMtx, spatialOpAffine, new CoordinateMapping(fields)) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, params DGField[] fields)
            : this(solver, AllTrue(fields.Length), spatialOpMtx, spatialOpAffine, new CoordinateMapping(fields)) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields)
            : this(solver, AllTrue(fields.Fields.Count), spatialOpMtx, spatialOpAffine, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, bool[] temporalOp, MsrMatrix spatialOpMtx, IList<double> spatialOpAffine, CoordinateMapping fields)
            : base(solver, temporalOp, spatialOpMtx, spatialOpAffine, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, SpatialOperator spatialOp, CoordinateMapping fields)
            : this(solver, AllTrue(fields.Fields.Count), spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, bool[] temporalOp, SpatialOperator spatialOp, CoordinateMapping fields)
            : base(solver, temporalOp, spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, SpatialOperator spatialOp, params DGField[] fields)
            : this(solver, AllTrue(fields.Length), spatialOp, fields) { }

        /// <summary>
        /// See <see cref="ImplicitTimeStepper"/>
        /// </summary>
        public CrankNicolson(ISparseSolverExt solver, bool[] temporalOp, SpatialOperator spatialOp, params DGField[] fields)
            : this(solver, temporalOp, spatialOp, new CoordinateMapping(fields)) { }

        /// <summary>
        /// Create an empty scheme.
        /// </summary>
        public CrankNicolson() {}

        /// <summary>
        /// The affine offset b from F(x) = Mx + b;
        /// can be null;
        /// </summary>
        /// <remarks>
        /// Compare to <see cref="ImplicitTimeStepper.m_AffineOffset1"/>:<br/>
        /// If boundary conditions are time-dependent, this vector may change over 
        /// time; Re-calculation can be implemented e.g. in <see cref="ImplicitTimeStepper.BeforeTimeStep"/>.
        /// At a given time <em>t</em>, for given initial values
        /// it is 
        /// the convention is that this vector represents the inhomogeneous b.c.
        /// at time <em>t</em>, i.e. at the current time.<br/>
        /// Can be null; If present, it will be overwritten during 
        /// the invocation of <see cref="PerformTimeStep"/>;
        /// </remarks>
        protected double[] m_AffineOffset0;

        /// <summary>
        /// Scaling between implicit and explicit Euler;
        /// <list type="bullet">
        ///   <item>0.0: explicit Euler </item>
        ///   <item>0.5: Crank-Nicholson (default)</item>
        ///   <item>1.0: implicit Euler</item>
        /// </list>
        /// </summary>
        protected double m_Theta = 0.5;

        /// <summary>
        /// Solves the linear system (diag(1/<paramref name="dt"/>) +
        /// <em>M</em> / 2) * x =
        /// <see cref="ImplicitTimeStepper.CurrentState"/> /
        /// <paramref name="dt"/> -
        /// <em>M</em> * 
        /// <see cref="ImplicitTimeStepper.CurrentState"/> / 2 -
        /// 0.5*(<see cref="ImplicitTimeStepper.m_AffineOffset1"/> +
        /// <see cref="m_AffineOffset0"/> )
        /// and writes the
        /// result do <see cref="ImplicitTimeStepper.CurrentState"/>.
        /// </summary>
        /// <param name="dt">The length of the timestep</param>
        protected override void PerformTimeStep(double dt) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {


                if (m_Theta <= 0 || m_Theta > 1)
                    throw new ApplicationException("m_Theta out of range; must be betwwen 0.0 (including) and 1.0 (including), but m_Theta = " + m_Theta);

                int n = Mapping.LocalLength;
                int np = m_diagVecOneSec.Length;

                double[] diag = new double[n];
                for (int i = 0; i < n; i++) {
                    diag[i] = m_diagVecOneSec[i % np];
                }
                BLAS.dscal(n, 1.0 / (m_Theta * dt), diag, 1);

                double[] rhs;
                if (m_AffineOffset0 != null) {
                    rhs = m_AffineOffset0;
                    BLAS.dscal(rhs.Length, -(1.0 - m_Theta) / m_Theta, rhs, 1);
                    BLAS.daxpy(rhs.Length, -1.0, m_AffineOffset1, 1, rhs, 1);
                } else {
                    rhs = (double[])m_AffineOffset1.Clone();
                }

                if (m_Theta != 1.0) {
                    // if m_Theta == 0.0, this has no effect and is a waste of comp. power
                    ISparseMatrix eM = m_Solver.GetMatrix();
                    eM.SpMV<CoordinateVector, double[]>(-(1.0 - m_Theta) / m_Theta, CurrentState, 1.0, rhs);
                }

                for (int i = 0; i < n; i++) {
                    rhs[i] += diag[i] * CurrentState[i];
                    // fraglich
                }

                tr.Info("Calling solver");
                LastSolverResult = m_Solver.Solve<double[], CoordinateVector, double[]>(1.0, diag, CurrentState, rhs);
                tr.Info("no. of iterations: " + LastSolverResult.NoOfIterations);
                tr.Info("converged? : " + LastSolverResult.Converged);
                tr.Info("Pure solver runtime: " + LastSolverResult.RunTime.TotalSeconds + " sec.");


                /*
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
                    rhs[i] += diag[i] * DGCoordinates[i];
                }

                m_Context.IOMaster.tracer.Message(ht, "Calling solver");
                LastSolverResult = m_Solver.Solve<double[], CoordinateVector, double[]>(1.0, diag, DGCoordinates, rhs);
                */

            }
        }

        /// <summary>
        /// contains the sparse solver statistics acquired in the last call to <see cref="ImplicitTimeStepper.Perform"/>
        /// </summary>
        public SolverResult LastSolverResult;
    }
}