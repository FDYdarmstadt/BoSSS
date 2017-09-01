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
using System.Linq;
using System.Text;

namespace ilPSP.LinSolvers.HYPRE {
    
    /// <summary>
    /// ilPSP wrapper for ParaSails.
    /// </summary>
    /// <remarks>
    /// From HYPRE manual:<br/>
    /// ParaSails is a parallel implementation of a sparse approximate inverse preconditioner, using a
    /// priori sparsity patterns and least-squares (Frobenius norm) minimization. Symmetric positive
    /// definite (SPD) problems are handled using a factored SPD sparse approximate inverse. General
    /// (nonsymmetric and/or indefinite) problems are handled with an unfactored sparse approximate
    /// inverse. It is also possible to precondition nonsymmetric but definite matrices with a factored,
    /// SPD preconditioner.<br/>
    /// ParaSails uses a priori sparsity patterns that are patterns of powers of sparsified matrices.
    /// ParaSails also uses a post-filtering technique to reduce the cost of applying the preconditioner. In
    /// advanced usage not described here, the pattern of the preconditioner can also be reused to generate
    /// preconditioners for different matrices in a sequence of linear solves.<br/>
    /// For more details about the ParaSails algorithm, see: E. Chow. A priori sparsity patterns for 
    /// parallel sparse approximate inverse preconditioners. SIAM J. Sci. Comput., 21:1804–1822, 2000.;
    /// </remarks>
    public class ParaSails : Solver, IImplicitPrecond {

        /// <summary>
        /// ctor
        /// </summary>
        public ParaSails() {
            base.m_NativeSetupFuncPtr = Wrappers.ParaSails.my.Delegate2FunctionPointer[
                Wrappers.ParaSails.my.HYPRE_ParaSailsSetup];
            base.m_NativeSolverFuncPtr = Wrappers.ParaSails.my.Delegate2FunctionPointer[
                Wrappers.ParaSails.my.HYPRE_ParaSailsSolve];
        }
        
        /// <summary>
        /// creates the native HYPRE object
        /// </summary>
        protected override void CreateSolver() {
            HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsCreate(MPI.Wrappers.csMPI.Raw._COMM.WORLD, out base.m_Solver));
        }

        /// <summary>
        /// see <see cref="Solver.CallSolver"/>
        /// </summary>
        /// <remarks>
        /// ParaSails should be used only as aa implicit preconditioner, this method
        /// exists only for reasons of compability and completeness.
        /// </remarks>
        protected override void CallSolver(out int NoOfIter, out bool Converged, IJVector Unknowns, IJVector Rhs) {
            if (m_Solver.p == IntPtr.Zero)
                throw new ApplicationException("solver not initialized");

            HypreException.Check(Wrappers.ParaSails.__HYPRE_ParaSailsSetup(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector));
            Wrappers.ParaSails.__HYPRE_ParaSailsSolve(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector);
            // We dont want to raise an exception for
            // a 'method did not converge'- or 'nomerical breakdown' - Error
            // that may occur in HYPRE
            Wrappers.Utilities.HYPRE_ClearAllErrors();                                            

            NoOfIter = 1;
            Converged = true;
        }

        /// <summary>
        /// destroys the native ParaSails object
        /// </summary>
        public override void Dispose() {
            base.Dispose();

            if (m_Solver.p != IntPtr.Zero) {
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsDestroy(m_Solver));
                m_Solver.p = IntPtr.Zero;
            }
        }

        /// <summary>
        /// From HYPRE manual: <br/>
        /// Value of filter parameter. The filter parameter
        /// isused to drop small nonzeros in the preconditioner, to
        /// reduce the cost of applying the preconditioner. Values
        /// from 0.05 to 0.1 are recommended. The default value
        /// is 0.1.
        /// </summary>
        public double Filter {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                double filter;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetFilter(m_Solver, out filter));
                return filter;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                if (value < 0.0)
                    throw new ArgumentOutOfRangeException("must be positive.");
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetFilter(m_Solver, value));
            }
        }

        /// <summary>
        /// From HYPRE manual: <br/>
        /// Value of the load balance parameter, 0 &#8804; <see cref="Loadbal"/> &#8804; 1.0.
        /// A zero value indicates thatno load balance is
        /// attempted; a value of unity indicatesthat perfect load
        /// balance will be attempted. The recommended value is
        /// 0.9 to balance the overhead ofdata exchanges for load
        /// balancing. No load balancingis needed if the preconditioner
        /// is very sparse andfast to construct. The default
        /// value when this parameter is not set is 0.
        /// </summary>
        public double Loadbal {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                double r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetFilter(m_Solver, out r));
                return r;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                if (value < 0.0 || value > 1.0)
                    throw new ArgumentOutOfRangeException("must be between 0 and 1;");
            }
        }

        /// <summary>
        ///  From HYPRE manual: <br/>
        ///  Value of the logging parameter. A nonzero valuesends
        ///  statistics of the setup procedure to stdout.The
        ///  default value when this parameter is not set is 0.
        /// </summary>
        public int Logging {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                int r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetLogging(m_Solver, out r));
                return r;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetLogging(m_Solver, value));
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public int NLevels {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                int r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetNlevels(m_Solver, out r));
                return r;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetNlevels(m_Solver, value));
            }
        }
        
        /// <summary>
        /// From HYPRE manual: <br/> 
        /// Value of the pattern reuse parameter. A nonzero
        /// valueindicates that the pattern of the preconditioner
        /// shouldbe reused for subsequent constructions of the
        /// preconditioner. A zero value indicates that the preconditioner
        /// should be constructed from scratch.The
        /// default value when this parameter is not set is 0.
        /// </summary>
        public bool Reuse {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                int r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetReuse(m_Solver, out r));
                return (r != 0);
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                int v = (value) ? 1 : 0;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetReuse(m_Solver, v));
            }
        }


        /// <summary>
        /// see <see cref="Sym"/>;
        /// </summary>
        public enum _Sym {
            /// <summary>
            /// from HYPRE manual: nonsymmetric and/or indefinite problem, and nonsymmetric 
            /// </summary>
            UnsymIndef = 0,

            /// <summary>
            /// from HYPRE manual: SPD problem, and SPD (factored) preconditioner
            /// </summary>
            SPD = 1,

            /// <summary>
            /// from HYPRE manual: nonsymmetric, definite problem, and SPD (factored) preconditioner
            /// </summary>
            UnsymDefSPD = 2
        }

        /// <summary>
        /// from HYPRE manual:<br/>
        /// Set the symmetry parameter for the ParaSails preconditioner.
        /// </summary>
        public _Sym Sym {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                int r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetSym(m_Solver, out r));
                return (_Sym)r;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetSym(m_Solver, (int)value));
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public double Thresh {
            get {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                double r;
                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsGetThresh(m_Solver, out r));
                return r;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                HypreException.Check(Wrappers.ParaSails.HYPRE_ParaSailsSetThresh(m_Solver, value));
            }
        }
    }
}
