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
using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE {
    
    /// <summary>
    /// Euclid preconditioner;
    /// </summary>
    /// <remarks>
    /// From HYPRE 2.4.0b Manual:<br/>
    /// The Euclid library is a scalable implementation of the Parallel ILU algorithm that was presented at
    /// SC99 [1], and published in expanded form in the SIAM Journal on Scientific Computing [2]. By
    /// scalable we mean that the factorization (setup) and application (triangular solve) timings remain
    /// nearly constant when the global problem size is scaled in proportion to the number of processors.
    /// As with all ILU preconditioning methods, the number of iterations is expected to increase with
    /// global problem size.
    /// <list type="bullet">
    /// <item>[1]: D. Hysom and A. Pothen. Efficient parallel computation of ILU(k) preconditioners. In Proceedings
    /// of Supercomputing ’99. ACM, November 1999. Published on CDROM, ISBN #1-58113-091-0, ACM Order #415990, IEEE Computer Society Press Order # RS00197.</item>
    /// <item>[2]: D. Hysom and A. Pothen. Efficient parallel computation of ILU(k) preconditioners. In Proceedings
    /// of Supercomputing ’99. ACM, November 1999. Published on CDROM, ISBN #1-58113-091-0, ACM Order #415990, IEEE Computer Society Press Order # RS00197.</item>
    /// </list>
    /// </remarks>
    public class Euclid : Solver, IImplicitPrecond {

        /// <summary>
        /// ctor
        /// </summary>
        public Euclid() {
            base.m_NativeSolverFuncPtr = Wrappers.Euclid.my.Delegate2FunctionPointer[Wrappers.Euclid.my.HYPRE_EuclidSolve];
            base.m_NativeSetupFuncPtr = Wrappers.Euclid.my.Delegate2FunctionPointer[Wrappers.Euclid.my.HYPRE_EuclidSetup];
        }

        /// <summary>
        /// creates solver for MPI_COMM_WORLD communicator
        /// </summary>
        protected override void CreateSolver() {
            HypreException.Check(Wrappers.Euclid.Create(csMPI.Raw._COMM.WORLD, out base.m_Solver));
        }


        /// <summary>
        /// uses ILUT and defines a drop tolerance relative to the largest absolute value of any entry in the row being
        /// factored
        /// </summary>
        public double ILUT {
            set {
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetILUT(m_Solver, value));
            }
        }

        bool m_BJ = false;

        /// <summary>
        /// From HYPRE manual:<br/>
        /// Use Block Jacobi ILU preconditioning instead of PILU. Default: 0 (false). Guidance: if subdomains
        /// contain relatively few nodes (less than 1,000), or the problem is not well partitioned,
        /// Block Jacobi ILU may give faster solution time than PILU        
        /// </summary>
        public bool BJ {
            set {
                m_BJ = value;
                int __BJ = value ? 1: 0;
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetBJ(m_Solver, __BJ));
            }
            get {
                return m_BJ;
            }
        }

        bool m_RowScale = false;

        /// <summary>
        /// From HYPRE manual:<br/>
        /// Scale values prior to factorization such that the largest value in any row is +1 or -1.
        /// Default: 0 (false). CAUTION: If the coefficient matrix A is symmetric, this setting is likely
        /// to cause the filled matrix, F = L + U − I, to be unsymmetric. Guidance: if the matrix is
        /// poorly scaled, turning on row scaling may help convergence.
        /// </summary>
        public bool RowScale {
            set {
                m_RowScale = value;
                int __RS = value ? 1: 0;
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetRowScale(m_Solver, __RS));
            }
            get {
                return m_RowScale;
            }
        }

        double m_SparseA = 0.0;

        /// <summary>
        /// From HYPRE manual:<br/>
        /// Drop-tolerance for ILU(k) factorization. Default: 0 (no dropping). Entries are
        /// treated as zero if their absolute value is less than (sparseA * max), where “max” is the largest
        /// absolute value of any entry in the row. Guidance: try this in conjunction with -rowScale.
        /// CAUTION: If the coefficient matrix A is symmetric, this setting is likely to cause the filled
        /// matrix, F = L+U−I, to be unsymmetric. This setting has no effect when ILUT factorization
        /// is selected.
        /// </summary>
        public double SparseA {
            get {
                return m_SparseA;
            }
            set {
                m_SparseA = value;
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetSparseA(m_Solver, value));
            }
        }

        bool m_Mem = false;

        /// <summary>
        /// From HYPRE manual:<br/>
        /// When Euclid’s destructor is called a summary of runtime settings and timing information
        /// is printed to stdout. Default: 0 (false). The timing marks in the report are the maximum
        /// over all processors in the MPI communicator.
        /// </summary>
        public bool Mem {
            get {
                return m_Mem;
            }
            set {
                m_Mem = value;
                int __RS = value ? 1: 0;
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetMem(m_Solver, __RS));
            }
        }

                
        /// <summary>
        /// usual solver call
        /// </summary>
        protected override void CallSolver(out int NoOfIter, out bool Converged, IJVector Unknowns, IJVector Rhs) {
            if (m_Solver.p == IntPtr.Zero)
                throw new ApplicationException("solver not initialized");

            HypreException.Check(Wrappers.Euclid.__HYPRE_EuclidSetup(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector));
            Wrappers.Euclid.__HYPRE_EuclidSolve(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector);
            // We dont want to raise an exception for
            // a 'method did not converge'- or 'nomerical breakdown' - Error
            // that may occur in HYPRE
            Wrappers.Utilities.HYPRE_ClearAllErrors();                                            


            NoOfIter = 1;
            Converged = true;
        }

        int m_Level = 1;

        /// <summary>
        /// from HYPRE manual:<br/>
        /// Factorization level for ILU(k). Default: 1. Guidance: for 2D convection-diffusion and
        /// similar problems, fastest solution time is typically obtained with levels 4 through 8. For 3D
        /// problems fastest solution time is typically obtained with level 1.
        /// </summary>
        public int Level {
            get {
                return m_Level;
            }
            set {
                if (value < 1)
                    throw new ArgumentOutOfRangeException("must be greater or equal to 1.");
                m_Level = value;
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidSetLevel(m_Solver, m_Level));
            }
        }


        /// <summary>
        /// disposal
        /// </summary>
        public override void Dispose() {
            base.Dispose();

            if (m_Solver.p != IntPtr.Zero) {
                HypreException.Check(Wrappers.Euclid.HYPRE_EuclidDestroy(m_Solver));
                m_Solver.p = IntPtr.Zero;
            }
        }
    }
}
