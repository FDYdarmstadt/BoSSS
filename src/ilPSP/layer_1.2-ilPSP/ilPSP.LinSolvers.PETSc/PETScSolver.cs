/* =======================================================================
Copyright 2026 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using System.Diagnostics;
using System.Runtime.InteropServices;
using MPI.Wrappers;
using ilPSP.Tracing;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// L0 implementation shared by PETSc solvers which call the PETSc ABI
    /// directly, rather than the historical BoSSS PETSc C wrapper.
    /// </summary>
    /// <remarks>
    /// This class owns the PETSc AIJ matrix, KSP, vectors, common solve
    /// protocol, error handling, and disposal. Process-lifetime PETSc ABI
    /// bindings are supplied by <see cref="PETScNativeLibraries"/>. Derived classes specify
    /// only their KSP/preconditioner policy. Future direct-ABI Krylov and
    /// multigrid solvers belong on this base as well.
    /// </remarks>
    public abstract unsafe class PETScSolver : ISparseSolver {

        /// <summary>Loaded direct PETSc ABI for the selected native variant.</summary>
        protected UnsafePETSc m_Native;
        /// <summary>Native KSP handle configured by the derived solver.</summary>
        protected IntPtr m_Ksp;
        /// <summary>Native PETSc solution vector.</summary>
        protected IntPtr m_Solution;
        /// <summary>Native PETSc right-hand-side vector.</summary>
        protected IntPtr m_RightHandSide;
        /// <summary>Native communicator used by the owned PETSc objects.</summary>
        protected IntPtr m_Communicator;
        PETScMPIAIJMatrix m_Matrix;
        UnsafePETSc.KspMonitor m_KspMonitor;
        int m_LocalRows = -1;
        bool m_Disposed;

        /// <summary>PETSc matrix ownership and transfer representation.</summary>
        public PETScMatrixInput MatrixInput { get; set; } = PETScMatrixInput.BorrowedCsr;

        /// <summary>
        /// Writes outer-KSP iteration residuals through a managed PETSc monitor.
        /// This avoids PETSc's viewer-based <c>-ksp_monitor</c> option path.
        /// </summary>
        public bool EnableKspMonitor { get; set; }

        /// <summary>Writes the final PETSc convergence reason to standard output.</summary>
        public bool EnableKspConvergedReason { get; set; }

        /// <summary>PETSc convergence reason from the most recent solve.</summary>
        public int ConvergenceReason { get; private set; }

        /// <summary>Number of KSP iterations from the most recent solve.</summary>
        public int IterationCount { get; private set; } = -1;

        /// <summary>Residual norm from the most recent solve.</summary>
        public double ResidualNorm { get; private set; } = double.NaN;

        /// <summary>Whether the most recent PETSc KSP solve did not converge.</summary>
        public bool DidNotConverge => ConvergenceReason <= 0;

        /// <summary>Native PETSc error description from the most recent failure.</summary>
        public string LastErrorMessage { get; private set; }

        /// <summary>Configures the derived solver's KSP and preconditioner.</summary>
        protected abstract void ConfigureKsp();

        /// <summary>Validates solver-specific MPI/runtime restrictions.</summary>
        protected virtual void ValidateParallelConfiguration(int mpiSize) {
        }

        /// <summary>Selects the matrix transfer mode for this KSP configuration.</summary>
        protected virtual PETScMatrixInput SelectMatrixInput(int mpiSize) {
            return MatrixInput;
        }

        /// <summary>Defines one distributed square matrix and creates its KSP.</summary>
        public void DefineMatrix(IMutableMatrixEx matrix) {
            using(var tr = new FuncTrace("PETSC2_DefineMatrix")) {
                if(matrix == null)
                    throw new ArgumentNullException(nameof(matrix));
                ThrowIfDisposed();
                if(m_Matrix != null)
                    throw new InvalidOperationException("DefineMatrix can only be called once in the lifetime of a PETSc solver.");
                if(matrix.NoOfRows != matrix.NoOfCols || matrix.NoOfRows <= 0)
                    throw new ArgumentException("PETSc requires a non-empty square matrix.", nameof(matrix));
                if(!matrix.MPI_Comm.Equals(csMPI.Raw._COMM.WORLD))
                    throw new NotSupportedException("The direct PETSc ABI currently supports MPI_COMM_WORLD only.");

                csMPI.Raw.Comm_Size(matrix.MPI_Comm, out int mpiSize);
                ValidateParallelConfiguration(mpiSize);

                MPICollectiveWatchDog.Watch(matrix.MPI_Comm);

                try {
                    m_Native = PETScNativeLibraries.Get();
                    EnsurePetscInitialized();
                    m_Communicator = GetNativeCommunicator(matrix.MPI_Comm);
                    using(new BlockTrace("PETSC2_NativeMatrixCreate", tr)) {
                        m_Matrix = PETScMPIAIJMatrix.Create(m_Native, m_Communicator, SelectMatrixInput(mpiSize), matrix);
                        m_LocalRows = m_Matrix.LocalRows;
                    }
                    using(new BlockTrace("PETSC2_NativeKspCreate", tr)) {
                        Check(m_Native.CreateKsp(m_Communicator, out m_Ksp), "PETSc KSPCreate");
                        Check(m_Native.SetKspOperators(m_Ksp, m_Matrix.Handle, m_Matrix.Handle), "PETSc KSPSetOperators");
                        ConfigureKsp();
                        Check(m_Native.CreateMatrixVectors(m_Matrix.Handle, out m_Solution, out m_RightHandSide), "PETSc MatCreateVecs");
                        Check(m_Native.SetUpKsp(m_Ksp), "PETSc KSPSetUp");
                        if(EnableKspMonitor) {
                            // Keep the delegate alive until KSPDestroy: PETSc holds its
                            // function pointer and may invoke it during every KSPSolve.
                            m_KspMonitor = WriteKspMonitor;
                            Check(m_Native.SetKspMonitor(m_Ksp, m_KspMonitor, IntPtr.Zero, IntPtr.Zero), "PETSc KSPMonitorSet");
                        }
                    }
                } catch {
                    DestroyNativeObjects();
                    m_Native = null;
                    m_LocalRows = -1;
                    throw;
                }
            }
        }

        /// <summary>Solves the matrix defined by <see cref="DefineMatrix"/>.</summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            using(var tr = new FuncTrace("PETSC2_Solve")) {
                ThrowIfDisposed();
                if(m_Ksp == IntPtr.Zero)
                    throw new InvalidOperationException("DefineMatrix must be called before Solve.");
                if(x == null)
                    throw new ArgumentNullException(nameof(x));
                if(rhs == null)
                    throw new ArgumentNullException(nameof(rhs));
                if(x.Count != m_LocalRows || rhs.Count != m_LocalRows)
                    throw new ArgumentException($"Expected local vectors of length {m_LocalRows}, got x={x.Count} and rhs={rhs.Count}.");

                MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                CopyToNativeVector(m_Solution, x);
                CopyToNativeVector(m_RightHandSide, rhs);
                Check(m_Native.SetKspInitialGuessNonzero(m_Ksp, ContainsNonzero(x) ? 1 : 0), "PETSc KSPSetInitialGuessNonzero");

                var stopwatch = Stopwatch.StartNew();
                Check(m_Native.SolveKsp(m_Ksp, m_RightHandSide, m_Solution), "PETSc KSPSolve");
                stopwatch.Stop();

                Check(m_Native.GetKspIterationNumber(m_Ksp, out int iterations), "PETSc KSPGetIterationNumber");
                Check(m_Native.GetKspResidualNorm(m_Ksp, out double residualNorm), "PETSc KSPGetResidualNorm");
                Check(m_Native.GetKspConvergedReason(m_Ksp, out int reason), "PETSc KSPGetConvergedReason");
                CopyFromNativeVector(m_Solution, x);

                ConvergenceReason = reason;
                IterationCount = iterations;
                ResidualNorm = residualNorm;
                if(EnableKspConvergedReason)
                    WriteConvergedReason(reason, iterations, residualNorm);
                ValidateResult(x, residualNorm, reason);
                return new SolverResult { Converged = reason > 0, NoOfIterations = iterations, RunTime = stopwatch.Elapsed };
            }
        }

        /// <summary>Applies PETSc command-line options using PETSc's process option database.</summary>
        protected void ApplyKspOptions(string options) {
            if(string.IsNullOrWhiteSpace(options))
                return;
            // Match the established C wrapper: PC and nested KSP objects also
            // consult PETSc's process option database during KSPSetFromOptions.
            Check(m_Native.OptionsClear(IntPtr.Zero), "PETSc PetscOptionsClear");
            Check(m_Native.OptionsInsertString(IntPtr.Zero, options.Trim()), "PETSc PetscOptionsInsertString");
            Check(m_Native.SetKspFromOptions(m_Ksp), "PETSc KSPSetFromOptions");
        }

        void EnsurePetscInitialized() {
            Check(m_Native.Initialized(out int initialized), "PETSc PetscInitialized");
            if(initialized == 0)
                Check(m_Native.InitializeNoArguments(), "PETSc PetscInitializeNoArguments");
        }

        void CopyToNativeVector(IntPtr vector, IList<double> source) {
            Check(m_Native.GetVectorArray(vector, out IntPtr values), "PETSc VecGetArray");
            try {
                double* destination = (double*)values;
                for(int i = 0; i < m_LocalRows; i++) destination[i] = source[i];
            } finally {
                Check(m_Native.RestoreVectorArray(vector, ref values), "PETSc VecRestoreArray");
            }
        }

        void CopyFromNativeVector(IntPtr vector, IList<double> destination) {
            Check(m_Native.GetVectorArrayRead(vector, out IntPtr values), "PETSc VecGetArrayRead");
            try {
                double* source = (double*)values;
                for(int i = 0; i < m_LocalRows; i++) destination[i] = source[i];
            } finally {
                Check(m_Native.RestoreVectorArrayRead(vector, ref values), "PETSc VecRestoreArrayRead");
            }
        }

        static bool ContainsNonzero(IList<double> values) {
            for(int i = 0; i < values.Count; i++) if(values[i] != 0.0) return true;
            return false;
        }

        static int WriteKspMonitor(IntPtr ksp, int iteration, double residualNorm, IntPtr context) {
            try {
                Console.WriteLine($"PETSC_KSP: iteration={iteration}, residual={residualNorm:R}");
            } catch {
                // A managed exception must never cross PETSc's native callback boundary.
            }
            return 0;
        }

        static void WriteConvergedReason(int convergenceReason, int iterations, double residualNorm) {
            try {
                Console.WriteLine($"PETSC_KSP: convergence_reason={convergenceReason}, iterations={iterations}, residual={residualNorm:R}");
            } catch {
                // Console failures must not change solver behavior.
            }
        }

        static void ValidateResult(IList<double> solution, double residualNorm, int reason) {
            if(reason <= 0) throw new ArithmeticException($"PETSc KSP did not converge; convergence reason {reason}, residual {residualNorm:R}.");
            if(double.IsNaN(residualNorm) || double.IsInfinity(residualNorm)) throw new ArithmeticException($"PETSc returned non-finite residual {residualNorm}.");
            for(int i = 0; i < solution.Count; i++)
                if(double.IsNaN(solution[i]) || double.IsInfinity(solution[i])) throw new ArithmeticException($"PETSc returned a non-finite solution at local index {i}.");
        }

        /// <summary>Converts a PETSc error code into a managed exception.</summary>
        protected void Check(int errorCode, string operation) {
            if(errorCode == 0) return;
            string message = null;
            if(m_Native != null && m_Native.ErrorMessage(errorCode, out IntPtr text, IntPtr.Zero) == 0 && text != IntPtr.Zero)
                message = Marshal.PtrToStringAnsi(text);
            LastErrorMessage = message;
            throw new ApplicationException($"{operation} failed with PETSc error code {errorCode}{(string.IsNullOrEmpty(message) ? string.Empty : ": " + message)}");
        }

        static IntPtr GetNativeCommunicator(MPI_Comm communicator) {
            int size = csMPI.Raw.MPI_Comm_f2c(communicator, out uint fourByteCommunicator, out ulong eightByteCommunicator);
            switch(size) {
                case 4: return new IntPtr(unchecked((int)fourByteCommunicator));
                case 8: return new IntPtr(unchecked((long)eightByteCommunicator));
                default: throw new NotSupportedException($"Unsupported native MPI_Comm size {size}.");
            }
        }

        /// <summary>Destroys all PETSc objects owned by this solver.</summary>
        public void Dispose() {
            if(m_Disposed) return;
            m_Disposed = true;
            Exception destroyError = null;
            try {
                if(m_Matrix != null) MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                DestroyNativeObjects();
            } catch(Exception error) {
                destroyError = error;
            } finally {
                m_Native = null;
                m_LocalRows = -1;
                m_Communicator = IntPtr.Zero;
            }
            if(destroyError != null) throw new ApplicationException("PETSc direct-ABI solver destruction failed.", destroyError);
        }

        /// <summary>Destroys common PETSc objects; derived classes may extend this cleanup.</summary>
        protected virtual void DestroyNativeObjects() {
            if(m_Native == null) return;
            if(m_Solution != IntPtr.Zero) Check(m_Native.DestroyVector(ref m_Solution), "PETSc VecDestroy(solution)");
            if(m_RightHandSide != IntPtr.Zero) Check(m_Native.DestroyVector(ref m_RightHandSide), "PETSc VecDestroy(rhs)");
            if(m_Ksp != IntPtr.Zero) Check(m_Native.DestroyKsp(ref m_Ksp), "PETSc KSPDestroy");
            m_KspMonitor = null;
            if(m_Matrix != null) { m_Matrix.Dispose(); m_Matrix = null; }
        }

        /// <summary>Throws when a public operation is attempted after disposal.</summary>
        protected void ThrowIfDisposed() {
            if(m_Disposed) throw new ObjectDisposedException(GetType().Name);
        }
    }
}
