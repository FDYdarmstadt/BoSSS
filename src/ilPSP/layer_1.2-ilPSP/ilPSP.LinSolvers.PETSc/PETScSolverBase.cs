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
    /// Common native matrix, solver, and lifecycle implementation for PETSc solvers.
    /// </summary>
    /// <remarks>
    /// The native implementation supports Linux and Windows and uses
    /// MPI_COMM_WORLD for the MPI variants. <see cref="PETScNativeVariant.OpenMp"/> uses PETSc MPIUNI
    /// and is therefore restricted to one BoSSS MPI rank. Derived classes only
    /// supply their PETSc options and optional field decomposition.
    /// </remarks>
    public abstract partial class PETScSolverBase : ISparseSolver {

        UnsafePETSc m_Native;
        IntPtr m_Solver;
        int m_LocalRows = -1;
        bool m_Disposed;
        UnsafePETSc._bosss_petsc_ksp_monitor m_KspMonitor;
        /// <summary>
        /// Selects the assembled CSR transfer path. The default is the split,
        /// borrowed CSR path.
        /// </summary>
        public PETScMatrixInput MatrixInput { get; set; } = PETScMatrixInput.BorrowedCsr;

        /// <summary>
        /// Matrix input path actually used by the most recent <see cref="DefineMatrix"/> call.
        /// </summary>
        public PETScMatrixInput? EffectiveMatrixInput { get; private set; }

        /// <summary>
        /// Diagnostic option which keeps the managed CSR buffers pinned until
        /// the native PETSc matrix has been destroyed.
        /// </summary>
        public bool PinManagedCsrBuffers { get; set; }

        /// <summary>
        /// Parallel runtime used by PETSc. This must be set before <see cref="DefineMatrix"/>.
        /// The default selects a variant from the active MPI and threading configuration.
        /// </summary>
        public PETScNativeVariant NativeVariant { get; set; } = PETScNativeVariant.Auto;

        /// <summary>
        /// Native variant selected when <see cref="DefineMatrix"/> was called.
        /// </summary>
        public PETScNativeVariant? EffectiveNativeVariant { get; private set; }

        /// <summary>
        /// Resolves an explicitly requested variant or selects one automatically.
        /// </summary>
        public static PETScNativeVariant ResolveNativeVariant(PETScNativeVariant requestedVariant, int mpiSize, int numThreads) {
            if(mpiSize < 1)
                throw new ArgumentOutOfRangeException(nameof(mpiSize), mpiSize, "The MPI size must be positive.");
            if(numThreads < 1)
                throw new ArgumentOutOfRangeException(nameof(numThreads), numThreads, "The thread count must be positive.");

            if(requestedVariant != PETScNativeVariant.Auto)
                return requestedVariant;
            if(mpiSize == 1)
                return PETScNativeVariant.OpenMp;
            if(numThreads == 1)
                return PETScNativeVariant.Mpi;
            return PETScNativeVariant.MpiOpenMp;
        }

        /// <summary>
        /// PETSc command-line-style options supplied by the concrete solver.
        /// </summary>
        protected abstract string SolverOptions { get; }

        /// <summary>
        /// Optional local portion of the first field's global row indices for PCFIELDSPLIT.
        /// Both field arrays must be supplied together.
        /// </summary>
        protected virtual int[] SolverFirstFieldRows => null;

        /// <summary>
        /// Optional local portion of the second field's global row indices for PCFIELDSPLIT.
        /// Both field arrays must be supplied together.
        /// </summary>
        protected virtual int[] SolverSecondFieldRows => null;

        /// <summary>
        /// Allows a concrete solver to add context to a native setup failure.
        /// </summary>
        protected virtual Exception CreateSetupException(Exception error) {
            return error;
        }

        /// <summary>
        /// Validates solver-specific restrictions of the active MPI and PETSc
        /// configuration before native PETSc resources are created.
        /// </summary>
        protected virtual void ValidateParallelConfiguration(int mpiSize, PETScNativeVariant effectiveNativeVariant) {
        }

        /// <summary>
        /// PETSc convergence reason returned by the most recent solve.
        /// Positive values indicate convergence; negative values indicate divergence.
        /// </summary>
        public int ConvergenceReason { get; private set; }

        /// <summary>
        /// Number of KSP iterations reported by the most recent solve, or -1
        /// when PETSc failed before it could report solver status.
        /// </summary>
        public int IterationCount { get; private set; } = -1;

        /// <summary>
        /// True if the most recent PETSc KSP solve did not converge.
        /// PETSc uses positive convergence reasons and zero or negative values
        /// for an unfinished or diverged solve.
        /// </summary>
        public bool DidNotConverge => ConvergenceReason <= 0;

        /// <summary>
        /// Residual norm returned by PETSc for the most recent solve.
        /// </summary>
        public double ResidualNorm { get; private set; } = double.NaN;

        /// <summary>
        /// PETSc's native error description from the most recent failed solve.
        /// This is null after a successful solve.
        /// </summary>
        public string LastErrorMessage { get; private set; }

        /// <summary>
        /// Defines and uploads a distributed square matrix to PETSc.
        /// </summary>
        public void DefineMatrix(IMutableMatrixEx matrix) {
            using(var tr = new FuncTrace("PETSC_DefineMatrix")) {
            // 1. Validate the BoSSS matrix and choose the matching PETSc DLL.
            if(matrix == null)
                throw new ArgumentNullException(nameof(matrix));
            ThrowIfDisposed();
            if(m_Solver != IntPtr.Zero)
                throw new InvalidOperationException("DefineMatrix can only be called once in the lifetime of a PETSc solver.");
            if(matrix.NoOfRows != matrix.NoOfCols || matrix.NoOfRows <= 0)
                throw new ArgumentException("PETSc requires a non-empty square matrix.", nameof(matrix));
            if(!matrix.MPI_Comm.Equals(csMPI.Raw._COMM.WORLD))
                throw new NotSupportedException("The current PETSc native ABI supports MPI_COMM_WORLD only.");

            csMPI.Raw.Comm_Size(matrix.MPI_Comm, out int mpiSize);
            PETScNativeVariant effectiveNativeVariant = ResolveNativeVariant(NativeVariant, mpiSize, ilPSP.Environment.NumThreads);
            EffectiveNativeVariant = effectiveNativeVariant;
            if(effectiveNativeVariant == PETScNativeVariant.OpenMp && mpiSize != 1)
                throw new NotSupportedException("The PETSc OpenMP variant uses MPIUNI and requires exactly one BoSSS MPI rank. Use MpiOpenMp for multiple ranks.");
            ValidateParallelConfiguration(mpiSize, effectiveNativeVariant);

            MPICollectiveWatchDog.Watch(matrix.MPI_Comm);

            // 2. Validate the PETSc 32-bit global partition and optional field split.
            long globalRows = matrix.RowPartitioning.TotalLength;
            long globalRowStart = matrix.RowPartitioning.i0;
            if(globalRows > int.MaxValue || globalRowStart > int.MaxValue)
                throw new NotSupportedException("The packaged PETSc uses 32-bit indices and this matrix exceeds that range.");

            ValidateFieldRows(SolverFirstFieldRows, SolverSecondFieldRows, globalRows);

            m_LocalRows = matrix.RowPartitioning.LocalLength;
            // 3. Convert once to either copied CSR or PETSc's split borrowed CSR.
            // MatCreateMPIAIJWithSplitArrays deliberately creates MATMPIAIJ even
            // for one rank. Top-level ILU/ICC and the UMFPACK factor package are
            // sequential and require MATSEQAIJ. MatCreateAIJ in the copied path
            // selects SEQAIJ for a one-rank communicator, so use it only for this
            // combination. MPI runs with more than one rank retain their requested
            // input path; PETScDirectSolver rejects top-level UMFPACK for such a run.
            bool requiresSequentialFactorMatrix = mpiSize == 1 && RequestsSequentialFactorMatrix(SolverOptions);
            bool splitCsr = MatrixInput == PETScMatrixInput.BorrowedCsr && !requiresSequentialFactorMatrix;
            CsrArrays copiedCsr = null;
            SplitCsrArrays borrowedCsr = null;
            using(new BlockTrace("PETSC_CsrExport", tr)) {
                if(splitCsr) {
                    borrowedCsr = ExportSplitCsr(matrix, globalRows, globalRowStart);
                    PinSplitCsrBuffers(borrowedCsr);
                } else {
                    copiedCsr = ExportCsr(matrix, globalRows, globalRowStart);
                    if(PinManagedCsrBuffers)
                        PinCsrBuffers(copiedCsr.RowPointers, copiedCsr.ColumnIndices, copiedCsr.Values);
                }
            }
            EffectiveMatrixInput = splitCsr ? PETScMatrixInput.BorrowedCsr : PETScMatrixInput.CopiedCsr;
            try {
                using(new BlockTrace("PETSC_NativeMatrixCreate", tr)) {
                    // 4. Load the selected native PETSc DLL and create its Mat object.
                    m_Native = new UnsafePETSc(effectiveNativeVariant);
                    // This records that BoSSS initialized PETSc, allowing application
                    // shutdown to finalize PETSc before MPI is finalized.
                    Check(m_Native.InitializeRuntime(), "PETSc initialization");
                    if(splitCsr) {
                        Check(m_Native.CreateSplitCsrMatrixSolver(
                            m_LocalRows,
                            (int)globalRows,
                            (int)globalRowStart,
                            borrowedCsr.DiagonalRowPointers,
                            borrowedCsr.DiagonalColumnIndices,
                            borrowedCsr.DiagonalValues,
                            borrowedCsr.OffdiagonalRowPointers,
                            borrowedCsr.OffdiagonalColumnIndices,
                            borrowedCsr.OffdiagonalValues,
                            out m_Solver), "PETSc split CSR matrix creation");
                    } else {
                        Check(m_Native.CreateMatrixSolver(
                            m_LocalRows,
                            (int)globalRows,
                            (int)globalRowStart,
                            copiedCsr.RowPointers,
                            copiedCsr.ColumnIndices,
                            copiedCsr.Values,
                            out m_Solver), "PETSc matrix creation");
                    }
                }
                using(new BlockTrace("PETSC_NativeKspCreate", tr)) {
                    // 5. Create the KSP object for that matrix.
                    Check(m_Native.CreateKspForMatrix(m_Solver), "PETSc KSP creation");
                }
                using(new BlockTrace("PETSC_NativeOptions", tr)) {
                    // 6. Apply KSP/PC options and optional field-split indices.
                    Check(m_Native.ConfigureKspSolver(
                        m_Solver,
                        SolverOptions ?? string.Empty,
                        SolverFirstFieldRows?.Length ?? 0,
                        SolverFirstFieldRows,
                        SolverSecondFieldRows?.Length ?? 0,
                        SolverSecondFieldRows), "PETSc KSP option configuration");
                }
                using(new BlockTrace("PETSC_KspSetup", tr)) {
                    // 7. Force PETSc to allocate/setup the selected solver now.
                    Check(m_Native.SetupKspSolver(m_Solver), "PETSc KSP setup");
                    if(RequestsKspMonitor(SolverOptions)) {
                        m_KspMonitor = WriteKspMonitor;
                        Check(m_Native.SetKspMonitor(m_Solver, m_KspMonitor, IntPtr.Zero), "PETSc KSP monitor setup");
                    }
                }
            } catch(Exception error) {
                if(m_Solver != IntPtr.Zero && m_Native != null) {
                    // The staged constructor can have created a native solver before
                    // KSPSetUp reports an error. Release it before unloading delegates.
                    m_Native.DestroySolver(m_Solver);
                    m_Solver = IntPtr.Zero;
                }
                m_Native?.Dispose();
                m_Native = null;
                m_LocalRows = -1;
                ReleasePinnedCsrBuffers();
                Exception setupError = CreateSetupException(error);
                if(ReferenceEquals(setupError, error))
                    throw;
                    throw setupError;
            }
            }
        }

        // This remains part of the assembled-CSR path: PETSc's sequential
        // ILU/ICC and UMFPACK require MATSEQAIJ for one MPI rank.
        static bool RequestsSequentialFactorMatrix(string options) {
            if(string.IsNullOrWhiteSpace(options))
                return false;

            string pcType = null;
            string factorBackend = null;
            string[] tokens = options.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            for(int index = 0; index < tokens.Length; index++) {
                string token = tokens[index];
                if(token.StartsWith("-pc_type=", StringComparison.Ordinal))
                    pcType = token.Substring("-pc_type=".Length);
                else if(token.Equals("-pc_type", StringComparison.Ordinal) && index + 1 < tokens.Length)
                    pcType = tokens[++index];
                else if(token.StartsWith("-pc_factor_mat_solver_type=", StringComparison.Ordinal))
                    factorBackend = token.Substring("-pc_factor_mat_solver_type=".Length);
                else if(token.Equals("-pc_factor_mat_solver_type", StringComparison.Ordinal) && index + 1 < tokens.Length)
                    factorBackend = tokens[++index];
            }
            return string.Equals(pcType, "ilu", StringComparison.OrdinalIgnoreCase)
                || string.Equals(pcType, "icc", StringComparison.OrdinalIgnoreCase)
                || string.Equals(factorBackend, "umfpack", StringComparison.OrdinalIgnoreCase);
        }

        /// <summary>
        /// Solves the matrix defined by <see cref="DefineMatrix"/>.
        /// </summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            using(var tr = new FuncTrace("PETSC_Solve")) {
            //tr.InfoToConsole = true;
            ThrowIfDisposed();
            if(m_Solver == IntPtr.Zero)
                throw new InvalidOperationException("DefineMatrix must be called before Solve.");
            if(x == null)
                throw new ArgumentNullException(nameof(x));
            if(rhs == null)
                throw new ArgumentNullException(nameof(rhs));
            if(x.Count != m_LocalRows || rhs.Count != m_LocalRows)
                throw new ArgumentException($"Expected local vectors of length {m_LocalRows}, got x={x.Count} and rhs={rhs.Count}.");

            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            double[] solution = x as double[] ?? CopyToArray(x);
            double[] rightHandSide = rhs as double[] ?? CopyToArray(rhs);
            bool useInitialGuess = false;
            for(int i = 0; i < solution.Length; i++)
                useInitialGuess |= solution[i] != 0.0;

            double[] nativeSolution = m_LocalRows > 0 ? solution : new double[1];
            double[] nativeRhs = m_LocalRows > 0 ? rightHandSide : new double[1];
            var stopwatch = Stopwatch.StartNew();
            int iterations;
            double residualNorm;
            int convergenceReason;
            LastErrorMessage = null;
            int solveError = m_Native.Solve(
                m_Solver,
                nativeSolution,
                nativeRhs,
                useInitialGuess ? 1 : 0,
                out iterations,
                out residualNorm,
                out convergenceReason);
            stopwatch.Stop();

            ConvergenceReason = convergenceReason;
            IterationCount = iterations;
            ResidualNorm = residualNorm;
            if(RequestsConvergedReason(SolverOptions))
                WriteConvergedReason(convergenceReason, iterations, residualNorm);
            if(solveError != 0) {
                LastErrorMessage = GetLastErrorMessage();
                throw new ApplicationException(
                    $"PETSc KSPSolve failed with PETSc error code {solveError}: {LastErrorMessage}; "
                    + $"KSP reason {convergenceReason}, iterations {iterations}, residual {residualNorm:R}.");
            }
            ValidateResult(solution, residualNorm, convergenceReason);

            if(!ReferenceEquals(solution, x)) {
                for(int i = 0; i < solution.Length; i++)
                    x[i] = solution[i];
            }

            return new SolverResult {
                Converged = convergenceReason > 0,
                NoOfIterations = iterations,
                RunTime = stopwatch.Elapsed
            };
            }
        }

        static double[] CopyToArray(IList<double> source) {
            var result = new double[source.Count];
            for(int i = 0; i < result.Length; i++)
                result[i] = source[i];
            return result;
        }

        static void ValidateFieldRows(int[] first, int[] second, long globalRows) {
            bool hasFirst = first != null && first.Length > 0;
            bool hasSecond = second != null && second.Length > 0;
            if(hasFirst != hasSecond)
                throw new ArgumentException("Both PETSc field-split row arrays must be supplied together.");
            if(!hasFirst)
                return;

            foreach(int row in first)
                if(row < 0 || row >= globalRows)
                throw new ArgumentOutOfRangeException("firstFieldRows", $"Field row {row} is outside [0,{globalRows}).");
            foreach(int row in second)
                if(row < 0 || row >= globalRows)
                    throw new ArgumentOutOfRangeException("secondFieldRows", $"Field row {row} is outside [0,{globalRows}).");
        }

        static void ValidateResult(double[] solution, double residualNorm, int convergenceReason) {
            if(convergenceReason <= 0)
                throw new ArithmeticException($"PETSc KSP did not converge; convergence reason {convergenceReason}, residual {residualNorm:R}.");
            if(double.IsNaN(residualNorm) || double.IsInfinity(residualNorm))
                throw new ArithmeticException($"PETSc returned non-finite residual {residualNorm}; convergence reason {convergenceReason}.");
            for(int i = 0; i < solution.Length; i++)
                if(double.IsNaN(solution[i]) || double.IsInfinity(solution[i]))
                    throw new ArithmeticException($"PETSc returned a non-finite solution at local index {i}; convergence reason {convergenceReason}.");
        }

        static void WriteKspMonitor(IntPtr context, int iteration, double residualNorm) {
            try {
                Console.WriteLine($"PETSC_KSP: iteration={iteration}, residual={residualNorm:R}");
            } catch {
                // Exceptions must never cross the native PETSc callback boundary.
            }
        }

        static bool RequestsKspMonitor(string options) {
            return ContainsPetscOption(options, "ksp_monitor", allowSuffix: true);
        }

        static bool RequestsConvergedReason(string options) {
            return ContainsPetscOption(options, "ksp_converged_reason", allowSuffix: false);
        }

        static bool ContainsPetscOption(string options, string optionName, bool allowSuffix) {
            if(string.IsNullOrWhiteSpace(options))
                return false;

            string prefix = "-" + optionName;
            foreach(string token in options.Split((char[])null, StringSplitOptions.RemoveEmptyEntries)) {
                if(token.Equals(prefix, StringComparison.Ordinal))
                    return true;
                if(allowSuffix
                    && token.StartsWith(prefix + "_", StringComparison.Ordinal)
                    && token.Length > prefix.Length + 1)
                    return true;
            }
            return false;
        }

        static void WriteConvergedReason(int convergenceReason, int iterations, double residualNorm) {
            try {
                Console.WriteLine($"PETSC_KSP: convergence_reason={convergenceReason}, iterations={iterations}, residual={residualNorm:R}");
            } catch {
                // Console failures must not change solver behavior.
            }
        }

        void Check(int errorCode, string operation) {
            if(errorCode == 0)
                return;
            string message = GetLastErrorMessage();
            throw new ApplicationException($"{operation} failed with PETSc error code {errorCode}: {message}");
        }

        string GetLastErrorMessage() {
            IntPtr messagePointer = m_Native?.LastError() ?? IntPtr.Zero;
            return messagePointer == IntPtr.Zero ? "No PETSc error message was provided." : Marshal.PtrToStringAnsi(messagePointer);
        }

        void ThrowIfDisposed() {
            if(m_Disposed)
                throw new ObjectDisposedException(GetType().Name);
        }

        /// <summary>
        /// Releases the native KSP, matrix, and dynamic-library handle.
        /// This is collective on MPI_COMM_WORLD when a matrix was defined.
        /// </summary>
        public void Dispose() {
            using(var tr = new FuncTrace("PETSC_Dispose")) {
            if(m_Disposed)
                return;

            Exception destroyError = null;
            if(m_Solver != IntPtr.Zero) {
                try {
                    MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
                    Check(m_Native.DestroySolver(m_Solver), "PETSc solver destruction");
                } catch(Exception error) {
                    destroyError = error;
                } finally {
                    m_Solver = IntPtr.Zero;
                }
            }

            ReleasePinnedCsrBuffers();
            m_Native?.Dispose();
            m_Native = null;
            m_KspMonitor = null;
            m_Disposed = true;
            if(destroyError != null)
                throw destroyError;
            }
        }

    }

    /// <summary>
    /// Options-driven PETSc Krylov solver for distributed BoSSS matrices.
    /// </summary>
    /// <remarks>
    /// PETSc KSP and preconditioner types are selected through <see cref="Options"/>.
    /// </remarks>
    public sealed class PETScKrylovSolver : PETScSolverBase {

        /// <summary>
        /// PETSc command-line-style options, for example
        /// <c>-ksp_type cg -pc_type jacobi -ksp_rtol 1e-10</c>.
        /// </summary>
        public string Options { get; set; } = "-ksp_type gmres -pc_type jacobi -ksp_rtol 1e-10 -ksp_atol 1e-50 -ksp_max_it 10000";

        /// <summary>
        /// Optional local portion of the first field's global row indices for PCFIELDSPLIT.
        /// Both field arrays must be supplied together.
        /// </summary>
        public int[] FirstFieldRows { get; set; }

        /// <summary>
        /// Optional local portion of the second field's global row indices for PCFIELDSPLIT.
        /// Both field arrays must be supplied together.
        /// </summary>
        public int[] SecondFieldRows { get; set; }

        /// <inheritdoc/>
        protected override string SolverOptions => Options;

        /// <inheritdoc/>
        protected override int[] SolverFirstFieldRows => FirstFieldRows;

        /// <inheritdoc/>
        protected override int[] SolverSecondFieldRows => SecondFieldRows;
    }
}
