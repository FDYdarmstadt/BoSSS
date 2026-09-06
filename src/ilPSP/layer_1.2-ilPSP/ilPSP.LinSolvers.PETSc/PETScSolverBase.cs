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
    public abstract class PETScSolverBase : ISparseSolver {

        UnsafePETSc m_Native;
        IntPtr m_Solver;
        int m_LocalRows = -1;
        bool m_Disposed;
        int[] m_PinnedRowPointers;
        int[] m_PinnedColumnIndices;
        double[] m_PinnedValues;
        GCHandle m_RowPointersHandle;
        GCHandle m_ColumnIndicesHandle;
        GCHandle m_ValuesHandle;
        int[] m_PinnedDiagonalRowPointers;
        int[] m_PinnedDiagonalColumnIndices;
        double[] m_PinnedDiagonalValues;
        int[] m_PinnedOffdiagonalRowPointers;
        int[] m_PinnedOffdiagonalColumnIndices;
        double[] m_PinnedOffdiagonalValues;
        GCHandle m_DiagonalRowPointersHandle;
        GCHandle m_DiagonalColumnIndicesHandle;
        GCHandle m_DiagonalValuesHandle;
        GCHandle m_OffdiagonalRowPointersHandle;
        GCHandle m_OffdiagonalColumnIndicesHandle;
        GCHandle m_OffdiagonalValuesHandle;
        UnsafePETSc._bosss_petsc_ksp_monitor m_KspMonitor;
        UnsafePETSc._bosss_petsc_apply m_MatrixFreeApply;
        IPETScMatrixFreeOperator m_MatrixFreeOperator;
        PETScVectorView m_OwnedInputView;
        PETScVectorView m_HaloInputView;
        PETScVectorView m_OutputView;
        Exception m_MatrixFreeCallbackError;

        /// <summary>
        /// Selects whether PETSc receives an owned CSR matrix or an iterative
        /// MatShell operator. The default is the split, borrowed CSR path.
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
        /// Whether this concrete solver can use a PETSc MatShell. Sparse direct
        /// factorization requires an assembled PETSc matrix and therefore returns false.
        /// </summary>
        protected virtual bool SupportsMatrixShell => true;

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
            PETScNativeVariant effectiveNativeVariant = ResolveNativeVariant(
                NativeVariant,
                mpiSize,
                ilPSP.Environment.NumThreads);
            EffectiveNativeVariant = effectiveNativeVariant;
            if(effectiveNativeVariant == PETScNativeVariant.OpenMp && mpiSize != 1)
                throw new NotSupportedException("The PETSc OpenMP variant uses MPIUNI and requires exactly one BoSSS MPI rank. Use MpiOpenMp for multiple ranks.");
            ValidateParallelConfiguration(mpiSize, effectiveNativeVariant);

            MPICollectiveWatchDog.Watch(matrix.MPI_Comm);

            long globalRows = matrix.RowPartitioning.TotalLength;
            long globalRowStart = matrix.RowPartitioning.i0;
            if(globalRows > int.MaxValue || globalRowStart > int.MaxValue)
                throw new NotSupportedException("The packaged PETSc uses 32-bit indices and this matrix exceeds that range.");

            ValidateFieldRows(SolverFirstFieldRows, SolverSecondFieldRows, globalRows);

            if(MatrixInput == PETScMatrixInput.MatrixFreeShell)
                throw new InvalidOperationException("MatrixFreeShell requires DefineMatrixFree(IPETScMatrixFreeOperator), not DefineMatrix.");
            if(MatrixInput == PETScMatrixInput.BlockMsrMatrixShell) {
                EffectiveMatrixInput = PETScMatrixInput.BlockMsrMatrixShell;
                DefineMatrixFreeCore(new AssembledMatrixShellOperator(matrix), effectiveNativeVariant, tr);
                return;
            }

            m_LocalRows = matrix.RowPartitioning.LocalLength;
            // MatCreateMPIAIJWithSplitArrays deliberately creates MATMPIAIJ even
            // for one rank. Top-level ILU/ICC and the UMFPACK factor package are
            // sequential and require MATSEQAIJ. MatCreateAIJ in the copied path
            // selects SEQAIJ for a one-rank communicator, so use it only for this
            // combination. MPI runs with more than one rank retain their requested
            // input path; PETScDirectSolver rejects top-level UMFPACK for such a run.
            bool requiresSequentialFactorMatrix = mpiSize == 1
                && RequestsSequentialFactorMatrix(SolverOptions);
            bool splitCsr = MatrixInput == PETScMatrixInput.BorrowedCsr
                && !requiresSequentialFactorMatrix;
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
                    Check(m_Native.CreateKspForMatrix(m_Solver), "PETSc KSP creation");
                }
                using(new BlockTrace("PETSC_NativeOptions", tr)) {
                    Check(m_Native.ConfigureKspSolver(
                        m_Solver,
                        SolverOptions ?? string.Empty,
                        SolverFirstFieldRows?.Length ?? 0,
                        SolverFirstFieldRows,
                        SolverSecondFieldRows?.Length ?? 0,
                        SolverSecondFieldRows), "PETSc KSP option configuration");
                }
                using(new BlockTrace("PETSC_KspSetup", tr)) {
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

        /// <summary>
        /// Defines a PETSc MatShell from a BoSSS matrix-free operator. This path is
        /// intended for iterative solvers; direct PETSc factorizations require CSR.
        /// </summary>
        public void DefineMatrixFree(IPETScMatrixFreeOperator matrixFreeOperator) {
            using(var tr = new FuncTrace("PETSC_DefineMatrix")) {
                if(MatrixInput != PETScMatrixInput.MatrixFreeShell)
                    throw new InvalidOperationException("Set MatrixInput to MatrixFreeShell before calling DefineMatrixFree.");
                if(matrixFreeOperator == null)
                    throw new ArgumentNullException(nameof(matrixFreeOperator));
                ThrowIfDisposed();
                if(m_Solver != IntPtr.Zero)
                    throw new InvalidOperationException("DefineMatrix can only be called once in the lifetime of a PETSc solver.");
                if(matrixFreeOperator.GlobalLength <= 0 || matrixFreeOperator.GlobalLength > int.MaxValue
                    || matrixFreeOperator.GlobalRowStart < 0 || matrixFreeOperator.GlobalRowStart > int.MaxValue
                    || matrixFreeOperator.LocalLength < 0
                    || matrixFreeOperator.GlobalRowStart + matrixFreeOperator.LocalLength > matrixFreeOperator.GlobalLength)
                    throw new ArgumentException("The matrix-free operator has an invalid 32-bit distributed row partition.", nameof(matrixFreeOperator));
                if(!matrixFreeOperator.MPIComm.Equals(csMPI.Raw._COMM.WORLD))
                    throw new NotSupportedException("The current PETSc native ABI supports MPI_COMM_WORLD only.");

                csMPI.Raw.Comm_Size(matrixFreeOperator.MPIComm, out int mpiSize);
                PETScNativeVariant effectiveNativeVariant = ResolveNativeVariant(
                    NativeVariant,
                    mpiSize,
                    ilPSP.Environment.NumThreads);
                EffectiveNativeVariant = effectiveNativeVariant;
                if(effectiveNativeVariant == PETScNativeVariant.OpenMp && mpiSize != 1)
                    throw new NotSupportedException("The PETSc OpenMP variant uses MPIUNI and requires exactly one BoSSS MPI rank. Use MpiOpenMp for multiple ranks.");
                ValidateParallelConfiguration(mpiSize, effectiveNativeVariant);
                ValidateFieldRows(SolverFirstFieldRows, SolverSecondFieldRows, matrixFreeOperator.GlobalLength);
                DefineMatrixFreeCore(matrixFreeOperator, effectiveNativeVariant, tr);
            }
        }

        void DefineMatrixFreeCore(IPETScMatrixFreeOperator matrixFreeOperator, PETScNativeVariant effectiveNativeVariant, Tmeas parentTrace) {
            if(!SupportsMatrixShell)
                throw new NotSupportedException($"{GetType().Name} requires an assembled CSR matrix and cannot use a PETSc MatShell.");
            if((SolverFirstFieldRows?.Length ?? 0) != 0 || (SolverSecondFieldRows?.Length ?? 0) != 0)
                throw new NotSupportedException("PETSc field splits are not yet implemented for MatShell operators.");
            ValidateMatrixShellOptions(SolverOptions);

            int[] ghostGlobalIds = matrixFreeOperator.GhostGlobalIds ?? Array.Empty<int>();
            long globalLength = matrixFreeOperator.GlobalLength;
            long globalRowStart = matrixFreeOperator.GlobalRowStart;
            foreach(int ghostGlobalId in ghostGlobalIds)
                if(ghostGlobalId < 0 || ghostGlobalId >= globalLength)
                    throw new ArgumentOutOfRangeException(nameof(matrixFreeOperator), $"Ghost index {ghostGlobalId} is outside [0,{globalLength}).");

            MPICollectiveWatchDog.Watch(matrixFreeOperator.MPIComm);
            m_LocalRows = matrixFreeOperator.LocalLength;
            m_MatrixFreeOperator = matrixFreeOperator;
            m_MatrixFreeApply = ApplyMatrixFree;
            try {
                using(new BlockTrace("PETSC_NativeMatrixCreate", parentTrace)) {
                    m_Native = new UnsafePETSc(effectiveNativeVariant);
                    Check(m_Native.InitializeRuntime(), "PETSc initialization");
                    Check(m_Native.CreateShellSolver(
                        m_LocalRows,
                        (int)globalLength,
                        (int)globalRowStart,
                        ghostGlobalIds.Length,
                        ghostGlobalIds,
                        SolverOptions ?? string.Empty,
                        m_MatrixFreeApply,
                        IntPtr.Zero,
                        out m_Solver), "PETSc MatShell creation");
                }
                using(new BlockTrace("PETSC_KspSetup", parentTrace)) {
                    // Force PETSc to validate the selected KSP/PC while DefineMatrix is
                    // still on the stack. Otherwise KSP lazily discovers an unsupported
                    // MatShell operation only from KSPSolve, which is much less useful.
                    Check(m_Native.SetupKspSolver(m_Solver), "PETSc MatShell KSP setup");
                    if(RequestsKspMonitor(SolverOptions)) {
                        m_KspMonitor = WriteKspMonitor;
                        Check(m_Native.SetKspMonitor(m_Solver, m_KspMonitor, IntPtr.Zero), "PETSc KSP monitor setup");
                    }
                }
            } catch {
                if(m_Solver != IntPtr.Zero && m_Native != null) {
                    m_Native.DestroySolver(m_Solver);
                    m_Solver = IntPtr.Zero;
                }
                m_Native?.Dispose();
                m_Native = null;
                m_MatrixFreeApply = null;
                m_MatrixFreeOperator = null;
                m_LocalRows = -1;
                throw;
            }
        }

        static void ValidateMatrixShellOptions(string options) {
            string pcType = GetPetscOptionValue(options, "-pc_type");
            if(pcType != null && !pcType.Equals("none", StringComparison.OrdinalIgnoreCase)) {
                throw new NotSupportedException(
                    $"PETSc MatShell currently supports only '-pc_type none', but the configuration requests '-pc_type {pcType}'. "
                    + "ASM, ILU, GAMG, Jacobi, and direct factor PCs require an assembled PETSc matrix. "
                    + "Use CopiedCsr or BorrowedCsr for those preconditioners, or select '-pc_type none' for BlockMsrMatrixShell/MatrixFreeShell.");
            }
        }

        static string GetPetscOptionValue(string options, string option) {
            if(string.IsNullOrWhiteSpace(options))
                return null;
            string[] tokens = options.Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            for(int index = 0; index < tokens.Length; index++) {
                if(tokens[index].StartsWith(option + "=", StringComparison.Ordinal))
                    return tokens[index].Substring(option.Length + 1);
                if(tokens[index].Equals(option, StringComparison.Ordinal) && index + 1 < tokens.Length)
                    return tokens[index + 1];
            }
            return null;
        }

        static bool RequestsSequentialFactorMatrix(string options) {
            string pcType = GetPetscOptionValue(options, "-pc_type");
            if(string.Equals(pcType, "ilu", StringComparison.OrdinalIgnoreCase)
                || string.Equals(pcType, "icc", StringComparison.OrdinalIgnoreCase))
                return true;

            return string.Equals(
                GetPetscOptionValue(options, "-pc_factor_mat_solver_type"),
                "umfpack",
                StringComparison.OrdinalIgnoreCase);
        }

        int ApplyMatrixFree(
            IntPtr context,
            IntPtr ownedInput,
            int ownedLength,
            IntPtr haloInput,
            int haloLength,
            IntPtr output) {
            try {
                if(m_MatrixFreeOperator == null || ownedLength != m_LocalRows || output == IntPtr.Zero)
                    return 1;
                if(haloLength < 0 || ownedInput == IntPtr.Zero && ownedLength != 0 || haloInput == IntPtr.Zero && haloLength != 0)
                    return 1;
                if(m_OwnedInputView == null) {
                    m_OwnedInputView = new PETScVectorView();
                    m_HaloInputView = new PETScVectorView();
                    m_OutputView = new PETScVectorView();
                }
                m_OwnedInputView.Reset(ownedInput, ownedLength);
                m_HaloInputView.Reset(haloInput, haloLength);
                m_OutputView.Reset(output, ownedLength);
                m_MatrixFreeOperator.Apply(m_OwnedInputView, m_HaloInputView, m_OutputView);
                return 0;
            } catch(Exception error) {
                m_MatrixFreeCallbackError = error;
                return 1;
            }
        }

        sealed class AssembledMatrixShellOperator : IPETScMatrixFreeOperator {
            readonly IMutableMatrixEx m_Matrix;
            readonly int[] m_GhostGlobalIds;
            readonly Dictionary<int, int> m_GhostMap;
            long[] m_Columns;
            double[] m_Values;

            public AssembledMatrixShellOperator(IMutableMatrixEx matrix) {
                m_Matrix = matrix ?? throw new ArgumentNullException(nameof(matrix));
                if(matrix.RowPartitioning.TotalLength != matrix.ColPartitioning.TotalLength
                    || matrix.RowPartitioning.i0 != matrix.ColPartitioning.i0
                    || matrix.RowPartitioning.LocalLength != matrix.ColPartitioning.LocalLength)
                    throw new NotSupportedException("An assembled PETSc MatShell requires identical row and column partitionings.");

                var ghosts = new HashSet<int>();
                long[] columns = null;
                long firstOwned = matrix.ColPartitioning.i0;
                long lastOwned = firstOwned + matrix.ColPartitioning.LocalLength;
                for(int localRow = 0; localRow < matrix.RowPartitioning.LocalLength; localRow++) {
                    int length = matrix.GetOccupiedColumnIndices(matrix.RowPartitioning.i0 + localRow, ref columns);
                    for(int entry = 0; entry < length; entry++) {
                        long column = columns[entry];
                        if(column < 0 || column >= matrix.ColPartitioning.TotalLength || column > int.MaxValue)
                            throw new NotSupportedException($"PETSc column index {column} is outside the supported 32-bit matrix range.");
                        if(column < firstOwned || column >= lastOwned)
                            ghosts.Add((int)column);
                    }
                }
                m_GhostGlobalIds = new List<int>(ghosts).ToArray();
                Array.Sort(m_GhostGlobalIds);
                m_GhostMap = new Dictionary<int, int>(m_GhostGlobalIds.Length);
                for(int i = 0; i < m_GhostGlobalIds.Length; i++)
                    m_GhostMap.Add(m_GhostGlobalIds[i], i);
            }

            public MPI_Comm MPIComm => m_Matrix.MPI_Comm;
            public long GlobalLength => m_Matrix.RowPartitioning.TotalLength;
            public long GlobalRowStart => m_Matrix.RowPartitioning.i0;
            public int LocalLength => m_Matrix.RowPartitioning.LocalLength;
            public int[] GhostGlobalIds => m_GhostGlobalIds;

            public void Apply(PETScVectorView ownedInput, PETScVectorView haloInput, PETScVectorView output) {
                long firstOwned = GlobalRowStart;
                long lastOwned = firstOwned + LocalLength;
                for(int localRow = 0; localRow < LocalLength; localRow++) {
                    int length = m_Matrix.GetRow(firstOwned + localRow, ref m_Columns, ref m_Values);
                    double accumulator = 0.0;
                    for(int entry = 0; entry < length; entry++) {
                        long column = m_Columns[entry];
                        double input;
                        if(column >= firstOwned && column < lastOwned) {
                            input = ownedInput[(int)(column - firstOwned)];
                        } else if(column >= int.MinValue && column <= int.MaxValue
                            && m_GhostMap.TryGetValue((int)column, out int ghostIndex)) {
                            input = haloInput[ghostIndex];
                        } else {
                            throw new InvalidOperationException($"Matrix entry ({firstOwned + localRow},{column}) is missing from the PETSc shell halo map.");
                        }
                        accumulator += m_Values[entry] * input;
                    }
                    output[localRow] = accumulator;
                }
            }
        }

        /// <summary>
        /// Solves the matrix defined by <see cref="DefineMatrix"/>.
        /// </summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
            using(var tr = new FuncTrace("PETSC_Solve")) {
            tr.InfoToConsole = true;
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
            m_MatrixFreeCallbackError = null;
            int solveError = m_Native.Solve(
                m_Solver,
                nativeSolution,
                nativeRhs,
                useInitialGuess ? 1 : 0,
                out iterations,
                out residualNorm,
                out convergenceReason);
            stopwatch.Stop();

            Exception matrixFreeCallbackError = m_MatrixFreeCallbackError;
            m_MatrixFreeCallbackError = null;
            if(matrixFreeCallbackError != null)
                throw new ApplicationException("The PETSc matrix-free callback failed.", matrixFreeCallbackError);

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
            m_MatrixFreeApply = null;
            m_MatrixFreeOperator = null;
            m_OwnedInputView = null;
            m_HaloInputView = null;
            m_OutputView = null;
            m_MatrixFreeCallbackError = null;
            m_Disposed = true;
            if(destroyError != null)
                throw destroyError;
            }
        }

        sealed class CsrArrays {
            internal int[] RowPointers;
            internal int[] ColumnIndices;
            internal double[] Values;
        }

        sealed class SplitCsrArrays {
            internal int[] DiagonalRowPointers;
            internal int[] DiagonalColumnIndices;
            internal double[] DiagonalValues;
            internal int[] OffdiagonalRowPointers;
            internal int[] OffdiagonalColumnIndices;
            internal double[] OffdiagonalValues;
        }

        CsrArrays ExportCsr(IMutableMatrixEx matrix, long globalRows, long globalRowStart) {
            var result = new CsrArrays {
                RowPointers = new int[m_LocalRows + 1]
            };
            var columns = new List<int>();
            var values = new List<double>();
            long[] rowColumns = null;
            double[] rowValues = null;
            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                int rowLength = matrix.GetRow(globalRowStart + localRow, ref rowColumns, ref rowValues);
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    ValidateCsrColumn(column, globalRows);
                    columns.Add((int)column);
                    values.Add(rowValues[entry]);
                }
                result.RowPointers[localRow + 1] = columns.Count;
            }
            result.ColumnIndices = columns.Count > 0 ? columns.ToArray() : new int[1];
            result.Values = values.Count > 0 ? values.ToArray() : new double[1];
            return result;
        }

        SplitCsrArrays ExportSplitCsr(IMutableMatrixEx matrix, long globalRows, long globalRowStart) {
            if(matrix.RowPartitioning.TotalLength != matrix.ColPartitioning.TotalLength
                || matrix.RowPartitioning.i0 != matrix.ColPartitioning.i0
                || matrix.RowPartitioning.LocalLength != matrix.ColPartitioning.LocalLength)
                throw new NotSupportedException("PETSc split CSR requires identical row and column partitionings.");

            long firstOwnedColumn = matrix.ColPartitioning.i0;
            long lastOwnedColumn = firstOwnedColumn + matrix.ColPartitioning.LocalLength;
            var result = new SplitCsrArrays {
                DiagonalRowPointers = new int[m_LocalRows + 1],
                OffdiagonalRowPointers = new int[m_LocalRows + 1]
            };
            long[] rowColumns = null;
            double[] rowValues = null;

            // Count first: PETSc requires one contiguous CSR array for its local and
            // remote portions. Add a structural zero diagonal when BoSSS has no entry
            // there; this permits ILU/ASM and direct factorization on saddle systems.
            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                long globalRow = globalRowStart + localRow;
                int rowLength = matrix.GetRow(globalRow, ref rowColumns, ref rowValues);
                int diagonalCount = 0;
                int offdiagonalCount = 0;
                bool hasStructuralDiagonal = false;
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    ValidateCsrColumn(column, globalRows);
                    if(column >= firstOwnedColumn && column < lastOwnedColumn)
                        diagonalCount++;
                    else
                        offdiagonalCount++;
                    hasStructuralDiagonal |= column == globalRow;
                }
                if(!hasStructuralDiagonal)
                    diagonalCount++;
                result.DiagonalRowPointers[localRow + 1] = checked(result.DiagonalRowPointers[localRow] + diagonalCount);
                result.OffdiagonalRowPointers[localRow + 1] = checked(result.OffdiagonalRowPointers[localRow] + offdiagonalCount);
            }

            int diagonalNonzeros = result.DiagonalRowPointers[m_LocalRows];
            int offdiagonalNonzeros = result.OffdiagonalRowPointers[m_LocalRows];
            result.DiagonalColumnIndices = new int[Math.Max(diagonalNonzeros, 1)];
            result.DiagonalValues = new double[Math.Max(diagonalNonzeros, 1)];
            result.OffdiagonalColumnIndices = new int[Math.Max(offdiagonalNonzeros, 1)];
            result.OffdiagonalValues = new double[Math.Max(offdiagonalNonzeros, 1)];
            int[] diagonalOffsets = (int[])result.DiagonalRowPointers.Clone();
            int[] offdiagonalOffsets = (int[])result.OffdiagonalRowPointers.Clone();

            // Populate the packed PETSc layout. Diagonal columns are local indices,
            // while off-diagonal columns remain global PETSc indices.
            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                long globalRow = globalRowStart + localRow;
                int rowLength = matrix.GetRow(globalRow, ref rowColumns, ref rowValues);
                bool hasStructuralDiagonal = false;
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    if(column >= firstOwnedColumn && column < lastOwnedColumn) {
                        int offset = diagonalOffsets[localRow]++;
                        result.DiagonalColumnIndices[offset] = checked((int)(column - firstOwnedColumn));
                        result.DiagonalValues[offset] = rowValues[entry];
                    } else {
                        int offset = offdiagonalOffsets[localRow]++;
                        result.OffdiagonalColumnIndices[offset] = checked((int)column);
                        result.OffdiagonalValues[offset] = rowValues[entry];
                    }
                    hasStructuralDiagonal |= column == globalRow;
                }
                if(!hasStructuralDiagonal) {
                    int offset = diagonalOffsets[localRow]++;
                    result.DiagonalColumnIndices[offset] = localRow;
                    result.DiagonalValues[offset] = 0.0;
                }
            }
            return result;
        }

        static void ValidateCsrColumn(long column, long globalRows) {
            if(column < 0 || column >= globalRows || column > int.MaxValue)
                throw new NotSupportedException($"PETSc column index {column} is outside the supported 32-bit matrix range.");
        }

        void PinCsrBuffers(int[] rowPointers, int[] columnIndices, double[] values) {
            m_PinnedRowPointers = rowPointers;
            m_PinnedColumnIndices = columnIndices;
            m_PinnedValues = values;
            m_RowPointersHandle = GCHandle.Alloc(m_PinnedRowPointers, GCHandleType.Pinned);
            m_ColumnIndicesHandle = GCHandle.Alloc(m_PinnedColumnIndices, GCHandleType.Pinned);
            m_ValuesHandle = GCHandle.Alloc(m_PinnedValues, GCHandleType.Pinned);
        }

        void PinSplitCsrBuffers(SplitCsrArrays csr) {
            m_PinnedDiagonalRowPointers = csr.DiagonalRowPointers;
            m_PinnedDiagonalColumnIndices = csr.DiagonalColumnIndices;
            m_PinnedDiagonalValues = csr.DiagonalValues;
            m_PinnedOffdiagonalRowPointers = csr.OffdiagonalRowPointers;
            m_PinnedOffdiagonalColumnIndices = csr.OffdiagonalColumnIndices;
            m_PinnedOffdiagonalValues = csr.OffdiagonalValues;
            m_DiagonalRowPointersHandle = GCHandle.Alloc(m_PinnedDiagonalRowPointers, GCHandleType.Pinned);
            m_DiagonalColumnIndicesHandle = GCHandle.Alloc(m_PinnedDiagonalColumnIndices, GCHandleType.Pinned);
            m_DiagonalValuesHandle = GCHandle.Alloc(m_PinnedDiagonalValues, GCHandleType.Pinned);
            m_OffdiagonalRowPointersHandle = GCHandle.Alloc(m_PinnedOffdiagonalRowPointers, GCHandleType.Pinned);
            m_OffdiagonalColumnIndicesHandle = GCHandle.Alloc(m_PinnedOffdiagonalColumnIndices, GCHandleType.Pinned);
            m_OffdiagonalValuesHandle = GCHandle.Alloc(m_PinnedOffdiagonalValues, GCHandleType.Pinned);
        }

        void ReleasePinnedCsrBuffers() {
            if(m_RowPointersHandle.IsAllocated)
                m_RowPointersHandle.Free();
            if(m_ColumnIndicesHandle.IsAllocated)
                m_ColumnIndicesHandle.Free();
            if(m_ValuesHandle.IsAllocated)
                m_ValuesHandle.Free();
            if(m_DiagonalRowPointersHandle.IsAllocated)
                m_DiagonalRowPointersHandle.Free();
            if(m_DiagonalColumnIndicesHandle.IsAllocated)
                m_DiagonalColumnIndicesHandle.Free();
            if(m_DiagonalValuesHandle.IsAllocated)
                m_DiagonalValuesHandle.Free();
            if(m_OffdiagonalRowPointersHandle.IsAllocated)
                m_OffdiagonalRowPointersHandle.Free();
            if(m_OffdiagonalColumnIndicesHandle.IsAllocated)
                m_OffdiagonalColumnIndicesHandle.Free();
            if(m_OffdiagonalValuesHandle.IsAllocated)
                m_OffdiagonalValuesHandle.Free();
            m_PinnedRowPointers = null;
            m_PinnedColumnIndices = null;
            m_PinnedValues = null;
            m_PinnedDiagonalRowPointers = null;
            m_PinnedDiagonalColumnIndices = null;
            m_PinnedDiagonalValues = null;
            m_PinnedOffdiagonalRowPointers = null;
            m_PinnedOffdiagonalColumnIndices = null;
            m_PinnedOffdiagonalValues = null;
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
