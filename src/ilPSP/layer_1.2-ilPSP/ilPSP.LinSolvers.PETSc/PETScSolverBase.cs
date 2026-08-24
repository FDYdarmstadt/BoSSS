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

        /// <summary>
        /// Parallel runtime used by PETSc. This must be set before <see cref="DefineMatrix"/>.
        /// </summary>
        public PETScNativeVariant NativeVariant { get; set; } = PETScNativeVariant.Mpi;

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
        /// PETSc convergence reason returned by the most recent solve.
        /// Positive values indicate convergence; negative values indicate divergence.
        /// </summary>
        public int ConvergenceReason { get; private set; }

        /// <summary>
        /// Residual norm returned by PETSc for the most recent solve.
        /// </summary>
        public double ResidualNorm { get; private set; } = double.NaN;

        /// <summary>
        /// Defines and uploads a distributed square matrix to PETSc.
        /// </summary>
        public void DefineMatrix(IMutableMatrixEx matrix) {
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
            if(NativeVariant == PETScNativeVariant.OpenMp && mpiSize != 1)
                throw new NotSupportedException("The PETSc OpenMP variant uses MPIUNI and requires exactly one BoSSS MPI rank. Use MpiOpenMp for multiple ranks.");

            MPICollectiveWatchDog.Watch(matrix.MPI_Comm);

            long globalRows = matrix.RowPartitioning.TotalLength;
            long globalRowStart = matrix.RowPartitioning.i0;
            if(globalRows > int.MaxValue || globalRowStart > int.MaxValue)
                throw new NotSupportedException("The packaged PETSc uses 32-bit indices and this matrix exceeds that range.");

            ValidateFieldRows(SolverFirstFieldRows, SolverSecondFieldRows, globalRows);

            m_LocalRows = matrix.RowPartitioning.LocalLength;
            int[] rowPointers = new int[m_LocalRows + 1];
            var columns = new List<int>();
            var values = new List<double>();
            long[] rowColumns = null;
            double[] rowValues = null;

            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                long globalRow = globalRowStart + localRow;
                int rowLength = matrix.GetRow(globalRow, ref rowColumns, ref rowValues);
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    if(column < 0 || column >= globalRows || column > int.MaxValue)
                        throw new NotSupportedException($"PETSc column index {column} is outside the supported 32-bit matrix range.");
                    columns.Add((int)column);
                    values.Add(rowValues[entry]);
                }
                rowPointers[localRow + 1] = columns.Count;
            }

            m_Native = new UnsafePETSc(NativeVariant);
            int[] columnArray = columns.Count > 0 ? columns.ToArray() : new int[1];
            double[] valueArray = values.Count > 0 ? values.ToArray() : new double[1];
            try {
                // This records that BoSSS initialized PETSc, allowing application
                // shutdown to finalize PETSc before MPI is finalized.
                Check(m_Native.InitializeRuntime(), "PETSc initialization");
                Check(
                    m_Native.CreateKspSolver(
                        m_LocalRows,
                        (int)globalRows,
                        (int)globalRowStart,
                        rowPointers,
                        columnArray,
                        valueArray,
                        SolverOptions ?? string.Empty,
                        SolverFirstFieldRows?.Length ?? 0,
                        SolverFirstFieldRows,
                        SolverSecondFieldRows?.Length ?? 0,
                        SolverSecondFieldRows,
                        out m_Solver),
                    "PETSc matrix/KSP setup");
            } catch(Exception error) {
                m_Native.Dispose();
                m_Native = null;
                m_LocalRows = -1;
                Exception setupError = CreateSetupException(error);
                if(ReferenceEquals(setupError, error))
                    throw;
                throw setupError;
            }
        }

        /// <summary>
        /// Solves the matrix defined by <see cref="DefineMatrix"/>.
        /// </summary>
        public SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double> {
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
            Check(
                m_Native.Solve(
                    m_Solver,
                    nativeSolution,
                    nativeRhs,
                    useInitialGuess ? 1 : 0,
                    out int iterations,
                    out double residualNorm,
                    out int convergenceReason),
                "PETSc KSPSolve");
            stopwatch.Stop();

            ConvergenceReason = convergenceReason;
            ResidualNorm = residualNorm;
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
            if(double.IsNaN(residualNorm) || double.IsInfinity(residualNorm))
                throw new ArithmeticException($"PETSc returned non-finite residual {residualNorm}; convergence reason {convergenceReason}.");
            for(int i = 0; i < solution.Length; i++)
                if(double.IsNaN(solution[i]) || double.IsInfinity(solution[i]))
                    throw new ArithmeticException($"PETSc returned a non-finite solution at local index {i}; convergence reason {convergenceReason}.");
        }

        void Check(int errorCode, string operation) {
            if(errorCode == 0)
                return;
            IntPtr messagePointer = m_Native?.LastError() ?? IntPtr.Zero;
            string message = messagePointer == IntPtr.Zero ? "No PETSc error message was provided." : Marshal.PtrToStringAnsi(messagePointer);
            throw new ApplicationException($"{operation} failed with PETSc error code {errorCode}: {message}");
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

            m_Native?.Dispose();
            m_Native = null;
            m_Disposed = true;
            if(destroyError != null)
                throw destroyError;
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
