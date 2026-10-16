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
using System.Globalization;
using System.Linq;
using System.Runtime.Serialization;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PETSc;
using MPI.Wrappers;

namespace BoSSS.Solution.AdvancedSolvers {


    /// <summary>
    /// Typed Level-3 configuration for a PETSc sparse direct solve.
    /// </summary>
    /// <remarks>
    /// The selected factor backend must be compiled into the native PETSc
    /// variant. The default backend supports the packaged single-rank
    /// factorization; distributed direct solves require MUMPS or SuperLU_DIST.
    /// </remarks>
    [Serializable]
    public sealed class PETScDirectSolverConfig : ISolverFactory {

        /// <summary>
        /// Native PETSc parallel runtime.
        /// </summary>
        [DataMember]
        public PETScNativeVariant NativeVariant = PETScNativeVariant.Auto;

        /// <summary>
        /// PETSc matrix representation. Borrowed split CSR retains the backing
        /// buffers for PETSc's lifetime and avoids a second native CSR copy.
        /// </summary>
        [DataMember]
        public PETScMatrixInput MatrixInput = PETScMatrixInput.BorrowedCsr;

        /// <summary>
        /// LU or Cholesky factorization.
        /// </summary>
        [DataMember]
        public PETScDirectFactorization Factorization = PETScDirectFactorization.LU;

        /// <summary>
        /// PETSc factorization package.
        /// </summary>
        [DataMember]
        public PETScDirectBackend Backend = PETScDirectBackend.Default;

        /// <summary>
        /// Optional package-specific PETSc factorization options.
        /// </summary>
        [DataMember]
        public string AdditionalOptions = string.Empty;

        /// <summary>
        /// Human-readable solver name.
        /// </summary>
        public string Name => $"PETSc direct {Factorization} ({Backend}; {PETScKrylovSolverConfig.DescribeNativeVariant(NativeVariant)}; {MatrixInput})";

        /// <summary>
        /// Compact solver name.
        /// </summary>
        public string Shortname => "PETSc-" + Factorization + "-" + Backend;

        /// <summary>
        /// Creates and initializes a direct PETSc solver for the supplied operator level.
        /// </summary>
        public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            ValidateConfiguration(level);

            var solver = new PETScDirectSolverAdapter {
                NativeVariant = NativeVariant,
                MatrixInput = MatrixInput,
                Factorization = Factorization,
                Backend = Backend,
                AdditionalOptions = AdditionalOptions
            };
            solver.Init(level);
            return solver;
        }

        /// <summary>
        /// Validates combinations whose support is determined by the selected PETSc
        /// package and parallel runtime, before PETSc starts KSP setup.
        /// </summary>
        void ValidateConfiguration(MultigridOperator level) {
            if(level == null)
                throw new ArgumentNullException(nameof(level));

            if(MatrixInput != PETScMatrixInput.BorrowedCsr && MatrixInput != PETScMatrixInput.CopiedCsr)
                throw new NotSupportedException($"Unknown PETSc direct matrix input mode '{MatrixInput}'.");

            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out int mpiSize);
            PETScNativeVariant effectiveVariant = PETScNativeVariantResolver.Resolve(
                NativeVariant, mpiSize, ilPSP.Environment.NumThreads, csMPI.Raw.Initialized());

            if(Backend == PETScDirectBackend.Umfpack && mpiSize > 1)
                throw new NotSupportedException(
                    "PETSc direct UMFPACK is a sequential factorization and requires exactly one MPI rank. "
                    + "Use PETSc direct MUMPS or SuperLU_DIST for a distributed direct solve, "
                    + "or use ASM with UMFPACK local subdomain factorizations.");

            if(Factorization == PETScDirectFactorization.Cholesky
                && (Backend == PETScDirectBackend.Umfpack || Backend == PETScDirectBackend.SuperLuDist))
                throw new NotSupportedException($"PETSc direct {Backend} supports LU factorization only.");

            if(Backend == PETScDirectBackend.SuperLuDist && effectiveVariant == PETScNativeVariant.OpenMp)
                throw new NotSupportedException(
                    "PETSc direct SuperLU_DIST requires an MPI-enabled PETSc variant. "
                    + "Select the MPI or MPI+OpenMP variant, or run with more than one MPI rank when using Auto.");
        }

        /// <summary>
        /// Compares two direct PETSc configurations.
        /// </summary>
        public bool Equals(ISolverFactory other) {
            return Equals(other as object);
        }

        /// <summary>
        /// Compares two direct PETSc configurations.
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as PETScDirectSolverConfig;
            return other != null
                && NativeVariant == other.NativeVariant
                && MatrixInput == other.MatrixInput
                && Factorization == other.Factorization
                && Backend == other.Backend
                && string.Equals(AdditionalOptions, other.AdditionalOptions, StringComparison.Ordinal);
        }

        /// <summary>
        /// Returns a hash code for this configuration.
        /// </summary>
        public override int GetHashCode() {
            unchecked {
                int hash = (int)NativeVariant;
                hash = hash * 397 ^ (int)MatrixInput;
                hash = hash * 397 ^ (int)Factorization;
                hash = hash * 397 ^ (int)Backend;
                hash = hash * 397 ^ (AdditionalOptions?.GetHashCode() ?? 0);
                return hash;
            }
        }
    }


    /// <summary>
    /// Level-3 adapter for one direct PETSc ABI solve.
    /// </summary>
    /// <remarks>
    /// This adapter retains the established Level-3 PETSc direct-solver API
    /// and delegates to the direct PETSc ABI implementation.
    /// </remarks>
    public sealed class PETScDirectSolverAdapter : ISolverSmootherTemplate {
        IPETScDirectSolverAbi m_Solver;
        MultigridOperator m_Operator;
        int m_Iterations;
        bool m_Converged;

        /// <summary>Native PETSc parallel runtime.</summary>
        public PETScNativeVariant NativeVariant { get; set; } = PETScNativeVariant.Auto;

        /// <summary>PETSc matrix ownership and transfer representation.</summary>
        public PETScMatrixInput MatrixInput { get; set; } = PETScMatrixInput.BorrowedCsr;

        /// <summary>LU or Cholesky factorization.</summary>
        public PETScDirectFactorization Factorization { get; set; } = PETScDirectFactorization.LU;

        /// <summary>PETSc factorization package.</summary>
        public PETScDirectBackend Backend { get; set; } = PETScDirectBackend.Default;

        /// <summary>Optional PETSc package-specific options.</summary>
        public string AdditionalOptions { get; set; } = string.Empty;

        /// <summary>Reason returned by PETSc for the most recent solve.</summary>
        public int ConvergenceReason => m_Solver?.ConvergenceReason ?? 0;

        /// <summary>Residual norm returned by PETSc for the most recent solve.</summary>
        public double ResidualNorm => m_Solver?.ResidualNorm ?? double.NaN;

        public int IterationsInNested => 0;
        public int ThisLevelIterations => m_Iterations;
        public bool Converged => m_Converged;

        public void Init(MultigridOperator op) {
            if(op == null)
                throw new ArgumentNullException(nameof(op));
            if(ReferenceEquals(op, m_Operator) && m_Solver != null)
                return;

            DisposeSolver();
            var solver = PETScDirectSolverAbi.Create();
            solver.NativeVariant = NativeVariant;
            solver.MatrixInput = MatrixInput;
            solver.Factorization = Factorization;
            solver.Backend = Backend;
            solver.AdditionalOptions = AdditionalOptions ?? string.Empty;
            try {
                solver.DefineMatrix(op.OperatorMatrix);
            } catch {
                solver.Dispose();
                throw;
            }
            m_Solver = solver;
            m_Operator = op;
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            if(m_Solver == null)
                throw new InvalidOperationException("PETSc direct solver is not initialized.");

            var result = m_Solver.Solve(X, B);
            m_Iterations += result.NoOfIterations;
            m_Converged = result.Converged;
        }

        public void ResetStat() {
            m_Iterations = 0;
            m_Converged = false;
        }

        public long UsedMemory() => 0;

        public object Clone() {
            return new PETScDirectSolverAdapter {
                NativeVariant = NativeVariant,
                MatrixInput = MatrixInput,
                Factorization = Factorization,
                Backend = Backend,
                AdditionalOptions = AdditionalOptions
            };
        }

        public void Dispose() {
            DisposeSolver();
            m_Operator = null;
        }

        void DisposeSolver() {
            m_Solver?.Dispose();
            m_Solver = null;
        }
    }



    /// <summary>
    /// Configuration for the options-driven PETSc Krylov solver adapter.
    /// </summary>
    [Serializable]
    public sealed class PETScKrylovSolverConfig : IterativeSolverConfig {

        /// <summary>
        /// PETSc command-line-style options. The KSP and preconditioner are not fixed by
        /// this adapter; for example, use <c>-ksp_type cg -pc_type jacobi</c>.
        /// The inherited <see cref="IterativeSolverConfig"/> fields are translated to
        /// PETSc options. Options supplied here take precedence; a warning is written
        /// when an option is also supplied by an inherited field.
        /// </summary>
        [DataMember]
        public string PETScOptions = "-ksp_type gmres -pc_type jacobi";

        /// <summary>
        /// Native PETSc parallel runtime. The default selects a variant automatically.
        /// </summary>
        [DataMember]
        public PETScNativeVariant NativeVariant = PETScNativeVariant.Auto;

        /// <summary>
        /// PETSc matrix representation used by this Krylov solver.
        /// </summary>
        [DataMember]
        public PETScMatrixInput MatrixInput = PETScMatrixInput.BorrowedCsr;

        /// <summary>
        /// Optional BoSSS variable indices that form PETSc field split 0. Both split arrays
        /// must be non-empty when <c>-pc_type fieldsplit</c> is used.
        /// </summary>
        [DataMember]
        public int[] FirstFieldVariables;

        /// <summary>
        /// Optional BoSSS variable indices that form PETSc field split 1.
        /// </summary>
        [DataMember]
        public int[] SecondFieldVariables;

        /// <summary>
        /// Uses BoSSS's explicit geometric hierarchy as a PETSc <c>PCMG</c>
        /// preconditioner instead of a single-level PETSc matrix.
        /// </summary>
        [DataMember]
        public bool UseGeometricMultigrid;

        /// <summary>
        /// Stops descending the BoSSS hierarchy once this many global unknowns
        /// or fewer remain; that level becomes PETSc's coarse level.
        /// </summary>
        [DataMember]
        public int CoarseKickIn = 90000;

        public override string Name => $"PETSc Krylov solver ({DescribeNativeVariant(NativeVariant)}; {MatrixInput}; {ComposeOptions()})";

        public override string Shortname => "PETSc-Krylov";

        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            var solver = new PETScKrylovSolverAdapter {
                Options = ComposeOptions(),
                NativeVariant = NativeVariant,
                MatrixInput = MatrixInput,
                FirstFieldVariables = FirstFieldVariables == null ? null : (int[])FirstFieldVariables.Clone(),
                SecondFieldVariables = SecondFieldVariables == null ? null : (int[])SecondFieldVariables.Clone(),
                UseGeometricMultigrid = UseGeometricMultigrid,
                CoarseKickIn = CoarseKickIn,
                EnableKspMonitor = WriteDebugInfoLevel >= 2,
                EnableKspConvergedReason = WriteDebugInfoLevel >= 1
            };
            solver.Init(level);
            return solver;
        }

        string ComposeOptions() {
            var generatedOptions = new List<(string Option, string Value, string Source)> {
                ("-ksp_rtol", ConvergenceCriterion.ToString("R", CultureInfo.InvariantCulture), nameof(ConvergenceCriterion)),
                ("-ksp_atol", Math.Max(ConvergenceCriterion, MachineEpsilon).ToString("R", CultureInfo.InvariantCulture), nameof(ConvergenceCriterion) + "/" + nameof(MachineEpsilon)),
                ("-ksp_min_it", MinSolverIterations.ToString(CultureInfo.InvariantCulture), nameof(MinSolverIterations)),
                ("-ksp_max_it", MaxSolverIterations.ToString(CultureInfo.InvariantCulture), nameof(MaxSolverIterations))
            };

            foreach(var generated in generatedOptions) {
                if(ContainsPetscOption(PETScOptions, generated.Option)) {
                    Console.Error.WriteLine(
                        $"PETSc configuration warning: {generated.Option} is set both by {generated.Source} "
                        + $"and PETScOptions; the PETScOptions value takes precedence.");
                }
            }

            string overlay = string.Join(" ", generatedOptions.Select(option =>
                option.Value == null ? option.Option : option.Option + " " + option.Value));
            return string.IsNullOrWhiteSpace(PETScOptions) ? overlay : overlay + " " + PETScOptions;
        }

        static bool ContainsPetscOption(string options, string option) {
            if(string.IsNullOrWhiteSpace(options))
                return false;

            foreach(string token in options.Split((char[])null, StringSplitOptions.RemoveEmptyEntries)) {
                if(token.Equals(option, StringComparison.Ordinal)
                    || token.StartsWith(option + "=", StringComparison.Ordinal))
                    return true;
                if(option.Equals("-ksp_monitor", StringComparison.Ordinal)
                    && token.StartsWith("-ksp_monitor_", StringComparison.Ordinal))
                    return true;
            }
            return false;
        }

        internal static string DescribeNativeVariant(PETScNativeVariant variant) {
            switch(variant) {
                case PETScNativeVariant.Auto:
                    return "native Auto";
                case PETScNativeVariant.OpenMp:
                    return "OMP";
                case PETScNativeVariant.Mpi:
                    return "MPI";
                case PETScNativeVariant.MpiOpenMp:
                    return "MPI+OMP";
                default:
                    return "native " + variant;
            }
        }
    }

    /// <summary>
    /// Level-3 adapter from a supplied operator level to a PETSc Krylov solver.
    /// </summary>
    public sealed class PETScKrylovSolverAdapter : ISolverSmootherTemplate {

        PETScSolver m_Solver;
        MultigridOperator m_Operator;
        int m_Iterations;
        bool m_Converged;

        /// <summary>
        /// PETSc command-line-style solver options.
        /// </summary>
        public string Options { get; set; }

        /// <summary>
        /// Native PETSc parallel runtime.
        /// </summary>
        public PETScNativeVariant NativeVariant { get; set; } = PETScNativeVariant.Auto;

        /// <summary>
        /// PETSc matrix representation used for the supplied operator matrix.
        /// </summary>
        public PETScMatrixInput MatrixInput { get; set; } = PETScMatrixInput.BorrowedCsr;

        /// <summary>
        /// Optional BoSSS variable indices for PETSc field split 0.
        /// </summary>
        public int[] FirstFieldVariables { get; set; }

        /// <summary>
        /// Optional BoSSS variable indices for PETSc field split 1.
        /// </summary>
        public int[] SecondFieldVariables { get; set; }

        /// <summary>Whether to install BoSSS's geometric hierarchy as PETSc PCMG.</summary>
        public bool UseGeometricMultigrid { get; set; }

        /// <summary>Global-DOF threshold at which the supplied BoSSS hierarchy stops.</summary>
        public int CoarseKickIn { get; set; } = 90000;

        /// <summary>Writes managed outer-KSP residual lines to job standard output.</summary>
        public bool EnableKspMonitor { get; set; }

        /// <summary>Writes the final PETSc KSP convergence reason to job standard output.</summary>
        public bool EnableKspConvergedReason { get; set; }

        /// <summary>
        /// PETSc convergence reason from the most recent solve.
        /// </summary>
        public int ConvergenceReason => m_Solver?.ConvergenceReason ?? 0;

        /// <summary>
        /// PETSc residual norm from the most recent solve.
        /// </summary>
        public double ResidualNorm => m_Solver?.ResidualNorm ?? double.NaN;

        /// <summary>
        /// Number of PETSc KSP iterations from the most recent solve, or -1 if
        /// the native solver failed before reporting a KSP status.
        /// </summary>
        public int LastSolveIterationCount => m_Solver?.IterationCount ?? -1;

        /// <summary>
        /// True if the most recent PETSc solve did not reach a positive PETSc
        /// convergence reason.
        /// </summary>
        public bool DidNotConverge => m_Solver?.DidNotConverge ?? false;

        /// <summary>
        /// Native PETSc error description from the most recent failed solve.
        /// </summary>
        public string LastErrorMessage => m_Solver?.LastErrorMessage;

        public int IterationsInNested => 0;

        public int ThisLevelIterations => m_Iterations;

        public bool Converged => m_Converged;

        public void Init(MultigridOperator op) {
            if(op == null)
                throw new ArgumentNullException(nameof(op));
            if(ReferenceEquals(op, m_Operator) && m_Solver != null)
                return;

            DisposeSolver();
            ValidateFieldVariables(op);

            PETScSolver solver;
            if(UseGeometricMultigrid) {
                if((FirstFieldVariables?.Length ?? 0) > 0 || (SecondFieldVariables?.Length ?? 0) > 0)
                    throw new NotSupportedException("PETSc PCMG cannot simultaneously use the Level-3 PCFIELDSPLIT variable configuration.");
                solver = new PETScGeometricMultigridSolver {
                    Options = ForcePcMultigrid(Options),
                    Hierarchy = BuildGeometricHierarchy(op),
                    NativeVariant = NativeVariant,
                    MatrixInput = MatrixInput,
                    EnableKspMonitor = EnableKspMonitor,
                    EnableKspConvergedReason = EnableKspConvergedReason
                };
            } else {
                solver = new PETScKrylovSolver {
                    Options = Options ?? string.Empty,
                    NativeVariant = NativeVariant,
                    MatrixInput = MatrixInput,
                    EnableKspMonitor = EnableKspMonitor,
                    EnableKspConvergedReason = EnableKspConvergedReason
                };
            }

            try {
                if(solver is PETScKrylovSolver krylovSolver && FirstFieldVariables != null && FirstFieldVariables.Length > 0) {
                    krylovSolver.FirstFieldRows = ConvertRows(op.Mapping.GetSubvectorIndices(FirstFieldVariables));
                    krylovSolver.SecondFieldRows = ConvertRows(op.Mapping.GetSubvectorIndices(SecondFieldVariables));
                }
                solver.DefineMatrix(op.OperatorMatrix);
            } catch {
                solver.Dispose();
                throw;
            }

            m_Solver = solver;
            m_Operator = op;
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> {
            if(m_Solver == null)
                throw new InvalidOperationException("PETSc solver is not initialized.");

            var result = m_Solver.Solve(X, B);
            m_Iterations += result.NoOfIterations;
            m_Converged = result.Converged;
        }

        public void ResetStat() {
            m_Iterations = 0;
            m_Converged = false;
        }

        public long UsedMemory() {
            return 0;
        }

        public object Clone() {
            return new PETScKrylovSolverAdapter {
                Options = Options,
                NativeVariant = NativeVariant,
                MatrixInput = MatrixInput,
                FirstFieldVariables = FirstFieldVariables == null ? null : (int[])FirstFieldVariables.Clone(),
                SecondFieldVariables = SecondFieldVariables == null ? null : (int[])SecondFieldVariables.Clone(),
                UseGeometricMultigrid = UseGeometricMultigrid,
                CoarseKickIn = CoarseKickIn,
                EnableKspMonitor = EnableKspMonitor,
                EnableKspConvergedReason = EnableKspConvergedReason
            };
        }

        public void Dispose() {
            DisposeSolver();
            m_Operator = null;
        }

        void DisposeSolver() {
            m_Solver?.Dispose();
            m_Solver = null;
        }

        void ValidateFieldVariables(MultigridOperator op) {
            bool hasFirst = FirstFieldVariables != null && FirstFieldVariables.Length > 0;
            bool hasSecond = SecondFieldVariables != null && SecondFieldVariables.Length > 0;
            if(hasFirst != hasSecond)
                throw new ArgumentException("Both PETSc field-split variable arrays must be supplied together.");
            if(!hasFirst)
                return;

            int fieldCount = op.Mapping.AggBasis.Length;
            ValidateVariables(FirstFieldVariables, fieldCount, nameof(FirstFieldVariables));
            ValidateVariables(SecondFieldVariables, fieldCount, nameof(SecondFieldVariables));
        }

        static void ValidateVariables(int[] variables, int fieldCount, string parameterName) {
            foreach(int variable in variables)
                if(variable < 0 || variable >= fieldCount)
                    throw new ArgumentOutOfRangeException(parameterName, $"Variable index {variable} is outside [0,{fieldCount}).");
        }

        static int[] ConvertRows(long[] rows) {
            var result = new int[rows.Length];
            for(int i = 0; i < rows.Length; i++) {
                if(rows[i] < 0 || rows[i] > int.MaxValue)
                    throw new NotSupportedException($"PETSc field row {rows[i]} exceeds the packaged 32-bit index range.");
                result[i] = (int)rows[i];
            }
            return result;
        }

        PETScGeometricMultigridHierarchy BuildGeometricHierarchy(MultigridOperator finest) {
            if(CoarseKickIn < 1)
                throw new ArgumentOutOfRangeException(nameof(CoarseKickIn), "The PETSc multigrid coarse-level threshold must be positive.");

            var levels = new List<MultigridOperator> { finest };
            MultigridOperator coarsest = finest;
            while(coarsest.CoarserLevel != null && coarsest.Mapping.TotalLength > CoarseKickIn) {
                coarsest = coarsest.CoarserLevel;
                levels.Add(coarsest);
            }
            _ = coarsest.OperatorMatrix;
            if(levels.Count < 2)
                throw new InvalidOperationException(
                    $"PETSc geometric multigrid requires at least two levels at or above CoarseKickIn={CoarseKickIn}.");

            var hierarchy = new PETScGeometricMultigridHierarchy();
            for(int index = 0; index < levels.Count - 1; index++) {
                MultigridOperator coarse = levels[index + 1];
                BlockMsrMatrix interpolation = coarse.RightChangeOfBasis == null
                    ? coarse.GetPrologonationOperator
                    : BlockMsrMatrix.Multiply(coarse.GetPrologonationOperator, coarse.RightChangeOfBasis);
                BlockMsrMatrix restriction = coarse.LeftChangeOfBasis == null
                    ? coarse.GetRestrictionOperator
                    : BlockMsrMatrix.Multiply(coarse.LeftChangeOfBasis, coarse.GetRestrictionOperator);
                hierarchy.Transitions.Add(new PETScMultigridTransition {
                    CoarseOperator = coarse.OperatorMatrix,
                    Interpolation = interpolation,
                    Restriction = restriction
                });
            }
            return hierarchy;
        }

        static string ForcePcMultigrid(string options) {
            var tokens = (options ?? string.Empty).Split((char[])null, StringSplitOptions.RemoveEmptyEntries);
            var retained = new List<string>();
            for(int index = 0; index < tokens.Length; index++) {
                if(tokens[index].Equals("-pc_type", StringComparison.Ordinal)) {
                    if(index + 1 < tokens.Length)
                        index++;
                    continue;
                }
                if(tokens[index].StartsWith("-pc_type=", StringComparison.Ordinal))
                    continue;
                retained.Add(tokens[index]);
            }
            retained.Add("-pc_type");
            retained.Add("mg");
            return string.Join(" ", retained);
        }
    }
}
