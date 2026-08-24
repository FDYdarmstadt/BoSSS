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
using System.Runtime.Serialization;
using ilPSP.LinSolvers.PETSc;

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
        public PETScNativeVariant NativeVariant = PETScNativeVariant.Mpi;

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
        public string Name => "PETSc direct " + Factorization + " (" + Backend + ")";

        /// <summary>
        /// Compact solver name.
        /// </summary>
        public string Shortname => "PETSc-" + Factorization + "-" + Backend;

        /// <summary>
        /// Creates and initializes a direct PETSc solver for a multigrid level.
        /// </summary>
        public ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            var solver = new PETScKSP {
                Options = PETScDirectSolver.BuildOptions(Factorization, Backend, AdditionalOptions),
                NativeVariant = NativeVariant
            };
            solver.Init(level);
            return solver;
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
                hash = hash * 397 ^ (int)Factorization;
                hash = hash * 397 ^ (int)Backend;
                hash = hash * 397 ^ (AdditionalOptions?.GetHashCode() ?? 0);
                return hash;
            }
        }
    }



    /// <summary>
    /// Configuration for the options-driven PETSc KSP adapter.
    /// </summary>
    [Serializable]
    public sealed class PETScKSPConfig : IterativeSolverConfig {

        /// <summary>
        /// PETSc command-line-style options. The KSP and preconditioner are not fixed by
        /// this adapter; for example, use <c>-ksp_type cg -pc_type jacobi</c>.
        /// Options supplied here take precedence over the tolerance and iteration defaults
        /// inherited from <see cref="IterativeSolverConfig"/>.
        /// </summary>
        [DataMember]
        public string PETScOptions = "-ksp_type gmres -pc_type jacobi";

        /// <summary>
        /// Native PETSc parallel runtime. MPI remains the default.
        /// </summary>
        [DataMember]
        public PETScNativeVariant NativeVariant = PETScNativeVariant.Mpi;

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

        public override string Name => "PETSc KSP (options-driven)";

        public override string Shortname => "PETSc-KSP";

        public override ISolverSmootherTemplate CreateInstance(MultigridOperator level) {
            var solver = new PETScKSP {
                Options = ComposeOptions(),
                NativeVariant = NativeVariant,
                FirstFieldVariables = FirstFieldVariables == null ? null : (int[])FirstFieldVariables.Clone(),
                SecondFieldVariables = SecondFieldVariables == null ? null : (int[])SecondFieldVariables.Clone()
            };
            solver.Init(level);
            return solver;
        }

        string ComposeOptions() {
            string defaults = string.Format(
                CultureInfo.InvariantCulture,
                "-ksp_rtol {0:R} -ksp_atol {0:R} -ksp_max_it {1}",
                ConvergenceCriterion,
                MaxSolverIterations);
            return string.IsNullOrWhiteSpace(PETScOptions) ? defaults : defaults + " " + PETScOptions;
        }
    }

    /// <summary>
    /// Level-3 adapter from a <see cref="MultigridOperator"/> to the ilPSP PETSc KSP solver.
    /// </summary>
    public sealed class PETScKSP : ISolverSmootherTemplate {

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
        public PETScNativeVariant NativeVariant { get; set; } = PETScNativeVariant.Mpi;

        /// <summary>
        /// Optional BoSSS variable indices for PETSc field split 0.
        /// </summary>
        public int[] FirstFieldVariables { get; set; }

        /// <summary>
        /// Optional BoSSS variable indices for PETSc field split 1.
        /// </summary>
        public int[] SecondFieldVariables { get; set; }

        /// <summary>
        /// PETSc convergence reason from the most recent solve.
        /// </summary>
        public int ConvergenceReason => m_Solver?.ConvergenceReason ?? 0;

        /// <summary>
        /// PETSc residual norm from the most recent solve.
        /// </summary>
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
            ValidateFieldVariables(op);

            var solver = new PETScSolver {
                Options = Options ?? string.Empty,
                NativeVariant = NativeVariant
            };

            try {
                if(FirstFieldVariables != null && FirstFieldVariables.Length > 0) {
                    solver.FirstFieldRows = ConvertRows(op.Mapping.GetSubvectorIndices(FirstFieldVariables));
                    solver.SecondFieldRows = ConvertRows(op.Mapping.GetSubvectorIndices(SecondFieldVariables));
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
                throw new InvalidOperationException("PETSc KSP is not initialized.");

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
            return new PETScKSP {
                Options = Options,
                NativeVariant = NativeVariant,
                FirstFieldVariables = FirstFieldVariables == null ? null : (int[])FirstFieldVariables.Clone(),
                SecondFieldVariables = SecondFieldVariables == null ? null : (int[])SecondFieldVariables.Clone()
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
    }
}
