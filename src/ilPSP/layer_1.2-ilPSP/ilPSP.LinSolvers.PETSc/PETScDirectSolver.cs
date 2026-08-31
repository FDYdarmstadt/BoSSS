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
using System.Text;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Factorization performed by a <see cref="PETScDirectSolver"/>.
    /// </summary>
    public enum PETScDirectFactorization {

        /// <summary>
        /// General sparse LU factorization.
        /// </summary>
        LU,

        /// <summary>
        /// Sparse Cholesky factorization for symmetric positive-definite matrices.
        /// </summary>
        Cholesky
    }

    /// <summary>
    /// PETSc matrix-factorization package requested by a
    /// <see cref="PETScDirectSolver"/>.
    /// </summary>
    /// <remarks>
    /// A backend must be present in the selected native PETSc build. The portable
    /// baseline currently contains only PETSc's built-in sequential factorization.
    /// </remarks>
    public enum PETScDirectBackend {

        /// <summary>
        /// Let PETSc select its default factorization package.
        /// </summary>
        Default,

        /// <summary>
        /// MUMPS sparse direct solver. A MUMPS-enabled PETSc build is required.
        /// </summary>
        Mumps,

        /// <summary>
        /// SuperLU_DIST sparse direct solver. A SuperLU_DIST-enabled PETSc build is required.
        /// </summary>
        SuperLuDist,

        /// <summary>
        /// SuiteSparse UMFPACK. This backend is intended for single-rank LU solves.
        /// </summary>
        Umfpack
    }

    /// <summary>
    /// Typed PETSc sparse direct solver.
    /// </summary>
    /// <remarks>
    /// This configures PETSc as <c>KSPPREONLY</c> with an LU or Cholesky
    /// factorization preconditioner. Distributed direct solves require a backend
    /// such as MUMPS or SuperLU_DIST to be compiled into PETSc.
    /// </remarks>
    public sealed class PETScDirectSolver : PETScSolverBase {

        /// <summary>
        /// Factorization to perform.
        /// </summary>
        public PETScDirectFactorization Factorization { get; set; } = PETScDirectFactorization.LU;

        /// <summary>
        /// Requested PETSc factorization backend.
        /// </summary>
        public PETScDirectBackend Backend { get; set; } = PETScDirectBackend.Default;

        /// <summary>
        /// Additional PETSc options appended after the typed direct-solver options,
        /// for example package-specific pivoting or ordering controls.
        /// </summary>
        public string AdditionalOptions { get; set; } = string.Empty;

        /// <summary>
        /// Produces the PETSc option string corresponding to a typed direct-solver
        /// configuration.
        /// </summary>
        public static string BuildOptions(
            PETScDirectFactorization factorization,
            PETScDirectBackend backend,
            string additionalOptions = null) {

            string pcType;
            switch(factorization) {
                case PETScDirectFactorization.LU:
                    pcType = "lu";
                    break;
                case PETScDirectFactorization.Cholesky:
                    pcType = "cholesky";
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(factorization), factorization, "Unknown PETSc direct factorization.");
            }

            string backendName;
            switch(backend) {
                case PETScDirectBackend.Default:
                    backendName = null;
                    break;
                case PETScDirectBackend.Mumps:
                    backendName = "mumps";
                    break;
                case PETScDirectBackend.SuperLuDist:
                    backendName = "superlu_dist";
                    break;
                case PETScDirectBackend.Umfpack:
                    backendName = "umfpack";
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(backend), backend, "Unknown PETSc direct backend.");
            }

            var options = new StringBuilder("-ksp_type preonly -pc_type ");
            options.Append(pcType);
            if(backendName != null) {
                options.Append(" -pc_factor_mat_solver_type ");
                options.Append(backendName);
            }
            if(!string.IsNullOrWhiteSpace(additionalOptions)) {
                options.Append(' ');
                options.Append(additionalOptions.Trim());
            }
            return options.ToString();
        }

        /// <inheritdoc/>
        protected override string SolverOptions => BuildOptions(Factorization, Backend, AdditionalOptions);

        /// <inheritdoc/>
        protected override void ValidateParallelConfiguration(int mpiSize, PETScNativeVariant effectiveNativeVariant) {
            if(Factorization == PETScDirectFactorization.Cholesky
                && (Backend == PETScDirectBackend.Umfpack || Backend == PETScDirectBackend.SuperLuDist)) {
                throw new NotSupportedException(
                    $"PETSc direct backend '{Backend}' supports LU factorization only; Cholesky is not available.");
            }
            if(Backend == PETScDirectBackend.SuperLuDist
                && effectiveNativeVariant == PETScNativeVariant.OpenMp) {
                throw new NotSupportedException(
                    "PETSc direct backend 'SuperLuDist' requires an MPI-capable PETSc native variant. "
                    + "Use the Mpi or MpiOpenMp variant.");
            }
            if(mpiSize > 1 && (Backend == PETScDirectBackend.Default || Backend == PETScDirectBackend.Umfpack)) {
                throw new NotSupportedException(
                    $"PETSc direct backend '{Backend}' is sequential and cannot run with {mpiSize} MPI ranks. "
                    + "Use the Mumps or SuperLuDist backend for distributed direct solves.");
            }
        }

        /// <inheritdoc/>
        protected override Exception CreateSetupException(Exception error) {
            if(Backend == PETScDirectBackend.Default)
                return error;
            return new InvalidOperationException(
                "PETSc could not initialize the requested direct-solver backend '"
                + Backend
                + "'. Ensure that this backend was compiled into the selected native PETSc variant.",
                error);
        }
    }
}
