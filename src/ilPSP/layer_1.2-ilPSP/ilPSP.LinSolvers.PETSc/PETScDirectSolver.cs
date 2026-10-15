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

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>Factorization performed by a <see cref="PETScDirectSolver"/>.</summary>
    public enum PETScDirectFactorization {
        /// <summary>General sparse LU factorization.</summary>
        LU,
        /// <summary>Sparse Cholesky factorization for symmetric positive-definite matrices.</summary>
        Cholesky
    }

    /// <summary>PETSc matrix-factorization package requested by a direct solver.</summary>
    public enum PETScDirectBackend {
        /// <summary>Let PETSc select its default factorization package.</summary>
        Default,
        /// <summary>MUMPS sparse direct solver.</summary>
        Mumps,
        /// <summary>SuperLU_DIST sparse direct solver.</summary>
        SuperLuDist,
        /// <summary>SuiteSparse UMFPACK.</summary>
        Umfpack
    }

    /// <summary>Stable managed surface for a direct PETSc ABI solve.</summary>
    public interface IPETScDirectSolverAbi : ISparseSolver {
        /// <summary>Factorization to perform.</summary>
        PETScDirectFactorization Factorization { get; set; }

        /// <summary>Requested PETSc factorization backend.</summary>
        PETScDirectBackend Backend { get; set; }

        /// <summary>Package-specific PETSc options.</summary>
        string AdditionalOptions { get; set; }

        /// <summary>Requested native PETSc runtime.</summary>
        PETScNativeVariant NativeVariant { get; set; }

        /// <summary>PETSc matrix ownership and transfer representation.</summary>
        PETScMatrixInput MatrixInput { get; set; }

        /// <summary>PETSc convergence reason from the latest solve.</summary>
        int ConvergenceReason { get; }

        /// <summary>PETSc residual norm from the latest solve.</summary>
        double ResidualNorm { get; }
    }

    /// <summary>Creates the direct PETSc ABI implementation for upper BoSSS layers.</summary>
    public static class PETScDirectSolverAbi {
        /// <summary>Creates a direct PETSc ABI solver.</summary>
        public static IPETScDirectSolverAbi Create() => new PETScDirectSolver();
    }

    /// <summary>
    /// Direct-solver policy built on the shared direct PETSc ABI L0 solver.
    /// </summary>
    /// <remarks>
    /// <see cref="PETScSolver"/> owns the native PETSc lifecycle, matrix and
    /// vector transfer, KSP creation, solve operation, and disposal. This class
    /// only chooses <c>KSPPREONLY</c>, the factorization PC, and its backend.
    /// </remarks>
    public sealed class PETScDirectSolver : PETScSolver, IPETScDirectSolverAbi {

        /// <summary>Factorization to perform.</summary>
        public PETScDirectFactorization Factorization { get; set; } = PETScDirectFactorization.LU;

        /// <summary>Requested PETSc factorization backend.</summary>
        public PETScDirectBackend Backend { get; set; } = PETScDirectBackend.Default;

        /// <summary>Additional PETSc package-specific options for this KSP.</summary>
        public string AdditionalOptions { get; set; } = string.Empty;

        /// <inheritdoc/>
        protected override void ConfigureKsp() {
            Check(m_Native.SetKspType(m_Ksp, "preonly"), "PETSc KSPSetType");
            Check(m_Native.GetKspPreconditioner(m_Ksp, out IntPtr pc), "PETSc KSPGetPC");
            string pcType = Factorization == PETScDirectFactorization.LU ? "lu"
                : Factorization == PETScDirectFactorization.Cholesky ? "cholesky"
                : throw new ArgumentOutOfRangeException(nameof(Factorization));
            Check(m_Native.SetPcType(pc, pcType), "PETSc PCSetType");

            string backend = GetBackendName(Backend);
            if(backend != null)
                Check(m_Native.SetPcFactorMatSolverType(pc, backend), "PETSc PCFactorSetMatSolverType");

            ApplyKspOptions(AdditionalOptions);
        }

        /// <inheritdoc/>
        protected override void ValidateParallelConfiguration(int mpiSize, PETScNativeVariant variant) {
            if(Factorization == PETScDirectFactorization.Cholesky
                && (Backend == PETScDirectBackend.Umfpack || Backend == PETScDirectBackend.SuperLuDist))
                throw new NotSupportedException($"PETSc direct backend '{Backend}' supports LU factorization only; Cholesky is not available.");
            if(Backend == PETScDirectBackend.SuperLuDist && variant == PETScNativeVariant.OpenMp)
                throw new NotSupportedException("PETSc direct backend 'SuperLuDist' requires an MPI-capable PETSc native variant.");
            if(mpiSize > 1 && (Backend == PETScDirectBackend.Default || Backend == PETScDirectBackend.Umfpack))
                throw new NotSupportedException($"PETSc direct backend '{Backend}' is sequential and cannot run with {mpiSize} MPI ranks.");
        }

        static string GetBackendName(PETScDirectBackend backend) {
            switch(backend) {
                case PETScDirectBackend.Default: return null;
                case PETScDirectBackend.Mumps: return "mumps";
                case PETScDirectBackend.SuperLuDist: return "superlu_dist";
                case PETScDirectBackend.Umfpack: return "umfpack";
                default: throw new ArgumentOutOfRangeException(nameof(backend), backend, "Unknown PETSc direct backend.");
            }
        }
    }
}
