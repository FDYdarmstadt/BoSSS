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

    /// <summary>
    /// Options-driven Krylov solver using the direct PETSc ABI L0 implementation.
    /// </summary>
    /// <remarks>
    /// KSP and PC behavior is supplied solely through <see cref="Options"/>.
    /// Optional field-index arrays configure PETSc PCFIELDSPLIT names <c>0</c>
    /// and <c>1</c> after the selected PC type is read from the options database.
    /// </remarks>
    public sealed unsafe class PETScKrylovSolver : PETScSolver {
        IntPtr m_FirstFieldIndexSet;
        IntPtr m_SecondFieldIndexSet;

        /// <summary>PETSc command-line-style KSP and PC options.</summary>
        public string Options { get; set; } = "-ksp_type gmres -pc_type jacobi -ksp_rtol 1e-10 -ksp_atol 1e-50 -ksp_max_it 10000";

        /// <summary>Optional global row indices for PCFIELDSPLIT field 0.</summary>
        public int[] FirstFieldRows { get; set; }

        /// <summary>Optional global row indices for PCFIELDSPLIT field 1.</summary>
        public int[] SecondFieldRows { get; set; }

        /// <inheritdoc/>
        protected override void ConfigureKsp() {
            bool hasFirst = FirstFieldRows != null && FirstFieldRows.Length > 0;
            bool hasSecond = SecondFieldRows != null && SecondFieldRows.Length > 0;
            if(hasFirst != hasSecond)
                throw new ArgumentException("Both PETSc field-split row arrays must be supplied together.");

            ApplyKspOptions(Options);
            Check(m_Native.GetKspPreconditioner(m_Ksp, out IntPtr pc), "PETSc KSPGetPC");
            if(hasFirst) {
                Check(m_Native.ObjectTypeCompare(pc, "fieldsplit", out int isFieldSplit), "PETSc PetscObjectTypeCompare");
                if(isFieldSplit != 0) {
                    CreateFieldSplitIndexSet(FirstFieldRows, "0", pc, out m_FirstFieldIndexSet);
                    CreateFieldSplitIndexSet(SecondFieldRows, "1", pc, out m_SecondFieldIndexSet);
                }
            }
        }

        /// <inheritdoc/>
        protected override PETScMatrixInput SelectMatrixInput(int mpiSize) {
            if(mpiSize == 1 && MatrixInput == PETScMatrixInput.BorrowedCsr && RequestsSequentialFactorMatrix(Options))
                return PETScMatrixInput.CopiedCsr;
            return base.SelectMatrixInput(mpiSize);
        }

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

        void CreateFieldSplitIndexSet(int[] rows, string name, IntPtr pc, out IntPtr indexSet) {
            fixed(int* indices = rows) {
                // PETSC_COPY_VALUES: PETSc copies the managed index array before this call returns.
                Check(m_Native.CreateGeneralIndexSet(m_Communicator, rows.Length, (IntPtr)indices, 0, out indexSet), "PETSc ISCreateGeneral");
            }
            Check(m_Native.SetPcFieldSplitIndexSet(pc, name, indexSet), "PETSc PCFieldSplitSetIS");
        }

        /// <inheritdoc/>
        protected override void DestroyNativeObjects() {
            try {
                base.DestroyNativeObjects();
            } finally {
                if(m_Native != null && m_FirstFieldIndexSet != IntPtr.Zero)
                    Check(m_Native.DestroyIndexSet(ref m_FirstFieldIndexSet), "PETSc ISDestroy(field 0)");
                if(m_Native != null && m_SecondFieldIndexSet != IntPtr.Zero)
                    Check(m_Native.DestroyIndexSet(ref m_SecondFieldIndexSet), "PETSc ISDestroy(field 1)");
            }
        }
    }
}
