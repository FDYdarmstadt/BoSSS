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

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Selects how a BoSSS operator is exposed to PETSc.
    /// </summary>
    public enum PETScMatrixInput {

        /// <summary>
        /// Exports BoSSS CSR data and lets PETSc copy it into an MPIAIJ matrix.
        /// This compatibility path supports direct and iterative solvers.
        /// </summary>
        CopiedCsr,

        /// <summary>
        /// Packs the BoSSS matrix once into PETSc's diagonal/off-diagonal MPIAIJ CSR
        /// layout and gives PETSc pinned references to the six resulting arrays.
        /// PETSc does not copy these arrays. They remain immutable and pinned until the
        /// solver is disposed. A structural zero diagonal is inserted where necessary,
        /// so this mode also supports PETSc factor preconditioners and direct solvers.
        /// The packing is necessary because a <c>BlockMsrMatrix</c> does not expose one
        /// contiguous scalar CSR representation.
        /// </summary>
        BorrowedCsr,

    }
}
