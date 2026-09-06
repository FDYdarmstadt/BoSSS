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
using MPI.Wrappers;

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

        /// <summary>
        /// Uses a PETSc MatShell which reads an existing BoSSS BlockMsrMatrix on
        /// every multiplication. This is an iterative-solver compatibility path;
        /// PETSc performs the halo exchange from the matrix sparsity pattern.
        /// </summary>
        BlockMsrMatrixShell,

        /// <summary>
        /// Uses a PETSc MatShell backed by an <see cref="IPETScMatrixFreeOperator"/>.
        /// PETSc exchanges only the declared ghost entries and no sparse matrix is
        /// constructed or stored by either side.
        /// </summary>
        MatrixFreeShell
    }

    /// <summary>
    /// A short-lived view of a PETSc-owned vector segment passed to a matrix-free
    /// callback. The view is valid only for the duration of <see cref="IPETScMatrixFreeOperator.Apply"/>.
    /// </summary>
    public sealed class PETScVectorView {
        IntPtr m_Data;
        int m_Length;

        internal void Reset(IntPtr data, int length) {
            m_Data = data;
            m_Length = length;
        }

        /// <summary>Number of entries in this local vector view.</summary>
        public int Count => m_Length;

        /// <summary>Gets or sets a vector entry.</summary>
        public unsafe double this[int index] {
            get {
                CheckIndex(index);
                return ((double*)m_Data)[index];
            }
            set {
                CheckIndex(index);
                ((double*)m_Data)[index] = value;
            }
        }

        void CheckIndex(int index) {
            if(index < 0 || index >= m_Length)
                throw new ArgumentOutOfRangeException(nameof(index));
        }
    }

    /// <summary>
    /// Contract for a distributed PETSc matrix-free operator.
    /// </summary>
    public interface IPETScMatrixFreeOperator {

        /// <summary>Communicator of the distributed operator; currently MPI_COMM_WORLD only.</summary>
        MPI_Comm MPIComm { get; }

        /// <summary>Global square-operator dimension.</summary>
        long GlobalLength { get; }

        /// <summary>First globally owned row and vector index on this rank.</summary>
        long GlobalRowStart { get; }

        /// <summary>Number of locally owned rows and vector entries.</summary>
        int LocalLength { get; }

        /// <summary>
        /// Global, zero-based vector indices required remotely by <see cref="Apply"/>.
        /// The returned array may be null or empty when no remote values are needed.
        /// </summary>
        int[] GhostGlobalIds { get; }

        /// <summary>
        /// Applies the local part of the operator. The implementation must write all
        /// output entries and must not retain any vector view after this call returns.
        /// GhostGlobalIds determines the order of haloInput.
        /// </summary>
        void Apply(PETScVectorView ownedInput, PETScVectorView haloInput, PETScVectorView output);
    }
}
