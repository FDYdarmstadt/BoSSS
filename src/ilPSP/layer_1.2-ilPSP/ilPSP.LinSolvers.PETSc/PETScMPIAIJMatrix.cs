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
    /// Base class and factory for BoSSS-owned PETSc AIJ matrices.
    /// </summary>
    /// <remarks>
    /// The subclasses differ in CSR-storage ownership:
    /// <see cref="PETScCopiedMPIAIJMatrix"/> lets PETSc copy the data, while
    /// <see cref="PETScBorrowedMPIAIJMatrix"/> pins and lends it to PETSc.
    /// </remarks>
    internal abstract class PETScMPIAIJMatrix : IDisposable {
        /// <summary>Validated PETSc-compatible dimensions of one BoSSS matrix.</summary>
        internal sealed class Dimensions {
            internal readonly int LocalRows;
            internal readonly int LocalColumns;
            internal readonly int GlobalRows;
            internal readonly int GlobalColumns;
            internal readonly int GlobalRowStart;
            internal readonly int MpiSize;

            internal Dimensions(IMutableMatrixEx matrix) {
                if(matrix == null)
                    throw new ArgumentNullException(nameof(matrix));
                LocalRows = matrix.RowPartitioning.LocalLength;
                LocalColumns = matrix.ColPartitioning.LocalLength;
                GlobalRows = CheckedIndex(matrix.RowPartitioning.TotalLength, "row count");
                GlobalColumns = CheckedIndex(matrix.ColPartitioning.TotalLength, "column count");
                GlobalRowStart = CheckedIndex(matrix.RowPartitioning.i0, "first row index");
                csMPI.Raw.Comm_Size(matrix.MPI_Comm, out MpiSize);
                if(MpiSize < 1)
                    throw new InvalidOperationException("The BoSSS matrix has an invalid MPI communicator.");
            }

            static int CheckedIndex(long value, string description) {
                if(value < 0 || value > int.MaxValue)
                    throw new NotSupportedException($"PETSc uses 32-bit indices and the matrix {description} {value} is outside that range.");
                return checked((int)value);
            }
        }

        /// <summary>
        /// Creates the requested PETSc matrix transfer implementation.
        /// </summary>
        /// <remarks>
        /// Keeping this selection here makes a solver independent of matrix
        /// transfer details; it owns exactly one <see cref="PETScMPIAIJMatrix"/>.
        /// </remarks>
        internal static PETScMPIAIJMatrix Create(UnsafePETSc native, IntPtr communicator,
            PETScMatrixInput matrixInput, IMutableMatrixEx matrix) {
            var dimensions = new Dimensions(matrix);
            switch(matrixInput) {
                case PETScMatrixInput.CopiedCsr:
                    return PETScCopiedMPIAIJMatrix.Create(native, communicator, matrix, dimensions);
                case PETScMatrixInput.BorrowedCsr:
                    return PETScBorrowedMPIAIJMatrix.Create(native, communicator, matrix, dimensions);
                default:
                    throw new ArgumentOutOfRangeException(nameof(matrixInput), matrixInput,
                        "Unknown PETSc matrix input mode.");
            }
        }

        /// <summary>Native PETSc matrix handle.</summary>
        public abstract IntPtr Handle { get; }

        /// <summary>Number of BoSSS rows owned by this MPI rank.</summary>
        public abstract int LocalRows { get; }

        /// <summary>Destroys the PETSc matrix and releases associated storage.</summary>
        public abstract void Dispose();
    }
}
