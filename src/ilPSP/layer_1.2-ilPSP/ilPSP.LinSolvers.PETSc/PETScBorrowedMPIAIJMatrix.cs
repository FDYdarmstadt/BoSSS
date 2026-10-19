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
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Owns a PETSc <c>MATMPIAIJ</c> matrix created from borrowed split CSR arrays.
    /// </summary>
    /// <remarks>
    /// <see cref="UnsafePETSc.CreateMpiAijWithSplitArrays"/> borrows all
    /// six arrays. This owner pins them until the native matrix is destroyed,
    /// preserving the zero-copy transfer used by the established BoSSS wrapper.
    /// </remarks>
    internal sealed class PETScBorrowedMPIAIJMatrix : PETScMPIAIJMatrix {

        sealed class SplitCsr {
            internal int[] DiagonalRowPointers;
            internal int[] DiagonalColumnIndices;
            internal double[] DiagonalValues;
            internal int[] OffdiagonalRowPointers;
            internal int[] OffdiagonalColumnIndices;
            internal double[] OffdiagonalValues;
        }

        readonly UnsafePETSc m_Native;
        readonly int m_LocalRows;
        IntPtr m_Handle;
        int[] m_DiagonalRowPointers;
        int[] m_DiagonalColumnIndices;
        double[] m_DiagonalValues;
        int[] m_OffdiagonalRowPointers;
        int[] m_OffdiagonalColumnIndices;
        double[] m_OffdiagonalValues;
        GCHandle m_DiagonalRowPointersHandle;
        GCHandle m_DiagonalColumnIndicesHandle;
        GCHandle m_DiagonalValuesHandle;
        GCHandle m_OffdiagonalRowPointersHandle;
        GCHandle m_OffdiagonalColumnIndicesHandle;
        GCHandle m_OffdiagonalValuesHandle;

        /// <summary>Native PETSc matrix handle.</summary>
        public override IntPtr Handle => m_Handle;

        /// <summary>Number of matrix rows owned by this rank.</summary>
        public override int LocalRows => m_LocalRows;

        PETScBorrowedMPIAIJMatrix(UnsafePETSc native, Dimensions dimensions) {
            m_Native = native ?? throw new ArgumentNullException(nameof(native));
            m_LocalRows = dimensions.LocalRows;
        }

        /// <summary>
        /// Exports and creates one PETSc MPIAIJ matrix. The source matrix's
        /// row and column partitionings must be identical.
        /// </summary>
        internal static PETScBorrowedMPIAIJMatrix Create(UnsafePETSc native, IntPtr communicator,
            IMutableMatrixEx matrix, Dimensions dimensions) {
            var result = new PETScBorrowedMPIAIJMatrix(native, dimensions);
            try {
                SplitCsr csr = ExportSplitCsr(matrix, dimensions.LocalRows, dimensions.GlobalColumns, dimensions.GlobalRowStart);
                result.Pin(csr);
                Check(native.CreateMpiAijWithSplitArrays(
                    communicator, dimensions.LocalRows, dimensions.LocalColumns,
                    dimensions.GlobalRows, dimensions.GlobalColumns,
                    result.m_DiagonalRowPointersHandle.AddrOfPinnedObject(),
                    result.m_DiagonalColumnIndicesHandle.AddrOfPinnedObject(),
                    result.m_DiagonalValuesHandle.AddrOfPinnedObject(),
                    result.m_OffdiagonalRowPointersHandle.AddrOfPinnedObject(),
                    result.m_OffdiagonalColumnIndicesHandle.AddrOfPinnedObject(),
                    result.m_OffdiagonalValuesHandle.AddrOfPinnedObject(), out result.m_Handle),
                    native, "PETSc MatCreateMPIAIJWithSplitArrays");
                return result;
            } catch {
                result.Dispose();
                throw;
            }
        }

        static SplitCsr ExportSplitCsr(IMutableMatrixEx matrix, int localRows, int globalRows, int globalRowStart) {
            if(matrix == null)
                throw new ArgumentNullException(nameof(matrix));
            if(matrix.RowPartitioning.TotalLength != matrix.ColPartitioning.TotalLength
                || matrix.RowPartitioning.i0 != matrix.ColPartitioning.i0
                || matrix.RowPartitioning.LocalLength != matrix.ColPartitioning.LocalLength)
                throw new NotSupportedException("PETSc split CSR requires identical row and column partitionings.");

            long firstOwnedColumn = matrix.ColPartitioning.i0;
            long lastOwnedColumn = firstOwnedColumn + matrix.ColPartitioning.LocalLength;
            var result = new SplitCsr {
                DiagonalRowPointers = new int[localRows + 1],
                OffdiagonalRowPointers = new int[localRows + 1]
            };
            long[] rowColumns = null;
            double[] rowValues = null;

            for(int localRow = 0; localRow < localRows; localRow++) {
                long globalRow = globalRowStart + localRow;
                int rowLength = matrix.GetRow(globalRow, ref rowColumns, ref rowValues);
                ValidateRow(rowLength, rowColumns, rowValues);
                int diagonalCount = 0;
                int offdiagonalCount = 0;
                bool hasStructuralDiagonal = false;
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    ValidateColumn(column, globalRows);
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

            int diagonalNonzeros = result.DiagonalRowPointers[localRows];
            int offdiagonalNonzeros = result.OffdiagonalRowPointers[localRows];
            result.DiagonalColumnIndices = new int[Math.Max(diagonalNonzeros, 1)];
            result.DiagonalValues = new double[Math.Max(diagonalNonzeros, 1)];
            result.OffdiagonalColumnIndices = new int[Math.Max(offdiagonalNonzeros, 1)];
            result.OffdiagonalValues = new double[Math.Max(offdiagonalNonzeros, 1)];
            int[] diagonalOffsets = (int[])result.DiagonalRowPointers.Clone();
            int[] offdiagonalOffsets = (int[])result.OffdiagonalRowPointers.Clone();

            for(int localRow = 0; localRow < localRows; localRow++) {
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

        void Pin(SplitCsr csr) {
            m_DiagonalRowPointers = csr.DiagonalRowPointers;
            m_DiagonalColumnIndices = csr.DiagonalColumnIndices;
            m_DiagonalValues = csr.DiagonalValues;
            m_OffdiagonalRowPointers = csr.OffdiagonalRowPointers;
            m_OffdiagonalColumnIndices = csr.OffdiagonalColumnIndices;
            m_OffdiagonalValues = csr.OffdiagonalValues;
            m_DiagonalRowPointersHandle = GCHandle.Alloc(m_DiagonalRowPointers, GCHandleType.Pinned);
            m_DiagonalColumnIndicesHandle = GCHandle.Alloc(m_DiagonalColumnIndices, GCHandleType.Pinned);
            m_DiagonalValuesHandle = GCHandle.Alloc(m_DiagonalValues, GCHandleType.Pinned);
            m_OffdiagonalRowPointersHandle = GCHandle.Alloc(m_OffdiagonalRowPointers, GCHandleType.Pinned);
            m_OffdiagonalColumnIndicesHandle = GCHandle.Alloc(m_OffdiagonalColumnIndices, GCHandleType.Pinned);
            m_OffdiagonalValuesHandle = GCHandle.Alloc(m_OffdiagonalValues, GCHandleType.Pinned);
        }

        public override void Dispose() {
            try {
                if(m_Handle != IntPtr.Zero)
                    Check(m_Native.DestroyMatrix(ref m_Handle), m_Native, "PETSc MatDestroy");
            } finally {
                if(m_DiagonalRowPointersHandle.IsAllocated) m_DiagonalRowPointersHandle.Free();
                if(m_DiagonalColumnIndicesHandle.IsAllocated) m_DiagonalColumnIndicesHandle.Free();
                if(m_DiagonalValuesHandle.IsAllocated) m_DiagonalValuesHandle.Free();
                if(m_OffdiagonalRowPointersHandle.IsAllocated) m_OffdiagonalRowPointersHandle.Free();
                if(m_OffdiagonalColumnIndicesHandle.IsAllocated) m_OffdiagonalColumnIndicesHandle.Free();
                if(m_OffdiagonalValuesHandle.IsAllocated) m_OffdiagonalValuesHandle.Free();
                m_DiagonalRowPointers = null;
                m_DiagonalColumnIndices = null;
                m_DiagonalValues = null;
                m_OffdiagonalRowPointers = null;
                m_OffdiagonalColumnIndices = null;
                m_OffdiagonalValues = null;
            }
        }

        static void ValidateRow(int rowLength, long[] columns, double[] values) {
            if(rowLength < 0 || (rowLength > 0 && (columns == null || values == null || columns.Length < rowLength || values.Length < rowLength)))
                throw new InvalidOperationException("The BoSSS matrix returned invalid row storage.");
        }

        static void ValidateColumn(long column, int globalRows) {
            if(column < 0 || column >= globalRows || column > int.MaxValue)
                throw new NotSupportedException($"PETSc column index {column} is outside the supported 32-bit matrix range.");
        }

        static void Check(int errorCode, UnsafePETSc native, string operation) {
            if(errorCode == 0)
                return;
            string message = null;
            if(native.ErrorMessage(errorCode, out IntPtr text, IntPtr.Zero) == 0 && text != IntPtr.Zero)
                message = Marshal.PtrToStringAnsi(text);
            throw new ApplicationException($"{operation} failed with PETSc error code {errorCode}{(string.IsNullOrEmpty(message) ? string.Empty : ": " + message)}");
        }
    }
}
