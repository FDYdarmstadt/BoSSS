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
    /// The fast MPI CSR path: pack PETSc's split layout and let PETSc borrow it.
    /// </summary>
    public abstract partial class PETScSolverBase {

        sealed class SplitCsrArrays {
            internal int[] DiagonalRowPointers;
            internal int[] DiagonalColumnIndices;
            internal double[] DiagonalValues;
            internal int[] OffdiagonalRowPointers;
            internal int[] OffdiagonalColumnIndices;
            internal double[] OffdiagonalValues;
        }

        // PETSc retains these arrays, so the GC handles must stay allocated for
        // exactly the lifetime of the associated native matrix.
        int[] m_PinnedDiagonalRowPointers;
        int[] m_PinnedDiagonalColumnIndices;
        double[] m_PinnedDiagonalValues;
        int[] m_PinnedOffdiagonalRowPointers;
        int[] m_PinnedOffdiagonalColumnIndices;
        double[] m_PinnedOffdiagonalValues;
        GCHandle m_DiagonalRowPointersHandle;
        GCHandle m_DiagonalColumnIndicesHandle;
        GCHandle m_DiagonalValuesHandle;
        GCHandle m_OffdiagonalRowPointersHandle;
        GCHandle m_OffdiagonalColumnIndicesHandle;
        GCHandle m_OffdiagonalValuesHandle;

        SplitCsrArrays ExportSplitCsr(IMutableMatrixEx matrix, long globalRows, long globalRowStart) {
            if(matrix.RowPartitioning.TotalLength != matrix.ColPartitioning.TotalLength
                || matrix.RowPartitioning.i0 != matrix.ColPartitioning.i0
                || matrix.RowPartitioning.LocalLength != matrix.ColPartitioning.LocalLength)
                throw new NotSupportedException("PETSc split CSR requires identical row and column partitionings.");

            long firstOwnedColumn = matrix.ColPartitioning.i0;
            long lastOwnedColumn = firstOwnedColumn + matrix.ColPartitioning.LocalLength;
            var result = new SplitCsrArrays {
                DiagonalRowPointers = new int[m_LocalRows + 1],
                OffdiagonalRowPointers = new int[m_LocalRows + 1]
            };
            long[] rowColumns = null;
            double[] rowValues = null;

            // 1. Count exact storage. A zero diagonal is inserted where needed:
            // PETSc factorization and ILU-like preconditioners require the slot.
            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                long globalRow = globalRowStart + localRow;
                int rowLength = matrix.GetRow(globalRow, ref rowColumns, ref rowValues);
                int diagonalCount = 0;
                int offdiagonalCount = 0;
                bool hasStructuralDiagonal = false;
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    ValidateCsrColumn(column, globalRows);
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

            // 2. Allocate the two PETSc CSR blocks and fill them. The diagonal
            // block has rank-local column indices; remote columns stay global.
            int diagonalNonzeros = result.DiagonalRowPointers[m_LocalRows];
            int offdiagonalNonzeros = result.OffdiagonalRowPointers[m_LocalRows];
            result.DiagonalColumnIndices = new int[Math.Max(diagonalNonzeros, 1)];
            result.DiagonalValues = new double[Math.Max(diagonalNonzeros, 1)];
            result.OffdiagonalColumnIndices = new int[Math.Max(offdiagonalNonzeros, 1)];
            result.OffdiagonalValues = new double[Math.Max(offdiagonalNonzeros, 1)];
            int[] diagonalOffsets = (int[])result.DiagonalRowPointers.Clone();
            int[] offdiagonalOffsets = (int[])result.OffdiagonalRowPointers.Clone();

            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
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

        void PinSplitCsrBuffers(SplitCsrArrays csr) {
            m_PinnedDiagonalRowPointers = csr.DiagonalRowPointers;
            m_PinnedDiagonalColumnIndices = csr.DiagonalColumnIndices;
            m_PinnedDiagonalValues = csr.DiagonalValues;
            m_PinnedOffdiagonalRowPointers = csr.OffdiagonalRowPointers;
            m_PinnedOffdiagonalColumnIndices = csr.OffdiagonalColumnIndices;
            m_PinnedOffdiagonalValues = csr.OffdiagonalValues;
            m_DiagonalRowPointersHandle = GCHandle.Alloc(m_PinnedDiagonalRowPointers, GCHandleType.Pinned);
            m_DiagonalColumnIndicesHandle = GCHandle.Alloc(m_PinnedDiagonalColumnIndices, GCHandleType.Pinned);
            m_DiagonalValuesHandle = GCHandle.Alloc(m_PinnedDiagonalValues, GCHandleType.Pinned);
            m_OffdiagonalRowPointersHandle = GCHandle.Alloc(m_PinnedOffdiagonalRowPointers, GCHandleType.Pinned);
            m_OffdiagonalColumnIndicesHandle = GCHandle.Alloc(m_PinnedOffdiagonalColumnIndices, GCHandleType.Pinned);
            m_OffdiagonalValuesHandle = GCHandle.Alloc(m_PinnedOffdiagonalValues, GCHandleType.Pinned);
        }

        void ReleaseBorrowedCsrBuffers() {
            if(m_DiagonalRowPointersHandle.IsAllocated)
                m_DiagonalRowPointersHandle.Free();
            if(m_DiagonalColumnIndicesHandle.IsAllocated)
                m_DiagonalColumnIndicesHandle.Free();
            if(m_DiagonalValuesHandle.IsAllocated)
                m_DiagonalValuesHandle.Free();
            if(m_OffdiagonalRowPointersHandle.IsAllocated)
                m_OffdiagonalRowPointersHandle.Free();
            if(m_OffdiagonalColumnIndicesHandle.IsAllocated)
                m_OffdiagonalColumnIndicesHandle.Free();
            if(m_OffdiagonalValuesHandle.IsAllocated)
                m_OffdiagonalValuesHandle.Free();

            m_PinnedDiagonalRowPointers = null;
            m_PinnedDiagonalColumnIndices = null;
            m_PinnedDiagonalValues = null;
            m_PinnedOffdiagonalRowPointers = null;
            m_PinnedOffdiagonalColumnIndices = null;
            m_PinnedOffdiagonalValues = null;
        }
    }
}
