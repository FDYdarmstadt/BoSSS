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

using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// The compatibility CSR path: export one local CSR block, then let PETSc copy it.
    /// </summary>
    public abstract partial class PETScSolverBase {

        sealed class CsrArrays {
            internal int[] RowPointers;
            internal int[] ColumnIndices;
            internal double[] Values;
        }

        // Normally PETSc copies this data immediately. Pinning is available only
        // for diagnostics and is released together with the PETSc matrix.
        int[] m_PinnedRowPointers;
        int[] m_PinnedColumnIndices;
        double[] m_PinnedValues;
        GCHandle m_RowPointersHandle;
        GCHandle m_ColumnIndicesHandle;
        GCHandle m_ValuesHandle;

        CsrArrays ExportCsr(IMutableMatrixEx matrix, long globalRows, long globalRowStart) {
            var result = new CsrArrays {
                RowPointers = new int[m_LocalRows + 1]
            };
            var columns = new List<int>();
            var values = new List<double>();
            long[] rowColumns = null;
            double[] rowValues = null;
            for(int localRow = 0; localRow < m_LocalRows; localRow++) {
                int rowLength = matrix.GetRow(globalRowStart + localRow, ref rowColumns, ref rowValues);
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    ValidateCsrColumn(column, globalRows);
                    columns.Add((int)column);
                    values.Add(rowValues[entry]);
                }
                result.RowPointers[localRow + 1] = columns.Count;
            }
            result.ColumnIndices = columns.Count > 0 ? columns.ToArray() : new int[1];
            result.Values = values.Count > 0 ? values.ToArray() : new double[1];
            return result;
        }

        void PinCsrBuffers(int[] rowPointers, int[] columnIndices, double[] values) {
            m_PinnedRowPointers = rowPointers;
            m_PinnedColumnIndices = columnIndices;
            m_PinnedValues = values;
            m_RowPointersHandle = GCHandle.Alloc(m_PinnedRowPointers, GCHandleType.Pinned);
            m_ColumnIndicesHandle = GCHandle.Alloc(m_PinnedColumnIndices, GCHandleType.Pinned);
            m_ValuesHandle = GCHandle.Alloc(m_PinnedValues, GCHandleType.Pinned);
        }

        void ReleaseCopiedCsrBuffers() {
            if(m_RowPointersHandle.IsAllocated)
                m_RowPointersHandle.Free();
            if(m_ColumnIndicesHandle.IsAllocated)
                m_ColumnIndicesHandle.Free();
            if(m_ValuesHandle.IsAllocated)
                m_ValuesHandle.Free();

            m_PinnedRowPointers = null;
            m_PinnedColumnIndices = null;
            m_PinnedValues = null;
        }
    }
}
