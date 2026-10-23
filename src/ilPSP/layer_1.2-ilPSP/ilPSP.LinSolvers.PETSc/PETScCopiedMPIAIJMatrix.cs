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
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.PETSc {

    /// <summary>
    /// Owns a PETSc AIJ matrix populated from copied BoSSS CSR data.
    /// </summary>
    /// <remarks>
    /// Uses PETSc's public <c>Mat*AIJSetPreallocationCSR</c> APIs. PETSc copies
    /// and assembles the supplied CSR arrays, so this owner has no pinned-buffer
    /// lifetime after <see cref="Create"/> returns.
    /// </remarks>
    internal sealed unsafe class PETScCopiedMPIAIJMatrix : PETScMPIAIJMatrix {

        sealed class Csr {
            internal int[] RowPointers;
            internal int[] ColumnIndices;
            internal double[] Values;
        }

        readonly UnsafePETSc m_Native;
        readonly int m_LocalRows;
        IntPtr m_Handle;

        /// <summary>Native PETSc matrix handle.</summary>
        public override IntPtr Handle => m_Handle;

        /// <summary>Number of matrix rows owned by this rank.</summary>
        public override int LocalRows => m_LocalRows;

        PETScCopiedMPIAIJMatrix(UnsafePETSc native, Dimensions dimensions) {
            m_Native = native ?? throw new ArgumentNullException(nameof(native));
            m_LocalRows = dimensions.LocalRows;
        }

        /// <summary>Creates an assembled sequential or distributed AIJ matrix.</summary>
        internal static PETScCopiedMPIAIJMatrix Create(UnsafePETSc native, IntPtr communicator,
            IMutableMatrixEx matrix, Dimensions dimensions) {
            var result = new PETScCopiedMPIAIJMatrix(native, dimensions);
            try {
                Csr csr = ExportCsr(matrix, dimensions.LocalRows, dimensions.GlobalColumns, dimensions.GlobalRowStart);
                Check(native.CreateAijMatrix(communicator, dimensions.LocalRows, dimensions.LocalColumns,
                    dimensions.GlobalRows, dimensions.GlobalColumns,
                    0, IntPtr.Zero, 0, IntPtr.Zero, out result.m_Handle), native, "PETSc MatCreateAIJ");
                fixed(int* pRowPointers = csr.RowPointers)
                fixed(int* pColumnIndices = csr.ColumnIndices)
                fixed(double* pValues = csr.Values) {
                    if(dimensions.MpiSize == 1)
                        Check(native.SetSequentialAijPreallocationCsr(result.m_Handle, (IntPtr)pRowPointers, (IntPtr)pColumnIndices, (IntPtr)pValues), native, "PETSc MatSeqAIJSetPreallocationCSR");
                    else
                        Check(native.SetMpiAijPreallocationCsr(result.m_Handle, (IntPtr)pRowPointers, (IntPtr)pColumnIndices, (IntPtr)pValues), native, "PETSc MatMPIAIJSetPreallocationCSR");
                }
                return result;
            } catch {
                result.Dispose();
                throw;
            }
        }

        static Csr ExportCsr(IMutableMatrixEx matrix, int localRows, int globalRows, int globalRowStart) {
            if(matrix == null)
                throw new ArgumentNullException(nameof(matrix));
            var columns = new List<int>();
            var values = new List<double>();
            var rowPointers = new int[localRows + 1];
            long[] rowColumns = null;
            double[] rowValues = null;
            for(int localRow = 0; localRow < localRows; localRow++) {
                int rowLength = matrix.GetRow(globalRowStart + localRow, ref rowColumns, ref rowValues);
                if(rowLength < 0 || (rowLength > 0 && (rowColumns == null || rowValues == null || rowColumns.Length < rowLength || rowValues.Length < rowLength)))
                    throw new InvalidOperationException("The BoSSS matrix returned invalid row storage.");
                int previousColumn = -1;
                for(int entry = 0; entry < rowLength; entry++) {
                    long column = rowColumns[entry];
                    if(column < 0 || column >= globalRows || column > int.MaxValue)
                        throw new NotSupportedException($"PETSc column index {column} is outside the supported 32-bit matrix range.");
                    int columnIndex = (int)column;
                    if(columnIndex <= previousColumn)
                        throw new InvalidOperationException("The PETSc copied CSR transfer requires strictly ascending, duplicate-free column indices in every BoSSS matrix row.");
                    columns.Add(columnIndex);
                    values.Add(rowValues[entry]);
                    previousColumn = columnIndex;
                }
                rowPointers[localRow + 1] = columns.Count;
            }
            return new Csr {
                RowPointers = rowPointers,
                ColumnIndices = columns.Count == 0 ? new int[1] : columns.ToArray(),
                Values = values.Count == 0 ? new double[1] : values.ToArray()
            };
        }

        public override void Dispose() {
            if(m_Handle != IntPtr.Zero)
                Check(m_Native.DestroyMatrix(ref m_Handle), m_Native, "PETSc MatDestroy");
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
