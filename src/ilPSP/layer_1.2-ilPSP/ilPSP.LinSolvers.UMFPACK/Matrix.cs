/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using System.Diagnostics;
using System.Runtime.InteropServices;
using MPI.Wrappers;
using ilPSP.Tracing;

namespace ilPSP.LinSolvers.UMFPACK {

    /// <summary>
    /// 0-based compressed sparse storage for UMFPACK.
    /// </summary>
    public class Matrix : IDisposable {

        /// <summary>
        /// row partition of the matrix.
        /// </summary>
        public IPartitioning RowPart;

        /// <summary>
        /// Matrix dimension.
        /// </summary>
        public long n;

        /// <summary>
        /// Column pointers for UMFPACK. This stores CSR row pointers of A, interpreted as CSC column pointers of A'.
        /// </summary>
        public long[] Ap;

        /// <summary>
        /// Row indices for UMFPACK. This stores CSR column indices of A, interpreted as CSC row indices of A'.
        /// </summary>
        public long[] Ai;

        /// <summary>
        /// Matrix entries. Only initialized on MPI processor rank 0.
        /// </summary>
        public IntPtr AxPtr;

        MPI_Comm m_comm;

        /// <summary>
        /// Initializes this matrix as a copy of <paramref name="M"/>.
        /// </summary>
        public Matrix(IMutableMatrixEx M) {
            using(var tr = new FuncTrace()) {
                if(M.RowPartitioning.IsMutable)
                    throw new NotSupportedException();
                if(M.ColPartitioning.IsMutable)
                    throw new NotSupportedException();
                if(M.NoOfCols != M.NoOfRows)
                    throw new ArgumentException("Matrix must be quadratic.", "M");

                RowPart = M.RowPartitioning;
                int size = M.RowPartitioning.MpiSize, rank = M.RowPartitioning.MpiRank;
                Debug.Assert(M.RowPartitioning.MpiRank == M.ColPartitioning.MpiRank);
                Debug.Assert(M.RowPartitioning.MpiSize == M.ColPartitioning.MpiSize);
                m_comm = M.MPI_Comm;

                long[] col = null;
                double[] val = null;

                if(size == 1) {
                    using(new BlockTrace("serial init", tr)) {
                        n = M.RowPartitioning.TotalLength;
                        long len = M.GetTotalNoOfNonZerosPerProcess();
                        int Nrows = M.RowPartitioning.LocalLength;

                        Ap = new long[n + 1];
                        Ai = new long[len];
                        AxPtr = Marshal.AllocHGlobal((IntPtr)(len * sizeof(double)));

                        unsafe {
                            double* Ax = (double*)AxPtr;
                            long cnt = 0;
                            for(int i = 0; i < Nrows; i++) {
                                Ap[i] = cnt;
                                long iRow = M.RowPartitioning.i0 + i;
                                int LR = M.GetRow(iRow, ref col, ref val);

                                for(int j = 0; j < LR; j++) {
                                    if(val[j] != 0.0) {
                                        Ai[cnt] = col[j];
                                        Ax[cnt] = val[j];
                                        cnt++;
                                    }
                                }
                            }
                            Ap[Nrows] = cnt;
                            if(len != cnt)
                                throw new ApplicationException("internal error.");
                        }
                    }
                } else {
                    using(new BlockTrace("Collect matrix on proc 0", tr)) {
                        int len_loc = checked((int)M.GetTotalNoOfNonZerosPerProcess());
                        Partitioning part = new Partitioning(len_loc, m_comm);

                        int n_loc = M.RowPartitioning.LocalLength;
                        long[] Ap_loc = new long[n_loc];
                        long[] Ai_loc = new long[len_loc];
                        double[] Ax_loc = new double[len_loc];

                        using(new BlockTrace("local matrix assembly", tr)) {
                            int cnt = 0;
                            long i0 = part.i0;

                            for(int i = 0; i < n_loc; i++) {
                                Ap_loc[i] = cnt + i0;
                                long iRow = i + M.RowPartitioning.i0;
                                int LR = M.GetRow(iRow, ref col, ref val);

                                for(int j = 0; j < LR; j++) {
                                    if(val[j] != 0.0) {
                                        Ai_loc[cnt] = col[j];
                                        Ax_loc[cnt] = val[j];
                                        cnt++;
                                    }
                                }
                            }

                            if(cnt != len_loc)
                                throw new ApplicationException("internal error.");
                        }

                        if(rank == 0) {
                            n = M.RowPartitioning.TotalLength;
                            Ap = new long[n + 1];
                            Ai = new long[part.TotalLength];
                            AxPtr = Marshal.AllocHGlobal((IntPtr)(part.TotalLength * sizeof(double)));
                        } else {
                            AxPtr = IntPtr.Zero;
                        }

                        unsafe {
                            int* displs = stackalloc int[size];
                            int* recvcounts = stackalloc int[size];
                            for(int i = 0; i < size; i++) {
                                recvcounts[i] = part.GetLocalLength(i);
                                displs[i] = checked((int)part.GetI0Offest(i));
                            }

                            fixed(double* pAx_loc = Ax_loc) {
                                double* pAx = (double*)AxPtr;
                                csMPI.Raw.Gatherv((IntPtr)pAx_loc, Ax_loc.Length, csMPI.Raw._DATATYPE.DOUBLE,
                                    (IntPtr)pAx, (IntPtr)recvcounts, (IntPtr)displs, csMPI.Raw._DATATYPE.DOUBLE, 0, m_comm);
                            }

                            fixed(long* pAi_loc = Ai_loc, pAi = Ai) {
                                csMPI.Raw.Gatherv((IntPtr)pAi_loc, Ai_loc.Length, csMPI.Raw._DATATYPE.LONG_LONG,
                                    (IntPtr)pAi, (IntPtr)recvcounts, (IntPtr)displs, csMPI.Raw._DATATYPE.LONG_LONG, 0, m_comm);
                            }

                            for(int i = 0; i < size; i++) {
                                displs[i] = checked((int)M.RowPartitioning.GetI0Offest(i));
                                recvcounts[i] = M.RowPartitioning.GetLocalLength(i);
                            }

                            fixed(long* pAp_loc = Ap_loc, pAp = Ap) {
                                csMPI.Raw.Gatherv((IntPtr)pAp_loc, Ap_loc.Length, csMPI.Raw._DATATYPE.LONG_LONG,
                                    (IntPtr)pAp, (IntPtr)recvcounts, (IntPtr)displs, csMPI.Raw._DATATYPE.LONG_LONG, 0, m_comm);
                            }
                        }

                        if(rank == 0)
                            Ap[M.RowPartitioning.TotalLength] = part.TotalLength;
                    }
                }
            }
        }

        /// <summary>
        /// Release unmanaged matrix values.
        /// </summary>
        public void Dispose() {
            if(AxPtr != IntPtr.Zero)
                Marshal.FreeHGlobal(AxPtr);
            AxPtr = IntPtr.Zero;
        }
    }
}
