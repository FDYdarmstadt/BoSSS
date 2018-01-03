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
using System.Linq;
using System.IO;
using System.Diagnostics;
using MPI.Wrappers;
using System.Collections.Generic;

namespace ilPSP.LinSolvers.PARDISO {
    class Matrix {

        /// <summary>
        /// row partiton of the matrix;
        /// </summary>
        public IPartitioning RowPart;

        /// <summary>
        /// PARDISO parameter: matrix dimension
        /// </summary>
        public int n;

        /// <summary>
        /// PARDISO parameter: index into <see cref="ja"/>, for each row;
        /// FORTRAN indexing: starts at 1! only initialized on MPI processor rank 0;
        /// </summary>
        public int[] ia;

        /// <summary>
        /// PARDISO parameter: column indices
        /// Attention: FORTRAN indexing: starts at 1! only initialized on MPI processor rank 0;
        /// </summary>
        public int[] ja;
        
        /// <summary>
        /// PARDISO parameter: matrix values
        /// Attention: FORTRAN indexing: starts at 1! only initialized on MPI processor rank 0.
        /// </summary>
        public double[] a;


        public bool Symmetric = false;

        MPI_Comm m_comm;

        /// <summary>
        /// initializes this matrix as a copy of the matrix <paramref name="M"/>.
        /// </summary>
        /// <param name="M"></param>
        public Matrix(IMutableMatrixEx M)  {
            if (M.RowPartitioning.IsMutable)
                throw new NotSupportedException();
            if (M.ColPartition.IsMutable)
                throw new NotSupportedException();

            if (M.NoOfCols != M.NoOfRows)
                throw new ArgumentException("Matrix must be quadratic.","M");
            this.Symmetric = (M is MsrMatrix) && ((MsrMatrix)M).AssumeSymmetric;
            RowPart = M.RowPartitioning;
            
            int size = M.RowPartitioning.MpiSize, rank = M.RowPartitioning.MpiRank;
            Debug.Assert(M.RowPartitioning.MpiRank == M.ColPartition.MpiRank);
            Debug.Assert(M.RowPartitioning.MpiSize == M.ColPartition.MpiSize);
            m_comm = M.MPI_Comm;

            int LR;
            int[] col = null;
            double[] val = null;

            if (size == 1) {
                // serial init on one processor
                // ++++++++++++++++++++++++++++

                n = (int)M.RowPartitioning.TotalLength;

                int len;
                if (Symmetric) {
                    // upper triangle + diagonal (diagonal entries are 
                    // always required, even if 0.0, for symmetric matrices in PARDISO)

                    len = M.GetGlobalNoOfUpperTriangularNonZeros() + n;
                } else {
                    len = M.GetTotalNoOfNonZerosPerProcess();
                }
                int Nrows = M.RowPartitioning.LocalLength;

                int cnt = 0;
                ia = new int[n + 1];
                ja = new int[len];
                a = new double[len];
                for (int i = 0; i < Nrows; i++) {
                    ia[i] = cnt + 1; // fortran indexing
                    int iRow = M.RowPartitioning.i0 + i;

                    LR = M.GetRow(iRow, ref col, ref val);

                    double diagelem = M[iRow, iRow];
                    if (Symmetric && diagelem == 0) {
                        // in the symmetric case, we always need to provide the diagonal element
                        ja[cnt] = iRow + 1; // fortran indexing
                        a[cnt] = 0.0;
                        cnt++;
                    }
                    for (int j = 0; j < LR; j++) {

                        if (val[j] != 0.0) {

                            if (Symmetric && col[j] < iRow) {
                                // entry is in lower triangular matrix -> ignore (for symmetric mtx.)
                                continue;
                            } else {
                                ja[cnt] = col[j] + 1; // fortran indexing
                                a[cnt] = val[j];
                                cnt++;
                            }
                        }
                    }

                    //if (M.GetTotalNoOfNonZeros() != len)
                    //    throw new Exception();
                }
                ia[Nrows] = cnt + 1; // fortran indexing

                if (len != cnt)
                    throw new ApplicationException("internal error.");

            } else {
                // collect matrix on processor 0
                // +++++++++++++++++++++++++++++

                // Number of elements, start indices for index pointers
                // ====================================================
                                
                int len_loc;
                if (Symmetric)
                    // number of entries is:
                    // upper triangle + diagonal (diagonal entries are 
                    // always required, even if 0.0, for symmetric matrices in PARDISO)
                    len_loc = M.GetLocalNoOfUpperTriangularNonZeros() + M.RowPartitioning.LocalLength;
                else
                    len_loc = M.GetTotalNoOfNonZerosPerProcess();

                Partitioning part = new Partitioning(len_loc, m_comm);
                if (part.TotalLength > int.MaxValue)
                    throw new ApplicationException("too many matrix entries for PARDISO - more than maximum 32-bit signed integer");

                
                // local matrix assembly
                // =====================
                int n_loc = M.RowPartitioning.LocalLength;
                int[] ia_loc = new int[n_loc];
                int[] ja_loc = new int[len_loc];
                double[] a_loc = new double[len_loc];
                {

                    int cnt = 0;
                    int i0 = (int)part.i0;
                    
                    for (int i = 0; i < n_loc; i++) {
                        ia_loc[i] = cnt + 1 + i0; // fortran indexing
                        int iRow = i + (int)M.RowPartitioning.i0;

                        LR = M.GetRow(iRow, ref col, ref val);
                        
                        double diagelem = M[iRow, iRow];
                        if (Symmetric && diagelem == 0) {
                            // in the symmetric case, we always need to provide the diagonal element
                            ja_loc[cnt] = iRow + 1; // fortran indexing
                            a_loc[cnt] = 0.0;
                            cnt++;
                        }
                        for (int j = 0; j < LR; j++) {

                            if (val[j] != 0.0) {

                                if (Symmetric && col[j] < iRow) {
                                    // entry is in lower triangular matrix -> ignore (for symmetric mtx.)
                                    continue;
                                } else {
                                    ja_loc[cnt] = col[j] + 1; // fortran indexing
                                    a_loc[cnt] = val[j];
                                    cnt++;
                                }
                            }
                        }

                    }

                    if (cnt != len_loc)
                        throw new ApplicationException("internal error.");
                }

                // assemble complete matrix on proc. 0
                // ===================================
                if (rank == 0) {

                    n = M.RowPartitioning.TotalLength;

                    // process 0: collect data from other processors
                    // +++++++++++++++++++++++++++++++++++++++++++++

                    this.ia = new int[M.RowPartitioning.TotalLength + 1];
                    this.ja = new int[part.TotalLength];
                    this.a = new double[part.TotalLength];
                    
                    Array.Copy(ia_loc, 0, this.ia, 0, ia_loc.Length);
                    Array.Copy(ja_loc, 0, this.ja, 0, ja_loc.Length);
                    Array.Copy(a_loc, 0, this.a, 0, a_loc.Length);


                    unsafe { 
                        fixed (int* pia = &this.ia[0], pja = &this.ja[0]) {
                            fixed (double* pa = &this.a[0]) {

                                for (int rcv_rank = 1; rcv_rank < size; rcv_rank++) {
                                    MPI_Status status;
                                    csMPI.Raw.Recv((IntPtr)(pa + part.GetI0Offest(rcv_rank)), part.GetLocalLength(rcv_rank), csMPI.Raw._DATATYPE.DOUBLE, rcv_rank, 321555 + rcv_rank, m_comm, out status);
                                    csMPI.Raw.Recv((IntPtr)(pja + part.GetI0Offest(rcv_rank)), part.GetLocalLength(rcv_rank), csMPI.Raw._DATATYPE.INT, rcv_rank, 32155 + rcv_rank, m_comm, out status);
                                    csMPI.Raw.Recv((IntPtr)(pia + M.RowPartitioning.GetI0Offest(rcv_rank)), M.RowPartitioning.GetLocalLength(rcv_rank), csMPI.Raw._DATATYPE.INT, rcv_rank, 3215 + rcv_rank, m_comm, out status);
                                }
                            }
                        }
                    }

                    this.ia[M.RowPartitioning.TotalLength] = (int)part.TotalLength + 1;
                    
                } else {

                    // send data to process 0
                    // ++++++++++++++++++++++

                    unsafe {
                        fixed (void* pia = &ia_loc[0], pja = &ja_loc[0], pa = &a_loc[0]) {

                            csMPI.Raw.Send((IntPtr)pa, a_loc.Length, csMPI.Raw._DATATYPE.DOUBLE, 0, 321555 + rank, m_comm);
                            csMPI.Raw.Send((IntPtr)pja, ja_loc.Length, csMPI.Raw._DATATYPE.INT, 0, 32155 + rank, m_comm);
                            csMPI.Raw.Send((IntPtr)pia, ia_loc.Length, csMPI.Raw._DATATYPE.INT, 0, 3215 + rank, m_comm);
                        }
                    }
                }

                ia_loc = null; ja_loc = null; a_loc = null;
                GC.Collect();
            }

            //if(rank == 0) {
            //    Console.WriteLine("PARDISO.Matrix: Number of nonzeros: " + this.a.Length);
            //    Console.WriteLine("PARDISO.Matrix: sum of entries: " + this.a.Sum());
            //}

            //SaveToTextFile("C:\\tmp\\pard.txt");
        }


        /// <summary>
        /// Writes the matrix into a text file (for debugging purposes).
        /// Basically, it uses a tabulator separated format (which can e.g. be
        /// imported into Matlab via <code>matrix = dlmread(path)</code>
        /// </summary>
        /// <remarks>
        /// The content of the file specified via <paramref name="path"/> will
        /// be overwritten.
        /// 
        /// If running on more than one process, multiple files are written 
        /// (they can be concatenated e.g. by 'cat' or directly in Matlab).
        /// </remarks>
        /// <param name="path">Path to the text file</param>
        public void SaveToTextFile(string path) {

            int Rank, Size;
            csMPI.Raw.Comm_Size(m_comm, out Size);
            csMPI.Raw.Comm_Rank(m_comm, out Rank);
            string append = "";
            if (Rank >= 1)
                return;

            int NoOfCols = (int)RowPart.TotalLength; // assume quadratic matrix
            StreamWriter writer = new StreamWriter(path + append);
            int cnt = 0;
            for (int i = 0; i < RowPart.TotalLength; i++) {

                string separator = "";
                int currentColumn = -1;
                int len = ia[i + 1] - ia[i];
                for (int j = 0; j < len; j++) {

                    // Add zeros for missing columns (the sparse format does
                    // not store zero values)
                    for (int k = currentColumn + 1; k < (ja[cnt]-1); k++) {
                        writer.Write(separator + "{0,14:F0}", 0.0);
                        separator = "\t";
                    }

                    // Enforce use of . as decimal separator in scientific format
                    writer.Write(separator + a[cnt].ToString(
                        "E", System.Globalization.CultureInfo.InvariantCulture).PadLeft(14));
                    currentColumn = ja[cnt]-1;
                    cnt++;
                    separator = "\t";
                }

                // Add zeros for columns following after the last entry
                for (int j = currentColumn + 1; j < NoOfCols; j++) {
                    writer.Write("\t{0,14:F0}", 0.0);
                }

                writer.Write("\n");
            }
            writer.Close();
        }
    }
}
