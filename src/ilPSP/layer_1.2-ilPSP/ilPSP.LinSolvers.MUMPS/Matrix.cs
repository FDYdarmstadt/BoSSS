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
using System.IO;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MPI.Wrappers;


namespace ilPSP.LinSolvers.MUMPS {
    class Matrix {

        /// <summary>
        /// row partiton of the matrix;
        /// </summary>
        public IPartitioning RowPart;

        /// <summary>
        /// Dimension of the MUMPS Matrix
        /// </summary>
        public int n;

        /// <summary>
        /// Number of nonzero entries in MUMPS Matrix
        /// </summary>
        public int nz;

        /// <summary>
        /// Array of row indices
        /// </summary>
        public int[] irn;

        /// <summary>
        /// Array of col indices
        /// </summary>
        public int[] jrn;

        /// <summary>
        /// Array of values
        /// </summary>
        public double[] a;

        public double[] a_loc;

        public int[] irn_loc;

        public int[] jrn_loc;

        public int n_loc;

        public int nz_loc;

        /// <summary>
        /// 0: unsymmetric
        /// 1: assumed to be symmetric positive definite so that pivots are taken from the diagonal
        /// without numerical pivoting during the factorization.With this option, non-positive definite
        /// matrices that do not require pivoting can also be treated in certain cases (see remark below).
        /// 2: general symmetric
        /// </summary>
        public int Symmetric = 0;

        public Matrix(IMutableMatrixEx M) {
            if (M.RowPartitioning.IsMutable)
                throw new NotSupportedException();
            if (M.ColPartition.IsMutable)
                throw new NotSupportedException();

            if (M.NoOfCols != M.NoOfRows)
                throw new ArgumentException("Matrix must be quadratic.", "M");
            if ((M is MsrMatrix) && ((MsrMatrix)M).AssumeSymmetric)
                this.Symmetric = 2;
            RowPart = M.RowPartitioning;

            int size = M.RowPartitioning.MpiSize, rank = M.RowPartitioning.MpiRank;
            Debug.Assert(M.RowPartitioning.MpiRank == M.ColPartition.MpiRank);
            Debug.Assert(M.RowPartitioning.MpiSize == M.ColPartition.MpiSize);
            var comm = M.MPI_Comm;

            int LR;
            int[] col = null;
            double[] val = null;

            if (size == 0) {
                // serial init on one processor
                // ++++++++++++++++++++++++++++

                n = (int)M.RowPartitioning.TotalLength;

                int len;
                if (Symmetric ==2) {
                    // upper triangle + diagonal (diagonal entries are 
                    // always required, even if 0.0, for symmetric matrices in PARDISO)

                    len = M.GetGlobalNoOfUpperTriangularNonZeros() + n;
                } else {
                    len = M.GetTotalNoOfNonZerosPerProcess();
                }
                int Nrows = M.RowPartitioning.LocalLength;
                nz = len;
                int cnt = 0;
                irn = new int[len];
                jrn = new int[len];
                a = new double[len];
                for (int i = 0; i < Nrows; i++) {
                    //irn[i] = cnt;
                    int iRow = M.RowPartitioning.i0 + i;

                    LR = M.GetRow(iRow, ref col, ref val);

                    double diagelem = M[iRow, iRow];
                    if (Symmetric==2 && diagelem == 0) {
                        // in the symmetric case, we always need to provide the diagonal element
                        jrn[cnt] = iRow; 
                        a[cnt] = 0.0;
                        cnt++;
                    }
                    for (int j = 0; j < LR; j++) {

                        if (val[j] != 0.0) {

                            if (Symmetric==2 && col[j] < iRow) {
                                // entry is in lower triangular matrix -> ignore (for symmetric mtx.)
                                continue;
                            } else {
                                irn[cnt] = iRow +1;
                                jrn[cnt] = col[j] +1; 
                                a[cnt] = val[j]; 
                                cnt++;
                            }
                        }
                    }

                    //if (M.GetTotalNoOfNonZeros() != len)
                    //    throw new Exception();
                }

                if (len != cnt)
                    throw new ApplicationException("internal error.");

            } else {
                // collect local matrices 
                // +++++++++++++++++++++++

                // Number of elements, start indices for index pointers
                // ====================================================

                int len_loc;
                if (Symmetric == 2)
                    // number of entries is:
                    // upper triangle + diagonal (diagonal entries are 
                    // always required, even if 0.0, for symmetric matrices in PARDISO)
                    len_loc = M.GetLocalNoOfUpperTriangularNonZeros() + M.RowPartitioning.LocalLength;
                else
                    len_loc = M.GetTotalNoOfNonZerosPerProcess();

                Partitioning part = new Partitioning(len_loc, comm);
                if (part.TotalLength > int.MaxValue)
                    throw new ApplicationException("too many matrix entries for MUMPS - more than maximum 32-bit signed integer");


                // local matrix assembly
                // =====================
                int n_loc = M.RowPartitioning.LocalLength;
                this.nz_loc = len_loc ;
                int[] ia_loc = new int[len_loc];
                int[] ja_loc = new int[len_loc];
                double[] a_loc = new double[len_loc];
                {

                    int cnt = 0;
                    int i0 = (int)part.i0;

                    for (int i = 0; i < n_loc; i++) {
                        //ia_loc[i] = cnt + 1 + i0; // fortran indexing
                        int iRow = i + (int)M.RowPartitioning.i0;

                        LR = M.GetRow(iRow, ref col, ref val);

                        double diagelem = M[iRow, iRow];
                        if (Symmetric == 2 && diagelem == 0) {
                            // in the symmetric case, we always need to provide the diagonal element
                            ja_loc[cnt] = iRow + 1; // fortran indexing
                            a_loc[cnt] = 0.0;
                            cnt++;
                        }
                        for (int j = 0; j < LR; j++) {

                            if (val[j] != 0.0) {

                                if (Symmetric ==2 && col[j] < iRow) {
                                    // entry is in lower triangular matrix -> ignore (for symmetric mtx.)
                                    continue;
                                } else {
                                    ia_loc[cnt] = iRow + 1;
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

                nz = part.TotalLength;

                n = (int)M.RowPartitioning.TotalLength;

                this.a_loc = a_loc;
                this.irn_loc = ia_loc;
                this.jrn_loc = ja_loc;
                this.n_loc = (int)M.RowPartitioning.TotalLength;

                GC.Collect();
            }
        }
    }
}
