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
using MPI.Wrappers;

namespace ilPSP.LinSolvers.HYPRE {
    
    /// <summary>
    /// object-oriented wrapper around HYPRE_IJMatrix objects
    /// </summary>
    public class IJMatrix : IDisposable, IMutableMatrix {

        /// <summary>
        /// calls dispose
        /// </summary>
        ~IJMatrix() {
            Dispose();
        }


        /// <summary>
        /// handle/pointer of the IJMatrix object;
        /// </summary>
        internal Wrappers.T_IJMatrix m_IJMatrix;

        /// <summary>
        /// handle/pointer of the ParCSR Matrix object which is used by  the
        /// IJMatrix object
        /// </summary>
        internal Wrappers.T_ParCSR_matrix m_ParCSR_matrix;
        
        IPartitioning m_RowPartition;
        IPartitioning m_ColPartition;

        /// <summary>
        /// partition of rows of this matrix;
        /// </summary>
        public IPartitioning RowPartitioning { get { return m_RowPartition; } }

        /// <summary>
        /// number of rows, over all mpi processors
        /// </summary>
        public int NoOfRows {
            get { return (int)RowPartitioning.TotalLength; }
        }

        /// <summary>
        /// number of columns, over all mpi processors
        /// </summary>
        public int NoOfCols {
            get { return (int)ColPartition.TotalLength; }
        }

        /// <summary>
        /// initializes this matrix to be a copy of <paramref name="mtx"/>;
        /// </summary>
        /// <param name="mtx"></param>
        public IJMatrix(IMutableMatrixEx mtx) {

            if (mtx.RowPartitioning.MPI_Comm != mtx.ColPartition.MPI_Comm)
                throw new ArgumentException();
            if (mtx.RowPartitioning.IsMutable)
                throw new ArgumentException();
            if (mtx.ColPartition.IsMutable)
                throw new ArgumentException();

            m_RowPartition = mtx.RowPartitioning;
            m_ColPartition = mtx.ColPartition;

          
            if (mtx.NoOfRows != mtx.NoOfCols)
                throw new ArgumentException("matrix must be quadratic.", "mtx");
            
            if (mtx.NoOfRows > (int.MaxValue - 2))
                throw new ApplicationException("unable to create HYPRE matrix: no. of matrix rows is larger than HYPRE index type (32 bit signed int);");
            if (mtx.NoOfCols > (int.MaxValue - 2))
                throw new ApplicationException("unable to create HYPRE matrix: no. of matrix columns is larger than HYPRE index type (32 bit signed int);");

            
            // matrix: init
            MPI_Comm comm = csMPI.Raw._COMM.WORLD;
            int ilower = (int) mtx.RowPartitioning.i0;
            int iupper = (int) ilower + mtx.RowPartitioning.LocalLength - 1;
            int jlower = (int) mtx.ColPartition.i0;
            int jupper = (int) jlower + mtx.ColPartition.LocalLength - 1;

            HypreException.Check(Wrappers.IJMatrix.Create(comm, ilower, iupper, jlower, jupper, out m_IJMatrix));
            HypreException.Check(Wrappers.IJMatrix.SetObjectType(m_IJMatrix, Wrappers.Constants.HYPRE_PARCSR));
            HypreException.Check(Wrappers.IJMatrix.Initialize(m_IJMatrix));

            // matrix: set values, row by row ...
            int nrows, lmax = mtx.GetMaxNoOfNonZerosPerRow();
            int[] rows = new int[1], cols = new int[lmax], ncols = new int[1];
            double[] values = new double[lmax];
            int LR;
            int[] col = null;
            double[] val = null;
            for (int i = 0; i < mtx.RowPartitioning.LocalLength; i++) {

                int iRowGlob = i + mtx.RowPartitioning.i0;
                LR = mtx.GetRow(iRowGlob, ref col, ref val);
                
                nrows = 1;
                rows[0] = iRowGlob;
                ncols[0] = LR;
                int cnt = 0;
                for (int j = 0; j < LR; j++) {
                    if (val[j] != 0.0) {
                        cols[cnt] = col[j];
                        values[cnt] = val[j];
                        cnt++;
                    }
                }
                if (cnt <= 0)
                    throw new ArgumentException(string.Format("Zero matrix row detected (local row index: {0}, global row index: {1}).",i,iRowGlob));

                HypreException.Check(Wrappers.IJMatrix.SetValues(m_IJMatrix, nrows, ncols, rows, cols, values));
            }

            // matrix: assembly
            HypreException.Check(Wrappers.IJMatrix.Assemble(m_IJMatrix));
            HypreException.Check(Wrappers.IJMatrix.GetObject(m_IJMatrix, out m_ParCSR_matrix));
        }

        /// <summary>
        /// returns the diagonal element in the <paramref name="row"/>-th row.
        /// </summary>
        /// <param name="row">global row/column index</param>
        /// <returns>value of diagonal element</returns>
        public double GetDiagonalElement(int row) {
            if (row < m_RowPartition.i0
                || row >= (m_RowPartition.i0 + m_RowPartition.LocalLength))
                throw new IndexOutOfRangeException("row index is not assigned to current processor.");

            int[] ncols = new int[] { 1 };
            int[] rows = new int[] { row };
            int[] cols = new int[] { row };
            double[] ret = new double[1];
            HypreException.Check(Wrappers.IJMatrix.GetValues(m_IJMatrix, 1, ncols, rows, cols, ret));

            return ret[0];
        }

        /// <summary>
        /// sets the diagonal element in the <paramref name="row"/>-th row
        /// to value <paramref name="val"/>
        /// </summary>
        /// <param name="row">global row/column index</param>
        /// <param name="val">new value of diagonal element</param>
        public void SetDiagonalElement(int row, double val) {
            if (row < m_RowPartition.i0
                || row >= (m_RowPartition.i0 + m_RowPartition.LocalLength))
                throw new IndexOutOfRangeException("row index is not assigned to current processor.");
            
            int[] ncols = new int[] { 1 };
            int[] rows = new int[] { row };
            int[] cols = new int[] { row };
            double[] _val = new double[] { val };
            HypreException.Check(Wrappers.IJMatrix.SetValues(m_IJMatrix, 1, ncols, rows, cols, _val));

            //// test code:
            //double vchk = GetDiagElement(row);
            //if (vchk != val)
            //    throw new ApplicationException();
        }



        /// <summary>
        /// destroys the HYPRE objects, if not already done;
        /// </summary>
        public void Dispose() {
            if (m_IJMatrix.p != IntPtr.Zero) {
                Wrappers.IJMatrix.Destroy(m_IJMatrix);
                
                m_IJMatrix.p = IntPtr.Zero;
                m_ParCSR_matrix.p = IntPtr.Zero;
            }
        }
        
        /// <summary>
        /// The column partition defines how a vector that should be
        /// multiplied with this matrix (matrix times column vector)
        /// must be partitioned.
        /// </summary>
        public IPartitioning ColPartition {
            get { return m_ColPartition; }
        }

        /// <summary>
        /// MPI Communicator on which this object lives on.
        /// </summary>
        public MPI_Comm MPI_Comm {
            get {
                if (RowPartitioning.MPI_Comm != ColPartition.MPI_Comm)
                    throw new ApplicationException("Internal error -- mismatch between row and column MPI communicator.");
                return RowPartitioning.MPI_Comm;
            }
        }

        /// <summary>
        /// general matrix/vector product; see <see cref="ISparseMatrix.SpMV"/>;
        /// </summary>
        /// <typeparam name="VectorType1"></typeparam>
        /// <typeparam name="VectorType2"></typeparam>
        /// <param name="alpha"></param>
        /// <param name="a"></param>
        /// <param name="beta"></param>
        /// <param name="acc"></param>
        public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : System.Collections.Generic.IList<double>
            where VectorType2 : System.Collections.Generic.IList<double> {

            if (a.Count < ColPartition.LocalLength)
                throw new ArgumentException("length/count of 'a' must be equal to local length of column partition", "a");
            if (acc.Count < RowPartitioning.LocalLength)
                throw new ArgumentException("length/count of 'acc' must be greater or equal to local length of row partition", "acc");
            if (object.ReferenceEquals(a, acc))
                throw new ArgumentException("in-place computation is not supported.", "a,acc");


            IJVector _acc = new IJVector(RowPartitioning);
            if( beta != 0.0)
                _acc.SetValues(acc);
            IJVector _a = new IJVector(ColPartition);
            _a.SetValues(a);
            
            HypreException.Check(Wrappers.ParCSRMatrix.ParCSRMatrixMatvec(alpha, m_ParCSR_matrix, _a.ParCRS_vector, beta, _acc.ParCRS_vector));

            _acc.GetValues(acc);
            
            _acc.Dispose();
            _a.Dispose();
        }


        #region IMutuableMatrix Members

        /// <summary>
        /// see <see cref="IMutableMatrix.GetValues"/>
        /// </summary>
        public double[] GetValues(int RowIndex, int[] ColumnIndices) {
            if( RowIndex < this.RowPartitioning.i0 || RowIndex >= (this.RowPartitioning.i0 + this.RowPartitioning.LocalLength))
                throw new ArgumentOutOfRangeException("RowIndex","row index not within local range");

            double[] ret = new double[ColumnIndices.Length];
            
            HypreException.Check(Wrappers.IJMatrix.GetValues(m_IJMatrix, 1, new int[] { ColumnIndices.Length }, new int[] { RowIndex }, ColumnIndices, ret));
            
            return ret;
        }

        /// <summary>
        /// see <see cref="IMutableMatrix.SetValues"/>
        /// </summary>
        public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
            if (RowIndex < this.RowPartitioning.i0 || RowIndex >= (this.RowPartitioning.i0 + this.RowPartitioning.LocalLength))
                throw new ArgumentOutOfRangeException("RowIndex", "row index not within local range");

            HypreException.Check(Wrappers.IJMatrix.SetValues(m_IJMatrix, 1, new int[] { ColumnIndices.Length }, new int[] { RowIndex }, ColumnIndices, newValues));
        }

        /// <summary>
        /// see <see cref="IMutableMatrix.this"/>
        /// </summary>
        public double this[int i, int j] {
            get {
                return GetValues(i, new int[] { j })[0];
            }
            set {
                SetValues(i, new int[] { j }, new double[] { value });
            }
        }

        /// <summary>
        /// Accumulates <paramref name="Block"/>*<paramref name="alpha"/> to this matrix,
        /// at the row/column offset <paramref name="i0"/> resp. <paramref name="j0"/>.
        /// </summary>
        /// <param name="i0">Row offset.</param>
        /// <param name="j0">Column offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation operation.</param>
        /// <param name="Block">Block to accumulate.</param>
        public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block) {
            this.AccBlock(i0, j0, alpha, Block, 1.0);
        }

        /// <summary>
        /// Accumulates a block of entries to this matrix.
        /// </summary>
        /// <param name="i0">Row index offset.</param>
        /// <param name="j0">Column index offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation.</param>
        /// <param name="Block">Block to add.</param>
        /// <param param name="beta">pre-scaling</param>
        public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block, double beta) {
            if (Block.Dimension != 2)
                throw new ArgumentException();
            int I = Block.NoOfRows;
            int J = Block.NoOfCols;

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    this[i0 + i, j0 + j] = this[i0 + i, j0 + j]*beta + alpha * Block[i, j];
        }

        /// <summary>
        /// see <see cref="IMutableMatrix.OccupationMutable"/>
        /// </summary>
        public bool OccupationMutable {
            get {
                return (m_ParCSR_matrix.p == IntPtr.Zero);
            }
        }

        #endregion
    }



}
