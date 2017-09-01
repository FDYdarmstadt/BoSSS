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
using System.Collections.Generic;
using System.IO;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

namespace BoSSS.Platform {

    /// <summary>
    /// Special, block-diagonal sparse matrix
    /// </summary>
    public class BlockDiagonalMatrix : IMutableMatrixEx, ICloneable, ISparseMatrix {

        /// <summary>
        /// ctor. The size of the block-diagonals is determined by the block
        /// size (see <see cref="Partitioning.BlockSize"/>) of the row and
        /// column mapping (<paramref name="_RowPartition"/>,
        /// <paramref name="ColPart"/>).
        /// </summary>
        public BlockDiagonalMatrix(IPartitioning _RowPartition, IPartitioning ColPart, int RowBlkSz, int ColBlkSz) {
            if (_RowPartition.MPI_Comm != ColPart.MPI_Comm)
                throw new ArgumentException();
            m_RowPart = _RowPartition;
            m_ColumnPart = ColPart;
            NoOfRowsPerBlock = RowBlkSz;
            NoOfColsPerBlock = ColBlkSz;

            if (_RowPartition.LocalLength % NoOfRowsPerBlock != 0)
                throw new ArgumentException();
            if (ColPart.LocalLength % NoOfRowsPerBlock != 0)
                throw new ArgumentException();


            Blox = new double[m_RowPart.LocalLength / NoOfRowsPerBlock, NoOfRowsPerBlock, NoOfColsPerBlock];
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
        /// used by <see cref="Clone"/>
        /// </summary>
        private BlockDiagonalMatrix() {
        }


        /// <summary>
        /// creates a non-shallow copy
        /// </summary>
        public object Clone() {
            var ret = new BlockDiagonalMatrix();
            ret.Blox = this.Blox.CloneAs();
            ret.m_RowPart = this.RowPartitioning; // immutable object, so shallow cloning is ok.
            ret.m_ColumnPart = this.ColPartition;
            return ret;
        }


        /// <summary>
        /// Number of rows in the the diagonal blocks; note that the blocks are
        /// not necessarily quadratic.
        /// </summary>
        public int NoOfRowsPerBlock {
            get;
            private set;
        }

        /// <summary>
        /// Number of columns in the the diagonal blocks; note that the blocks
        /// are not necessarily quadratic.
        /// </summary>
        public int NoOfColsPerBlock {
            get;
            private set;
        }

        private IPartitioning m_RowPart;

        /// <summary>
        /// Main storage of this object
        /// <list type="bullet">
        ///     <item>1st index: Block index</item>
        ///     <item>2nd index: Row within block</item>
        ///     <item>3rd index: Column within block</item>
        /// </list>
        /// </summary>
        private double[, ,] Blox;

        /// <summary>
        /// ctor, creates a matrix with quadratic blocks.
        /// </summary>
        public BlockDiagonalMatrix(IPartitioning _RowPartition, int BlockSz)
            : this(_RowPartition, _RowPartition, BlockSz, BlockSz) {
        }

        /// <summary>
        /// ctor, creates a matrix with quadratic blocks.
        /// </summary>
        public BlockDiagonalMatrix(int totalNoOfRows, int BlockSize)
            : this(new Partitioning(totalNoOfRows, MPI.Wrappers.csMPI.Raw._COMM.WORLD), BlockSize) {
        }

        /// <summary>
        /// Ctor, creates a matrix initialized with the diagonal blocks of
        /// <paramref name="Src"/>.
        /// </summary>        
        /// <param name="Src"></param>
        public BlockDiagonalMatrix(MsrMatrix Src, int RowBlkSz, int ColBlkSz)
            : this(Src.RowPartitioning, Src.ColPartition, RowBlkSz, ColBlkSz) //
        {
            int i0Row = Src.RowPartitioning.i0;
            int i0Col = Src.ColPartition.i0;

            for (int BlockNo = 0; BlockNo < NoOfBlocks; BlockNo++) {
                for (int RowBlock = 0; RowBlock < this.NoOfRowsPerBlock; RowBlock++) {
                    for (int ColBlock = 0; ColBlock < this.NoOfColsPerBlock; ColBlock++) {
                        int i, j;
                        GetIndex(BlockNo, RowBlock, ColBlock, out i, out j);
                        Blox[BlockNo, RowBlock, ColBlock] = Src[i0Row + i, i0Col + j];
                    }
                }
            }
        }

        /// <summary>
        /// Determines the indices (<paramref name="i"/>, <paramref name="j"/>)
        /// into the matrix for a given position (<paramref name="block"/>,
        /// <paramref name="rowInBlock"/>, <paramref name="columnInBlock"/>).
        /// </summary>
        /// <param name="block"></param>
        /// <param name="rowInBlock"></param>
        /// <param name="columnInBlock"></param>
        /// <param name="i"></param>
        /// <param name="j"></param>
        private void GetIndex(int block, int rowInBlock, int columnInBlock, out int i, out int j) {
            i = block * this.NoOfRowsPerBlock + rowInBlock;
            j = block * this.NoOfColsPerBlock + columnInBlock;
        }

        /// <summary>
        /// The total number of diagonal blocks.
        /// </summary>
        public int NoOfBlocks {
            get {
                return Blox.GetLength(0);
            }
        }

#pragma warning disable 1591

        #region IMutableMatrixEx Members

        public int GetOccupiedColumnIndices(int RowIndex, ref int[] ColumnIndices) {
        //public int[] GetOccupiedColumnIndices(int RowIndex) {
            int _NoOfColsPerBlock = NoOfRowsPerBlock;

            if (ColumnIndices == null || ColumnIndices.Length < _NoOfColsPerBlock)
                ColumnIndices = new int[_NoOfColsPerBlock];

            TransformRow(RowIndex); // tests the range of 'RowIndex'
            int iRowGlob = RowIndex / NoOfRowsPerBlock;

            int j0 = iRowGlob * _NoOfColsPerBlock;
            for (int j = 0; j < _NoOfColsPerBlock; j++)
                ColumnIndices[j] = j0 + j;

            return _NoOfColsPerBlock;
        }

        public int GetRow(int RowIndex, ref int[] ColumnIndices, ref double[] Values) {
            //public MsrMatrix.MatrixEntry[] GetRow(int RowIndex) {
            //MsrMatrix.MatrixEntry[] ret = new MsrMatrix.MatrixEntry[NoOfColsPerBlock];

            TransformRow(RowIndex); // tests the range of 'RowIndex'
            int iRowGlob = RowIndex / NoOfRowsPerBlock;
            int j0 = iRowGlob * NoOfColsPerBlock;

            if (ColumnIndices == null || ColumnIndices.Length < NoOfColsPerBlock)
                ColumnIndices = new int[NoOfColsPerBlock];
            if (Values == null || Values.Length < NoOfColsPerBlock)
                Values = new double[NoOfColsPerBlock];

            int iBlk, sub_i0, sub_j0;
            TransformIndexToBlockIndex(RowIndex, j0, out iBlk, out sub_i0, out sub_j0);

            for (int j = 0; j < NoOfColsPerBlock; j++) {
                Values[j] = Blox[iBlk, sub_i0, sub_j0 + j];
                ColumnIndices[j] = j0 + j;
            }

            return NoOfColsPerBlock;
        }

        #endregion

        #region IMutableMatrix Members

        public double[] GetValues(int RowIndex, int[] ColumnIndices) {
            double[] ret = new double[ColumnIndices.Length];
            for (int i = 0; i < ret.Length; i++)
                ret[i] = this[RowIndex, ColumnIndices[i]];
            return ret;
        }

        public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
            if (ColumnIndices.Length != newValues.Length)
                throw new ArgumentException();
            for (int i = 0; i < ColumnIndices.Length; i++)
                this[RowIndex, ColumnIndices[i]] = newValues[i];
        }

        private int TransformRow(int i) {
            int r = i;
            r -= (int)m_RowPart.i0;
            if (r < 0 || r >= m_RowPart.LocalLength)
                throw new ArgumentException("row index is out of range of current mpi process.");
            return r;
        }

        private bool TransformIndexToBlockIndex(int i, int j, out int BlockInd, out int RowWithinBlk, out int ColWithinBlk) {
            int iLoc = TransformRow(i);

            int BlockRowGlob = i / NoOfRowsPerBlock;
            int BlockColGlob = j / NoOfRowsPerBlock;
            if (BlockColGlob != BlockRowGlob) {
                // (i,j) - pair OUTSIDE of a diagonal block
                // ++++++++++++++++++++++++++++++++++++++++

                BlockInd = int.MinValue;
                RowWithinBlk = int.MinValue;
                ColWithinBlk = int.MinValue;
                return false;
            } else {
                // (i,j) - pair INSIDE of a diagonal block
                // +++++++++++++++++++++++++++++++++++++++

                BlockInd = iLoc / NoOfRowsPerBlock;
                RowWithinBlk = iLoc % NoOfRowsPerBlock;
                ColWithinBlk = j % NoOfRowsPerBlock;
                return true;
            }
        }

        public double this[int i, int j] {
            get {
                int a, b, c;
                bool res = TransformIndexToBlockIndex(i, j, out a, out b, out c);
                if (!res)
                    return 0.0;
                else
                    return Blox[a, b, c];
            }
            set {
                int a, b, c;
                bool res = TransformIndexToBlockIndex(i, j, out a, out b, out c);
                if (!res)
                    throw new IndexOutOfRangeException("Index outside diagonal block: unable to set element, no memory allocated for given row/column index.");
                Blox[a, b, c] = value;
            }
        }

        /// <summary>
        /// Accumulates a block of entries to this matrix.
        /// </summary>
        /// <param name="i0">Row index offset.</param>
        /// <param name="j0">Column index offset.</param>
        /// <param name="alpha">Scaling factor for the accumulation.</param>
        /// <param name="Block">Block to add.</param>
        public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block) {
            if (Block.Dimension != 2)
                throw new ArgumentException();
            int I = Block.NoOfRows;
            int J = Block.NoOfCols;

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    this[i0 + i, j0 + j] += alpha * Block[i, j];
        }

        /// <summary>
        /// Provides direct access to the matrix using the block index
        /// <paramref name="block"/> and the row and column index <b>within</b>
        /// this block.
        /// </summary>
        /// <param name="block">
        /// The selected block
        /// </param>
        /// <param name="blockRow">
        /// A row within the selected <paramref name="block"/>
        /// </param>
        /// <param name="blockCol">
        /// A column within the selected <paramref name="block"/>
        /// </param>
        /// <returns>
        /// The stored value
        /// </returns>
        public double this[int block, int blockRow, int blockCol] {
            get {
                return Blox[block, blockRow, blockCol];
            }
            set {
                Blox[block, blockRow, blockCol] = value;
            }
        }

        public bool OccupationMutable {
            get {
                return false;
            }
        }

        #endregion

        #region ISparseMatrix Members

        public double GetDiagonalElement(int row) {
            return this[row, row];
        }

        public void SetDiagonalElement(int row, double val) {
            this[row, row] = val;
        }

        public IPartitioning RowPartitioning {
            get {
                return m_RowPart;
            }
        }

        /// <summary>
        /// matrix/vector product:
        /// <paramref name="acc"/>=<paramref name="acc"/>*<paramref name="beta"/>
        /// + <paramref name="alpha"/>*this*<paramref name="a"/>
        /// </summary>
        public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {

            if (acc.Count != m_RowPart.LocalLength)
                throw new ArgumentException("mismatch between accumulator and number of rows, i.e. row partition.", "acc");
            if (a.Count != m_ColumnPart.LocalLength)
                throw new ArgumentException("mismatch between input vector 'a' and number of columns, i.e. column partition.", "a");


            int M1 = NoOfColsPerBlock;
            int N1 = NoOfRowsPerBlock;
            double[] aBuf = new double[M1];

            for (int cnt = 0; cnt < NoOfBlocks; cnt++) {
                int i0 = cnt * N1;
                int j0 = cnt * M1;

                for (int j = 0; j < M1; j++)
                    aBuf[j] = a[j + j0];

                for (int i = 0; i < N1; i++) {
                    double _acc = 0;

                    for (int j = 0; j < M1; j++)
                        _acc += Blox[cnt, i, j] * aBuf[j];


                    acc[i + i0] = acc[i + i0] * beta + _acc * alpha;
                }
            }
        }

        #endregion

        private IPartitioning m_ColumnPart;

        /// <summary>
        /// column partition that is induced by the row partition
        /// </summary>
        public IPartitioning ColPartition {
            get {
                return m_ColumnPart;
            }
        }

        /// <summary>
        /// matrix-matrix-multiplication, see <see cref="Multiply"/>.
        /// </summary>
        static public MsrMatrix operator *(BlockDiagonalMatrix left, MsrMatrix right) {
            return Multiply(left, right);
        }

        /// <summary>
        /// multiplies an <see cref="BlockDiagonalMatrix"/> by an <see cref="MsrMatrix"/>
        /// </summary>
        public static MsrMatrix Multiply(BlockDiagonalMatrix left, MsrMatrix right) {
            using (new FuncTrace()) {
                if (left.NoOfCols != right.NoOfRows)
                    throw new ArgumentException("matrix size mismatch");

                int _N = left.NoOfRowsPerBlock;   // left blocks: no. of rows
                int _M = left.NoOfColsPerBlock;  // left blocks: no. of columns = right blocks: no. of rows
                int _L = _M; // right blocks: no. of columns. Hope that fits, otherwise inefficient

                MsrMatrix result = new MsrMatrix(left.RowPartitioning, right.ColPartition);

                MultidimensionalArray lBlock = MultidimensionalArray.Create(_N, _M);
                //MsrMatrix.MSREntry[][] BlockRow = new MsrMatrix.MSREntry[_M][];
                //List<int> OccupiedBlocks = new List<int>();
                var OccupiedBlocks = new SortedDictionary<int, MultidimensionalArray>();
                List<MultidimensionalArray> MatrixPool = new List<MultidimensionalArray>(); // avoid realloc
                int poolCnt = 0;

                MultidimensionalArray resBlock = MultidimensionalArray.Create(_N, _L);

                // loop over all blocks ...
                int NoBlks = left.m_RowPart.LocalLength / left.NoOfRowsPerBlock;
                for (int iBlk = 0; iBlk < NoBlks; iBlk++) {

                    // reset:
                    OccupiedBlocks.Clear();
                    poolCnt = 0;
                    foreach (var M in MatrixPool)
                        M.Clear();

                    // upper left corner of block
                    int i0 = iBlk * _N + (int)left.RowPartitioning.i0;
                    int j0 = (i0 / _N) * _M;

                    // load block
                    for (int n = 0; n < _N; n++)
                        for (int m = 0; m < _M; m++)
                            lBlock[n, m] = left.Blox[iBlk, n, m];

                    // read Msr rows
                    int Last_rBlkCol = int.MinValue;
                    MultidimensionalArray Block = null;
                    for (int m = 0; m < _N; m++) {
                        var BlockRow = right.GetRowShallow(j0 + m);
                        foreach (var e in BlockRow) {
                            if (e.m_ColIndex >= 0) {
                                int rBlkCol = e.m_ColIndex / _L;

                                if (rBlkCol != Last_rBlkCol) {

                                    if (!OccupiedBlocks.TryGetValue(rBlkCol, out Block)) {
                                        if (MatrixPool.Count <= poolCnt)
                                            MatrixPool.Add(MultidimensionalArray.Create(_M, _L));

                                        Block = MatrixPool[poolCnt];
                                        OccupiedBlocks.Add(rBlkCol, Block);
                                        poolCnt++;
                                    }

                                    Last_rBlkCol = rBlkCol;
                                }

                                int jj = e.m_ColIndex % _L;

                                Block[m, jj] = e.Value;
                            }
                        }
                    }
                    Block = null;

                    // execute multiplys
                    foreach (KeyValuePair<int, MultidimensionalArray> rBlock in OccupiedBlocks) {
                        resBlock.GEMM(1.0, lBlock, rBlock.Value, 0.0);

                        // upper left edge of resBlock
                        int _i0 = i0;
                        int _j0 = rBlock.Key * _L;

                        // write block
                        for (int i = 0; i < _M; i++)
                            for (int j = 0; j < _L; j++)
                                result[i + _i0, j + _j0] = resBlock[i, j];
                    }

                }


                // return
                return result;
            }
        }

        /// <summary>
        /// Calculates the inverse of this matrix and returns the result.
        /// </summary>
        /// <returns>Inverse of this matrix</returns>
        public BlockDiagonalMatrix Invert() {
            // Some checks
            if (this.NoOfRows != this.NoOfCols)
                throw new NotSupportedException("Can't invert non-square matrix.");

            if (this.NoOfRowsPerBlock != this.NoOfColsPerBlock)
                throw new NotSupportedException("Can't invert matrix with non-square blocks.");

            // Matrix for the result
            BlockDiagonalMatrix res = new BlockDiagonalMatrix(RowPartitioning, ColPartition, NoOfRowsPerBlock, NoOfColsPerBlock);

            // FullMatrix to calculate inverse
            for (int BlockNo = 0; BlockNo < NoOfBlocks; BlockNo++) {

                var tmp = MultidimensionalArray.Create(this.NoOfRowsPerBlock, NoOfColsPerBlock);
                // Copy values from current block to FullMatrix
                for (int i = 0; i < this.NoOfRowsPerBlock; i++)
                    for (int j = 0; j < this.NoOfColsPerBlock; j++)
                        tmp[i, j] = this.Blox[BlockNo, i, j];

                // Calculate inverse of FullMatrix
                var BlockInverse = tmp.GetInverse();

                // Write inverse of FullMatrix to current block
                for (int i = 0; i < this.NoOfRowsPerBlock; i++)
                    for (int j = 0; j < this.NoOfColsPerBlock; j++)
                        res.Blox[BlockNo, i, j] = BlockInverse[i, j];
            }

            return res;
        }

        public BlockDiagonalMatrix InvertSymmetrical(bool ignoreEmptyBlocks = false) {
            if (this.NoOfRows != this.NoOfCols) {
                throw new NotSupportedException("Can't invert non-square matrix.");
            }

            if (this.NoOfRowsPerBlock != this.NoOfColsPerBlock) {
                throw new NotSupportedException("Can't invert matrix with non-square blocks.");
            }

            // Matrix for the result
            BlockDiagonalMatrix res = new BlockDiagonalMatrix(RowPartitioning, ColPartition, NoOfRowsPerBlock, NoOfColsPerBlock);

            // FullMatrix to calculate inverse
            var tmp = MultidimensionalArray.Create(this.NoOfRowsPerBlock, NoOfColsPerBlock);
            for (int BlockNo = 0; BlockNo < NoOfBlocks; BlockNo++) {

                // Copy values from current block to FullMatrix
                for (int i = 0; i < this.NoOfRowsPerBlock; i++) {
                    for (int j = 0; j < this.NoOfColsPerBlock; j++) {
                        tmp[i, j] = this.Blox[BlockNo, i, j];
                    }
                }

                if (ignoreEmptyBlocks && tmp.AbsSum() < 1e-15) {
                    // No inversion of empty blocks requested
                    continue;
                }

                // Calculate inverse of FullMatrix
                tmp.InvertSymmetrical();

                // Write inverse of FullMatrix to current block
                for (int i = 0; i < this.NoOfRowsPerBlock; i++) {
                    for (int j = 0; j < this.NoOfColsPerBlock; j++) {
                        res.Blox[BlockNo, i, j] = tmp[i, j];
                    }
                }
            }

            return res;
        }

#pragma warning restore 1591



        /// <summary>
        /// number of matrix rows
        /// </summary>
        public int NoOfRows {
            get {
                return (int)m_RowPart.TotalLength;
            }
        }

        /// <summary>
        /// number of matrix columns
        /// </summary>
        public int NoOfCols {
            get {
                return (int)ColPartition.TotalLength;
            }
        }

        /// <summary>
        /// Sets all entries to zero
        /// </summary>
        public void Clear() {
            Array.Clear(Blox, 0, Blox.Length);
        }

        /// <summary>
        /// Scales each entry by <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        public void Scale(double a) {
            for (int i = 0; i < Blox.GetLength(0); i++) {
                for (int j = 0; j < Blox.GetLength(1); j++) {
                    for (int k = 0; k < Blox.GetLength(2); k++) {
                        Blox[i, j, k] *= a;
                    }
                }
            }
        }

        /// <summary>
        /// Accumulates <paramref name="a"/>*<paramref name="M"/> to this
        /// matrix. Note that 
        /// </summary>
        /// <param name="a"></param>
        /// <param name="M">
        /// The matrix to be added. Note that the given matrix must not contain
        /// entries outside of the diagonal blocks represented by this matrix.
        /// </param>
        public void Acc(double a, IMatrix M) {
            int block;
            int rowInBlock;
            int columnIndBlock;
            for (int i = 0; i < M.NoOfRows; i++) {
                for (int j = 0; j < M.NoOfCols; j++) {
                    bool blockExists = this.TransformIndexToBlockIndex(
                        i, j, out block, out rowInBlock, out columnIndBlock);
                    if (blockExists) {
                        throw new ArgumentException(
                            "Cannot add matrix with off-block-diagonal entries to a block-diagonal matrix",
                            "M");
                    }

                    Blox[block, rowInBlock, columnIndBlock] += a * M[i, j];
                }
            }
        }

        /// <summary>
        /// Copies the entries of this matrix to <paramref name="dest"/>.
        /// </summary>
        /// <param name="dest"></param>
        /// <param name="i0"></param>
        /// <param name="i1"></param>
        public void CopyTo(double[,] dest, int i0, int i1) {
            for (int block = 0; block < NoOfBlocks; block++) {
                for (int rowInBlock = 0; rowInBlock < NoOfRowsPerBlock; rowInBlock++) {
                    for (int columnInBlock = 0; columnInBlock < NoOfColsPerBlock; columnInBlock++) {
                        int i, j;
                        GetIndex(block, rowInBlock, columnInBlock, out i, out j);
                        dest[i0 + i, i1 + j] = Blox[block, rowInBlock, columnInBlock];
                    }
                }
            }
        }

        /// <summary>
        /// Not implemented (not used in the code anyway?)
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="array"></param>
        /// <param name="RowWise"></param>
        /// <param name="arrayoffset"></param>
        public void CopyTo<T>(T array, bool RowWise, int arrayoffset) where T : IList<double> {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Creates a copy of block <paramref name="blockIndex"/>
        /// </summary>
        /// <param name="blockIndex">
        /// Index of the block
        /// </param>
        /// <returns>
        /// A copy of the selected block as a matrix
        /// </returns>
        public IMatrix GetBlock(int blockIndex) {
            MultidimensionalArray block = MultidimensionalArray.Create(NoOfRowsPerBlock, NoOfColsPerBlock);
            for (int i = 0; i < NoOfRowsPerBlock; i++) {
                for (int j = 0; j < NoOfColsPerBlock; j++) {
                    block[i, j] = Blox[blockIndex, i, j];
                }
            }
            return block;
        }

        /// <summary>
        /// Saves the contents of this matrix into a text file for debugging
        /// purposes.
        /// </summary>
        /// <param name="fileName"></param>
        /// <param name="fm"></param>
        public void SaveToTextFile(string fileName, FileMode fm = FileMode.Create) {
            MultidimensionalArray temp = MultidimensionalArray.Create(NoOfRows, NoOfCols);
            for (int block = 0; block < NoOfBlocks; block++) {
                for (int i = 0; i < NoOfRowsPerBlock; i++) {
                    for (int j = 0; j < NoOfColsPerBlock; j++) {
                        int row, col;
                        GetIndex(block, i, j, out row, out col);
                        temp[row, col] = this[block, i, j];
                    }
                }
            }
            temp.SaveToTextFile(fileName, fm);
        }
    }
    
}
