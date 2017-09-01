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
using System.Text;
using System.Runtime.InteropServices;
using System.Collections;
using ilPSP.Utils;
using ilPSP.Tracing;

namespace ilPSP.LinSolvers.monkey {
    partial class MatrixBase {


        /// <summary>
        /// used to collect all info from the <see cref="MsrMatrix"/> into a temporary CRS structure
        /// </summary>
        public class TempCSR {
            /// <summary>
            /// ctor
            /// </summary>
            public TempCSR() {
                RowStart.Add(0);
            }
            
            /// <summary>
            /// matrix values
            /// </summary>
            public List<double> Vals = new List<double>();
            
            /// <summary>
            /// (local, i.e. within the current MPI process) column indices of matrix
            /// columns
            /// </summary>
            public List<int> ColInd = new List<int>();
            
            /// <summary>
            /// index: (local, i.e. within the current MPI process) row index <em>r</em> <br/>
            /// content: index into <see cref="Vals"/> and <see cref="ColInd"/>,
            /// at which row <em>r</em> starts.
            /// </summary>
            public List<int> RowStart = new List<int>();
            
            /// <summary>
            /// index: (local, i.e. within the current MPI process) row index <br/>
            /// content: an index into <see cref="Vals"/>, denoting the position of the 
            /// diagonal element; if negative, the diagonal element is not contained in the internal
            /// part of the matrix (but maybe in the external one).
            /// </summary>
            public List<int> DiagIndex = new List<int>();

            int m_RowInd = 0;
            int m_ColInd = -1;

            int m_Cnt = 0;

            bool diagAssigned = false;

            /// <summary>
            /// total number of locally (i.e. within the current MPI process)
            /// stored rows.
            /// </summary>
            public int NoOfRows { get { return RowStart.Count - 1; } }

            /// <summary>
            /// adds the next nonzero entry to the CSR structure
            /// </summary>
            /// <param name="__ColInd"></param>
            /// <param name="val"></param>
            /// <param name="bIsDiag">
            /// true, if this entry is a diagonal entry;
            /// When running on more than one MPI process, this cannot be determined;
            /// </param>
            public void AddEntry(int __ColInd, double val, bool bIsDiag) {
                if (__ColInd <= m_ColInd)
                    throw new ArgumentException("rows must be specified from left to right (ascending column index).");
                m_ColInd = __ColInd;
                
                Vals.Add(val);
                ColInd.Add(__ColInd);

                if (bIsDiag) {
                    if (diagAssigned)
                        throw new ApplicationException("internal error: two diag elements.");
                    DiagIndex.Add(m_Cnt);
                    diagAssigned = true;
                }

                m_Cnt++;
            }

            /// <summary>
            /// must be called before moving to the next row
            /// </summary>
            public void NextRow() {
                if (!diagAssigned)
                    DiagIndex.Add(int.MinValue);
                m_RowInd++;
                RowStart.Add(m_Cnt);
                m_ColInd = -1;
                diagAssigned = false;
            }
        }
        
        /// <summary>
        /// Common baseclass for all sparse matix data formats
        /// </summary>
        abstract public class FormatBase {

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="_csr">source object</param>
            protected FormatBase(TempCSR _csr) {
                csr = _csr;
            }

            /// <summary>
            /// the source data object
            /// </summary>
            protected TempCSR csr;
            
            /// <summary>
            /// matrix entries.
            /// </summary>
            public double[] Val;
            
            /// <summary>
            /// computes the storage position (i.e. index into <see cref="Val"/>)
            /// for the matrix entry (<paramref name="row"/>,<paramref name="col"/>);
            /// </summary>
            /// <param name="row"></param>
            /// <param name="col"></param>
            /// <returns>
            /// An index into <see cref="Val"/>, that corresponds with matrix entry (<paramref name="row"/>,<paramref name="col"/>);
            /// If this entry is not set (i.e. it is zero), a negative value is returned;
            /// </returns>
            /// <remarks>
            /// the purpose of this function is to support the <see cref="IMutableMatrix"/>-methods;
            /// </remarks>
            abstract public int GetEntryIndex(int row, int col);

            /// <summary>
            /// provides information about columns
            /// </summary>
            /// <param name="row"></param>
            /// <param name="Out_MtxColIndices"></param>
            /// <param name="Out_PointersIntoVal"></param>
            /// <param name="Out_Values"></param>
            /// <param name="UsedLen"></param>
            abstract public void GetAllOccupiedColumns(int row, ref int[] Out_MtxColIndices, ref int[] Out_PointersIntoVal, ref double[] Out_Values, out int UsedLen);

            /// <summary>
            /// performs, for reference and testing purposes, the operation <br/>
            /// <paramref name="y"/> = <paramref name="y"/>*<paramref name="beta"/> + <paramref name="alpha"/>*this*<paramref name="x"/>
            /// </summary>
            abstract public void RefSpMv(double alpha, double[] x, double beta, double[] y);
        }

        /// <summary>
        /// Compressed Sparse Row format;
        /// </summary>
        public class CSR : FormatBase {

            /// <summary>
            /// ctor
            /// </summary>
            public CSR(TempCSR _csr)
                : base(_csr) {

                this.RowStart = _csr.RowStart.ToArray();
                //this.DiagIndex = _csr.DiagIndex.ToArray();
                this.ColInd = _csr.ColInd.ToArray();
                this.Val = _csr.Vals.ToArray();
            }
                        
            /// <summary>
            /// local column indices
            /// </summary>
            public int[] ColInd;
            
            /// <summary>
            /// indices/points into <see cref="FormatBase.Val"/> and <see cref="ColInd"/>, 
            /// denoting where matrix rows are starting.<br/>
            /// index: row index
            /// </summary>
            public int[] RowStart;

            /// <summary>
            /// total number of rows
            /// </summary>
            public int NoOfRows {
                get { return RowStart.Length - 1; }
            }

            /// <summary>
            /// <see cref="FormatBase.RefSpMv"/>
            /// </summary>
            public override void RefSpMv(double alpha, double[] x, double beta, double[] y) {
                if (y.Length != NoOfRows)
                    throw new ArgumentException("wrong length", "y");


                for (int Row = 0; Row < NoOfRows; Row++) {
                    int stRow = this.RowStart[Row];
                    int enRow = this.RowStart[Row + 1];

                    double rowacc = 0;
                    for (int i = stRow; i < enRow; i++) {
                        int Col = this.ColInd[i];
                        double MtxEntry = this.Val[i];

                        rowacc += MtxEntry * x[Col];
                    }

                    y[Row] = y[Row] * beta + rowacc * alpha;
                }
            }

            /// <summary>
            /// see <see cref="FormatBase.GetEntryIndex"/>
            /// </summary>
            public override int GetEntryIndex(int row, int col) {
                int stRow = this.RowStart[row];
                int enRow = this.RowStart[row + 1];

                for (int c = stRow; c < enRow; c++) {
                    if (ColInd[c] == col)
                        return c; // found entry
                }
                return int.MinValue; // did not found entry
            }

            /// <summary>
            /// see <see cref="FormatBase.GetAllOccupiedColumns"/>;
            /// </summary>
            public override void GetAllOccupiedColumns(int row, ref int[] Out_MtxColIndices, ref int[] Out_PointersIntoVal, ref double[] Out_Values, out int UsedLen) {
                int stRow = this.RowStart[row];
                int enRow = this.RowStart[row + 1];
                int L = enRow-stRow;
                UsedLen = L;

                if (Out_MtxColIndices == null || Out_MtxColIndices.Length < L) { Out_MtxColIndices = new int[L]; }
                if (Out_PointersIntoVal == null || Out_PointersIntoVal.Length < L) { Out_PointersIntoVal = new int[L]; }
                if (Out_Values == null || Out_Values.Length < L) { Out_Values = new double[L]; };
                double[] Values = Out_Values;
                int[] MtxColIndices = Out_MtxColIndices;
                int[] PointersIntoVal = Out_PointersIntoVal;
                

                int i = 0;
                for (int c = stRow; c < enRow; c++) {
                    PointersIntoVal[i] = c;
                    MtxColIndices[i] = ColInd[c];
                    Values[i] = Val[PointersIntoVal[i]];
                    i++;
                }
            }

            /*
            /// <summary>
            /// see <see cref="FormatBase.GetRowStat"/>
            /// </summary>
            public override void GetRowStat(int row, out int NoOfNonZeros, out int MinCol, out int MaxCol) {
                int stRow = this.RowStart[row];
                int enRow = this.RowStart[row + 1];

                NoOfNonZeros = enRow - stRow;

                MinCol = int.MaxValue;
                MaxCol = int.MinValue;
                for (int c = stRow; c < enRow; c++) {
                    MinCol = Math.Min(ColInd[c], MinCol);
                    MaxCol = Math.Max(ColInd[c], MaxCol);
                }
            }
             */
        }

        /// <summary>
        /// Baseclass for cell-based formats; By a cell, we refer to a 
        /// <see cref="CellSize"/>x<see cref="CellSize"/> - submatrix (of the complete matrix);
        /// these 'blocks' are refered to as cells to avoid confusion with CUDA-(thread)-blocks.
        /// </summary>
        /// <remarks>
        /// Index computation: the (i,j)-th entry of the k-th cell
        /// is computed by k*<see cref="CellStride"/> + j*<see cref="ColStride"/> + i.
        /// </remarks>
        abstract public class CelledFormats : FormatBase {

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="_csr"></param>
            /// <param name="AlignmentInBytes">
            /// memory aligmnet; the first entry of each cell will lie on an address (relative tor the address of the first entry)
            /// that is divideable by this value;
            /// </param>
            /// <param name="_CellSize"></param>
            protected CelledFormats(TempCSR _csr, int AlignmentInBytes, int _CellSize)
                : base(_csr) {

                CellSize = _CellSize;
                if (csr.NoOfRows % _CellSize != 0)
                    throw new ArgumentException("blocksize must be a divider of number of rows");
                NoOfCellRows = csr.NoOfRows / _CellSize;

                ColStride = _CellSize;
                CellStride = _CellSize * _CellSize;
                int CellStrideInBytes = CellStride * 8;
                while (CellStrideInBytes % AlignmentInBytes != 0) {
                    CellStrideInBytes += 8;
                    CellStride++;
                }
                Console.Write("");
            }

            /// <summary>
            /// finds all occuopied cell-columns in a given cell-row
            /// </summary>
            /// <param name="CellRow"></param>
            /// <returns>
            /// collection of cell - columns, which contain the non-empty bocks
            /// </returns>
            protected int[] FindCells(int CellRow) {

                if (CellRow < 0 || CellRow >= csr.NoOfRows / CellSize)
                    throw new IndexOutOfRangeException("matrix cell index out of range.");

                int RowSt = CellRow * CellSize;
                int RowEnd = RowSt + CellSize;

                List<int> CellColColl = new List<int>();
                for (int i = RowSt; i < RowEnd; i++) {

                    int c0 = csr.RowStart[i];
                    int c1 = csr.RowStart[i + 1];

                    for (int j = c0; j < c1; j++) {
                        int col = csr.ColInd[j];

                        int CellCol = col / CellSize;
                        if (!CellColColl.Contains(CellCol))
                            CellColColl.Add(CellCol);
                    }
                }

                return CellColColl.ToArray();
            }

            /// <summary>
            /// computes the total number of occupied cells in cell row <paramref name="CellRow"/>;
            /// </summary>
            protected int Compute_NoOfCellsPerCellRow(int CellRow) {
                return FindCells(CellRow).Length;
            }

            /// <summary>
            /// computes the total number of occupied cells in the whole matrix
            /// </summary>
            protected int Compute_TotalNoOfCells() {
                int Sum = 0;
                int L = csr.NoOfRows/CellSize;
                for (int i = 0; i < L; i++)
                    Sum += Compute_NoOfCellsPerCellRow(i);
                return Sum;
            }

            /// <summary>
            /// guess what?
            /// </summary>
            protected int Compute_MaxNoOfCellsPerRow() {
                int MAx = 0;
                int L = csr.NoOfRows / CellSize;
                for (int i = 0; i < L; i++)
                    MAx = Math.Max(MAx,Compute_NoOfCellsPerCellRow(i));
                return MAx;
            }

            /// <summary>
            /// defines the indexing of the cells.
            /// </summary>
            /// <param name="GlobalCellIndex">global cell index (cells are counted row-wise)</param>
            /// <param name="ColWithinCell">column index, locally within cell</param>
            /// <param name="RowWithinCell">row index, locally within cell</param>
            /// <returns>an index into <see cref="FormatBase.Val"/></returns>
            protected int GetIndexIntoVal(int GlobalCellIndex, int RowWithinCell, int ColWithinCell) {
                if( ColWithinCell < 0 || ColWithinCell > CellSize)
                    throw new IndexOutOfRangeException();
                if( RowWithinCell < 0 || RowWithinCell > CellSize)
                    throw new IndexOutOfRangeException();

                return GlobalCellIndex * CellStride + ColWithinCell * ColStride + RowWithinCell;
            }

            /// <summary>
            /// loads the content of one cell from the <see cref="TempCSR"/>-structure (see <see cref="FormatBase.csr"/>),
            /// into the <see cref="FormatBase.Val"/>-array.
            /// </summary>
            /// <param name="CellRow"></param>
            /// <param name="CellCol"></param>
            /// <param name="GlobalCellIdx"></param>
            protected void LoadCell(int CellRow, int CellCol, int GlobalCellIdx) {
                if ((GlobalCellIdx + 1) * CellStride > Val.Length)
                    throw new IndexOutOfRangeException("global block index");
                if (CellRow < 0 || CellRow >= csr.NoOfRows / CellSize)
                    throw new IndexOutOfRangeException("CellRow out of range");
                if (CellCol < 0)
                    throw new IndexOutOfRangeException("CellCol out of range");
                
                //int i0 = GlobalCellIdx * CellStride;

                int rowSt = CellRow * CellSize;
                int rowEn = rowSt + CellSize;

                int colSt = CellCol * CellSize;
                int colEn = colSt + CellSize;

                for (int i = rowSt; i < rowEn; i++) { // loop over row
                    int c0 = csr.RowStart[i];
                    int c1 = csr.RowStart[i + 1];

                    int DiagCol = csr.DiagIndex[i];

                    for (int j = c0; j < c1; j++) { // loop over column
                        if (csr.ColInd[j] >= colSt && csr.ColInd[j] < colEn) {
                            int iCell = i % CellSize;
                            int jCell = csr.ColInd[j] % CellSize;

                            int iii = GetIndexIntoVal(GlobalCellIdx, iCell, jCell);

                            Val[iii] = csr.Vals[j];

                            if (DiagCol == csr.ColInd[j])
                                csr.DiagIndex[i] = iii;
                        }
                    }
                }
            }

            /// <summary>
            /// defines the global cell index for the firs cell in cell-row <paramref name="CellRowInd"/>
            /// </summary>
            abstract protected int GetCellRowStart(int CellRowInd);

            /// <summary>
            /// <see cref="FormatBase.RefSpMv"/>
            /// </summary>
            public override void RefSpMv(double alpha, double[] x, double beta, double[] y) {
                if (y.Length != NoOfCellRows * CellSize)
                    throw new ArgumentException("wrong length of y", "y");

                for (int i = 0; i < y.Length; i++)
                    y[i] *= beta;

                int CellCnt = 0; // global cell counter

                // loop over all cell rows
                for (int i = 0; i < NoOfCellRows; i++) {
                    int cellSt = GetCellRowStart(i);   // first cell in this row
                    int cellEn = GetCellRowStart(i+1); // first cell in next row

                    int i0_y = i * CellSize; // offset into y-vector for the current cell-row

                    // loop over cells columns
                    for (int j = cellSt; j < cellEn; j++) {
                        int CelCol = CellColumn[j];

                        int j0_x = CelCol * CellSize; // offset into x-vector for the current cell-column

                        // loop over cell rows
                        for (int ii = 0; ii < CellSize; ii++) {
                            double acc = 0;
                            // loop over cell columns
                            for (int jj = 0; jj < CellSize; jj++) {
                                acc += x[j0_x + jj] * Val[GetIndexIntoVal(CellCnt, ii, jj)];
                            }
                            y[i0_y + ii] += acc * alpha;
                        }

                        CellCnt++;
                    }
                }
            }
            
            /// <summary>
            /// size of the cell (a cell is a <see cref="CellSize"/>x<see cref="CellSize"/> - submatrix);
            /// </summary>
            public int CellSize;

            /// <summary>
            /// stride between columns in cell, in 8-byte units (length of double)
            /// </summary>
            public int ColStride;

            /// <summary>
            /// stride between cells, in 8-byte units (length of double)
            /// </summary>
            public int CellStride; 

            /// <summary>
            /// for each cell, the cell column
            /// </summary>
            public int[] CellColumn;

            /// <summary>
            /// total number of cell columns, in all rows
            /// </summary>
            public int TotalNumberOfCells { get { return CellColumn.Length; } }

            /// <summary>
            /// number of cell-rows in the matrix
            /// </summary>
            public int NoOfCellRows;


            /// <summary>
            /// see <see cref="FormatBase.GetEntryIndex"/>
            /// </summary>
            public override int GetEntryIndex(int row, int col) {
                int cellRow = row / CellSize;

                int cellSt = GetCellRowStart(cellRow);     // first cell in cellRow
                int cellEn = GetCellRowStart(cellRow + 1); // first cell in next cell-Row
                                
                // loop over cells columns
                for (int jCell = cellSt; jCell < cellEn; jCell++) {
                    int CelCol = CellColumn[jCell];

                    // column range of cell (cellRow,CelCol)
                    int j0 = CelCol * CellSize;
                    int jE = (CelCol + 1) * CellSize;

                    if (col >= j0 && col < jE) {
                        // Das ist ein Bingo!
                        // Sagt man das so, 'das ist ein Bingo'?

                        // row&col indices within cell:
                        int ii = row - cellRow * CellSize;
                        int jj = col - j0;

                        return GetIndexIntoVal(jCell,ii,jj);
                    }
                }

                // not found
                return int.MinValue;
            }
            
            /// <summary>
            /// see <see cref="FormatBase.GetAllOccupiedColumns"/>;
            /// </summary>
            public override void GetAllOccupiedColumns(int row, ref int[] Out_MtxColIndices, ref int[] Out_PointersIntoVal, ref double[] Out_Values, out int UsedLen) {
                int cellRow = row / CellSize;
                int cellSt = GetCellRowStart(cellRow);     // first cell in cellRow
                int cellEn = GetCellRowStart(cellRow + 1); // first cell in next cell-Row

                int L = (cellEn - cellSt) * CellSize;
                UsedLen = L;

                if (Out_MtxColIndices == null || Out_MtxColIndices.Length < L) Out_MtxColIndices = new int[L];
                if (Out_PointersIntoVal == null || Out_PointersIntoVal.Length < L) Out_PointersIntoVal = new int[L]; 
                if (Out_Values == null || Out_Values.Length < L) { Out_Values = new double[L]; };
                int[] MtxColIndices = Out_MtxColIndices;
                int[] PointersIntoVal = Out_PointersIntoVal;

                // loop over cells columns
                int i = 0;
                for (int jCell = cellSt; jCell < cellEn; jCell++) {
                    int CelCol = CellColumn[jCell];

                    // column range of cell (cellRow,CelCol)
                    int j0 = CelCol * CellSize;
                    int jE = (CelCol + 1) * CellSize;

                    
                    for( int j = j0; j < jE; j++) { // loop over columns within cell
                        // row&col indices within cell:
                        int ii = row - cellRow * CellSize;
                        int jj = j - j0;

                        PointersIntoVal[i] = GetIndexIntoVal(jCell, ii, jj);
                        MtxColIndices[i] = j;
                        Out_Values[i] = Val[PointersIntoVal[i]];
                        i++;
                    }
                }
            }

            /*
            public override void GetRowStat(int row, out int NoOfNonZeros, out int MinCol, out int MaxCol) {
                int cellRow = row / CellSize;
                int cellSt = GetCellRowStart(cellRow);     // first cell in cellRow
                int cellEn = GetCellRowStart(cellRow + 1); // first cell in next cell-Row


                MinCol = int.MaxValue;
                MaxCol = int.MinValue;
                for (int jCell = cellSt; jCell < cellEn; jCell++) {
                    int CelCol = CellColumn[jCell];

                    // column range of cell (cellRow,CelCol)
                    int j0 = CelCol * CellSize;
                    int jE = (CelCol + 1) * CellSize;

                    // min, max
                    MinCol = Math.Min(j0, MinCol);
                    MaxCol = Math.Max(j0, MaxCol);
                }
            }
             */
        }


        /// <summary>
        /// Block CSR format
        /// </summary>
        public class BCSR : CelledFormats {
            
            /// <summary>
            /// index: a cell row index <em>r</em>; <br/>
            /// content: the global cell index of the first cell in cell-row <em>r</em>.
            /// </summary>
            public int[] CellRowStart;

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="csr"></param>
            /// <param name="AlignmentInBytes"></param>
            /// <param name="_BlockSize"></param>
            public BCSR(TempCSR csr, int AlignmentInBytes, int _BlockSize)
                : base(csr, AlignmentInBytes, _BlockSize) {

                int NumberOfCells = base.Compute_TotalNoOfCells();
                Val = new double[NumberOfCells * CellStride];

                CellRowStart = new int[NoOfCellRows+1];
                CellColumn = new int[NumberOfCells];

                int Cnt = 0;
                for (int i = 0; i < NoOfCellRows; i++) {
                    CellRowStart[i] = Cnt;
                    int[] CellCols = base.FindCells(i);
                    foreach (int cc in CellCols) {
                        LoadCell(i, cc, Cnt);
                        CellColumn[Cnt] = cc;
                        Cnt++;
                    }
                }
                CellRowStart[NoOfCellRows] = Cnt;

            }
            
            /// <summary>
            /// see <see cref="CelledFormats.GetCellRowStart"/>, see <see cref="CellRowStart"/>;
            /// </summary>
            protected override int GetCellRowStart(int CellRowInd) {
                return CellRowStart[CellRowInd];
            }
        }

        /// <summary>
        /// A variant of the Block CSR format (see <see cref="BCSR"/>) with
        /// constant number of cells per cell row, at the cost of memory-waste due to zero-blocks.
        /// (I.e. ELLPACK with Cells (Blocks));
        /// </summary>
        public class CCBCSR : CelledFormats {

            /// <summary>
            /// number of cells in each cell-row 
            /// </summary>
            public int NoOfCellsPerRow;

            /// <summary>
            /// ctor
            /// </summary>
            public CCBCSR(TempCSR csr, int AlignmentInBytes, int _CellSize)
                : base(csr, AlignmentInBytes, _CellSize) {

                if(_CellSize <= 1) {
                    Console.WriteLine("WARNING FROM MONKEY: internal cell (aka. block matrix) size is 1 -- no decent performance should be expected.");
                }

                NoOfCellsPerRow = base.Compute_MaxNoOfCellsPerRow();
                int NoOfCells = NoOfCellsPerRow * NoOfCellRows;
                Val = new double[NoOfCells * CellStride];
                CellColumn = new int[NoOfCells];

                int Cnt = 0;
                for (int i = 0; i < NoOfCellRows; i++) {
                    int[] CellCols = base.FindCells(i);
                    int j;
                    for (j = 0; j < CellCols.Length; j++) {
                        LoadCell(i, CellCols[j], Cnt);
                        CellColumn[Cnt] = CellCols[j];
                        Cnt++;
                    }

                    for (; j < NoOfCellsPerRow; j++) {
                        if (CellCols.Length > 0)
                            CellColumn[Cnt] = CellCols[0];
                        else
                            CellColumn[Cnt] = 0;
                        Cnt++;
                    }

                }
            }

            /// <summary>
            /// see <see cref="CelledFormats.GetCellRowStart"/>, see <see cref="NoOfCellsPerRow"/>;
            /// </summary>
            protected override int GetCellRowStart(int CellRowInd) {
                return NoOfCellsPerRow * CellRowInd;
            }
        }


        /// <summary>
        /// formats which are similar to ELLPACK
        /// </summary>
        abstract public class ELLPACKlike : FormatBase {

            /// <summary>
            /// maximum number (over all rows) of nonzero entries per row
            /// </summary>
            protected int Compute_MaxNoOfEntriesPerRow() {
                int Max = 0;
                for (int i = 0; i < NoOfRows; i++) {
                    int RowLen = csr.RowStart[i + 1] - csr.RowStart[i];
                    Max = Math.Max(Max, RowLen);
                }
                return Max;
            }

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="_csr">
            /// Source data
            /// </param>
            /// <param name="AlignmentInBytes">
            /// the first element of each column in each block
            /// will be at an (relative) address (with respect to the base address) that is a multiple
            /// of this number
            /// </param>
            /// <param name="_RowsPerBlock">
            /// number of rows that will be packed into one block
            /// </param>
            public ELLPACKlike(TempCSR _csr, int AlignmentInBytes, int _RowsPerBlock)
                : base(_csr) {
                // basic init, check arguments
                // ---------------------------
                NoOfRows = _csr.NoOfRows;
                if (_RowsPerBlock < 1)
                    throw new ArgumentException("must have at least one row per block.", "_RowsPerBlock");
                RowsPerBlock = _RowsPerBlock;
                NoOfPackedCols = Compute_MaxNoOfEntriesPerRow();
                if (AlignmentInBytes < 1)
                    throw new ArgumentException("alignement cannot be smaller than one byte.", "AlignmentInBytes");
                m_AlignemetInBytes = AlignmentInBytes;

                // compute number of blocks
                // ------------------------
                NoOfBlocks = (NoOfRows / _RowsPerBlock);
                if (NoOfBlocks * _RowsPerBlock < NoOfRows)
                    NoOfBlocks++;

                // allocate storage for matrix entries
                // -----------------------------------
                MtxEntries = new PacK<double>(this, AlignmentInBytes);
                Val = MtxEntries.Values;
                ArrayTools.SetAll(Val, double.NaN);

                // pack matrix
                // -----------
                // loop over rows of 'csr' ...
                for (int i = 0; i < NoOfRows; i++) {
                    int c1 = csr.RowStart[i];
                    int c2 = csr.RowStart[i + 1];

                    // loop over columns of 'csr' ...
                    int PackedCol = 0;
                    for (int j = c1; j < c2; j++) {
                        MtxEntries[i, PackedCol] = csr.Vals[j];
                        //ColInd[iii] = csr.ColInd[j];
                        SetColInd(i, PackedCol, csr.ColInd[j]);

                        PackedCol++;
                    }
                    // fill unused entries ...
                    for (; PackedCol < NoOfPackedCols; PackedCol++) {
                        MtxEntries[i, PackedCol] = 0;
                        //ColInd[iii] = csr.ColInd[c1];
                        SetColInd(i, PackedCol, csr.ColInd[c1]);
                    }
                }
            }

            /// <summary>
            /// records the column index;
            /// </summary>
            internal protected abstract void SetColInd(int iRow, int PackedCol, int ColInex);

            /// <summary>
            /// ELLPACK-like storage for matrix entries
            /// </summary>
            internal PacK<double> MtxEntries;

            /// <summary>
            /// numbero of rows that is packed into one (CUDA) block
            /// </summary>
            protected int RowsPerBlock;

            /// <summary>
            /// number of packed columns, i.e. the maximum number (over all rows) of nonzero entries per row
            /// </summary>
            internal int NoOfPackedCols;

            /// <summary>
            /// total number of blocks
            /// </summary>
            protected int NoOfBlocks;

            /// <summary>
            /// local number of rows in the matrix
            /// </summary>
            internal int NoOfRows;

            /// <summary>
            /// the first element of each column in each block
            /// will be at an (relative) address (with respect to the base address) that is a multiple
            /// of this number
            /// </summary>
            protected int m_AlignemetInBytes;

            /// <summary>
            /// 
            /// </summary>
            /// <typeparam name="T"></typeparam>
            public class PacK<T>
                where T : struct {

                internal PacK(ELLPACKlike o, int AlignmentInBytes) {
                    owner = o;
                    

                    this.ColStride = owner.RowsPerBlock;
                    int ColStrideInBytes = ColStride * Marshal.SizeOf(typeof(T));
                    while (ColStrideInBytes % AlignmentInBytes != 0) {
                        this.ColStride++;
                        ColStrideInBytes += Marshal.SizeOf(typeof(T));
                    }
                    
                    Values = new T[owner.NoOfPackedCols * owner.NoOfBlocks * this.ColStride];
                }

                /// <summary>
                /// linear memory for the ELLPACK-like formats, indices into this array are computed by
                /// <see cref="ComputeIndex"/>.
                /// </summary>
                public T[] Values;

                /// <summary>
                /// the stride, in number of items, for the column; <br/>
                /// Refer to external doc. and the implementation of <see cref="ComputeIndex"/> for further information.
                /// </summary>
                public int ColStride;

                /// <summary>
                /// The column stride (see also <see cref="ColStride"/>), in bytes.<br/>
                /// I.e., if the 0-th entry of the 0-th column (in some block) is at address <em>p0</em>, than the
                /// 0-th entry of the <em>k</em>-th column is at <see cref="ColStrideInBytes"/>*<em>k</em>;
                /// </summary>
                public int ColStrideInBytes {
                    get { return Marshal.SizeOf(typeof(T)) * ColStride; }
                }
                
                /// <summary>
                /// 
                /// </summary>
                private ELLPACKlike owner;
                
                /// <summary>
                /// computes an index into the linear memory <see cref="Values"/>.
                /// </summary>
                /// <param name="iRow">local matrix row index</param>
                /// <param name="iPackedCol">packed column index ( smaller thatn <see cref="ELLPACKlike.NoOfPackedCols"/></param>
                /// <returns>an index into <see cref="Values"/></returns>
                public int ComputeIndex(int iRow, int iPackedCol) {
                    if (iRow < 0 || iRow >= owner.NoOfRows)
                        throw new IndexOutOfRangeException();
                    if (iPackedCol < 0 || iPackedCol >= owner.NoOfPackedCols)
                        throw new IndexOutOfRangeException();

                    int blkRow = iRow / owner.RowsPerBlock;
                    int rstRow = iRow - blkRow * owner.RowsPerBlock;

                    int ind = blkRow * owner.NoOfPackedCols * ColStride;
                    ind += iPackedCol * ColStride;
                    ind += rstRow;

                    return ind;
                }

                /// <summary>
                /// gets/sets an element, indices are computed via <see cref="ComputeIndex"/>;
                /// </summary>
                public T this[int iRow, int iPackedCol] {
                    get {
                        int iii = ComputeIndex(iRow, iPackedCol);
                        return Values[iii];
                    }
                    set {
                        int iii = ComputeIndex(iRow, iPackedCol);
                        Values[iii] = value;
                    }
                }
            }

        }
        
        /// <summary>
        /// GPU - optimized ELLPACK format
        /// </summary>
        public class ELLPACKmod : ELLPACKlike {

            /// <summary>
            /// ctor; see <see cref="ELLPACKlike.ELLPACKlike"/>
            /// </summary>
            public ELLPACKmod(TempCSR _csr, int AlignmentInBytes, int _RowsPerBlock)
                : base(_csr, AlignmentInBytes, _RowsPerBlock) {
            }
            
            /// <summary>
            /// column indices for each matrix entry
            /// </summary>
            public PacK<int> ColInd;
            
            /// <summary>
            /// see <see cref="FormatBase.RefSpMv"/>
            /// </summary>
            public override void RefSpMv(double alpha, double[] x, double beta, double[] y) {
                if (y.Length != NoOfRows)
                    throw new ArgumentException("wrong length", "y");
                
                for (int i = 0; i < NoOfRows; i++) {
                    double acc = 0;
                    for (int j = 0; j < NoOfPackedCols; j++) {
                        double MtxVal = base.MtxEntries[i, j];
                        int MtxCol = this.ColInd[i,j];

                        acc += MtxVal * x[MtxCol];
                    }

                    y[i] = y[i] * beta + alpha * acc;
                }
            }
            
            /// <summary>
            /// see <see cref="ELLPACKlike.SetColInd"/>
            /// </summary>
            internal protected override void SetColInd(int iRow, int PackedCol, int ColInex) {
                if (ColInd == null) {
                    ColInd = new PacK<int>(this, base.m_AlignemetInBytes);
                    ArrayTools.SetAll(ColInd.Values, int.MinValue);
                }

                this.ColInd[iRow,PackedCol] = ColInex;
            }

            /// <summary>
            /// see <see cref="FormatBase.GetEntryIndex"/>
            /// </summary>
            public override int GetEntryIndex(int row, int col) {

                for (int j = 0; j < NoOfPackedCols; j++) {
                    //double MtxVal = base.MtxEntries[row, j];
                    int MtxCol = this.ColInd[row, j];

                    if (MtxCol == col)
                        // Bingo!
                        return base.MtxEntries.ComputeIndex(row, j);

                    // remark: In ELLPACK-like formats, unused columns are filled by 0.
                    //         The column index of the unused columns is one that is allready present in the used columns.
                    //         Because we return always the first column with 'MtxCol == col', the unused entries should not 
                    //         cause any problems.
                }

                return int.MinValue;
            }

            /// <summary>
            /// see <see cref="FormatBase.GetAllOccupiedColumns"/>;
            /// </summary>
            public override void GetAllOccupiedColumns(int row, ref int[] Out_MtxColIndices, ref int[] Out_PointersIntoVal, ref double[] Out_Values, out int UsedLen) {
                int L = NoOfPackedCols;

                if (Out_MtxColIndices == null || Out_MtxColIndices.Length < L) { Out_MtxColIndices = new int[L];  }
                if (Out_PointersIntoVal == null || Out_PointersIntoVal.Length < L) { Out_PointersIntoVal = new int[L];  }
                if (Out_Values == null || Out_Values.Length < L) { Out_Values = new double[L]; };
                double[] Values = Out_Values;
                int[] MtxColIndices = Out_MtxColIndices;
                int[] PointersIntoVal = Out_PointersIntoVal;


                //int iBlock = row / base.RowsPerBlock;


                int i = 0;
                for (int j = 0; j < NoOfPackedCols; j++) {
                    int MtxCol = this.ColInd[row, j];

                    MtxColIndices[i] = MtxCol;
                    PointersIntoVal[i] = base.MtxEntries.ComputeIndex(row, j);
                    Values[i] = Val[PointersIntoVal[i]];
                    i++;
                }
                UsedLen = i;
            }

        }

        /// <summary>
        /// an ELLPACK variant for devicec where caching needs to be controlled manually
        /// </summary>
        public class ManualCacheELLPACK : ELLPACKlike {

            /// <summary>
            /// ctor
            /// </summary>
            /// <param name="_csr"></param>
            /// <param name="AlignmentInBytes"></param>
            /// <param name="_RowsPerBlock">
            /// </param>
            /// <param name="_WarpSize">
            /// number of concurrent threads in one NVIDIA CUDA block;
            /// currently not used; 
            /// </param>
            public ManualCacheELLPACK(TempCSR _csr, int AlignmentInBytes, int _RowsPerBlock, int _WarpSize)
                : base(_csr, AlignmentInBytes, _RowsPerBlock) {
                using (new FuncTrace()) {
                    if (_WarpSize < 1)
                        throw new ArgumentException("WrapSize must be at least 1.", "_WrapSize");
                    if (_WarpSize > _RowsPerBlock)
                        throw new ArgumentException("WrapSize must smaller or equal to _RowsPerBlock.", "_WrapSize");
                    if (_RowsPerBlock % _WarpSize != 0)
                        throw new ArgumentException("WrapSize must be a divider of _RowsPerBlock", "_WrapSize");
                    RowsPerWrap = _WarpSize;

                    CollectBlocks();
                    ColIndTmp = null;
                }
            }

            /// <summary>
            /// number of matrix rows per wrap
            /// </summary>
            public int RowsPerWrap;

            /// <summary>
            /// column indices, into the sub-vector of the block;
            /// we use 16-bit numbers to save memory bandwith.
            /// </summary>
            public PacK<ushort> ColIndBlock;


            private PacK<int> ColIndTmp;

            /// <summary>
            /// Defines, for each block, a sub-vector of elements that can be loaded to
            /// shared memory, and is sufficent for the execution of the block.
            /// 1st index: block index <br/>
            /// 2nd index: sub-vector index
            /// </summary>
            public int[,] BlockSubVector;

            /// <summary>
            /// Length of each block - subvector, <see cref="BlockSubVector"/><br/>
            /// index: block index;
            /// </summary>
            public int[] BlockSubVectorLength;


            private void CollectBlocks() {
                using (new FuncTrace()) {
                    // collect the indices of all entries in x ( for SpMV M*x ) that we need in each block ...
                    // =======================================================================================
                    List<int>[] LoadIndicesTmp = new List<int>[NoOfBlocks];
                    int max = 0;
                    HashSet<int> LoadIndicesBlk = new HashSet<int>();
                    // for all blocks ...
                    for (int blk = 0; blk < NoOfBlocks; blk++) {
                        LoadIndicesBlk.Clear();

                        int i0 = base.RowsPerBlock * blk;                      // start row of block 'blk'
                        int IE = Math.Min(i0 + base.RowsPerBlock, NoOfRows);   // start of next block

                        // loop over rows of the block...
                        for (int i = i0; i < IE; i++) {
                            // loop over columns...
                            for (int j = 0; j < NoOfPackedCols; j++) {
                                int ColIdx = ColIndTmp[i, j];
                                LoadIndicesBlk.Add(ColIdx);
                            }
                        }

                        // convert to list and sort
                        List<int> LoadIndicesBlkList = new List<int>(LoadIndicesBlk);
                        LoadIndicesBlkList.Sort();
                        max = Math.Max(max, LoadIndicesBlkList.Count);
                        LoadIndicesTmp[blk] = LoadIndicesBlkList;
                    }

                    if (max > ushort.MaxValue)
                        throw new ApplicationException("too much different x-Values per block needed -> use smaller blocks;");

                    // convert temporary list
                    BlockSubVector = new int[NoOfBlocks, max];
                    ArrayTools.SetAll(BlockSubVector, int.MinValue);
                    BlockSubVector = new int[NoOfBlocks, max];
                    BlockSubVectorLength = new int[NoOfBlocks];
                    for (int blk = 0; blk < NoOfBlocks; blk++) {
                        IList<int> LoadIndicesBlkList = LoadIndicesTmp[blk];
                        BlockSubVectorLength[blk] = LoadIndicesBlkList.Count;
                        for (int j = 0; j < BlockSubVectorLength[blk]; j++) {
                            BlockSubVector[blk, j] = LoadIndicesBlkList[j];
                        }
                    }

                    // assemble 'ColIndBlock'
                    // ======================
                    ColIndBlock = new PacK<ushort>(this, base.m_AlignemetInBytes);
                    ArrayTools.SetAll(ColIndBlock.Values, ushort.MaxValue);
                    for (int blk = 0; blk < NoOfBlocks; blk++) {
                        List<int> LoadIndicesBlkList = LoadIndicesTmp[blk];

                        int i0 = base.RowsPerBlock * blk;                      // start row of block 'blk'
                        int IE = Math.Min(i0 + base.RowsPerBlock, NoOfRows);   // start of next block

                        for (int i = i0; i < IE; i++) {
                            for (int j = 0; j < NoOfPackedCols; j++) {
                                int iii = LoadIndicesBlkList.BinarySearch(ColIndTmp[i, j]);
                                if (iii < 0)
                                    throw new ApplicationException("should not happen - internal error.");
                                ColIndBlock[i, j] = (ushort)iii;
                            }
                        }
                    }
                }
            }

            void Swap(int iRow, int PackedCol_A, int PackedCol_B) {
                {
                    double buf = base.MtxEntries[iRow, PackedCol_A];
                    base.MtxEntries[iRow, PackedCol_A] = this.ColIndBlock[iRow, PackedCol_B];
                    base.MtxEntries[iRow, PackedCol_B] = buf;
                }
                {
                    ushort buf = this.ColIndBlock[iRow, PackedCol_A];
                    this.ColIndBlock[iRow, PackedCol_A] = this.ColIndBlock[iRow, PackedCol_B];
                    this.ColIndBlock[iRow, PackedCol_B] = buf;
                }
            }


            /// <summary>
            /// experimental resorting, specially for NVIDIA GPU's to reduce bank conflicts;
            /// </summary>
            /// <param name="NoOfSharedMemBanks"></param>
            /// <remarks>
            /// see NVIDIA CUDA Programming Guids, Version 3.0, Appendix G.3.3.4, page 146
            /// </remarks>
            private void ResortForUnblockingWraps(int NoOfSharedMemBanks) {
                // loop over all blocks ...
                for (int blk = 0; blk < NoOfBlocks; blk++) {
                    BitArray BankUsedMarker = new BitArray(NoOfSharedMemBanks);

                    // loop over all wraps ...
                    for (int wrp = 0; wrp < (RowsPerBlock / RowsPerWrap); wrp++) {
                        
                        int i0 = base.RowsPerBlock * blk + wrp*RowsPerWrap;  // start row of wrap
                        int IE = i0 + RowsPerWrap;                           // start of next wrap

                        // for all packed columns ...
                        for (int j = 0; j < NoOfPackedCols; j++) {
                            BankUsedMarker.SetAll(false);
                            int BcastBank = -1;
                                                        
                            for (int i = i0; i < IE; i++) {
                                int ColInd = this.ColIndBlock[i, j];

                                int Bank = ColInd % NoOfSharedMemBanks; // the success of the algorithm depends on 
                                //                                         wehter this formula really computes the right shared mem. bank

                                if (BankUsedMarker[Bank]) {
                                    if (BcastBank < 0) {
                                        // one broadcast is OK!
                                        BcastBank = Bank; // we will use broadcasting on this bank!
                                    } else if (Bank == BcastBank) {
                                        // one broadcast is OK!
                                    } else {

                                        throw new NotImplementedException("toDo");

                                    }
                                } else {
                                    BankUsedMarker[Bank] = true;
                                }
                            }
                        }
                    }
                }
            }

            /// <summary>
            /// see <see cref="ELLPACKlike.SetColInd"/>
            /// </summary>
            internal protected override void SetColInd(int iRow, int PackedCol, int ColInex) {
                if (ColIndTmp == null) {
                    ColIndTmp = new PacK<int>(this, base.m_AlignemetInBytes);
                    ArrayTools.SetAll(ColIndTmp.Values, int.MinValue);
                }

                this.ColIndTmp[iRow, PackedCol] = ColInex;
            }

            /// <summary>
            /// see <see cref="FormatBase.RefSpMv"/>
            /// </summary>
            public override void RefSpMv(double alpha, double[] x, double beta, double[] y) {
                // loop over all blocks ...
                for (int blk = 0; blk < NoOfBlocks; blk++) {

                    // load sub-vector for block
                    double[] xSub = new double[this.BlockSubVectorLength[blk]];
                    for (int i = 0; i < xSub.Length; i++) {
                        xSub[i] = x[this.BlockSubVector[blk,i]];
                    }

                    int i0 = base.RowsPerBlock * blk;                      // start row of block 'blk'
                    int IE = Math.Min(i0 + base.RowsPerBlock, NoOfRows);   // start of next block
                    for (int i = i0; i < IE; i++) {

                        double acc = 0;
                        for (int j = 0; j < NoOfPackedCols; j++) {
                            double MtxVal = Val[base.MtxEntries.ComputeIndex(i, j)];
                            int MtxCol = this.ColIndBlock.Values[this.ColIndBlock.ComputeIndex(i, j)];

                            acc += MtxVal * xSub[MtxCol];
                        }

                        y[i] = y[i] * beta + alpha * acc;
                    }
                }
            }
            
            /// <summary>
            /// see <see cref="FormatBase.GetEntryIndex"/>
            /// </summary>
            public override int GetEntryIndex(int row, int col) {
                
                int iBlock = row / base.RowsPerBlock;

                for (int j = 0; j < NoOfPackedCols; j++) {
                    //double MtxVal = Val[base.MtxEntries.ComputeIndex(i, j)];
                    ushort MtxColBlk = this.ColIndBlock.Values[this.ColIndBlock.ComputeIndex(row, j)];
                    int MtxCol = this.BlockSubVector[iBlock, MtxColBlk];
                    
                    if (MtxCol == col)
                        // Bingo!
                        return base.MtxEntries.ComputeIndex(row, j);

                    // remark: In ELLPACK-like formats, unused columns are filled by 0.
                    //         The column index of the unused columns is one that is already present in the used columns.
                    //         Because we return always the first column with 'MtxCol == col', the unused entries should not 
                    //         cause any problems.
                }

                return int.MinValue;
            }

            /// <summary>
            /// see <see cref="FormatBase.GetAllOccupiedColumns"/>;
            /// </summary>
            public override void GetAllOccupiedColumns(int row, ref int[] Out_MtxColIndices, ref int[] Out_PointersIntoVal, ref double[] Out_Values, out int UsedLen) {
                int L = NoOfPackedCols;
                
                if (Out_MtxColIndices == null || Out_MtxColIndices.Length < L) { Out_MtxColIndices = new int[L]; }
                if (Out_PointersIntoVal == null || Out_PointersIntoVal.Length < L) { Out_PointersIntoVal = new int[L]; }
                if (Out_Values == null || Out_Values.Length < L) Out_Values = new double[L];
                double[] Values = Out_Values;
                int[] MtxColIndices = Out_MtxColIndices;
                int[] PointersIntoVal = Out_PointersIntoVal;
                
                int iBlock = row / base.RowsPerBlock;
                
                int i = 0;
                for (int j = 0; j < NoOfPackedCols; j++) {
                    ushort MtxColBlk = this.ColIndBlock.Values[this.ColIndBlock.ComputeIndex(row, j)];
                    int MtxCol = this.BlockSubVector[iBlock, MtxColBlk];

                    MtxColIndices[i] = MtxCol;
                    PointersIntoVal[i] = base.MtxEntries.ComputeIndex(row, j);
                    Values[i] = Val[PointersIntoVal[i]];
                    i++;
                }

                UsedLen = i;
            }
        }
    }
}
