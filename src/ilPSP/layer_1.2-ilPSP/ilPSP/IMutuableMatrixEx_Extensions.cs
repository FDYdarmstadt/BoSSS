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

using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;

namespace ilPSP.LinSolvers {

    /// <summary>
    /// extension methods for <see cref="IMutableMatrix"/> and <see cref="IMutableMatrixEx"/>
    /// </summary>
    public static class IMutuableMatrixEx_Extensions {

        /// <summary>
        /// Returns a collection of all occupied columns in a the row <paramref name="RowIndex"/>;
        /// </summary>
        /// <returns>
        /// The column indices of non-zero/allocated entries.
        /// </returns>
        static public int[] GetOccupiedColumnIndices(this IMutableMatrixEx M, int iRow) {
            int[] idx = null;
            int L = M.GetOccupiedColumnIndices(iRow, ref idx);
            if (L != idx.Length)
                Array.Resize(ref idx, L);
            return idx;
        }

        /// <summary>
        /// Sets all entries in a row to 0
        /// </summary>
        /// <param name="i">row index in global indices</param>
        static public void ClearRow(this IMutableMatrixEx M, int i) {
            int[] ColIdx = null;
            double[] Values = null;
            int L = M.GetRow(i, ref ColIdx, ref Values);
            Array.Clear(Values, 0, L);
            Debug.Assert(ColIdx.Length >= L);
            Debug.Assert(Values.Length >= L);
            if (ColIdx.Length > L)
                Array.Resize(ref ColIdx, L);
            if (Values.Length > L)
                Array.Resize(ref Values, L);
            M.SetValues(i, ColIdx, Values);
        }


        /// <summary>
        /// returns the number of non-zero entries in Row <paramref name="iRow"/>;
        /// </summary>
        /// <param name="iRow">local row index</param>
        /// <param name="M">matrix to operate on</param>
        /// <returns></returns>
        static public int GetNoOfNonZerosPerRow(this IMutableMatrixEx M, int iRow) {

            int cnt = 0;
            int[] col = null;
            double[] val = null;
            int L = M.GetRow(iRow, ref col, ref val);
            for (int l = 0; l < L; l++) {
                if (val[l] != 0.0)
                    cnt++;
            }

            return cnt;
        }

        /// <summary>
        /// returns the number of non-zero entries
        /// in the upper triangle 
        /// in Row <paramref name="iRow"/>;
        /// </summary>
        /// <param name="iRow">local row index</param>
        /// <param name="M">matrix to operate on</param>
        /// <returns></returns>
        static public int GetNoOfUpperTriangualNonZerosPerRow(this IMutableMatrixEx M, int iRow) {
            int cnt = 0;
            int[] col = null;
            double[] val = null;
            int L = M.GetRow(iRow, ref col, ref val);
            for(int l = 0; l < L; l++) {
                if (col[l] > iRow && val[l] != 0.0)
                    cnt++;
            }
            return cnt;
        }

        /// <summary>
        /// returns the number of off-diagonal non-zero entries in Row <paramref name="iRow"/>;
        /// </summary>
        /// <param name="iRow">local row index</param>
        /// <param name="M">matrix to operate on</param>
        /// <returns></returns>
        static public int GetNoOfOffDiagonalNonZerosPerRow(this IMutableMatrixEx M, int iRow) {

            int cnt = 0;
            int[] col = null;
            double[] val = null;
            int L = M.GetRow(iRow, ref col, ref val);
            for (int l = 0; l < L; l++) {
                if (col[l] != iRow && val[l] != 0.0)
                    cnt++;
            }
            return cnt;
        }

        /// <summary>
        /// returns the maximum (over all locally stored rows) number of off-diagonal non-zero entries per row
        /// </summary>
        /// <returns></returns>
        static public int GetMaxNoOfOffDiagonalNonZerosPerRow(this IMutableMatrixEx M) {

            int r = 0;
            var _RowPartiton = M.RowPartitioning;

            int i0 = (int)_RowPartiton.i0;
            for (int i = (int)(_RowPartiton.i0 + _RowPartiton.LocalLength - 1); i >= i0; i--) {
                r = Math.Max(M.GetNoOfOffDiagonalNonZerosPerRow(i), r);
            }
            return r;
        }

        /// <summary>
        /// returns the maximum (over all locally stored rows) number of off-diagonal non-zero entries per row
        /// </summary>
        /// <returns></returns>
        static public int GetMaxNoOfNonZerosPerRow(this IMutableMatrixEx M) {
            int r = 0;
            var _RowPartiton = M.RowPartitioning;

            int i0 = (int)_RowPartiton.i0;
            for (int i = (int)(_RowPartiton.i0 + _RowPartiton.LocalLength - 1); i >= i0; i--) {
                r = Math.Max(M.GetNoOfNonZerosPerRow(i), r);
            }
            return r;
        }

        /// <summary>
        /// returns the number of non-zero elements outside the diagonal in all rows
        /// (on the current MPI process);
        /// </summary>
        /// <returns></returns>
        static public int GetTotalNoOfOffDiagonalNonZeros(this IMutableMatrixEx M) {
            
            int odnz = 0;
            var rowPart = M.RowPartitioning;

            int i0 = (int)rowPart.i0;
            for (int i = (int)(rowPart.i0 + rowPart.LocalLength - 1); i >= i0; i--) {
                odnz += M.GetNoOfOffDiagonalNonZerosPerRow(i);
            }

            return odnz;
        }

        /// <summary>
        /// Sum of All Processors:
        /// Returns the number of non-zero elements  
        /// in the upper triangle.
        /// </summary>
        /// <returns></returns>
        static public int GetGlobalNoOfUpperTriangularNonZeros(this IMutableMatrixEx M) {

            int odnz = M.GetLocalNoOfUpperTriangularNonZeros();

            odnz = odnz.MPISum(M.MPI_Comm);

            return odnz;
        }


        /// <summary>
        /// Only for this processor:
        /// </summary>
        /// Number of non-zero elements  
        /// in the upper triangle.
        /// <returns></returns>
        static public int GetLocalNoOfUpperTriangularNonZeros(this IMutableMatrixEx M) {

            int odnz = 0;
            var rowPart = M.RowPartitioning;
            int i0 = (int)rowPart.i0;
            for (int i = (int)(rowPart.i0 + rowPart.LocalLength - 1); i >= i0; i--) {
                odnz += M.GetNoOfUpperTriangualNonZerosPerRow(i);
            }

            return odnz;

        }

        /// <summary>
        /// returns the number of non-zero elements in the matrix on the current MPI process.
        /// </summary>
        /// <returns></returns>
        static public int GetTotalNoOfNonZerosPerProcess(this IMutableMatrixEx M) {
          

            int odnz = 0;
            var rowPart = M.RowPartitioning;
            int i0 = (int)rowPart.i0;
            for (int i = (int)(rowPart.i0 + rowPart.LocalLength - 1); i >= i0; i--) {
                odnz += M.GetNoOfNonZerosPerRow(i);
            }
           
            return odnz;
        }

        /// <summary>
        /// returns the number of non-zero elements in the matrix
        /// (on the current MPI process);
        /// </summary>
        /// <returns></returns>
        static public int GetTotalNoOfNonZeros(this IMutableMatrixEx M) {
            MPICollectiveWatchDog.Watch(M.RowPartitioning.MPI_Comm);
            int odnz = M.GetTotalNoOfNonZerosPerProcess();
            int odnz_global = odnz.MPISum(M.MPI_Comm);
            return odnz_global;
        }

        /// <summary>
        /// checks all (nonzero) matrix entries for infinity or NAN - values, and
        /// throws an <see cref="ArithmeticException"/> if found;
        /// </summary>
        static public void CheckForNanOrInfM<T>(this T M)
            where T : IMutableMatrixEx //
        {
            int[] col = null;
            double[] val = null;
            for (int i = M.RowPartitioning.i0; i < M.RowPartitioning.iE; i++) {
                int L = M.GetRow(i, ref col, ref val);
                for (int l = 0; l < L; l++) {
                    if (double.IsInfinity(val[l]))
                        throw new ArithmeticException("element (" + i + "," + col[l] + ") is infinity.");
                    if (double.IsNaN(val[l]))
                        throw new ArithmeticException("element (" + i + "," + col[l] + ") is NAN.");
                }
            }
        }

        /// <summary>
        /// converts an arbitrary mutable matrix to an <see cref="MsrMatrix"/>.
        /// </summary>
        /// <param name="M"></param>
        /// <returns></returns>
        static public MsrMatrix ToMsrMatrix(this IMutableMatrixEx M) {
            using (new FuncTrace()) {

                MsrMatrix R = new MsrMatrix(M.RowPartitioning, M.ColPartition);

                int[] col = null;
                double[] val = null;
                int i0 = (int)R.RowPartitioning.i0, L = R.RowPartitioning.LocalLength;
                for (int i = 0; i < L; i++) {
                    int iRow = i0 + i;

                    int Lr = M.GetRow(iRow, ref col, ref val);
                    R.SetRow(iRow, col, val, Lr); 
                }

                return R;
            }
        }

        /// <summary>
        /// performs the operation: <paramref name="Acc"/> = <paramref name="Acc"/> + <paramref name="alpha"/>*<paramref name="M"/>
        /// </summary>
        /// <param name="Acc">
        /// Input/Output: the accumulator
        /// </param>
        /// <param name="alpha">
        /// scaling for accumulation
        /// </param>
        /// <param name="M">
        /// Input: the matrix that is accumulated; unchanged on exit.
        /// </param>
        public static void Acc(this IMutableMatrix Acc, double alpha, IMutableMatrixEx M) {

            if (Acc.NoOfCols != M.NoOfCols)
                throw new ArgumentException("mismatch in number of columns");
            if (Acc.NoOfRows != M.NoOfRows)
                throw new ArgumentException("mismatch in number of rows");

            if (!Acc.RowPartitioning.EqualsPartition(M.RowPartitioning))
                throw new ArgumentException("unable to perform Acc - operation: matrices must have equal row partition.");

            MsrMatrix _M = M as MsrMatrix;

            int I = Acc.RowPartitioning.LocalLength;
            int i0 = (int)Acc.RowPartitioning.i0;

            double[] val = null;
            int[] col = null;
            int L;
            for (int i = 0; i < I; i++) {

                int iRow = i + i0;

                L = M.GetRow(iRow, ref col, ref val);

                for(int l = 0; l < L; l++) {
                    Acc[iRow, col[l]] += alpha * val[l];
                }
            }
        }

        /// <summary>
        /// finds maximum and minimum entry -- within the part that is stored on the local MPI process --
        /// of some matrix.
        /// </summary>
        /// <param name="M">input; the matrix to work on</param>
        /// <param name="Min">output: the minimum entry of <paramref name="M"/> within the local MPI process</param>
        /// <param name="MinRow">output: the row index, where <paramref name="Min"/> is located</param>
        /// <param name="MinCol">output: the column index, where <paramref name="Min"/> is located</param>
        /// <param name="Max">output: the maximum entry of <paramref name="M"/> within the local MPI process</param>
        /// <param name="MaxRow">output: the row index, where <paramref name="Max"/> is located</param>
        /// <param name="MaxCol">output: the column index, where <paramref name="Max"/> is located</param>
        public static void GetMinimumAndMaximum_MPILocal(this IMutableMatrixEx M,
            out double Min, out int MinRow, out int MinCol,
            out double Max, out int MaxRow, out int MaxCol) {

            Min = double.MaxValue;
            Max = double.MinValue;
            MinRow = int.MinValue;
            MinCol = int.MinValue;
            MaxRow = int.MinValue;
            MaxCol = int.MinValue;

            MsrMatrix _M = M as MsrMatrix;

            bool t = false;
            int I = M.RowPartitioning.LocalLength;
            int i0 = (int)M.RowPartitioning.i0;
            double[] val = null;
            int[] col = null;
            int L;
            for (int i = 0; i < I; i++) {

                int iRow = i + i0;

                L = M.GetRow(iRow, ref col, ref val);

                for(int l = 0; l < L; l++) {
                    t = true;

                    int ColIndex = col[l];
                    double Value = val[l];

                    if (Value > Max) {
                        Max = Value;
                        MaxRow = iRow;
                        MaxCol = ColIndex;
                    }

                    if (Value < Min) {
                        Min = Value;
                        MinRow = iRow;
                        MinCol = ColIndex;
                    }
                }
            }

            if (!t) {
                // matrix is completely empty -> min and max is 0.0, 1st occurrence per def. @ (0,0)
                Min = 0;
                MinCol = 0;
                MinRow = 0;
                Max = 0;
                MaxCol = 0;
                MaxRow = 0;
            }
        }

        /// <summary>
        /// collects all locally stored rows of matrix <paramref name="M"/>
        /// </summary>
        static public MsrMatrix.MatrixEntry[][] GetAllEntries(this IMutableMatrixEx M) {
            int i0 = (int)(M.RowPartitioning.i0), L = M.RowPartitioning.LocalLength;
            MsrMatrix.MatrixEntry[][] ret = new MsrMatrix.MatrixEntry[L][];

            double[] val = null;
            int[] col = null;
            int Lr;

            for (int i = 0; i < L; i++) {
                //ret[i] = M.GetRow(i + i0);
                Lr = M.GetRow(i + i0, ref col, ref val);
                var row = new MsrMatrix.MatrixEntry[Lr];
                for(int lr = 0; lr < Lr; lr++) {
                    row[lr].m_ColIndex = col[lr];
                    row[lr].Value = val[lr];
                }
                ret[i] = row;
            }
            return ret;
        }

        // workaround for some mono bug in BinaryFormatter
        [Serializable()]
        class Helper {
            public MsrMatrix.MatrixEntry[][] entries;
        }

        /// <summary>
        /// Saves the matrix in a custom sparse text format;
        /// Mainly for importing into MATLAB;
        /// </summary>
        /// <param name="path">Path to the text file</param>
        /// <remarks>
        /// MATLAB code for importing the matrix (save as file 'ReadMsr.m'):
        /// <code>
        /// function Mtx = ReadMsr(filename)
        /// 
        /// fid = fopen(filename);
        /// % matrix dimensions
        /// % -----------------
        /// NoOfRows = fscanf(fid,'%d',1);
        /// NoOfCols = fscanf(fid,'%d',1);
        /// NonZeros = fscanf(fid,'%d',1);
        /// cnt = 1;
        /// % read row and column array
        /// % -------------------------
        /// iCol = zeros(NonZeros,1);
        /// iRow = zeros(NonZeros,1);
        /// entries = zeros(NonZeros,1);
        /// l0 = 0;
        /// str = char(zeros(1,6));
        /// for i = 1:NoOfRows
        ///     NonZerosInRow = fscanf(fid,'%d',1);
        ///     if(l0 ~= NonZerosInRow)
        ///         str = char(zeros(1,NonZerosInRow*6));
        ///         for j = 1:NonZerosInRow
        ///             i0 = 1+(j-1)*6;
        ///             str(i0:i0+1) = '%f';
        ///             str(i0+3:i0+4) = '%f';
        ///         end
        ///     end
        ///     R = fscanf(fid,str,2*NonZerosInRow);
        ///     R2 = reshape(R',2,NonZerosInRow);
        ///     ind = cnt:(cnt+NonZerosInRow-1);
        ///     iCol(ind) = R2(1,:);
        ///     iRow(ind) = i;
        ///     entries(ind) = R2(2,:);
        /// 
        ///     cnt = cnt + NonZerosInRow;
        /// end
        /// fclose(fid);
        /// 
        /// if (cnt-1) &lt; NonZeros
        ///     iCol = iCol(1:(cnt-1),1);
        ///     iRow = iRow(1:(cnt-1),1);
        ///     entries = entries(1:(cnt-1),1);
        /// end
        /// 
        /// % create sparse matrix
        /// % --------------------
        /// Mtx = sparse(iRow,iCol+1,entries,NoOfRows,NoOfCols,NonZeros);
        /// 
        /// </code>
        /// </remarks>
        /// <param name="M">
        /// this pointer of extension method
        /// </param>
        static public void SaveToTextFileSparse(this IMutableMatrixEx M, string path) {
            using (new FuncTrace()) {
                int rank, size;

                csMPI.Raw.Comm_Rank(M.MPI_Comm, out rank);
                csMPI.Raw.Comm_Size(M.MPI_Comm, out size);

                SerialisationMessenger sms = new SerialisationMessenger(M.MPI_Comm);

                int NoOfNonZeros = M.GetTotalNoOfNonZeros();
                
                if (rank == 0) {
                    sms.CommitCommPaths();

                    // receive data from other processors
                    MsrMatrix.MatrixEntry[][] entries = M.GetAllEntries();
                    if (size > 1) {
                        Array.Resize(ref entries, (int)(M.RowPartitioning.TotalLength));
                    }

                    Helper rcvdata; int rcvRank;
                    while (sms.GetNext(out rcvRank, out rcvdata)) {
                        Array.Copy(rcvdata.entries, 0, entries, (int)M.RowPartitioning.GetI0Offest(rcvRank), rcvdata.entries.Length);
                    }

                    // open file
                    StreamWriter stw = new StreamWriter(path);

                    // serialize matrix data
                    stw.WriteLine(M.RowPartitioning.TotalLength); // number of rows
                    stw.WriteLine(M.NoOfCols);           // number of columns
                    stw.WriteLine(NoOfNonZeros);                 // number of non-zero entries in Matrix (over all MPI-processors)

                    for (int i = 0; i < entries.Length; i++) {
                        MsrMatrix.MatrixEntry[] row = entries[i];

                        int NonZPRow = 0;
                        foreach (MsrMatrix.MatrixEntry e in row) {
                            if (e.ColIndex >= 0 && e.Value != 0.0) NonZPRow++;
                        }
                        stw.Write(NonZPRow);
                        stw.Write(" ");

                        foreach (MsrMatrix.MatrixEntry e in row) {
                            if (e.ColIndex >= 0  && e.Value != 0.0 ) {
                                stw.Write(e.ColIndex);
                                stw.Write(" ");
                                stw.Write(e.Value.ToString("E16", NumberFormatInfo.InvariantInfo));
                                stw.Write(" ");
                            }
                        }

                        stw.WriteLine();
                    }

                    // finalize
                    stw.Flush();
                    stw.Close();
                } else {
                    sms.SetCommPath(0);
                    sms.CommitCommPaths();


                    var entries = M.GetAllEntries();
                    var c = new Helper();
                    c.entries = entries;
                    sms.Transmitt(0, c);


                    MsrMatrix.MatrixEntry[][] dummy; int dummy_;
                    if (sms.GetNext<MsrMatrix.MatrixEntry[][]>(out dummy_, out dummy))
                        throw new ApplicationException("error in app");

                }

                sms.Dispose();
            }
        }

        /// <summary>
        /// Adds <paramref name="factor"/> to all diagonal entries of <paramref name="M"/>.
        /// </summary>
        static public void AccEyeSp(this ISparseMatrix M, double factor = 1.0) {
            using (new FuncTrace()) {
                if (M.RowPartitioning.LocalLength != M.ColPartition.LocalLength)
                    throw new ArgumentException("supported only for quadratic matrices");

                var rm = M.RowPartitioning;

                int i0 = (int)rm.i0;
                int L = rm.LocalLength;

                for (int i = 0; i < L; i++) {
                    M.SetDiagonalElement(i + i0, M.GetDiagonalElement(i + i0) + factor);
                }

            }
        }
        
    
        /// <summary>
        /// accumulates a dense matrix <paramref name="FullMtx"/> to a sparse matrix
        /// -- certainly, only adviseable for small matrices.
        /// </summary>
        static public void AccDenseMatrix(this IMutableMatrixEx tis, double alpha, IMatrix FullMtx) {
            if(tis.RowPartitioning.LocalLength != FullMtx.NoOfRows)
                throw new ArgumentException("Mismatch in number of rows.");
            if(tis.ColPartition.TotalLength != FullMtx.NoOfCols)
                throw new ArgumentException("Mismatch in number of columns.");

            int i0 = tis.RowPartitioning.i0;
            int I = tis.RowPartitioning.LocalLength;
            int J = tis.ColPartition.TotalLength;
            int[] col = null;
            double[] val = null;

            for(int i = 0; i < I; i++) {
                
                int Lr = tis.GetRow(i + i0, ref col, ref val);
                var oldRow = new MsrMatrix.MatrixEntry[Lr];
                for (int lr = 0; lr < Lr; lr++) {
                    oldRow[i].m_ColIndex = col[lr];
                    oldRow[i].Value = val[lr];
                }

                List<int> NewColIdx = new List<int>(J);
                List<double> NewVals = new List<double>(J);
                for(int j = 0; j < J; j++) {
                    double FMij = FullMtx[i, j];
                    if(FMij != 0.0) {
                        NewVals.Add( alpha * FMij);
                        NewColIdx.Add(j);
                    }
                }

                Array.Sort<MsrMatrix.MatrixEntry>(oldRow);
                int k1 = 0, k2 = 0, K1 = oldRow.Length, K2 = NewVals.Count;
                while(k1 < K1 && k2 < K2) {
                    int j1 = oldRow[k1].m_ColIndex;
                    int j2 = NewColIdx[k2];

                    if(j1 < 0)
                        // should also chrash in RELEASE, therefor -> Exception.
                        throw new ApplicationException("expecting a row without un-allocated entries.");

                    if(j1 > j2) {
                        // 
                        k2++; // new row neds to catch up
                    } else if(j1 < j2) {
                        k1++;
                    } else {
                        NewVals[k2] += oldRow[k1].Value;
                        k1++;
                        k2++;
                    }
                }

                tis.SetValues(i + i0, NewColIdx.ToArray(), NewVals.ToArray());
            }
        }

        /// <summary>
        /// only for debugging and testing; converts the matrix to a full matrix;
        /// when running in parallel, the matrix is collected on process 0.
        /// </summary>
        /// <returns>
        /// null on all MPI processes, except on rank 0;
        /// </returns>
        static public MultidimensionalArray ToFullMatrixOnProc0(this IMutableMatrixEx tis) {
            int Rank, Size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out Size);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);


            double[,] ret = null;
            if (Rank == 0)
                ret = new double[tis.NoOfRows, tis.NoOfCols];

            SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);

            if (Rank > 0)
                sms.SetCommPath(0);
            sms.CommitCommPaths();

            Tuple<int,int[],double[]>[] data;
            {
                int L = tis.RowPartitioning.LocalLength;
                data = new Tuple<int, int[], double[]>[L];

                int i0 = (int)tis.RowPartitioning.i0;
                for (int i = 0; i < L; i++) {
                    double[] val = null; // this mem. must be inside the loop/allocated for every i and cannot be reused because we need all rows later.
                    int[] col = null;
                    int Lr = tis.GetRow(i + i0, ref col, ref val);
                    data[i] = new Tuple<int, int[], double[]>(Lr, col, val);
                }
            }

            if (Rank > 0)
                sms.Transmitt(0, data);

            int rcvProc = 0;
            if (Rank == 0) {
                do {
                    int i0 = (int)tis.RowPartitioning.GetI0Offest(rcvProc);
                    if (data.Length != tis.RowPartitioning.GetLocalLength(rcvProc))
                        throw new ApplicationException("internal error");

                    for (int i = 0; i < data.Length; i++) {
                        //foreach (MsrMatrix.MatrixEntry entry in data[i]) {
                        //    if (entry.m_ColIndex >= 0)
                        //        ret[i + i0, entry.m_ColIndex] = entry.Value;
                        //}
                        int Lr = data[i].Item1;
                        int[] col = data[i].Item2;
                        double[] val = data[i].Item3;
                        for (int lr = 0; lr < Lr; lr++)
                            ret[i + i0, col[lr]] = val[lr];
                    }

                } while (sms.GetNext(out rcvProc, out data));
            } else {
                if (sms.GetNext(out rcvProc, out data))
                    throw new ApplicationException("internal error");
            }

            if (Rank == 0) {
                var _ret = MultidimensionalArray.Create(ret.GetLength(0), ret.GetLength(1));
                for (int i = 0; i < _ret.NoOfRows; i++)
                    for (int j = 0; j < _ret.NoOfCols; j++)
                        _ret[i, j] = ret[i, j];
                return _ret;
            } else {
                return null;
            }
        }

        /// <summary>
        /// Writes the matrix (including all zeros) into a string (for debugging purposes).
        /// Basically, it uses a tabulator separated format (which can e.g. be
        /// imported into Matlab via <code>matrix = dlmread(path)</code>
        /// </summary>
        /// <param name="tis">
        /// the matrix to save
        /// </param>
        static public string ToStringDense(this IMutableMatrixEx tis) {

            string OutputString = "";

            double[] val = null;
            int[] col = null;
            int L;

            for (int i = 0; i < tis.RowPartitioning.LocalLength; i++) {

                L = tis.GetRow(i + tis.RowPartitioning.i0, ref col, ref val);
                int currentColumn = -1;

                string separator = "";
                for (int j = 0; j < L; j++) {
                    // Beware of undefined entries (see MSREntry)
                    

                    // Add zeros for missing columns (the sparse format does
                    // not store zero values)
                    for (int k = currentColumn + 1; k < col[j]; k++) {
                        OutputString += String.Format(separator + "{0,14:F0}", 0.0);
                        separator = "\t";
                    }

                    // Enforce use of . as decimal separator in scientific format
                    OutputString += (separator + val[j].ToString(
                        "E", System.Globalization.CultureInfo.InvariantCulture).PadLeft(14));
                    currentColumn = col[j];
                    separator = "\t";
                }

                // Add zeros for columns following after the last entry
                for (int j = currentColumn + 1; j < tis.NoOfCols; j++) {
                    OutputString += String.Format(separator + "{0,14:F0}", 0.0);
                }

                OutputString +=("\n");
            }
            return OutputString;
        }

        /// <summary>
        /// Writes the matrix (including all zeros) into a text file (for debugging purposes).
        /// Basically, it uses a tabulator separated format (which can e.g. be
        /// imported into Matlab via <code>matrix = dlmread(path)</code>
        /// </summary>
        /// <remarks>
        /// The content of the file specified via <paramref name="path"/> will
        /// be overwritten.<br/>
        /// If running on more than one process, multiple files are written 
        /// (they can be concatenated e.g. by 'cat' or directly in Matlab).
        /// </remarks>
        /// <param name="path">Path to the text file</param>
        /// <param name="tis">
        /// the matrix to save
        /// </param>
        static public void SaveToTextFile(this IMutableMatrixEx tis, string path) {

            int Rank, Size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out Size);
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);
            string append = "";
            if (Size > 1)
                append = "_" + (Rank + 1) + "of" + Size;

            StreamWriter writer = new StreamWriter(path + append);
            writer.Write(tis.ToStringDense());
            writer.Close();

        }

        /// <summary>
        /// The Infinity-Norm (maximum absolute row sum norm) of this matrix;
        /// </summary>
        /// <returns></returns>
        static public double InfNorm(this IMutableMatrixEx M) {
            using (var tr = new FuncTrace()) {
                double normLoc = 0;

                int L = M.RowPartitioning.LocalLength;
                int[] col = null;
                double[] val = null;
                int Lr;
                for (int i = 0; i < L; i++) {
                    double rownrm = 0;
                    Lr = M.GetRow(i + M.RowPartitioning.i0, ref col, ref val);
                    for (int j = 0; j < Lr; j++) {
                        rownrm += Math.Abs(val[j]);
                    }

                    normLoc = Math.Max(normLoc, rownrm);
                }

                //Console.WriteLine("local norm (R=" + M.RowPartitioning.Rank + ") = " + normLoc);
                //tr.Info("local norm " + normLoc);

                double normGlob = double.NaN;
                unsafe {
                    csMPI.Raw.Allreduce((IntPtr)(&normLoc), (IntPtr)(&normGlob), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                }
                return normGlob;
            }
        }

        /// <summary>
        /// extracts the diagonal vector from a matrix.
        /// </summary>
        static public double[] GetDiagVector(this IMutableMatrix M) {
            int i0 = M.RowPartitioning.i0;
            int L = M.RowPartitioning.LocalLength;
            double[] diag = new double[L];
            for(int i = 0; i < L; i++) {
                diag[i] = M[i + i0, i + i0];
            }
            return diag;
        }




        /// <summary>
        /// Extracts a sub-matrix from this one, see also
        /// <see cref="AccSubMatrixTo{V1,V2,V3,V4}"/>.
        /// </summary>
        static public MsrMatrix GetSubMatrix<V1, V3>(this IMutableMatrixEx OrgMtx, V1 RowIndicesSource, V3 ColumnIndiceSource)
            where V1 : IList<int>
            where V3 : IList<int> //
        {
            throw new NotImplementedException();
        }


        /// <summary>
        /// Similar to <see cref="AccSubMatrixTo{V1,V2,V3,V4}"/>, 
        /// but the destination (<paramref name="target"/>) is cleared before the accumulation.
        /// </summary>
        static public void WriteSubMatrixTo<V1, V2, V3, V4>(this IMutableMatrixEx OrgMtx, MsrMatrix target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColumnInidcesTarget)
            where V1 : IList<int>
            where V2 : IList<int>
            where V3 : IList<int>
            where V4 : IList<int> {

            throw new NotImplementedException();
        }
    }
}
