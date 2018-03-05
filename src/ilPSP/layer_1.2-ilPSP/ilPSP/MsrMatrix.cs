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
using System.Diagnostics;
using System.IO;
using System.Runtime.InteropServices;
using System.Runtime.Serialization.Formatters.Binary;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;

namespace ilPSP.LinSolvers {

    /// <summary>
    /// Class for sparse matrices - it is mutable and used during equation
    /// assembly. When the system of equations is complete, this object can be
    /// handed over to an <see cref="ISparseSolver"/> - implementation that
    /// solves the system. For performance reasons, this solver typically will
    /// convert this matrix into his internal format.<br/>
    /// This matrix is addressed by global row/column indices; 
    /// </summary>
    /// <remarks>
    /// MSR stands for 'M'utuble 'S'parse 'R'ow;
    /// </remarks>
    public class MsrMatrix : ICloneable, IMutableMatrixEx {

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
        /// <see cref="RowPartitioning"/>;
        /// </summary>
        IPartitioning m_RowPartitioning;

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IPartitioning RowPartitioning {
            get {
                return m_RowPartitioning;
            }
        }

        /// <summary>
        /// <see cref="RowPartitioning"/>;
        /// </summary>
        IPartitioning m_ColPartitioning;

        /// <summary>
        /// distribution of matrix rows over MPI processors
        /// </summary>
        public IPartitioning ColPartition {
            get {
                return m_ColPartitioning;
            }
        }

        


        /// <summary>
        /// counterpart to <see cref="SaveToFile"/>
        /// </summary>
        /// <param name="path"></param>
        /// <param name="mpi_comm">
        /// mpi communicator; currently, only MPI_COMM_WORLD is supported; 
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// Does not work when running in parallel;
        /// </remarks>
        /// <param name="RowPart">
        /// Optional row partition of the output matrix; if null, a row partition is chosen automatically.
        /// </param>
        /// <param name="ColPart">
        /// Optional column partition of the output matrix; if null, a column partition is chosen automatically.
        /// </param>
        public static MsrMatrix LoadFromFile(string path, MPI_Comm mpi_comm, IPartitioning RowPart = null, IPartitioning ColPart = null) {

            int rank, size;
            csMPI.Raw.Comm_Rank(mpi_comm, out rank);
            csMPI.Raw.Comm_Size(mpi_comm, out size);

            SerialisationMessenger sms = new SerialisationMessenger(mpi_comm);

            MsrMatrix M;
            if (rank == 0) {

                FileStream fs = new FileStream(path, FileMode.Open);
                BinaryFormatter bf = new BinaryFormatter();

                // deserialize No of Cols, Rows, create Matrix
                long[] Dim = new long[2];
                Dim[0] = (long)bf.Deserialize(fs);
                Dim[1] = (long)bf.Deserialize(fs);
                unsafe {
                    fixed (long* pDim = &Dim[0]) {
                        csMPI.Raw.Bcast((IntPtr)pDim, 2, csMPI.Raw._DATATYPE.LONG_LONG_INT, 0, mpi_comm);
                    }
                }
                RowPart = FindPartitioning(RowPart, mpi_comm, rank, size, Dim[0]);
                ColPart = FindPartitioning(ColPart, mpi_comm, rank, size, Dim[1]);


                //Partition par = new Partition((int)NoOfRows);
                M = new MsrMatrix(RowPart, ColPart);

                // load matrix
                MatrixEntry[][] entries = (MatrixEntry[][])bf.Deserialize(fs);
                //foreach (var row in entries) {
                //    foreach (var entri in row) {
                //        if (entri.ColIndex >= 0 && entri.val != 0.0) {
                //            Console.WriteLine("fuck me.");
                //        }
                //    }
                //}
                Array.Copy(entries, 0, M.m_Entries, 0, (int)RowPart.GetLocalLength(0));

                // close file
                fs.Close();

                // distribute data
                for (int r = 1; r < size; r++)
                    sms.SetCommPath(r);
                sms.CommitCommPaths();

                for (int r = 1; r < size; r++) {
                    MatrixEntry[][] entriesSend = new MatrixEntry[RowPart.GetLocalLength(r)][];
                    Array.Copy(entries, (int)RowPart.GetI0Offest(r), entriesSend, 0, entriesSend.Length);
                    sms.Transmitt(r, entriesSend);
                }

                MatrixEntry[][] dummy;
                int _dummy;
                if (sms.GetNext(out _dummy, out dummy))
                    throw new ApplicationException("internal error.");

            } else {
                sms.CommitCommPaths();

                // receive No of Cols, Rows, create Matrix
                long[] Dim = new long[2];
                unsafe {
                    fixed (long* pDim = &Dim[0]) {
                        csMPI.Raw.Bcast((IntPtr)pDim, 2, csMPI.Raw._DATATYPE.LONG_LONG_INT, 0, csMPI.Raw._COMM.WORLD);
                    }
                }
                RowPart = FindPartitioning(RowPart, mpi_comm, rank, size, Dim[0]);
                ColPart = FindPartitioning(ColPart, mpi_comm, rank, size, Dim[1]);
                M = new MsrMatrix(RowPart, RowPart);

                // receive data
                int rcvRank;
                sms.GetNext(out rcvRank, out M.m_Entries);
                if (M.m_Entries.Length != M.RowPartitioning.LocalLength)
                    throw new ApplicationException("internal error.");

                if (rcvRank != 0)
                    throw new ApplicationException("internal error.");
                MatrixEntry[][] dummy;
                int _dummy;
                if (sms.GetNext(out _dummy, out dummy))
                    throw new ApplicationException("internal error.");
            }

            sms.Dispose();
            return M;
        }

        private static IPartitioning FindPartitioning(IPartitioning Part, MPI_Comm mpi_comm, int rank, int size, long Dim) {
            if (Part != null) {
                if (Dim != Part.TotalLength)
                    throw new ArgumentException("mismatch between provided row partition and no of rows in file.", "RowPart");
            } else {
                int locLen = (int)(Dim * (rank + 1) / size - Dim * rank / size);
                Part = new Partitioning(locLen, mpi_comm);
            }
            return Part;
        }

        /// <summary>
        /// writes this matrix to a file (binary format, using 
        /// serialization and
        /// an <see cref="BinaryFormatter"/>);
        /// Mainly for debugging purposes;
        /// </summary>
        /// <param name="path"></param>
        /// <param name="fm">
        /// file mode; accepts only <see cref="FileMode.Create"/> 
        /// or <see cref="FileMode.CreateNew"/>;
        /// </param>
        /// <remarks>
        /// In this method  IO is only done on process 0, it does not scale;
        /// This method fails on mono, see https://bugzilla.xamarin.com/show_bug.cgi?id=10699.
        /// </remarks>
        public void SaveToFile(string path, FileMode fm = FileMode.Create) {
            if (fm != FileMode.Create && fm != FileMode.CreateNew)
                throw new ArgumentException("only 'Create' or 'CreateNew' are accepted;", "fm");

            int rank, size;
            MPI_Comm comm = RowPartitioning.MPI_Comm;
            csMPI.Raw.Comm_Rank(comm, out rank);
            csMPI.Raw.Comm_Size(comm, out size);

            SerialisationMessenger sms = new SerialisationMessenger(comm);

            if (rank == 0) {
                sms.CommitCommPaths();

                // receive data from other processors
                MatrixEntry[][] entries = m_Entries;
                if (size > 1) {
                    entries = new MatrixEntry[m_RowPartitioning.TotalLength][];
                    Array.Copy(m_Entries, entries, m_Entries.Length);
                }

                MatrixEntry[][] rcvdata;
                int rcvRank;
                while (sms.GetNext(out rcvRank, out rcvdata)) {
                    Array.Copy(rcvdata, 0, entries, (int)this.RowPartitioning.GetI0Offest(rcvRank), rcvdata.Length);
                }

                // open file
                FileStream fs = new FileStream(path, fm);
                BinaryFormatter bf = new BinaryFormatter();

                // serialize matrix data
                bf.Serialize(fs, (long)(this.RowPartitioning.TotalLength));
                bf.Serialize(fs, (long)(this.ColPartition.TotalLength));

                bf.Serialize(fs, entries);

                // finalize
                fs.Flush();
                fs.Close();
            } else {
                sms.SetCommPath(0);
                sms.CommitCommPaths();

                sms.Transmitt(0, m_Entries);

                MatrixEntry[][] dummy;
                int dummy_;
                if (sms.GetNext<MatrixEntry[][]>(out dummy_, out dummy))
                    throw new ApplicationException("error in app");
            }

            sms.Dispose();
        }

        /// <summary>
        /// constructs a new rectangular MSR matrix, i.e. number of columns and rows can be different. 
        /// </summary>
        /// <param name="LocalNoOfRows">
        /// number of rows that is assigned to the current MPI process
        /// </param>
        /// <param name="LocalNoOfColumns">
        /// number of columns that this matrix should have:
        /// the total number of columns is the sum over the local number, over all MPI processes.
        /// in total over all MPI processes.
        /// This argument must be EQUAL on all MPI processors (not checked for performance
        /// reasons;
        /// </param>
        /// <param name="_ColPerBlock">Block length of the column partition; if in doubt, use 1.</param>
        /// <param name="_RowsPerBlock">Block length of the row partition; if in Doubt, use 1</param>
        public MsrMatrix(int LocalNoOfRows, int LocalNoOfColumns, int _RowsPerBlock, int _ColPerBlock)
            : this(new Partitioning(LocalNoOfRows, csMPI.Raw._COMM.WORLD),
                   new Partitioning(LocalNoOfColumns, csMPI.Raw._COMM.WORLD)) {
        }

        /// <summary>
        /// constructs a new quadratic MSR matrix;
        /// </summary>
        /// <param name="LocalNumberOfRows">
        /// number of rows that is assigned to the current MPI process
        /// </param>
        /// <param name="BlockSize">
        /// blocking size.
        /// </param>
        public MsrMatrix(int LocalNumberOfRows, int BlockSize = 1)
            : this(new Partitioning(LocalNumberOfRows, csMPI.Raw._COMM.WORLD)) {
        }

        /// <summary>
        /// constructs a new quadratic MSR matrix;
        /// </summary>
        /// <param name="RowColPart">
        /// row and column partition
        /// </param>
        public MsrMatrix(IPartitioning RowColPart)
            : this(RowColPart, RowColPart) {
        }

        /// <summary>
        /// constructs a new MSR matrix; 
        /// </summary>
        /// <param name="RowPart">
        /// See <see cref="ISparseMatrix.RowPartitioning"/>.
        /// </param>
        /// <param name="ColPart">
        /// See <see cref="ISparseMatrix.ColPartition"/>.
        /// </param>
        /// <param name="AllocEntriesPerRow">
        /// the number of entries that are allocated on construction, for each row
        /// </param>
        public MsrMatrix(IPartitioning RowPart, IPartitioning ColPart, int AllocEntriesPerRow = 50) {
            //if (_RowMap.Count != _ColMap.Count)
            //    throw new ArgumentException("size mismatch between RowMap and ColMap.");

            if (RowPart.MPI_Comm != ColPart.MPI_Comm)
                throw new ArgumentException("row and column partition must share the same MPI communicator.");

            m_ColPartitioning = ColPart.IsMutable ? ColPart.GetImmutablePartition() : ColPart;
            m_RowPartitioning = RowPart.IsMutable ? RowPart.GetImmutablePartition() : RowPart;

            int _NoOfRows = m_RowPartitioning.LocalLength;
            int l = AllocEntriesPerRow;
            if (l < 0)
                throw new ArgumentOutOfRangeException();
            m_Entries = new MatrixEntry[_NoOfRows][];
            for (int i = 0; i < _NoOfRows; i++) {
                if (l > 0) {
                    m_Entries[i] = new MatrixEntry[l];
                    unsafe
                    {
                        fixed (MatrixEntry* pRow = m_Entries[i])
                        {
                            MatrixEntry* pEntry = pRow;
                            for (int k = 0; k < l; k++) {
                                pEntry->m_ColIndex = -1;
                                pEntry++;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// initializes this matrix to be a clone of <paramref name="otherMtx"/>
        /// </summary>
        public MsrMatrix(IMutableMatrixEx otherMtx) :
            this(otherMtx.RowPartitioning.IsMutable ? otherMtx.RowPartitioning.GetImmutablePartition() : otherMtx.RowPartitioning,
                otherMtx.ColPartition.IsMutable ? otherMtx.ColPartition.GetImmutablePartition() : otherMtx.ColPartition) {
            this.Acc(1.0, otherMtx);
        }

        /// <summary>
        /// a flag which indicates that a solver can assume this matrix to be symmetric;<br/>
        /// For performance reasons, we rely on the user to set this correctly (checking
        /// for symmetry is very expensive even on one pro
        /// </summary>
        public bool AssumeSymmetric {
            get {
                return m_AssumeSymmetric;
            }
            set {
                m_AssumeSymmetric = value;
            }
        }

        bool m_AssumeSymmetric;

        /// <summary>
        /// creates a non-shallow copy of this object
        /// </summary>
        public MsrMatrix CloneAs() {
            MsrMatrix ret = new MsrMatrix(this.RowPartitioning, this.ColPartition);

            int L = this.m_Entries.Length;
            for (int i = 0; i < L; i++) {
                ret.m_Entries[i] = this.m_Entries[i] != null ? ((MatrixEntry[])(this.m_Entries[i].Clone())) : null;
            }

            return ret;
        }

        /// <summary>
        /// creates a non-shallow copy of this object
        /// </summary>
        public object Clone() {
            return CloneAs();
        }

        /// <summary>
        /// sets all entries to 0
        /// </summary>
        public void Clear() {
            for (int i = m_Entries.Length - 1; i >= 0; i--) {
                MatrixEntry[] row = m_Entries[i];
                if (row != null) {
                    int l = row.Length;
                    for (int k = 0; k < l; k++) {
                        row[k].m_ColIndex = -1;
                        row[k].Value = double.NaN;
                    }
                }
            }
        }

        /// <summary>
        /// deletes the <paramref name="ind"/>-th entry of some row,
        /// by swapping the <paramref name="ind"/>-th element with the last
        /// valid one;
        /// </summary>
        /// <param name="row"></param>
        /// <param name="ind">
        /// index into <paramref name="row"/>;
        /// </param>
        private void ClearEntry(MatrixEntry[] row, int ind) {
            int L = row.Length;
            int len = 0;
            while (len < L && row[len].ColIndex >= 0) {
                len++;
            }
            len--;

            if (ind != len) {
                row[ind] = row[len];
            }
            row[len].m_ColIndex = -12345;

        }

        /// <summary>
        /// Sparse Matrix/Vector multiplication;
        /// the calculation 
        /// <paramref name="acc"/> = <paramref name="acc"/>*<paramref name="beta"/>
        /// + this*<paramref name="a"/>*<paramref name="alpha"/>
        /// is performed;
        /// </summary>
        /// <typeparam name="VectorType1"></typeparam>
        /// <typeparam name="VectorType2"></typeparam>
        /// <param name="alpha"></param>
        /// <param name="a">
        /// length must be at least least the count
        /// </param>
        /// <param name="beta"></param>
        /// <param name="acc">
        /// length of accumulator must be at least the update length 
        /// </param>
        public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {

            int size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            if (size != 1)
                throw new NotImplementedException("SpMV(...) does not work for more than 1 MPI process.");

            if (acc.Count != this.m_RowPartitioning.LocalLength)
                throw new ArgumentException("length of accumulator must be at least the update length of the row map.");
            if (a.Count != this.ColPartition.LocalLength)
                throw new ArgumentException("length of a must be at least the count (updated + external) of the column map.");
            if (object.ReferenceEquals(a, acc))
                throw new ArgumentException("in-place computation is not supported.", "a,acc");


            int noOfRows = m_RowPartitioning.LocalLength;

            for (int i = 0; i < noOfRows; i++) {
                MatrixEntry[] row = m_Entries[i];
                double s = 0.0;
                if (row != null) {
                    for (int j = 0; j < row.Length; j++) {
                        if (row[j].m_ColIndex < 0)
                            break;

                        s += row[j].Value * a[row[j].m_ColIndex];
                    }
                    s *= alpha;
                }
                if (beta != 0.0)
                    s += beta * acc[i];
                acc[i] = s;
                Debug.Assert(acc[i] == s);
            }

        }


        /// <summary>
        /// Extracts a sub-matrix from this one, see also
        /// <see cref="AccSubMatrixTo{V1,V2,V3,V4}"/>.
        /// </summary>
        public MsrMatrix GetSubMatrix<V1, V3>(V1 RowIndicesSource, V3 ColumnIndiceSource)
            where V1 : IList<int>
            where V3 : IList<int> //
        {
            MsrMatrix res = new MsrMatrix(RowIndicesSource.Count, ColumnIndiceSource.Count, 1, 1);
            this.AccSubMatrixTo(
                1.0,
                res,
                RowIndicesSource,
                default(int[]),
                ColumnIndiceSource,
                default(int[]));
            return res;
        }


        /// <summary>
        /// Similar to <see cref="AccSubMatrixTo{V1,V2,V3,V4}"/>, 
        /// but the destination (<paramref name="target"/>) is cleared before the accumulation.
        /// </summary>
        public void WriteSubMatrixTo<V1, V2, V3, V4>(MsrMatrix target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColumnInidcesTarget)
            where V1 : IList<int>
            where V2 : IList<int>
            where V3 : IList<int>
            where V4 : IList<int> {

            target.Clear();
            this.AccSubMatrixTo(
                1.0,
                target,
                RowIndicesSource,
                RowIndicesTarget,
                ColumnIndicesSource,
                ColumnInidcesTarget);
        }

        /// <summary>
        /// Extracts a submatrix from this matrix and accumulates it to
        /// another matrix.
        /// </summary>
        /// <param name="alpha">
        /// scaling factor
        /// </param>
        /// <param name="RowIndicesSource">row indices into this matrix</param>
        /// <param name="RowIndicesTarget">
        /// if null is specified, this array is assumed to be
        /// {0,1, ... , <paramref name="RowIndicesSource"/>.Length-1]};
        /// </param>
        /// <param name="ColumnIndicesSource">
        /// column indices into this matrix
        /// </param>
        /// <param name="ColIndicesTarget">
        /// if null is specified, this array is assumed to be
        /// {0,1, ... , <paramref name="ColumnIndicesSource"/>.Length-1]};
        /// </param>
        /// <param name="Target">
        /// Output:
        /// The [<paramref name="RowIndicesTarget"/>[i],<paramref name="ColIndicesTarget"/>[j]]-th
        /// entry of the returned matrix is equal to the 
        /// [<paramref name="RowIndicesSource"/>[i],<paramref name="ColumnIndicesSource"/>[j]]-th
        /// entry of this matrix.
        /// </param>
        public void AccSubMatrixTo<V1, V2, V3, V4>(
            double alpha, MsrMatrix Target,
            V1 RowIndicesSource, V2 RowIndicesTarget,
            V3 ColumnIndicesSource, V4 ColIndicesTarget)
            where V1 : IList<int>
            where V2 : IList<int>
            where V3 : IList<int>
            where V4 : IList<int> {
            using(new FuncTrace()) {

                //if (RowIndicesTarget != null)
                //    if (RowIndicesSource.Count != RowIndicesTarget.Count)
                //        throw new ArgumentException("RowIndicesSource- and RowIndicesTarget - arrays must match in length;");
                //if (ColumnInidcesTarget != null)
                //    if (ColumnIndicesSource.Count != ColumnInidcesTarget.Count)
                //        throw new ArgumentException("ColumnIndicesSource- and ColumnInidcesTarget - arrays must match in length;");

                int L;
                if(RowIndicesSource == null && RowIndicesTarget == null) {
                    if(this.RowPartitioning.LocalLength != Target.RowPartitioning.LocalLength)
                        throw new ArgumentException("Mismatch in local number of rows.");

                    L = this.RowPartitioning.LocalLength;
                } else if(RowIndicesSource != null && RowIndicesTarget != null) {
                    if(RowIndicesSource.Count != RowIndicesTarget.Count)
                        throw new ArgumentException("Mismatch in length of source and target row index vector.");
                    L = RowIndicesSource.Count;
                } else if(RowIndicesSource != null) {
                    L = RowIndicesSource.Count;
                } else if(RowIndicesTarget != null) {
                    L = RowIndicesTarget.Count;
                } else {
                    throw new ApplicationException("should never occur.");
                }

                int Q;
                if(ColumnIndicesSource == null && ColIndicesTarget == null) {
                    if(this.ColPartition.LocalLength != Target.ColPartition.LocalLength)
                        throw new ArgumentException("Mismatch in local number of columns.");

                    Q = this.ColPartition.LocalLength;
                } else if(ColumnIndicesSource != null && ColIndicesTarget != null) {
                    if(ColumnIndicesSource.Count != ColIndicesTarget.Count)
                        throw new ArgumentException("Mismatch in length of source and target column index vector.");
                    Q = ColumnIndicesSource.Count;
                } else if(ColumnIndicesSource != null) {
                    Q = ColumnIndicesSource.Count;
                } else if(ColIndicesTarget != null) {
                    Q = ColIndicesTarget.Count;
                } else {
                    throw new ApplicationException("should never occur.");
                }


                int j0Src = this.m_ColPartitioning.i0;
                int j0Dst = Target.m_ColPartitioning.i0;
                int JlocSrc = this.m_ColPartitioning.LocalLength;
                int jESrc = this.m_ColPartitioning.iE;
                Debug.Assert(jESrc - j0Src == JlocSrc);
                int[] ColTrf  = null;

                if(ColumnIndicesSource != null || ColIndicesTarget != null) {
                    ColTrf = new int[this.m_ColPartitioning.TotalLength];

                    int[] LocColTrf = new int[this.m_ColPartitioning.LocalLength];
                    LocColTrf.SetAll(int.MinValue);

                    for(int j = 0; j < Q; j++) {
                        int jColSrc = ColumnIndicesSource != null ? ColumnIndicesSource[j] : (j + j0Src); // source row index
                        this.m_ColPartitioning.TestIfInLocalRange(jColSrc);
                        int jColDst = (ColIndicesTarget != null) ? ColIndicesTarget[j] : (j + j0Dst); // destination row index

                        LocColTrf[jColSrc - j0Src] = jColDst;
                    }

                    unsafe {
                        int[] i0s = this.m_ColPartitioning.GetI0s();
                        Debug.Assert(i0s.Length == this.m_ColPartitioning.MpiSize + 1);
                        Debug.Assert(i0s[0] == 0);
                        int[] LL = new int[i0s.Length - 1];
                        for(int i = 0; i < LL.Length; i++) {
                            LL[i] = i0s[i+1] - i0s[i];
                        }
                        
                        fixed(int* pColTrf = ColTrf, pLocColTrf = LocColTrf, pi0s = i0s, pLL = LL) {
                            csMPI.Raw.Allgatherv((IntPtr)pLocColTrf, LocColTrf.Length,
                                csMPI.Raw._DATATYPE.INT,
                                (IntPtr)pColTrf, (IntPtr)pLL, (IntPtr)pi0s,
                                csMPI.Raw._DATATYPE.INT,
                                this.m_ColPartitioning.MPI_Comm);
                        }

#if DEBUG
                        for(int j = j0Src; j < jESrc; j++) {
                            Debug.Assert(ColTrf[j] == LocColTrf[j - j0Src]);
                        }
#endif
                    }

                    // this is a poor-man solution,
                    // should be replaced when the Matrix data structure is renewed.

                }



                int i0Src = this.m_RowPartitioning.i0;
                int i0Dst = Target.m_RowPartitioning.i0;

                //int JE = ColumnIndicesSource != null ? ColumnIndicesSource.Count : this.ColPartition.TotalLength;

                for(int i = 0; i < L; i++) { // loop over rows ...
                    int iRowSrc = RowIndicesSource != null ? RowIndicesSource[i] : (i + i0Src); // source row index
                    this.m_RowPartitioning.TestIfInLocalRange(iRowSrc);
                    int irowLoc = iRowSrc - i0Src;
                    TrimAndSortRow(irowLoc);

                    //int jMin, jMax;
                    //GetRowRange(irowLoc, out jMin, out jMax);

                    int iRowDst = (RowIndicesTarget != null) ? RowIndicesTarget[i] : (i + i0Dst); // destination row index
                    

                    MsrMatrix.MatrixEntry[] Row = this.GetRowShallow(iRowSrc);
                    foreach(var Entry in Row) {
                        if(ColTrf == null) {
                            // no column transformation required
                            // +++++++++++++++++++++++++++++++++
                            Target[iRowDst, Entry.ColIndex] += alpha * Entry.Value;
                        } else {

                            int iColDest = ColTrf[Entry.m_ColIndex];
                            if(iColDest >= 0)
                                Target[iRowDst, iColDest] += alpha * Entry.Value;
                           
                        }

                    }



                    /*
                    for(int j = 0; j < JE; j++) {
                        int jColSrc = ColumnIndicesSource != null ? ColumnIndicesSource[j] : j;

                        if(jColSrc < jMin || jColSrc > jMax)
                            continue; // nothing to do for this row

                        double entry = this[iRowSrc, jColSrc];
                        if(entry == 0.0)
                            continue; // no need to set this entry

                        int jColDest;
                        if(ColumnInidcesTarget != null)
                            jColDest = ColumnInidcesTarget[j];
                        else
                            jColDest = j;


                        if(entry != 0.0)
                            Target[iRowDst, jColDest] += alpha * entry;
                    }
                     */
                }
            }
        }


        /// <summary>
        /// structure for inter-process data exchange, when <see cref="Transpose"/>
        /// is executed parallel.
        /// </summary>
        [Serializable]
        struct TransposeHelper {
            public int Row;
            public int Col;
            public double Val;
        }

        /// <summary>
        /// computes the transpose of this matrix;
        /// may involve a lot of MPI communication;
        /// </summary>
        public MsrMatrix Transpose() {

            MsrMatrix T = new MsrMatrix(this.ColPartition, this.RowPartitioning);

            int Size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out Size);


            SerialisationMessenger sms = null;
            SortedDictionary<int, List<TransposeHelper>> exchangeData = null;
            if (Size > 1) {
                sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                exchangeData = new SortedDictionary<int, List<TransposeHelper>>();
            }


            // local transpose
            // ---------------
            IPartitioning rowPart = this.ColPartition;

            int i0targ = (int)rowPart.i0;
            int iEndtarg = i0targ + rowPart.LocalLength;

            for (int i = 0; i < m_Entries.Length; i++) {
                MatrixEntry[] row = m_Entries[i];
                int iGlob = (int)(i + this.m_RowPartitioning.i0);

                if (row != null) {
                    foreach (MatrixEntry e in row) {
                        if (e.m_ColIndex >= 0) {

                            if (e.m_ColIndex >= i0targ && e.m_ColIndex < iEndtarg) {
                                // local transpose
                                T[e.ColIndex, iGlob] = e.Value;
                            } else {
                                // must send to another processor

                                int proc = rowPart.FindProcess(e.m_ColIndex);
                                if (!exchangeData.ContainsKey(proc))
                                    exchangeData.Add(proc, new List<TransposeHelper>());
                                List<TransposeHelper> data_proc = exchangeData[proc];

                                TransposeHelper th;
                                th.Col = iGlob;
                                th.Row = e.ColIndex;
                                th.Val = e.Value;
                                data_proc.Add(th);
                            }
                        }
                    }
                }
            }

            // communicate
            // ===========
            if (sms != null) {

                sms.SetCommPathsAndCommit(exchangeData.Keys);
                foreach (int proc in exchangeData.Keys)
                    sms.Transmitt(proc, exchangeData[proc].ToArray());

                TransposeHelper[] recvdat;
                int dummy;
                while (sms.GetNext(out dummy, out recvdat)) {
                    for (int k = 0; k < recvdat.Length; k++) {
                        T[recvdat[k].Row, recvdat[k].Col] = recvdat[k].Val;
                    }
                }
                sms.Dispose();
            }

            return T;
        }

        /*
        /// <summary>
        /// The Infinity-Norm (maximum absolute row sum norm) of this matrix;
        /// </summary>
        /// <returns></returns>
        public double InfNorm() {
            using (var tr = new FuncTrace()) {
                double normLoc = 0;

                for (int i = 0; i < m_Entries.Length; i++) {
                    double rownrm = 0;
                    MatrixEntry[] row = m_Entries[i];
                    if (row != null) {
                        for (int j = 0; j < row.Length; j++) {
                            if (row[j].m_ColIndex >= 0) {
                                rownrm += Math.Abs(row[j].Value);
                            }
                        }
                    }

                    normLoc = Math.Max(normLoc, rownrm);
                }

                //Console.WriteLine("local norm" + normLoc);
                tr.Info("local norm " + normLoc);

                double normGlob = double.NaN;
                unsafe {
                    csMPI.Raw.Allreduce((IntPtr)(&normLoc), (IntPtr)(&normGlob), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.MAX, csMPI.Raw._COMM.WORLD);
                }
                return normGlob;
            }
        }
        */
        /*
        /// <summary>
        /// Determines the minimum and maximum column index for a specific row.
        /// </summary>
        /// <param name="iRow">
        /// the row index (local index (!)) which should be investigated
        /// </param>
        /// <param name="jMin">
        /// the minimum of all column indices (global index) of nonzero
        /// elements in row <paramref name="iRow"/>
        /// </param>
        /// <param name="jMax">
        /// the maximum of all column indices (global index)of nonzero elements
        /// in row <paramref name="iRow"/>
        /// </param>
        private void GetRowRange(int iRow, out int jMin, out int jMax) {
            MatrixEntry[] row = m_Entries[iRow];

            jMin = int.MaxValue;
            jMax = int.MinValue;

            if (row != null) {
                foreach (MatrixEntry e in row) {
                    int j = e.m_ColIndex;
                    if (j < 0)
                        continue;

                    jMin = Math.Min(jMin, j);
                    jMax = Math.Max(jMax, j);
                }
            }
        }
        */

        /// <summary>
        /// Computes the deviation of this matrix from symmetry
        /// </summary>
        /// <returns>
        /// The accumulated sum of the differences between corresponding
        /// off-diagonal entries.
        /// </returns>
        public double SymmetryDeviation() {
            int size;
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
            if (size != 1)
                throw new NotImplementedException("currently, this method does not work mpi parallel.");

            if (NoOfRows != NoOfCols)
                throw new ApplicationException("symmetry cannot be determined for non-quadratic matrices");

            // test for symmetry:
            double delta = 0;
            int i0 = (int)m_RowPartitioning.i0;
            //List<double> delta_vec = new List<double>();
            for (int i = 0; i < m_Entries.Length; i++) {
                MatrixEntry[] row = m_Entries[i];

                int iRow = i0 + i;

                if (row != null) {
                    foreach (MatrixEntry e in row) {
                        if (e.ColIndex == iRow)
                            continue;
                        if (e.ColIndex < 0)
                            continue;

                        double diffe = e.Value - this[e.ColIndex, i];

                        delta += diffe * diffe;
                    }
                }
            }

            return delta;
        }

        /// <summary>
        /// Performs the operation: this = this + <paramref name="Ascale"/>*<paramref name="A"/>;
        /// </summary>
        /// <param name="A">
        /// another matrix with same size and equal <see cref="RowPartitioning"/>;
        /// </param>
        /// <param name="Ascale">
        /// scaling
        /// </param>
        public void Acc(MsrMatrix A, double Ascale) {

            if (A.NoOfCols != this.NoOfCols)
                throw new ArgumentException("other and this matrix must have the same number of columns;", "A");
            if (!A.RowPartitioning.Equals(this.RowPartitioning))
                throw new ArgumentException("other and this matrix must have an equal row partition;", "A");

            if (A.m_Entries.Length != this.m_Entries.Length)
                // should never happen
                throw new ApplicationException("internal error; inconsistency between row partition and m_Entries");

            int i0 = (int)m_RowPartitioning.i0;
            for (int i = m_Entries.Length - 1; i >= 0; i--) {
                A.TrimAndSortRow(i);
                MatrixEntry[] ARow = A.m_Entries[i];
                int irow = i0 + i;
                if (ARow != null) {
                    foreach (MatrixEntry e in ARow) {
                        this[irow, e.m_ColIndex] = this[irow, e.m_ColIndex] + Ascale * e.Value;
                    }
                }
            }
        }

        /// <summary>
        /// multiplies this matrix by the factor <paramref name="a"/>
        /// </summary>
        /// <param name="a"></param>
        public void Scale(double a) {
            if (a == 0.0) {
                this.Clear();
            } else {

                foreach (MatrixEntry[] row in m_Entries) {
                    if (row != null) {
                        for (int j = row.Length - 1; j >= 0; j--) {

                            if (row[j].m_ColIndex >= 0)
                                row[j].Value *= a;
                        }
                    }
                }
            }
        }

        /*
        /// <summary>
        /// greatest common divider
        /// </summary>
        private static int GCD(int a, int b) {
            if (b == 0)
                return a;
            else
                return GCD(b, a % b);
        }
        */

        /// <summary>
        /// matrix-matrix-multiplication, see <see cref="Multiply"/>.
        /// </summary>
        static public MsrMatrix operator *(MsrMatrix left, MsrMatrix right) {
            return Multiply(left, right);
        }

        /// <summary>
        /// Multiply matrix <paramref name="M"/> by a scalar <paramref name="a"/>.
        /// </summary>
        static public MsrMatrix operator *(double a, MsrMatrix M) {
            MsrMatrix R = M.CloneAs();
            R.Scale(a);
            return R;
        }

        /// <summary>
        /// Multiply matrix <paramref name="M"/> by a scalar <paramref name="a"/>.
        /// </summary>
        static public MsrMatrix operator *(MsrMatrix M, double a) {
            MsrMatrix R = M.CloneAs();
            R.Scale(a);
            return R;
        }

        /// <summary>
        /// summation of matrices, see <see cref="Acc"/>.
        /// </summary>
        static public MsrMatrix operator +(MsrMatrix left, MsrMatrix right) {
            var ReturnMatrix = left.CloneAs();
            ReturnMatrix.Acc(1.0, right);
            return ReturnMatrix;
        }

        /// <summary>
        /// Difference of matrices, see <see cref="Acc"/>.
        /// </summary>
        static public MsrMatrix operator -(MsrMatrix left, MsrMatrix right) {
            var ReturnMatrix = left.CloneAs();
            ReturnMatrix.Acc(-1.0, right);
            return ReturnMatrix;
        }



        /// <summary>
        /// multiplies two matrices; this method is both, memory and
        /// computationally intensive;
        /// </summary>
        /// <param name="left">1st operand</param>
        /// <param name="right">2nd operand</param>
        /// <returns>the matrix <paramref name="left"/>*<paramref name="right"/></returns>
        static public MsrMatrix Multiply(MsrMatrix left, MsrMatrix right) {
            using (new FuncTrace()) {
                if (left.NoOfCols != right.NoOfRows)
                    throw new ArgumentException("number of columns of left matrix must match number of rows of right matrix");

                int size, rank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                //MsrMatrix Target = null;


                int i0 = (int)left.RowPartitioning.i0;
                int L = left.RowPartitioning.LocalLength;
                int j0 = (int)right.RowPartitioning.i0;
                int K = right.RowPartitioning.LocalLength;

                //int resC0 = (int)Target.ColPartition.i0;
                //int resCL = (int)Target.ColPartition.LocalLength;

                // local multiply, create send data
                // ================================

                // data packet sent to each processor
                List<MultiplyHelper>[] sendData = new List<MultiplyHelper>[size]; 


                int MaxLen = 0;
                for (int l = 0; l < K; l++) {
                    int ll = right.TrimAndSortRow(l);
                    MaxLen = Math.Max(MaxLen, ll);
                }

                MsrMatrix res = new MsrMatrix(
                    new Partitioning(left.RowPartitioning.LocalLength),
                    new Partitioning(right.ColPartition.LocalLength),
                    (int)Math.Ceiling(MaxLen * 4.0));

                unsafe {
                    int CurBufLen = MaxLen * 100;
                    MatrixEntry* ResltBuffer = (MatrixEntry*)Marshal.AllocHGlobal(sizeof(MatrixEntry) * CurBufLen);
                    MatrixEntry* MergeBuffer = (MatrixEntry*)Marshal.AllocHGlobal(sizeof(MatrixEntry) * CurBufLen);


                    MatrixEntry* pResltBuffer = ResltBuffer, pMergeBuffer = MergeBuffer;
                    for (int k = CurBufLen; k > 0; k--) {
                        pResltBuffer->m_ColIndex = -1;
                        pMergeBuffer->m_ColIndex = -1;

                        pMergeBuffer++;
                        pResltBuffer++;
                    }


                    {
                        // We want to compute 
                        //   res = left*right
                        // we follow the same strategy as Tim Davis in his ssmul code (part of suitesparse).
                        // Unfortunately, his code is not MPI parallel, so we need our own implementation.
                        //
                        // We compute the i-th row of 'res' as 
                        //   res[i,-]  =  sum_{over all 'k' for which left[i,k] != 0}  left[i,k]*right[k,-]
                        //
                        // loop 'i', over all rows of 'left'...
                        //   C[i,-] = 0
                        //   loop over all non-zero columns 'k' in row 'left[i,-]'...
                        //     val = left[i,k]
                        //     C[i,-] += val*right[k,-]
                        //
                        // The performance-critical part is indeed the '+=' operation. Therefor, we have
                        // to use a merge sort. ('USAXPY' Unsorted? Sparse AXPY (alpha times X plus Y))



                        for (int i = 0; i < L; i++) { // loop over left matrix rows...
                            int ResultItems = 0;

                            MatrixEntry[] rowA = left.m_Entries[i];
                            int R = rowA != null ? rowA.Length : 0;
                            if (R <= 0)
                                continue;
                            fixed (MatrixEntry* pRowA = rowA) {

                                for (MatrixEntry* pe = pRowA; R > 0 && pe->m_ColIndex >= 0; pe++) { // loop over all entries of row left[i,:] (row i of left matrix)
                                    R--; // elements to do: decrement 

                                    double val = pe->Value;
                                    
                                    int iColLoc = pe->ColIndex - j0;
                                    if (iColLoc >= 0 && iColLoc < K) {
                                        // row of Right matrix stored on this proc.
                                        // ++++++++++++++++++++++++++++++++++++++++
                                        MatrixEntry[] rowB = right.m_Entries[iColLoc];
                                        int C = rowB != null ? rowB.Length : 0;
                                        if (ResultItems + C > CurBufLen) {
                                            CurBufLen += MaxLen * 2;

                                            // re-alloc memory
                                            MergeBuffer = (MatrixEntry*)Marshal.ReAllocHGlobal((IntPtr)MergeBuffer, (IntPtr)(sizeof(MatrixEntry) * CurBufLen));
                                            ResltBuffer = (MatrixEntry*)Marshal.ReAllocHGlobal((IntPtr)ResltBuffer, (IntPtr)(sizeof(MatrixEntry) * CurBufLen));
                                        }

                                        // let's do USAXPY!!
                                        pResltBuffer = ResltBuffer; // the 'Reslt' buffer contains all entries of res[i,-] computed so far...
                                        pMergeBuffer = MergeBuffer; // the 'Merge' buffer will be used to merge and accumulate res[i,-] += val*right[iColLoc,-] 
                                        //                             later, the buffers will be swapped.
                                        int MergeItems = 0;
                                        Debug.Assert(pMergeBuffer - MergeBuffer == MergeItems);

                                        if (C > 0) {
                                            fixed (MatrixEntry* pRowB = rowB) {
                                                for (MatrixEntry* p_ee = pRowB; C > 0 && p_ee->m_ColIndex >= 0; p_ee++) { // loop over all columns in row 'right[iColLoc,-]'...
                                                    C--;

                                                    //
                                                    // this is the merge loop
                                                    //

                                                    int p_ee_col = p_ee->m_ColIndex;

                                                    while (pResltBuffer->m_ColIndex < p_ee_col && ResultItems > 0) {
                                                        // 
                                                        *pMergeBuffer = *pResltBuffer;
                                                        pMergeBuffer++;
                                                        pResltBuffer++;
                                                        ResultItems--;
                                                        MergeItems++;

                                                        Debug.Assert(pMergeBuffer - MergeBuffer == MergeItems);
                                                        Debug.Assert(MergeItems <= CurBufLen);
                                                    }

                                                    if (p_ee_col == pResltBuffer->m_ColIndex && ResultItems > 0) {
                                                        pMergeBuffer->m_ColIndex = p_ee_col;
                                                        pMergeBuffer->Value = pResltBuffer->Value + p_ee->Value * val;
                                                        ResultItems--;
                                                        MergeItems++;
                                                        pMergeBuffer++;
                                                        pResltBuffer++;

                                                        Debug.Assert(pMergeBuffer - MergeBuffer == MergeItems);
                                                        Debug.Assert(MergeItems <= CurBufLen);
                                                    } else {
                                                        *pMergeBuffer = *p_ee;
                                                        pMergeBuffer->Value *= val;
                                                        pMergeBuffer++;
                                                        MergeItems++;

                                                        Debug.Assert(pMergeBuffer - MergeBuffer == MergeItems);
                                                        Debug.Assert(MergeItems <= CurBufLen);
                                                    }
                                                }
                                            }
                                        }

                                        while (ResultItems > 0) {
                                            *pMergeBuffer = *pResltBuffer;
                                            pMergeBuffer++;
                                            pResltBuffer++;
                                            ResultItems--;
                                            MergeItems++;

                                            Debug.Assert(MergeItems <= CurBufLen);
                                        }

                                        Debug.Assert(pMergeBuffer - MergeBuffer == MergeItems);
                                        Debug.Assert(MergeItems <= CurBufLen);

                                        // buffer swap
                                        {
                                            ResultItems = MergeItems;
                                            MergeItems = 0;
                                            MatrixEntry* bbb = MergeBuffer;
                                            MergeBuffer = ResltBuffer;
                                            ResltBuffer = bbb;
                                        }


                                    } else {
                                        
                                        int owner_rank = right.RowPartitioning.FindProcess(pe->ColIndex);
                                        // row of Right matrix stored on proc. 'owner_rank'
                                        // ++++++++++++++++++++++++++++++++++++++++++++++++

                                        if (sendData[owner_rank] == null) {
                                            sendData[owner_rank] = new List<MultiplyHelper>();
                                        }

                                        MultiplyHelper pkt;
                                        pkt.i = i;
                                        pkt.j = pe->ColIndex;
                                        pkt.val = pe->Value;

                                        sendData[owner_rank].Add(pkt);
                                         
                                    }
                                }

                            }

                            var Res_iRow = res.m_Entries[i];
                            if (ResultItems > 0 && (Res_iRow == null || Res_iRow.Length < ResultItems)) {
                                //Array.Resize(ref Target.m_Entries[i], ResultItems);
                                res.m_Entries[i] = new MatrixEntry[ResultItems];
                                Res_iRow = res.m_Entries[i];
                            }

                            if (ResultItems > 0) {
                                fixed (MatrixEntry* pRes_iRow = Res_iRow)
                                {
                                    pResltBuffer = ResltBuffer;
                                    MatrixEntry* pdest = pRes_iRow;
                                    for (int j = 0; j < ResultItems; j++) {
                                        *pdest = *pResltBuffer;
                                        pResltBuffer++;
                                        pdest++;
                                    }
                                }
                            }
                        }
                    }


                    Marshal.FreeHGlobal((IntPtr)ResltBuffer);
                    Marshal.FreeHGlobal((IntPtr)MergeBuffer);
                }

#if DEBUG
                res.VerifyDataStructure();
#endif


                // send/receive data (question) ...
                // ================================

                MultiplyHelper[][] data_received = new MultiplyHelper[size][];
                {

                    SerialisationMessenger sms = new SerialisationMessenger(res.RowPartitioning.MPI_Comm);
                    for (int p = 0; p < size; p++) {
                        if (sendData[p] != null)
                            sms.SetCommPath(p);
                    }
                    sms.CommitCommPaths();


                    for (int p = 0; p < size; p++) {
                        if (sendData[p] != null) {
                            sms.Transmitt(p, sendData[p].ToArray());
                            sendData[p] = null;
                        }
                    }

                    int tp;
                    MultiplyHelper[] r;
                    while (sms.GetNext(out tp, out r)) {
                        data_received[tp] = r;
                    }
                }

                // external multiply
                // =================

                // loop over processors where data was received from ...
                for (int p = 0; p < size; p++) {
                    MultiplyHelper[] r = data_received[p];
                    if (r == null)
                        // no data received from process 'p' -- nothing to do for process 'p'
                        continue;

                    if (sendData[p] == null)
                        sendData[p] = new List<MultiplyHelper>();

                    List<MultiplyHelper> s = sendData[p];

                    foreach (var pkt in r) {
                        MatrixEntry[] row = right.m_Entries[pkt.j - j0];

                        if (row != null) {
                            foreach (var ee in row) {
                                if (ee.ColIndex < 0)
                                    break;

                                MultiplyHelper pkt_ret;
                                pkt_ret.i = pkt.i;
                                pkt_ret.j = ee.ColIndex;
                                pkt_ret.val = pkt.val * ee.Value;
                                s.Add(pkt_ret);
                            }
                        }
                    }
                }

                // send/receive data (answer) ...
                // ==============================

                data_received = new MultiplyHelper[size][];
                {

                    SerialisationMessenger sms = new SerialisationMessenger(res.RowPartitioning.MPI_Comm);
                    for (int p = 0; p < size; p++) {
                        if (sendData[p] != null)
                            sms.SetCommPath(p);
                    }
                    sms.CommitCommPaths();

                    for (int p = 0; p < size; p++) {
                        if (sendData[p] != null) {
                            sms.Transmitt(p, sendData[p].ToArray());
                            sendData[p] = null;
                        }
                    }

                    int tp;
                    MultiplyHelper[] r;
                    while (sms.GetNext(out tp, out r)) {
                        data_received[tp] = r;
                    }
                }

                // accumulate external 
                // ===================
                for (int p = 0; p < size; p++) {
                    MultiplyHelper[] r = data_received[p];
                    if (r == null)
                        // no data received from process 'p' -- nothing to do for process 'p'
                        continue;

                    foreach (var v in r) {
                        res[v.i + i0, v.j] += v.val;
                    }
                }


                // return result
                // =============

                return res;
            }
        }

        /// <summary>
        /// auxiliary data struct used in <see cref="Multiply"/>;
        /// </summary>
        [Serializable]
        struct MultiplyHelper {
            public int i;
            public int j;
            public double val;
        }

        /// <summary>
        /// total number of columns (over all MPI processors)
        /// </summary>
        public int NoOfCols {
            get {
                return (int)m_ColPartitioning.TotalLength;
                //return m_Size; 
            }
        }

        /// <summary>
        /// total number of rows (over all MPI processors)
        /// </summary>
        public int NoOfRows {
            get {
                long no = m_RowPartitioning.TotalLength;
                if (no > int.MaxValue)
                    throw new ApplicationException("64/32 bit index overflow");
                return (int)no;
            }
        }

        /// <summary>
        /// structure to store matrix during matrix assembly;
        /// <see cref="m_Entries"/>.
        /// </summary>
        /// <remarks>
        /// Note that this structure does <b>not</b> store the row it is
        /// associated with. 
        /// </remarks>
        [Serializable]
        public struct MatrixEntry {

            /// <summary>
            /// <see cref="ColIndex"/>;
            /// Attention: Setting this value may put the  
            /// owning <see cref="MsrMatrix"/> into a undefined state.
            /// </summary>
            public int m_ColIndex;

            /// <summary>
            /// column index (global coordinates) of this entry; negative value
            /// indicates a not-assigned entry;
            /// </summary>
            public int ColIndex {
                get {
                    return m_ColIndex;
                }
            }

            /// <summary>
            /// value of this entry
            /// </summary>
            public double Value;

            /// <summary>
            /// for sorting according to column index <see cref="m_ColIndex"/>;
            /// </summary>
            /// <param name="a"></param>
            /// <param name="b"></param>
            /// <returns></returns>
            static internal int Compare(MatrixEntry a, MatrixEntry b) {
                return a.m_ColIndex - b.m_ColIndex;
            }
        }

        /// <summary>
        /// matrix entries;
        /// 1st index: local row index; 2nd index: no interpretation;
        /// A row may contain unset elements which is indicated by a negative
        /// column index (<see cref="MatrixEntry.ColIndex"/>); All entries with
        /// valid column index are at the beginning of a row  array, i.e. after
        /// the first invalid entry in a row array there are no other valid
        /// entries.
        /// </summary>
        public MatrixEntry[][] m_Entries;

        /// <summary>
        /// returns row number <paramref name="i"/>;
        /// the row contains only valid entries and it is sorted in ascending order
        /// according to column index;
        /// </summary>
        /// <param name="i">row index in global indices</param>
        /// <returns>
        /// the return value is a reference to an internal data structure, not a copy;
        /// So, modifying the matrix (by operations like <see cref="ClearRow"/>) will also affect the 
        /// returned array; If this behavior is not desired, the returned array must be cloned;
        /// </returns>
        public MatrixEntry[] GetRowShallow(int i) {

            //int idiagcol = i;
            i -= (int)m_RowPartitioning.i0;

            TrimAndSortRow(i);
            if (m_Entries[i] == null)
                return new MatrixEntry[0];
            else
                return m_Entries[i];
        }

        /// <summary>
        /// sets row number <paramref name="i"/>;
        /// All previous entries in this row are overwritten;
        /// </summary>
        /// <param name="i"> row index in global indices</param>
        /// <param name="row">
        /// </param>
        public void SetRow(int i, params MatrixEntry[] row) {
            // check arguments
            i -= (int)m_RowPartitioning.i0;
            if (i < 0 || i >= m_RowPartitioning.LocalLength)
                throw new IndexOutOfRangeException("invalid row index");

            int J = this.NoOfCols;
            foreach (MatrixEntry e in row)
                if (e.m_ColIndex < 0 || e.m_ColIndex >= J)
                    throw new IndexOutOfRangeException("invalid column index");


            m_Entries[i] = row;
            TrimAndSortRow(i);
        }

        /// <summary>
        /// sets row number <paramref name="i"/>;
        /// All previous entries in this row are overwritten;
        /// </summary>
        /// <param name="i"> row index in global indices</param>
        /// <param name="col">
        /// Column indices.
        /// </param>
        /// <param name="val">
        /// Matrix entries.
        /// </param>
        /// <param name="L">
        /// Number of entries which are actually used in <paramref name="col"/> rep. <paramref name="val"/>.
        /// If negative, the length of <paramref name="col"/>.
        /// </param>
        public void SetRow(int i, int[] col, double[] val, int L = -1) {
            // check arguments
            i -= (int)m_RowPartitioning.i0;
            if (i < 0 || i >= m_RowPartitioning.LocalLength)
                throw new IndexOutOfRangeException("invalid row index");

            if(L < 0) {
                if (col.Length != val.Length)
                    throw new ArgumentException();
                L = col.Length;
            }

            if(m_Entries[i].Length < L) {
                m_Entries[i] = new MatrixEntry[L];
            }
            var row = m_Entries[i];

            int J = this.NoOfCols;
            for (int l = 0; l < L; l++) {
                if (col[l] < 0 || col[l] >= J)
                    throw new IndexOutOfRangeException("invalid column index");
                row[l].m_ColIndex = col[l];
                row[l].Value = val[l];
            }
            for (int l = L; l < row.Length; l++) {
                row[l].m_ColIndex = -1;
            }

            TrimAndSortRow(i);
        }

        /// <summary>
        /// Sets all entries in a row to 0
        /// </summary>
        /// <param name="i">row index in global indices</param>
        public void ClearRow(int i) {
            i -= (int)m_RowPartitioning.i0;

            MatrixEntry[] row = m_Entries[i];
            if (row != null) {
                for (int j = row.Length - 1; j >= 0; j--) {
                    row[j].m_ColIndex = -1;
                }
            }
        }

        /// <summary>
        /// releases unused entries in row <paramref name="i"/> and sorts all
        /// entries in this row in ascending order (with respect to their
        /// column index).
        /// </summary>
        /// <param name="i">row index in local indices (!)</param>
        private int TrimAndSortRow(int i) {
            MatrixEntry[] row = m_Entries[i];

            if (row == null)
                return 0;

            // trim the row
            int len = 0;
            while (len < row.Length && row[len].ColIndex >= 0) {
                len++;
            }
            if (len != row.Length) {
                Array.Resize<MatrixEntry>(ref row, len);
                m_Entries[i] = row;
            }

            // sort the row
            Array.Sort<MatrixEntry>(row, MatrixEntry.Compare);

            // return length;
            return len;
        }
        
        /// <summary>
        /// 
        /// </summary>
        /// <param name="i">local row index</param>
        /// <param name="j">local column index</param>
        /// <param name="allocate">if true, a new entry is allocated if there is
        /// actually no entry set for the specified column;</param>
        /// <returns>
        /// the 2nd index into the <see cref="m_Entries"/>-field;
        /// negative value if no entry exists at position <paramref name="i"/>,<paramref name="j"/>;
        /// </returns>
        int GetNewEntries2ndIndex(int i, int j, bool allocate) {
            i -= (int)m_RowPartitioning.i0;

            if (i < 0 || i >= m_Entries.Length)
                throw new IndexOutOfRangeException("row index out of range");
            if (j < 0 || j >= this.NoOfCols)
                throw new IndexOutOfRangeException("column index out of range");

            int r = 0;
            MatrixEntry[] row = m_Entries[i];

            
            if (row == null) {
                if (!allocate) {
                    // no entry found
                    return -1;
                } else {
                    row = new MatrixEntry[5];
                    m_Entries[i] = row;
                    for (int l = 1; l < row.Length; l++)
                        row[l].m_ColIndex = -1;
                    row[0].m_ColIndex = j;
                    return 0;
                }
            }


            for (r = 0; r < row.Length; r++) {
                if (row[r].m_ColIndex == j)
                    return r;
                if (row[r].m_ColIndex < 0)
                    break;
            }

            // no entry found
            if (!allocate)
                return -1;

            if (r >= row.Length) {
                int lold = row.Length;
                Array.Resize<MatrixEntry>(ref m_Entries[i], row.Length + 5);
                row = m_Entries[i];
                for (; lold < row.Length; lold++)
                    row[lold].m_ColIndex = -1;
            }
            row[r].m_ColIndex = j;
            return r;
        }

        //int last_i = -1;
        //int last_ind = -1;

        /// <summary>
        /// get/set an entry; Setting a zero element (to a value unequal to zero)
        /// allocates a new entry; Setting an entry to 0.0 releases the
        /// corresponding memory;
        /// </summary>
        /// <param name="i">global (over all MPI processes) row index</param>
        /// <param name="j">
        /// global (over all MPI processes) column index
        /// </param>
        /// <returns></returns>
        public double this[int i, int j] {
            get {
                long iLoc = i - m_RowPartitioning.i0;
                //if(i == last_i && m_Entries[i][last_ind].m_ColIndex == j) {
                //    return m_Entries[iLoc][last_ind].Value;
                //}
                
                // matrix contents are stored in m_NewEntries
                int ind = GetNewEntries2ndIndex(i, j, false);
                if(ind < 0) {
                    return 0.0;
                } else {
                    //last_i = i;
                    //last_ind = ind;
                    return m_Entries[iLoc][ind].Value;
                }
            }
            set {
                long iLoc = i - m_RowPartitioning.i0;
                if (value == 0.0) {
                    int ind = GetNewEntries2ndIndex(i, j, false);
                    if (ind >= 0) {
                        ClearEntry(m_Entries[iLoc], ind);
                    }
                    //last_i = -1;
                    //last_ind = -1;
                } else {
                    //if(i == last_i && m_Entries[i][last_ind].m_ColIndex == j) {
                    //    m_Entries[iLoc][last_ind].Value = value;
                    //}

                    int ind = GetNewEntries2ndIndex(i, j, true);
                    m_Entries[iLoc][ind].Value = value;
                    //last_i = i;
                    //last_ind = ind;
                }
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
        /// <param name="beta">Scaling applied to this matrix before accumulation</param>
        public void AccBlock(int i0, int j0, double alpha, MultidimensionalArray Block, double beta) {
            if (Block.Dimension != 2)
                throw new ArgumentException();
            int I = Block.NoOfRows;
            int J = Block.NoOfCols;

            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                    this[i0 + i, j0 + j] = this[i0 + i, j0 + j] * beta + alpha * Block[i, j];
        }

        /// <summary>
        /// Allocates memory for some entry, even if the entry is zero;
        /// </summary>
        /// <param name="i">global (over all MPI processes) row index</param>
        /// <param name="j">global (over all MPI processes) column index</param>
        public void AllocateEvenWhenZero(int i, int j) {
            GetNewEntries2ndIndex(i, j, true);
        }

        /// <summary>
        /// true, if this matrix and <paramref name="other"/> matrix have the 
        /// zeros/nonzeros at the same indices;
        /// </summary>
        /// <param name="other"></param>
        /// <returns> </returns>
        /// <remarks>
        /// the two matrices must have the same partition to be comparable
        /// </remarks>
        public bool FillingEquals(MsrMatrix other) {
            if (!this.RowPartitioning.Equals(other.RowPartitioning))
                throw new ArgumentException("unequal row partition - cannot compare;", "other");

            int locRes = 1;

            int L = this.m_Entries.Length;
            int i0 = (int)this.RowPartitioning.i0;
            for (int i = 0; i < L; i++) {
                MatrixEntry[] rowthis = this.GetRowShallow(i0 + i);
                MatrixEntry[] rowother = this.GetRowShallow(i0 + i);

                if ((rowthis == null) && (rowother == null)) {
                    continue;
                }

                if ((rowthis != null) != (rowother != null)) {
                    locRes = 0;
                    break;
                }
                
                if (rowthis.Length != rowother.Length) {
                    locRes = 0;
                    break;
                }

                int K = rowother.Length;
                for (int k = 0; k < K; k++) {
                    if (rowthis[k].ColIndex != rowother[k].ColIndex) {
                        locRes = 0;
                        break;
                    }
                }

                if (locRes == 0)
                    break;
            }


            unsafe {
                int globRes = -1;
                csMPI.Raw.Allreduce((IntPtr)(&locRes), (IntPtr)(&globRes), 1, csMPI.Raw._DATATYPE.INT, csMPI.Raw._OP.MIN, csMPI.Raw._COMM.WORLD);

                return (locRes > 0);
            }
        }

        #region IMutuableMatrix Members

        /// <summary>
        /// see <see cref="IMutableMatrix.GetValues"/>;
        /// </summary>
        public double[] GetValues(int RowIndex, int[] ColumnIndices) {
            double[] ret = new double[ColumnIndices.Length];
            for (int i = 0; i < ret.Length; i++) {
                ret[i] = this[RowIndex, ColumnIndices[i]];
            }
            return ret;
        }

         

        /// <summary>
        /// see <see cref="IMutableMatrix.SetValues"/>;
        /// </summary>
        public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
            if(newValues.Length != ColumnIndices.Length) {
                throw new ArgumentException("Mismatch in array length; column index array and value array must have the same length.");
            }
            int L = ColumnIndices.Length;
            MatrixEntry[] newRow = new MatrixEntry[L];
            this.RowPartitioning.TestIfInLocalRange(RowIndex);
            int NoCols = this.ColPartition.TotalLength;
            for(int i = 0; i < L; i++) {
                if(ColumnIndices[i] < 0 || ColumnIndices[i] >= NoCols)
                    throw new ArgumentOutOfRangeException("Column index out of range.");
                newRow[i].m_ColIndex = ColumnIndices[i];
                newRow[i].Value = newValues[i];

            }
            this.m_Entries[RowIndex - this.RowPartitioning.i0] = newRow;
        }

        /// <summary>
        /// always true, see <see cref="IMutableMatrix.OccupationMutable"/>;
        /// </summary>
        public bool OccupationMutable {
            get {
                return true;
            }
        }

        #endregion

        #region ISparseMatrix Members

        /// <summary>
        /// see <see cref="ISparseMatrix.GetDiagonalElement"/>;
        /// </summary>
        public double GetDiagonalElement(int row) {
            return this[row, row];
        }

        /// <summary>
        /// see <see cref="ISparseMatrix.SetDiagonalElement"/>;
        /// </summary>
        public void SetDiagonalElement(int row, double val) {
            this[row, row] = val;
        }

        #endregion

        #region IMutuableMatrixEx Members

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetOccupiedColumnIndices"/>
        /// </summary>
        public int GetOccupiedColumnIndices(int RowIndex, ref int[] ColumnIndices) {
            TrimAndSortRow(m_RowPartitioning.TransformIndexToLocal(RowIndex));

            MatrixEntry[] row = m_Entries[m_RowPartitioning.TransformIndexToLocal(RowIndex)];
            int nnz = 0;
            if (row != null) {
                while (nnz < row.Length && row[nnz].ColIndex >= 0)
                    nnz++;
            }

            if (ColumnIndices == null || ColumnIndices.Length < nnz)
                ColumnIndices = new int[nnz];

            int i;
            for (i = 0; i < nnz; i++) {
                Debug.Assert(row[i].ColIndex >= 0);
                ColumnIndices[i] = row[i].ColIndex;
            }
#if DEBUG
            if (row != null) {
                for (; i < row.Length; i++) {
                    Debug.Assert(row[i].ColIndex < 0);
                }
            }
#endif
            return nnz;
        }

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetOccupiedColumnIndices"/>;
        /// in contrast to <see cref="GetRowShallow"/>, the returned array is a non-shallow copy of the row;
        /// </summary>
        public int GetRow(int RowIndex, ref int[] ColumnIndices, ref double[] Values) {
            //MatrixEntry[] row = GetRowShallow(RowIndex);
            //return (MatrixEntry[])row.Clone();
            TrimAndSortRow(m_RowPartitioning.TransformIndexToLocal(RowIndex));

            MatrixEntry[] row = m_Entries[m_RowPartitioning.TransformIndexToLocal(RowIndex)];
            int nnz = 0;
            if (row != null) {
                while (nnz < row.Length && row[nnz].ColIndex >= 0)
                    nnz++;
            }


            if (ColumnIndices == null || ColumnIndices.Length < nnz)
                ColumnIndices = new int[nnz];
            if (Values == null || Values.Length < nnz)
                Values = new double[nnz];

            int i;
            for (i = 0; i < nnz; i++) {
                Debug.Assert(row[i].ColIndex >= 0);
                ColumnIndices[i] = row[i].ColIndex;
                Values[i] = row[i].Value;
            }
#if DEBUG
            if (row != null) {
                for (; i < row.Length; i++) {
                    Debug.Assert(row[i].ColIndex < 0);
                }
            }
#endif
            return nnz;

        }

        #endregion


        /// <summary>
        /// checks the integrity of the data structure
        /// </summary>
        public void VerifyDataStructure() {
            using (new FuncTrace()) {
                int L = this.m_Entries.Length;

                for (int i = 0; i < L; i++) {
                    this.TrimAndSortRow(i);
                    var row = this.m_Entries[i];
                    if (row != null) {
                        for (int j = 0; j < row.Length - 1; j++) {
                            if (row[j].ColIndex == row[j + 1].ColIndex)
                                throw new ApplicationException("Error in MSR data structure.");

                        }
                    }
                }
            }
        }

        /*
        /// <summary>
        /// Sets Diagonal Element to 1 
        /// if the whole line of the matrix only contains zeros
        /// This is only valid if 
        /// the corresponding vector <paramref name="rhs"/>
        /// is zero in this line.
        /// </summary>
        /// <param name="rhs">
        /// Vector used for checks
        /// </param>
        public void FillDiagonalWithOnes(IList<double> rhs = null) {
            int i0 = this.RowPartitioning.i0;
            int L = this.RowPartitioning.LocalLength;
            for (int l = 0; l < L; l++) {
                if (this.GetNoOfNonZeros(l + i0) == 0) {
                    if (rhs!=null && rhs[l] != 0.0)
                        throw new ArithmeticException("No solution for given RHS.");
                    this.SetDiagonalElement(l + i0, 1.0);
                }
            }
        }
        */
    }
}
