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
using System.Collections;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Utils;
using System.Linq;

namespace ilPSP.LinSolvers {
        
    /// <summary>
    /// A helper to execute Row and Col operations (scaling, addition, swapping) 
    /// in parallel on an <see cref="MsrMatrix"/>.
    /// </summary>
    public class MatrixOp {
        
        /// <summary>
        /// constructor.
        /// </summary>
        /// <param name="M">matrix to manipulate</param>
        /// <param name="ExtCol">
        /// key: processor rank 'p' <br/>
        /// value: a list of column indices (within the local range of columns of <paramref name="M"/>),
        /// which should be editable at rank 'p'.
        /// </param>
        public MatrixOp(IMutableMatrixEx M, IDictionary<int,int[]> ExtCol = null) {
            m_Matrix = new MsrExtMatrix(M, ExtCol);
        }
        
        MsrExtMatrix m_Matrix;

        /// <summary>
        /// the matrix on which we fummel herum on
        /// </summary>
        public IMutableMatrixEx Matrix {
            get {
                return m_Matrix.Mtx;
            }
        }

        /// <summary>
        /// sets an matrix entry; using this method (instead of directly modifying the matrix) ensures 
        /// that the additional list for manipulating columns are updated.
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <param name="val"></param>
        public void SetMatrixValue(int row, int col, double val) {
            m_Matrix[row, col] = val;
        }
        
        /// <summary>
        /// internally used encapsulation of <see cref="MsrMatrix"/>, which allows 'fast' access 
        /// of columns
        /// </summary>
        class MsrExtMatrix  {

            /// <summary>
            /// ctor.
            /// </summary>
            /// <param name="M"></param>
            /// <param name="ExtCol">
            /// key: processor rank 'p' <br/>
            /// value: a list of column indices (within the local range of columns of <paramref name="M"/>),
            /// which should be editable at rank 'p'.
            /// </param>
            public MsrExtMatrix(IMutableMatrixEx M, IDictionary<int,int[]> ExtCol) {

                this.ColPart = M.ColPartition;
                int i0Row = (int)M.RowPartitioning.i0, I = M.RowPartitioning.LocalLength,
                    i0Col = (int) this.ColPart.i0, J = this.ColPart.LocalLength;
                Mtx = M;


                // init
                // ====
                ColToRowLocal = new List<int>[ColPart.LocalLength];
                for (int j = 0; j < ColToRowLocal.Length; j++)
                    ColToRowLocal[j] = new List<int>();
                
                // build Column to row - mapping
                // =============================
                ColToRowExternal = new Dictionary<int, List<int>>();
                // key: global column index j, within the range of processor 'p'
                // values: global row indices 


                SortedDictionary<int, List<int>> ColForProc = new SortedDictionary<int, List<int>>();
                // key: MPI processor index 'p'
                // values: a set of global column indices, within the range of processor 'p',
                //         that contain nonzero entries on this processor

                int[] col = null;
                int L;

                // loop over all rows...
                for (int i = 0; i < I; i++) {

                    L = M.GetOccupiedColumnIndices(i + i0Row, ref col);

                    // loop over all nonzero entries in the row...
                    for(int l = 0; l < L; l++) {
                        int ColIndex = col[l];

                        int localColInd = ColIndex - i0Col;
                        if (localColInd >= 0 && localColInd < J) {
                            // column of 'entry' is within the local range of this processor
                            // + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                            ColToRowLocal[localColInd].Add(i + i0Row);
                        } else {
                            // column of 'entry' belongs to external processor 'proc'
                            // + + + + + + + + + + + + + + + + + + + + + + + + + + + +

                            int proc = this.ColPart.FindProcess(ColIndex);

                            {
                                //SortedDictionary<int, List<int>> ColToRowExt_proc;
                                //if (!ColToRowExternal.ContainsKey(proc)) {
                                //    ColToRowExt_proc = new SortedDictionary<int, List<int>>();
                                //    ColToRowExternal.Add(proc, ColToRowExt_proc);
                                //} else {
                                //    ColToRowExt_proc = ColToRowExternal[proc];
                                //}

                                int j = ColIndex;

                                List<int> Rows4Col;
                                if (!ColToRowExternal.ContainsKey(j)) {
                                    Rows4Col = new List<int>();
                                    ColToRowExternal.Add(j, Rows4Col);
                                } else {
                                    Rows4Col = ColToRowExternal[j];
                                }

                                Rows4Col.Add(i + i0Row);
                            }

                            {
                                List<int> ColForProc_proc;
                                if (!ColForProc.ContainsKey(proc)) {
                                    ColForProc_proc = new List<int>();
                                    ColForProc.Add(proc, ColForProc_proc);
                                } else {
                                    ColForProc_proc = ColForProc[proc];
                                }

                                if(!ColForProc_proc.Contains(ColIndex))
                                    ColForProc_proc.Add(ColIndex);
                            }
                        }
                    }
                }

                // communicate
                // ===========

                //SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw.MPI_COMM_WORLD);
                //{
                //    foreach (int proc in ColToRowExternal.Keys) {
                //        sms.SetCommPath(proc);
                //    }
                //    sms.CommitCommPaths();

                //    // send 
                //    foreach (int proc in ColToRowExternal.Keys) {
                //        SortedDictionary<int, List<int>> ColToRowExt_proc = ColToRowExternal[proc];
                //        SortedList _ColToRowExt_proc = new SortedList();
                //        foreach (int iCol in ColToRowExt_proc.Keys)
                //            _ColToRowExt_proc.Add(iCol, ColToRowExt_proc[iCol].ToArray());
                //    }

                //    // receive
                //    int p; SortedList rcv; 
                //    sms.GetNext(out p, out rcv);
                //    while (rcv != null) {

                //        foreach (int col in rcv.Keys) {
                //            int[] rowList = (int[])rcv[col];
                //            ColToRowLocal[col].AddRange(rowList);
                //        }

                //        sms.GetNext(out p, out rcv);
                //    }
                //}

                //sms.Dispose();

                // build 'ColProcessors'
                // =====================

                //ColProcessors = new List<int>[ColToRowLocal.Length];

                //for (int j = 0; j < ColToRowLocal.Length; j++) {
                //    List<int> mpiRank = null;
                //    foreach (int rowind in ColToRowLocal[j]) {
                //        int riloc = rowind - i0Row;

                //        if (riloc < 0 || riloc >= I) {
                //            if (mpiRank == null)
                //                mpiRank = new List<int>();
                //            mpiRank.Add(M.RowPartiton.FindProcess(rowind));
                //        }

                //        ColProcessors[j] = mpiRank;
                //    }
                //}

                // communicate: build 'ColProcessors'
                // ==================================
                {
                    ColProcessors = new List<int>[ColPart.LocalLength];

                    SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);

                    sms.SetCommPathsAndCommit(ColForProc.Keys);

                    foreach (int proc in ColForProc.Keys)
                        sms.Transmitt(proc, ColForProc[proc].ToArray());

                    int i0Loc = (int) ColPart.i0;

                    int rcvproc; int[] ColIndices;
                    while(sms.GetNext(out rcvproc, out ColIndices)) {

                        int Rank;                        
                        {

                            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD,out Rank);
                            //Console.WriteLine("P# " + Rank + ": receiving from P# " + rcvproc);
                        }

                        foreach (int ColInd in ColIndices) {
                            int localColInd = ColInd - i0Loc;
                            if (localColInd < 0 || localColInd >= ColPart.LocalLength)
                                throw new IndexOutOfRangeException("internal error");


                            if (ColProcessors[localColInd] == null)
                                ColProcessors[localColInd] = new List<int>();

                            if (ColProcessors[localColInd].Contains(rcvproc))
                                throw new ApplicationException("internal error.");
                            ColProcessors[localColInd].Add(rcvproc);
                        }
                    }
                    
                    sms.Dispose();
                }

                if (ExtCol != null) {
                    var send = new Dictionary<int,List<Tuple<int,List<int>>>>();
                    int myRank = M.RowPartitioning.MpiRank;

                    foreach(var kv in ExtCol) {
                        int rank = kv.Key;
                        int[] ColIdx = kv.Value;

                        var sendToRank = new List<Tuple<int, List<int>>>();

                        foreach(int iCol in ColIdx) {
                            List<int> c2p = ColProcessors[iCol - i0Col];
                            var t = new Tuple<int, List<int>>(iCol, c2p != null ? new List<int>(c2p) : new List<int>());
                            t.Item2.Add(myRank);
                            sendToRank.Add(t);
                        }

                        send.Add(rank, sendToRank);
                    }

                    var receive = SerialisationMessenger.ExchangeData(send,csMPI.Raw._COMM.WORLD);

                    ColProcessorsExternal = new Dictionary<int,List<int>>();
                    foreach (var kv in receive) {
                        var val = kv.Value;

                        foreach (var t in val) {
                            int iCol = t.Item1;
                            List<int> ranks = t.Item2;
                            int iMyRank = ranks.IndexOf(myRank);
                            if (iMyRank >= 0)
                                ranks.RemoveAt(iMyRank);

                            ColProcessorsExternal.Add(t.Item1, t.Item2);

                            Debug.Assert(this.ColPart.FindProcess(t.Item1) == kv.Key);
                        }
                    }

#if DEBUG
                    foreach (var procList in ColProcessorsExternal.Values) {
                        Debug.Assert(procList.Contains(myRank) == false);
                    }
#endif

                }
            }

            /// <summary>
            /// defines which columns an MPI process is allowed to manipulate with;
            /// </summary>
            public IPartitioning ColPart;

            /// <summary>
            /// the aggregated matrix
            /// </summary>
            public IMutableMatrixEx Mtx;
            

            /// <summary>
            /// for each column (within the local range of column partition <see cref="ColPart"/>),
            /// the indices of those locally stored rows that contain non-zero entries;<br/>
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>1st index: local column index (as defined by column partition <see cref="ColPart"/>)</item>
            ///   <item>2nd index: only enumeration</item>
            /// </list>
            /// </remarks>
            public List<int>[] ColToRowLocal;

            /// <summary>
            /// for each column (NOT within the local range of column partition <see cref="ColPart"/>),
            /// the indices of those locally stored rows that contain non-zero entries;<br/>
            /// key: global column index j, within the range of processor 'p'<br/>
            /// values: enumeration of global row indices 
            /// </summary>
            public Dictionary<int,List<int>> ColToRowExternal = new Dictionary<int, List<int>>();
            
            /// <summary>
            /// for each column, the 
            /// list of MPI processor ranks, that hold nonzero entries of this column.<br/>
            /// index: local column index (as defined by column partition <see cref="ColPart"/>);<br/>
            /// the current processor ('my rank') is not included.
            /// </summary>
            public List<int>[] ColProcessors;

            /// <summary>
            /// for each certain extern columns, the 
            /// list of MPI processor ranks, that hold nonzero entries of this column.<br/>
            /// index: local column index (as defined by column partition <see cref="ColPart"/>);<br/>
            /// the current processor ('my rank') is not included.
            /// </summary>
            public Dictionary<int,List<int>> ColProcessorsExternal;


            /// <summary>
            /// gets/sets the values of <see cref="Mtx"/>;
            /// using this property (for setting) ensures that the column-to-row - mapping (<see cref="ColToRowLocal"/>, <see cref="ColToRowExternal"/>)
            /// gets updated when the matrix population pattern (non-zero entries may become zero and vice-versa)
            /// changes.
            /// </summary>
            /// <param name="row">row index</param>
            /// <param name="col">column index</param>
            /// <returns></returns>
            public double this[int row, int col] {
                get {
                    return Mtx[row, col];
                }
                set {
                    int jLocal = col - (int)ColPart.i0;
                    List<int> _col;
                    if (jLocal >= 0 && jLocal < ColPart.LocalLength) {
                        _col = ColToRowLocal[col - ColPart.i0];
                    } else {
                        //int proc = ColPart.FindProcess(col);

                        if (!ColToRowExternal.TryGetValue(col, out _col)) {
                            if (value != 0.0) {
                                _col = new List<int>();
                                _col.Add(row);
                                ColToRowExternal.Add(col, _col);
                            }

                        }
                    }

                    if (_col != null && _col.Contains(row)) {
                        if (value == 0.0)
                            _col.Remove(row);
                    } else {
                        if (value != 0.0) {
                            _col.Add(row);
                        }
                    }

                    Mtx[row, col] = value;
                }
            }
        }
        
        /// <summary>
        /// checks whether <paramref name="i"/> is within the local range of the
        /// partition <paramref name="p"/> - if not, an exception is thrown.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="p"></param>
        void TestIndex( int i, IPartitioning p) {
            i -= (int)p.i0;
            if(i < 0 || i >= p.LocalLength)
                throw new IndexOutOfRangeException("row/col index out of range.");
        }

         /// <summary>
        /// adds row <paramref name="iSrc"/> times <paramref name="alpha"/> to row <paramref name="iDst"/>
        /// </summary>
        /// <param name="iSrc"></param>
        /// <param name="iDst">row index of the "accumulator" row</param>
        /// <param name="alpha">scaling for the source row</param>
        public void RowAddition(int iSrc, int iDst, double alpha) {
            TestIndex(iSrc,m_Matrix.Mtx.RowPartitioning);
            TestIndex(iDst,m_Matrix.Mtx.RowPartitioning);

            //MsrMatrix.MatrixEntry[] row = m_Matrix.Mtx.GetRow(iSrc);
            double[] val = null;
            int[] col = null;
            int L;
            L = m_Matrix.Mtx.GetRow(iSrc, ref col, ref val);

            //double rowSum = 0;
            for (int j = 0; j < L; j++) {
                int ColIdx = col[L];
                m_Matrix[iDst, ColIdx] += alpha * val[j];
                //rowSum += row[j].val;
            }
        }


        List<ColOp> DeferredColOpList = new List<ColOp>();


        /// <summary>
        /// baseclass for all kinds of column operations
        /// </summary>
        [Serializable]
        abstract class ColOp {

            /// <summary>
            /// global column index of the altered row
            /// </summary>
            public int jCol;
        }

        /// <summary>
        /// represents <see cref="ColClearDeferred"/> in the <see cref="DeferredColOpList"/>-list.
        /// </summary>
        [Serializable]
        class ColClear : ColOp {

        }

        /// <summary>
        /// represents <see cref="ColMulDeferred"/> in the <see cref="DeferredColOpList"/>-list.
        /// </summary>
        [Serializable]
        class ColMul : ColOp {

            /// <summary>
            /// multiplication factor
            /// </summary>
            public double alpha;
        }

        /// <summary>
        /// represents <see cref="ColAdditionDeferred"/> in the <see cref="DeferredColOpList"/>-list.
        /// </summary>
        [Serializable]
        class ColAddition : ColOp {

            /// <summary>
            /// global index of the row to add
            /// </summary>
            public int iSrc;

            /// <summary>
            /// scaling for row <see cref="iSrc"/>
            /// </summary>
            public double alpha;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="iSrc"></param>
        /// <param name="iDst"></param>
        /// <param name="alpha"></param>
        public void ColAdditionDeferred(int iSrc, int iDst, double alpha) {
            TestIndex(iDst, m_Matrix.Mtx.ColPartition);
#if DEBUG
            if (!m_Matrix.ColPart.IsInLocalRange(iSrc)) {
                if (!m_Matrix.ColProcessorsExternal.ContainsKey(iSrc))
                    throw new IndexOutOfRangeException();
            }
#endif
            ColAddition ca = new ColAddition();
            ca.alpha = alpha;
            ca.iSrc = iSrc;
            ca.jCol = iDst;
            DeferredColOpList.Add(ca);
        }

        /// <summary>
        /// multiplication of row <paramref name="iRow"/> by scalar <paramref name="alpha"/>
        /// </summary>
        /// <param name="iRow"></param>
        /// <param name="alpha"></param>
        /// <remarks>
        /// Row operations (unlike column operation) can be always executed
        /// on the local MPI processor.
        /// </remarks>
        public void RowMul(int iRow, double alpha) {
            if (alpha == 0.0) {
                RowClear(iRow);
            } else {

                TestIndex(iRow, m_Matrix.Mtx.RowPartitioning);

                int[] icol = null;
                int L = m_Matrix.Mtx.GetOccupiedColumnIndices(iRow, ref icol);
                for (int j = 0; j < L; j++)
                    m_Matrix.Mtx[iRow, icol[j]] *= alpha;
            }
        }

        /// <summary>
        /// sets all entries of row <paramref name="iRow"/> to 0.0;
        /// </summary>
        /// <param name="iRow">
        /// global row index
        /// </param>
        public void RowClear(int iRow) {
            TestIndex(iRow, m_Matrix.Mtx.RowPartitioning);
            int[] icol = null;
            int L = m_Matrix.Mtx.GetOccupiedColumnIndices(iRow, ref icol);

            for (int j = 0; j < L; j++)
                m_Matrix[iRow, icol[j]] = 0;
        }

        /// <summary>
        /// sets all entries of column <paramref name="iCol"/> to 0.0;
        /// </summary>
        /// <param name="iCol">
        /// global column index
        /// </param>        
        public void ColClearDeferred(int iCol) {
#if DEBUG
            if (!m_Matrix.ColPart.IsInLocalRange(iCol)) {
                if (!m_Matrix.ColProcessorsExternal.ContainsKey(iCol))
                    throw new IndexOutOfRangeException();
            }
#endif
            ColClear cm = new ColClear();
            cm.jCol = iCol;
            DeferredColOpList.Add(cm);
        }


        /// <summary>
        /// Multiplication of column <paramref name="iCol"/> by scalar <paramref name="alpha"/>;
        /// Note that column operations are deferred operation and require a call to <see cref="CompleteColOperation"/>
        /// for being executed.
        /// </summary>
        /// <param name="iCol"></param>
        /// <param name="alpha"></param>
        /// <remarks>
        /// Note that column operations, unlike row operations, can invoke data changes on more than one,
        /// but hopefully not too many MPI processors.
        /// Because MPI communication for each operation would be too slow,
        /// its better to collect a bunch of operations in a list and then execute them 'at once'.
        /// Therefore, column operations are designed as deferred calls, i.e. they need a call
        /// to <see cref="CompleteColOperation"/> to be executed.
        /// </remarks>
        public void ColMulDeferred(int iCol, double alpha) {
            ColMul cm = new ColMul();
            cm.alpha = alpha;
            cm.jCol = iCol;
            DeferredColOpList.Add(cm);
        }

        /// <summary>
        /// An MPI-collective call, which executes all column operations.
        /// </summary>
        public void CompleteColOperation() {

            int j0Loc = (int)m_Matrix.ColPart.i0;
            int LenLoc = m_Matrix.ColPart.LocalLength;


            // sort operations according to processor
            // ======================================

            // keys: MPI processor rank p
            // values: list of operations to execute on p
            SortedDictionary<int, List<ColOp>> OperationsPerProcessor = new SortedDictionary<int, List<ColOp>>();

            List<int> InvokedProc;
            for (int i = 0; i < DeferredColOpList.Count; i++) {
                ColOp op = DeferredColOpList[i];
                
                bool skip = false;
                ColAddition ca = op as ColAddition;
                List<int> InvokesProcSrc = null;
                if (ca != null) {
                    // we have a column addition - this requires some special treatments
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                    // problem 1: if 
                    //          * the destination col. (aka. accumulator column), i.e. column no. 'op.jCol'
                    //            is zero,
                    //          * and the source row is nonzero
                    //
                    //  then we need to send the command not only to processors which contain
                    //  nonzeros in destination row, but also the source row.
                    //

                    // problem 2: for subsequent operations, we may expect that some column which has
                    // originally been zero now contains nonzero elements.
                    // So, therefore, we have to add the processor set of the source row
                    // (i.e. 'm_Matrix.ColProcessors[ca.iSrc]' to the processor set of the 
                    // destination row.


                    if (m_Matrix.ColPart.IsInLocalRange(ca.iSrc)) {
                        InvokesProcSrc = m_Matrix.ColProcessors[ca.iSrc - j0Loc];
                    } else {
                        if (!this.m_Matrix.ColProcessorsExternal.TryGetValue(ca.iSrc, out InvokesProcSrc))
                            throw new IndexOutOfRangeException("manipulation operation not available on current processor");
                    }
                    
                    if (InvokesProcSrc != null) {
                        if (m_Matrix.ColProcessors[op.jCol - j0Loc] == null)
                            m_Matrix.ColProcessors[op.jCol - j0Loc] = InvokesProcSrc;
                        else {
                            InvokedProc = m_Matrix.ColProcessors[op.jCol - j0Loc];

                            foreach (int ps in InvokesProcSrc) {
                                if (!InvokedProc.Contains(ps))
                                    InvokedProc.Add(ps);
                            }
                        }
                    } else {
                        // optimization: source column is zero (on other processors) -> nothing to do
                        skip = true;
                    }
                }

                if (m_Matrix.ColPart.IsInLocalRange(op.jCol)) {
                    InvokedProc = m_Matrix.ColProcessors[op.jCol - j0Loc];
                } else {
                    if (!this.m_Matrix.ColProcessorsExternal.TryGetValue(op.jCol, out InvokedProc))
                        throw new IndexOutOfRangeException("manipulation operation not available on current processor");
                }

                if (InvokedProc != null && !skip) {
                    foreach (int proc in InvokedProc) {

                        bool skip2 = false;
                        if (ca != null) {
                            // optimization: don't need to send column addition if
                            //               the source row is zero
                            if (!InvokesProcSrc.Contains(proc))
                                skip2 = true;
                        }

                        if (!skip2) {
                            List<ColOp> DeferredOp_proc;
                            if (OperationsPerProcessor.ContainsKey(proc)) {
                                DeferredOp_proc = OperationsPerProcessor[proc];
                            } else {
                                DeferredOp_proc = new List<ColOp>();
                                OperationsPerProcessor.Add(proc, DeferredOp_proc);
                            }
                            DeferredOp_proc.Add(op);
                        }
                    }
                }
            }

            // transmit to other processors
            // ============================

            SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
            foreach (int proc in OperationsPerProcessor.Keys)
                sms.SetCommPath(proc);
            sms.CommitCommPaths();

            foreach (int proc in OperationsPerProcessor.Keys)
                sms.Transmitt(proc, OperationsPerProcessor[proc].ToArray());

            int rcvp; ColOp[] rcv;
            while (sms.GetNext(out rcvp, out rcv)) {

                DeferredColOpList.AddRange(rcv); // operations from different processors
                                                 // commute (because they are bound to the column partition)
                                                 // therefore, it doesn't matter how they are added

//#if DEBUG
//                {
//                    int Rank;
//                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out Rank);
//                    //Console.WriteLine("P# " + Rank + ": " + rcv.Length + " operation(s) received from P# " + rcvp);

//                    foreach (var op in rcv) {
//                        if ((op.jCol >= m_Matrix.ColPart.i0 && op.jCol < (m_Matrix.ColPart.i0 + m_Matrix.ColPart.LocalLength)))
//                            throw new ApplicationException("internal error");
//                    }
//                }
//#endif
            }

            // execute operations
            // ==================

            int L = DeferredColOpList.Count;
            for (int l = 0; l < L; l++) {
                ColOp op = DeferredColOpList[l];

                int jlocal = op.jCol - j0Loc;

                if (op is ColAddition) {
                    double alpha = ((ColAddition)op).alpha;
                    int iSrc = ((ColAddition)op).iSrc;

                    int[] col;
                    if (iSrc >= j0Loc && iSrc < (j0Loc + LenLoc))
                        // operation in local column range
                        col = m_Matrix.ColToRowLocal[iSrc - j0Loc].ToArray();
                    else {
                        // operation comes from other processor
                        List<int> _col;
                        if (m_Matrix.ColToRowExternal.TryGetValue(iSrc, out _col)) {
                            col = _col.ToArray();
                        } else {
                            col = new int[0];
                        }
                        //col = m_Matrix.ColToRowExternal[iSrc].ToArray();
                    }

                    foreach (int irow in col) {
                        m_Matrix[irow, op.jCol] += m_Matrix[irow, iSrc] * alpha;
                    }
                } else {
                    int[] col;
                    if (jlocal >= 0 && jlocal < LenLoc)
                        // operation in local column range
                        col = m_Matrix.ColToRowLocal[jlocal].ToArray();
                    else {
                        // operation comes from other processor
                        List<int> _col;
                        if (m_Matrix.ColToRowExternal.TryGetValue(op.jCol, out _col)) {
                            col = _col.ToArray();
                        } else {
                            col = new int[0];
                        }
                        //col = m_Matrix.ColToRowExternal[op.jCol].ToArray();
                    }

                    if (op is ColMul) {
                        double alpha = ((ColMul)op).alpha;
                        foreach (int irow in col) {
                            m_Matrix[irow, op.jCol] *= alpha;
                        }
                    } else if (op is ColClear) {
                        foreach (int irow in col) {
                            m_Matrix[irow, op.jCol] = 0;
                        }
                    } else {
                        throw new NotImplementedException();
                    }
                }
            }

            // finish & return
            // ===============
            DeferredColOpList.Clear();
        }
    }
}
