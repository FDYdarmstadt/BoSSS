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
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Utils;

namespace ilPSP.LinSolvers.monkey {


    /// <summary>
    /// Base class for all MSR sparse matrices
    /// </summary>
    /// <remarks>
    /// Example of the parallel MSR format:
    /// <code>
    /// e.g. a matrix on two processors:
    ///
    ///       Col.Ind. -> 0   1   2   3   4  |  5   6   7
    /// Row Ind.         (0) (1) (2) (3) (4)   (0) (1) (2)
    ///   |   0 (0)     [ a   0   j   0   0  |  0   0   0 ]
    ///   V   1 (1)     [ 0   b   0   0   0     0   0   0 ]
    ///       2 (2)     [ 0   0   c   0   k  |  l   m   0 ]
    ///       3 (3)     [ i   0   0   d   0     0   n   0 ]       Processor 0
    ///        - - - - -[- - - - - - - - - - |- - - - - - - - - - - - - - - - - -
    ///       4 (0)     [ 0   o   0   0   e     0   0   0 ]       Processor 1
    ///       5 (1)     [ 0   p   0   0   0  |  f   0   0 ]
    ///       6 (2)     [ 0   0   q   0   0     0   h   0 ]
    ///                                                         () Braces denotes local indices
    ///
    ///                    [ internal part   |  external part ]
    /// General structure  [  of Proc. 0     |   of Proc. 0   ]
    /// of the matrix:     [ -------------------------------- ]
    ///                    [ external part   |  internal part ]
    ///                    [  of Proc. 1     |   of Proc. 1   ]
    ///
    ///
    ///       Row Partition: Processor 0: i0 = 0, LocalLength = 4
    ///                      Processor 1: i0 = 4, LocalLength = 3
    ///       Column Partition: Processor 0: i0 = 0, LocalLength = 5
    ///                         Processor 1: i0 = 5, LocalLength = 7
    ///           => every vector that should be multiplied by the matrix needs to have the
    ///              same partition as the column partition.
    ///
    ///       Comm List(s) on Processor 0:   P1: { 1, 2, 4 }, i.e. Proc. 0 needs to send 
    ///                                                       rows 1,2,4 (local indices) to Proc. 1.
    ///       Comm List(s) on Processor 1:   P1: { 0, 1 },    i.e. Proc. 0 needs to send 
    ///                                                       rows 0,1 (local indices) to Proc. 0.
    ///       Internal part on Processor 0: (CSR - Format)
    ///                           RowStart = {0, 2, 4, 6, 7}       row 0 is from 0 to 3 (excl.), row 1 from 3 to 4 (excl.), ...
    ///                           Val      = {a, j, b, c, k, i, d}
    ///                           ColInd   = {0, 2, 1, 2, 4, 0, 3}   local column indices of the Val's
    ///
    ///       Internal part on Processor 1: (CSR - Format)
    ///                           RowStart = {0, 0, 1, 2}       
    ///                           Val      = {f, h}
    ///                           ColInd   = {0, 1}
    ///
    ///
    ///        External part on Processor 0:
    ///                           RowInd   = {2,3}, i.e. contains entries in rows 2 and 3  (local row indices)
    ///                           RowStart = {0, 2, 3}
    ///                           Val =      {l, m, n}
    ///                           ColInd =   {0, 1, 0}, indices into the comm list !!!
    ///
    ///        External part on Processor 1:
    ///                           RowInd   = {0,1,2}, i.e. contains entries in rows 0,1,2 (local row indices)
    ///                           RowStart = {0, 2, 3}
    ///                           Val =      {o, e, p, q}
    ///                           ColInd =   {0, 2, 0, 1}, indices into the comm list !!!
    ///
    /// </code>
    /// </remarks>
    public abstract partial class MatrixBase : LockAbleObject, IMutableMatrixEx {

        /// <summary>
        /// Not Implemented.
        /// </summary>
        public void Clear() {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Column partition
        /// </summary>
        protected IPartitioning m_ColPart;

        /// <summary>
        /// Row partition
        /// </summary>
        protected IPartitioning m_RowPart;

        /// <summary>
        /// gets the column partition (the partition that a vector must habe to be multiplied by this matrix);
        /// </summary>
        virtual public IPartitioning ColPartition { get { return m_ColPart; } }

        /// <summary>
        /// gets the row partition (how the rows of the matrix are distributed among the different 
        /// processes in the current MPI communicator)
        /// </summary>
        virtual public IPartitioning RowPartitioning { get { return m_RowPart; } }

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
        /// constructor
        /// </summary>
        protected MatrixBase(MsrMatrix M) {
            ilPSP.MPICollectiveWatchDog.Watch();
            if (M.RowPartitioning.IsMutable)
                throw new NotSupportedException();
            if (M.ColPartition.IsMutable)
                throw new NotSupportedException();
            m_CellSize = 1; // ggT(M.RowPartitioning.BlockSize, M.ColPartition.BlockSize);

            
            m_ColPart = M.ColPartition;
            m_RowPart = M.RowPartitioning;

            //if (M.ColPerBlock != M.RowsPerBlock)
            //    throw new ApplicationException("not supported.");
            
        }

        


        static int ggT(int n, int m) {
            return n == m ? n : n < m ? ggT(n, m - n) : ggT(n - m, m);
        }

        private int m_CellSize;

        /// <summary>
        /// size of sub - matrices (aka. 'Cells');
        /// </summary>
        public int CellSize {
            get { return m_CellSize; }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <remarks>
        /// the key of the dictionary are the processor ranks of those processors which
        /// send data to this processor.
        /// </remarks>
        public IDictionary<int, External> ExtMatrix;

        /// <summary>
        /// used to collect all info in an  <see cref="External"/>-object while 'parsing' the matrix;
        /// </summary>
        class ExternalTmp {

            protected List<double> Vals = new List<double>();
            protected List<int> ColInd = new List<int>();
            protected List<int> RowStart = new List<int>();
            protected List<int> RowIndices = new List<int>();
            protected List<int> GlobalColInd = new List<int>();

            int m_RowInd = 0;
            int m_ColInd = -1;

            int m_Cnt = 0;

            public void AddEntry(int __ColInd, int __GloablColInd, double val) {
                if (__ColInd <= m_ColInd)
                    throw new ArgumentException("rows must be specified from left to right (ascending column index).");

                if (m_ColInd < 0) {
                    // a non- emtpy row starts
                    RowIndices.Add(m_RowInd);
                    RowStart.Add(m_Cnt);
                }

                m_ColInd = __ColInd;
                Vals.Add(val);
                ColInd.Add(__ColInd);
                GlobalColInd.Add(__GloablColInd);
                m_Cnt++;
            }

            /// <summary>
            /// must be called before moving to the next row
            /// </summary>
            public void NextRow() {
                m_RowInd++;
                m_ColInd = -1;
            }

            public External GetFinalObj() {
                RowStart.Add(m_Cnt);

                External r;
                r.RowStart = this.RowStart.ToArray();
                this.RowStart = null;

                r.ColInd = this.ColInd.ToArray();
                this.ColInd = null;

                r.Val = this.Vals.ToArray();
                this.Vals = null;

                r.GlobalColInd = this.GlobalColInd.ToArray();
                this.GlobalColInd = null;

                r.rowInd = this.RowIndices.ToArray();
                this.RowIndices = null;

                return r;
            }

        }

        /// <summary>
        /// External part of the matrix, i.e. those entries for which communication is necessary
        /// (when doing a GEMV).
        /// </summary>
        public struct External {

            /// <summary>
            /// index: parallel with <see cref="RowStart"/>;<br/>
            /// content: row indices (local indices), i.e. all rows that are non-empty;
            /// </summary>
            public int[] rowInd;

            /// <summary>
            /// index: parallel with index of <see cref="ColInd"/>;<br/>
            /// content: matrix entries
            /// </summary>
            public double[] Val;

            /// <summary>
            /// index: parallel with index of <see cref="Val"/>;<br/>
            /// content: local (on this MPI process) column indices 
            /// </summary>
            public int[] ColInd;

            /// <summary>
            /// index: parallel with index of <see cref="Val"/>;<br/>
            /// content: global (over all MPI processors) column indices
            /// </summary>
            public int[] GlobalColInd;

            /// <summary>
            /// index: parallel with <see cref="RowStart"/>;<br/>
            /// content: entry [i] is a pointer into <see cref="Val"/> and <see cref="ColInd"/>
            /// that denotes where row <see cref="rowInd"/>[i] starts;
            /// </summary>
            public int[] RowStart;

        }

        /// <summary>
        /// the local part of the matrix;
        /// </summary>
        protected FormatBase m_LocalMtx;

        /// <summary>
        /// matrix assembly; must be called by each implementation, 
        /// </summary>
        /// <param name="M"></param>
        protected void PackMatrix(IMutableMatrixEx M) {
            ilPSP.MPICollectiveWatchDog.Watch();
            
           
            IPartitioning rp = M.RowPartitioning;
            IPartitioning cp = m_ColPart;

            // define Comm List
            // ================
            
            SortedDictionary<int, List<int>> CommLists = new SortedDictionary<int, List<int>>();
            // keys: processor rank p
            // values: List of global indices, which processor p needs to send to this processor

            int Lr;
            double[] val = null;
            int[] col = null;

            int L = rp.LocalLength;
            int i0 = (int)rp.i0;
            for (int iLoc = 0; iLoc < L; iLoc++) { // loop over all matrix rows...
                int iGlob = i0 + iLoc;

                //MsrMatrix.MatrixEntry[] row = (asMsr==null) ? M.GetRow(iGlob) : asMsr.GetRowShallow(iGlob);
                Lr = M.GetOccupiedColumnIndices(iGlob, ref col);

                for (int j = 0; j < Lr; j++) { // loop over all nonzero entries in row 'iGlob'
                    int jGlob = col[j];

                    if (cp.i0 <= jGlob && jGlob < (cp.i0 + cp.LocalLength)) {
                        // Entry on current processor

                    } else {
                        int proc = cp.FindProcess(jGlob);
                        // Entry on Processor proc

                        if (!CommLists.ContainsKey(proc))
                            CommLists.Add(proc, new List<int>());

                        List<int> CommList_proc = CommLists[proc];
                        if (!CommList_proc.Contains(jGlob)) // a lot of room for optimization
                            CommList_proc.Add(jGlob);
                    }
                }
            }

            // sort com list
            // =============
            {
                foreach (List<int> cl in CommLists.Values)
                    cl.Sort();
            }

            // define matrix
            // =============
            {
                TempCSR intTmp = new TempCSR();

                SortedDictionary<int, ExternalTmp> extTmp = new SortedDictionary<int, ExternalTmp>();
                foreach (int proc in CommLists.Keys)
                    extTmp.Add(proc, new ExternalTmp());

                for (int iLoc = 0; iLoc < L; iLoc++) {
                    int iGlob = i0 + iLoc;

                    Lr = M.GetRow(iGlob, ref col, ref val);


                    for (int j = 0; j < Lr; j++) {
                        int jGlob = col[j];
                        double Value = val[j];

                        bool bIsDiag = (iGlob == jGlob);

                        if (cp.i0 <= jGlob && jGlob < (cp.i0 + cp.LocalLength)) {
                            // Entry on current processor

                            intTmp.AddEntry(jGlob - (int)cp.i0, Value, bIsDiag);

                        } else {
                            int proc = cp.FindProcess(jGlob);
                            // Entry on Processor proc

                            List<int> CommList_proc = CommLists[proc];
                            int jloc = CommList_proc.IndexOf(jGlob);

                            ExternalTmp et = extTmp[proc];
                            et.AddEntry(jloc, jGlob, Value);
                        }
                    }

                    intTmp.NextRow();
                    foreach (ExternalTmp et in extTmp.Values)
                        et.NextRow();

                    
                }

                m_LocalMtx = AssembleFinalFormat(intTmp);

                ExtMatrix = new Dictionary<int, External>();
                foreach (int proc in extTmp.Keys)
                    ExtMatrix.Add(proc, extTmp[proc].GetFinalObj());
            }

            // send/receive & transform  Comm lists
            // ====================================
            {
                SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                SortedDictionary<int, int[]> CommListsTo = new SortedDictionary<int, int[]>();

                foreach (int proc in CommLists.Keys)
                    sms.SetCommPath(proc);
                sms.CommitCommPaths();

                foreach (int proc in CommLists.Keys)
                    sms.Transmitt(proc, CommLists[proc].ToArray());

                int _proc;
                int[] CommListReceived;
                sms.GetNext(out _proc, out CommListReceived);
                int Lcol = m_ColPart.LocalLength;
                int i0col = (int)m_ColPart.i0;
                while (CommListReceived != null) {

                    // convert indices to local coordinates
                    for (int i = 0; i < CommListReceived.Length; i++) {
                        CommListReceived[i] -= i0col;

                        // check:
                        if (CommListReceived[i] < 0 || CommListReceived[i] >= Lcol)
                            throw new ApplicationException("internal error: something wrong with received Comm List.");
                    }

                    CommListsTo.Add(_proc, CommListReceived);

                    sms.GetNext(out _proc, out CommListReceived);
                }

                sms.Dispose();

                m_SpmvCommPattern = new SpmvCommPattern();
                m_SpmvCommPattern.ComLists = CommListsTo;

            }

            // record the number of elements which we receive
            // ==============================================
            {
                m_SpmvCommPattern.NoOfReceivedEntries = new Dictionary<int, int>();
                foreach (int p in CommLists.Keys) {
                    m_SpmvCommPattern.NoOfReceivedEntries.Add(p, CommLists[p].Count);
                }
            }
        }




        /// <summary>
        /// converts the temporary CSR format into the final matrix format
        /// </summary>
        abstract protected FormatBase AssembleFinalFormat(TempCSR tmp);

        /// <summary>
        /// 
        /// </summary>
        protected SpmvCommPattern m_SpmvCommPattern;

        /// <summary>
        /// the communication pattern for multiplying a vector with this matrix, see <see cref="SpmvCommPattern.ComLists"/>
        /// </summary>
        public SpmvCommPattern _SpmvCommPattern {
            get { return m_SpmvCommPattern; }
        }

        /// <summary>
        /// the communication pattern for multiplying a vector with a sparse matrix, see <see cref="ComLists"/>
        /// </summary>
        public class SpmvCommPattern {

            /// <summary>
            /// the Communication lists for each processor which receives data from this processor.
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>keys <em>p</em>: ranks of processors to send data to</item>
            ///   <item>values: the comm list for the processor <em>p</em></item>
            /// </list>
            /// If a distributed MSR matrix is multiplied by a distributed vector <em>v</em>
            /// (the matrix rows are distributed according to its row partition
            /// <see cref="MatrixBase.RowPartitioning"/>, while <em>v</em> is distributed
            /// according to the column partition <see cref="VectorBase.Part"/>)
            /// specific entries of <em>v</em> must be exchanged between MPI processes.<br/>
            /// The indices of these entries are stored in the comm lists;
            /// A comm list for processor <em>p</em> is an sorted list (in ascending order) of local indices
            /// which must be sent from this process to process <em>p</em>.
            /// </remarks>
            public IDictionary<int, int[]> ComLists;

            /// <summary>
            /// the number of entries which this processor receives from other processors.
            /// </summary>
            /// <remarks>
            /// <list type="bullet">
            ///   <item>keys <em>p</em>: ranks of processors form which this process receives data</item>
            ///   <item>values: the number of entries (i.e. length of the comm list) 
            ///   which are received from the processor <em>p</em></item>
            /// </list>
            /// </remarks>
            public IDictionary<int, int> NoOfReceivedEntries;
        }

        #region ISparseMatrix Members


        /// <summary>
        /// used by <see cref="SpMV"/> to convert its input arguments into <see cref="VectorBase"/>-objects
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="a"></param>
        /// <param name="len">
        /// partition of the vector that should be created;
        /// local length of the partition (<see cref="Partitioning.LocalLength"/>)
        /// must be equal to the length/count of <paramref name="a"/>
        /// </param>
        /// <param name="CopyIsShallow">
        /// true on exit, if the returned vale is a shallow copy of <paramref name="a"/>, i.e.
        /// changing the returned value would alter also <paramref name="a"/>;
        /// </param>
        /// <returns></returns>
        protected abstract VectorBase CreateVec<T>(T a, IPartitioning len, out bool CopyIsShallow)
            where T : IList<double>;


        abstract internal void SpMV_Local_Start(double alpha, VectorBase a, double beta, VectorBase acc);
        abstract internal void SpMV_Local_Middle(double alpha, VectorBase a, double beta, VectorBase acc);
        abstract internal void SpMV_Local_End(double alpha, VectorBase a, double beta, VectorBase acc);

        abstract internal void SpMV_External_Begin(double alpha, double beta, VectorBase acc);
        abstract internal void SpMV_External_RecvCallBack(int procRank, IntPtr values);
        abstract internal void SpMV_External_Finalize();


        /// <summary>
        /// sparse matrix/vector product; requires a locked object (<see cref="LockAbleObject.Lock"/>); This function should be
        /// used by the iterative solvers.
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="_a_comm"></param>
        /// <param name="beta"></param>
        /// <param name="_acc"></param>
        public void SpMV_Expert(double alpha, VectorBase.CommVector _a_comm, double beta, VectorBase _acc) {
            if (!this.IsLocked || !_a_comm.Owner.IsLocked || !_acc.IsLocked)
                throw new ApplicationException("objects must be locked.");

            if (!Object.ReferenceEquals(this, _a_comm.Mtx))
                throw new ArgumentException("input vector was not specified for this matrix.", "_a_comm");

            VectorBase _a = _a_comm.Owner;
            
            // MPI: fill send buffer
            // ---------------------
            //if (Enviroment.MPIEnv.MPI_Rank == 0)
            //    Debugger.Break();
            _a_comm.FillSendBuffer();

            

            // Start multiplying inner Part
            SpMV_Local_Start(alpha, _a, beta, _acc); //   A GPU implementation would start the 
            //                                            computation at this point;
            //                                            A CPU implementation, which blocks the main processor
            //                                            will leave this method empty, because the 
            //                                            MPI communication (transmission) should be started first

            // MPI communication: start transmission
            // -------------------------------------
            _a_comm.StartTransmissionImReturn();
            _a_comm.InitReceiveImReturn();

            SpMV_Local_Middle(alpha, _a, beta, _acc); // On a CPU implementation, this is an ideal point for implementing 
            //                                           the real workload; The communication threads are started and are
            //                                           waiting for the external data to arrive

            // wait for transmission to finish/do external parts
            // -------------------------------------------------
            SpMV_External_Begin(alpha, beta, _acc);

            _a_comm.WaitCommFinish(SpMV_External_RecvCallBack); // calls the method 'GEMV_External_RecvCallBack' every time an
            //                                                     other processor delivers data.

            SpMV_Local_End(alpha, _a, beta, _acc); // Blocks, until the computation of the local part is finished;
            //                                        A CPU - implementation will typically leave this method 
            //                                        empty. 

            SpMV_External_Finalize(); // Here, a GPU implementation will combine the 
            //                           locally computed part (done on GPU) with the external
            //                           parts (done on CPU).
        }


        /// <summary>
        /// performs <paramref name="acc"/> = <paramref name="acc"/>*<paramref name="beta"/> + <paramref name="alpha"/>*this*<paramref name="a"/>;
        /// </summary>
        /// <typeparam name="VectorType1"></typeparam>
        /// <typeparam name="VectorType2"></typeparam>
        /// <param name="alpha"></param>
        /// <param name="a"></param>
        /// <param name="beta"></param>
        /// <param name="acc"></param>
        /// <remarks>
        /// works only in unlocked matrix state (see <see cref="LockAbleObject.Lock"/>, <see cref="LockAbleObject.Unlock"/>);
        /// </remarks>
        virtual public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {

            if (acc.Count < this.RowPartitioning.LocalLength)
                throw new ArgumentException("array is too short - must be as least as big as the local length of the row partition.", "acc");
            if (a.Count < this.ColPartition.LocalLength)
                throw new ArgumentException("array is too short - must be as least as big as the local length of the column partition.", "a");
            if (object.ReferenceEquals(a, acc))
                throw new ArgumentException("in-place computation is not supported.", "a,acc");

            // create vector objects
            bool dummy;
            VectorBase _a = CreateVec(a, this.ColPartition, out dummy);
            using (VectorBase.CommVector _a_comm = _a.CreateCommVector(this)) {
                bool notWriteBackReq;
                VectorBase _acc = CreateVec(acc, this.RowPartitioning, out notWriteBackReq);

                // lock objects
                this.Lock();
                _a.Lock();
                _acc.Lock();

                // check args
                if (!_a.Part.Equals(this.ColPartition))
                    throw new ArgumentException("mismatch between column partition and partition of a.", "a");
                if (!_acc.Part.Equals(this.RowPartitioning))
                    throw new ArgumentException("mismatch between row partition and partition of acc.", "acc");

                // real work:
                SpMV_Expert(alpha, _a_comm, beta, _acc);

                // unlock
                _a.Unlock();
                _acc.Unlock();
                this.Unlock();

                // copy back result (if required)
                if (!notWriteBackReq) {
                    _acc.GetValues(acc, 0, 0, this.RowPartitioning.LocalLength);
                }
            }
        }
        
        /// <summary>
        /// used by <see cref="GetDiagonalElement"/> and <see cref="SetDiagonalElement"/>;
        /// </summary>
        int[] m_DiagElementPointer;

        /// <summary>
        /// gets the diagonal element, if present, for row number <paramref name="row"/>
        /// </summary>
        /// <param name="row"> global (''over all MPI processes'') index</param>
        /// <remarks>
        /// works only in unlocked matrix state (see <see cref="LockAbleObject.Lock"/>, <see cref="LockAbleObject.Unlock"/>);
        /// </remarks>
        virtual public double GetDiagonalElement(int row) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_RowPart.TestIfInLocalRange(row);
            int rloc = m_RowPart.TransformIndexToLocal(row);

            if (m_DiagElementPointer == null) {
                m_DiagElementPointer = new int[m_RowPart.LocalLength];
                ArrayTools.SetAll(m_DiagElementPointer, int.MinValue);
            }

            if (m_DiagElementPointer[rloc] < 0) {
                m_DiagElementPointer[rloc] = m_LocalMtx.GetEntryIndex(rloc, row - (int) m_ColPart.i0); // cache the pointer;
                                                                                         // altering diagonal elments is used very often (e.g. implicit euler),
                                                                                         // so a special optimization seems to be justified.
            }

            if (m_DiagElementPointer[rloc] < 0)
                // diag entry is eiter not present, or in some external matrix
                // -- we assume this is a seldome case => let someone else do the job
                return GetValues(row, new int[] { row })[0];
                
            else
                return m_LocalMtx.Val[m_DiagElementPointer[rloc]];
        }

        /// <summary>
        /// gets the diagonal element, if present, for row number <paramref name="row"/>
        /// </summary>
        /// <param name="row">global (''over all MPI processes'') index</param>
        /// <param name="val">new value</param>
        /// <remarks>
        /// works only in unlocked matrix state (see <see cref="LockAbleObject.Lock"/>, <see cref="LockAbleObject.Unlock"/>);
        /// </remarks>
        virtual public void SetDiagonalElement(int row, double val) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_RowPart.TestIfInLocalRange(row);
            int rloc = m_RowPart.TransformIndexToLocal(row);

            if (m_DiagElementPointer == null) {
                m_DiagElementPointer = new int[m_RowPart.LocalLength];
                ArrayTools.SetAll(m_DiagElementPointer, int.MinValue);
            }

            if (m_DiagElementPointer[rloc] < 0) {
                m_DiagElementPointer[rloc] = m_LocalMtx.GetEntryIndex(rloc, row - (int)m_ColPart.i0); // cache the pointer;
                                                                                           // altering diagonal elments is used very often (e.g. implicit euler),
                                                                                           // so a special optimization seems to be justified.
            }

            if (m_DiagElementPointer[rloc] < 0)
                // diag entry is eiter not present, or in some external matrix
                // -- we assume this is a seldome case => let someone else do the job
                SetValues(row, new int[] { row }, new double[] { val });
            else
                m_LocalMtx.Val[m_DiagElementPointer[rloc]] = val;
        }

        #endregion

        #region IMutuableMatrixEx Members

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetOccupiedColumnIndices"/>
        /// </summary>
        public int GetOccupiedColumnIndices(int RowIndex, ref int[] ret) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_helper.Update(RowIndex, m_LocalMtx, ExtMatrix, m_RowPart, (int)m_ColPart.i0);

            if(ret == null || ret.Length < m_helper.UsedLen)
                ret = new int[m_helper.UsedLen];
            Array.Copy(m_helper.ColIndices, 0, ret, 0, m_helper.UsedLen);
            return m_helper.UsedLen;
        }

        /// <summary>
        /// see <see cref="IMutableMatrixEx.GetRow"/>
        /// </summary>
        public int GetRow(int RowIndex, ref int[] retCol, ref double[] retVal) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_helper.Update(RowIndex, m_LocalMtx, ExtMatrix, m_RowPart, (int)m_ColPart.i0);

            if (retCol == null || retCol.Length < m_helper.UsedLen)
                retCol = new int[m_helper.UsedLen];
            if (retVal == null || retVal.Length < m_helper.UsedLen)
                retVal = new double[m_helper.UsedLen];

            for (int i = 0; i < m_helper.UsedLen; i++) {
                retCol[i] = m_helper.ColIndices[i];
                retVal[i] = m_helper.Values[i];
            }

            return m_helper.UsedLen;
        }

        #endregion

        #region IMutuableMatrix Members

        /// <summary>
        /// Helper sutructure, used mainly by the methods that are defined by the
        /// <see cref="IMutableMatrixEx"/>-interface.
        /// </summary>
        struct Helper {

            public int[] ColIndices;
            public int[] PointersIntoVal;
            public double[] Values;
            public int row;
            public int UsedLen;
            public bool _1stUse; // if false, the structure is not initialized and 
            public int MinCol;
            public int MaxCol;
            public int[] OwnerProcRank;

            public void Update(int RowIndex, FormatBase m_LocalMtx, IDictionary<int, External> ExtMatrix, IPartitioning rowPart, int Col0) {

                rowPart.TestIfInLocalRange(RowIndex);
                int _row = rowPart.TransformIndexToLocal(RowIndex);
            
                if( !_1stUse || row != _row) {
                    m_LocalMtx.GetAllOccupiedColumns(_row, ref ColIndices, ref PointersIntoVal, ref Values, out UsedLen);
                    row = _row;
                    _1stUse = true;

                    // find Min And Max Col;
                    MinCol = int.MinValue;
                    MaxCol = int.MaxValue;
                    for (int i = 0; i < UsedLen; i++) {
                        MinCol = Math.Max(MinCol, ColIndices[i]);
                        MaxCol = Math.Max(MaxCol, ColIndices[i]);
                    }

                    // filter double values;
                    if (m_LocalMtx is CCBCSR || m_LocalMtx is ELLPACKlike) {
                        // these matrices may have unused 'dummy' entries (to fill up memory)

                        for (int i = 0; i < UsedLen; i++) {
                            int MtxCol = ColIndices[i];
                            if (Array.IndexOf<int>(ColIndices, MtxCol, 0, i) >= 0) {
                                // found a dummy

                                UsedLen++;
                                for (int ii = i; ii < UsedLen; ii++) {
                                    ColIndices[ii] = ColIndices[ii + 1];
                                    PointersIntoVal[ii] = PointersIntoVal[ii + 1];
                                    Values[ii] = Values[ii + 1];
                                }
                            }
                        }
                    }

                    if (OwnerProcRank == null || OwnerProcRank.Length < UsedLen)
                        OwnerProcRank = new int[UsedLen];
                    int rnk = rowPart.MpiRank;
                    for (int i = 0; i < UsedLen; i++)
                        OwnerProcRank[i] = rnk;

                    // transform column indices to global indices (for local matrix)...
                    for (int i = 0; i < UsedLen; i++) {
                        ColIndices[i] += Col0;
                    }
                                    
                    // find entries in external cells
                    foreach (KeyValuePair<int,External> sdfjbh in ExtMatrix) {
                        rnk = sdfjbh.Key;
                        External extm = sdfjbh.Value;

                        int k = Array.IndexOf<int>(extm.rowInd, row);
                        if (k >= 0) {
                            // need to look into this external matrix

                            int rSt = extm.RowStart[k];
                            int rEn = extm.RowStart[k + 1];

                            for (int i = rSt; i < rEn; i++) {

                                // resize arrays ...
                                if (UsedLen > PointersIntoVal.Length - 10)
                                    Array.Resize(ref PointersIntoVal, PointersIntoVal.Length + 10);
                                if (UsedLen > ColIndices.Length - 10)
                                    Array.Resize(ref ColIndices, ColIndices.Length + 10);
                                if (UsedLen > Values.Length - 10)
                                    Array.Resize(ref Values, Values.Length + 10);
                                if (UsedLen > OwnerProcRank.Length - 10)
                                    Array.Resize(ref Values, OwnerProcRank.Length + 10);

                                // record entry of external matrix
                                OwnerProcRank[UsedLen] = rnk;
                                PointersIntoVal[UsedLen] = i;
                                Values[UsedLen] = extm.Val[i];
                                ColIndices[UsedLen] = extm.GlobalColInd[i];
                            }
                        }
                    }
                }
            }
        }

        Helper m_helper;

        /// <summary>
        /// see <see cref="IMutableMatrix.GetValues"/>
        /// </summary>
        public double[] GetValues(int RowIndex, int[] ColumnIndices) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_helper.Update(RowIndex, this.m_LocalMtx, this.ExtMatrix, this.RowPartitioning, (int)m_ColPart.i0);
            double[] ret = new double[ColumnIndices.Length];

            for (int i = 0; i < ret.Length; i++) {
                int ptr = Array.IndexOf<int>(m_helper.ColIndices, ColumnIndices[i], 0, m_helper.UsedLen);
                if (ptr < 0) {
                    ret[i] = 0;
                } else {
                    int proc = m_helper.OwnerProcRank[ptr];
                    if (proc == this.m_RowPart.MpiRank) {
                        ret[i] = m_LocalMtx.Val[ptr];
                    } else {
                        ret[i] = this.ExtMatrix[proc].Val[ptr];
                    }
                }
            }

            return ret;
        }

        /// <summary>
        /// see <see cref="IMutableMatrix.SetValues"/>
        /// </summary>
        public void SetValues(int RowIndex, int[] ColumnIndices, double[] newValues) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked, call illegal");

            m_helper.Update(RowIndex, this.m_LocalMtx, this.ExtMatrix, this.RowPartitioning, (int)m_ColPart.i0);


            for (int i = 0; i < newValues.Length; i++) {
                int ptr = Array.IndexOf<int>(m_helper.ColIndices, ColumnIndices[i], 0, m_helper.UsedLen);
                if (ptr < 0) {
                    throw new ArgumentException("entry (" + RowIndex + "," + ColumnIndices[i] + ") is not present in matrix (matrix does not support allocation of new entries);");
                } else {
                    int proc = m_helper.OwnerProcRank[ptr];
                    if (proc == this.m_RowPart.MpiRank) {
                        m_LocalMtx.Val[ptr] = newValues[i];
                    } else {
                        this.ExtMatrix[proc].Val[ptr] = newValues[i];
                    }
                }
            }
        }

        /// <summary>
        /// set/get an arbitrary entry; setting an entry where no space is allocated will produce an excepion.
        /// </summary>
        public double this[int i, int j] {
            get {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked, call illegal");
                if (i == j)
                    return GetDiagonalElement(i); // may be faster


                int ptr = this.m_LocalMtx.GetEntryIndex(i, j);
                if (ptr < 0) {
                    return GetValues(i, new int[] { j })[0];
                } else {
                    return this.m_LocalMtx.Val[ptr];
                }
            }
            set {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked, call illegal");

                if (i == j)
                    SetDiagonalElement(i, value); // may be faster

                int ptr = this.m_LocalMtx.GetEntryIndex(i, j);
                if (ptr < 0) {
                    SetValues(i, new int[] { j }, new double[] { value });
                } else {
                    this.m_LocalMtx.Val[ptr] = value;
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
        /// for this class and subclasses, always false; see <see cref="IMutableMatrix.OccupationMutable"/>;
        /// </summary>
        public bool OccupationMutable {
            get { return false; }
        }

        /// <summary>
        /// Pseudo-implementation.
        /// </summary>
        public object Clone() {
            return new MsrMatrix(this);
        }

        #endregion
    }
}
