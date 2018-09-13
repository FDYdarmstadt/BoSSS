using BoSSS.Foundation.Grid;
using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Multigrid {


    /// <summary>
    /// Parallelization/MPI data exchange for a vector.
    /// </summary>
    /// <typeparam name="T">vector/array</typeparam>
    public class MPIexchange<T>
        where T : IList<double> {

        IGridData m_master;
        T m_vector;
        MultigridMapping m_map;
       
        /// <summary>
        /// Locally stored data, handed over by the constructor. 
        /// </summary>
        public T Vector {
            get {
                return m_vector;
            }
        }

        double[] m_Vector_Ext;

        /// <summary>
        /// Data associated with external cells.
        /// </summary>
        public double[] Vector_Ext {
            get {
                return m_Vector_Ext;
            }
        }


        /// <summary>
        /// ctor.
        /// </summary>
        public MPIexchange(MultigridMapping map, T vector) {

            // misc init
            // =========

            IGridData master = map.AggGrid;
            int J = master.iLogicalCells.Count;
            if (vector.Count != map.LocalLength)
                throw new ArgumentException("wrong length of input vector.");

            
            m_vector = vector;
            m_master = master;
            m_map = map;

            var Para = m_master.iParallel;
            var rvcProc = Para.ProcessesToReceiveFrom;
            var sndProc = Para.ProcessesToSendTo;
            rqst = new MPI_Request[sndProc.Length + rvcProc.Length];
            
            // allocate send buffers
            // =====================
            {
                SendBuffers = new double[sndProc.Length][];
                for (int i = 0; i < SendBuffers.Length; i++) {
                    int p = sndProc[i];

                    // compute length of send list
                    int L = 0;
                    foreach (int jCell in Para.SendCommLists[p]) {
                        Debug.Assert(map.IsLocalBlock(jCell + map.FirstBlock));
                        Debug.Assert(map.GetLength(jCell) == map.GetBlockLen(jCell + map.FirstBlock));
                        L += map.GetLength(jCell);
                    }

                    // alloc send buffer
                    SendBuffers[i] = new double[L];
                }
                SendBufferPin = new GCHandle[sndProc.Length];
            }

            // allocate receive buffers
            // ========================

            {
                int totL = 0;
                RcvBuffer = new double[rvcProc.Length][];
                for (int i = 0; i < RcvBuffer.Length; i++) {
                    int p = rvcProc[i];

                    // compute length of receive list
                    int L = 0;
                    int J0 = Para.RcvCommListsInsertIndex[p];
                    int JE = Para.RcvCommListsNoOfItems[p] + J0;
                    for(int jCell = J0; jCell < JE; jCell++) {
                        Debug.Assert(jCell >= map.LocalNoOfBlocks);
                        L += map.GetLength(jCell);
                    }
                    totL += L;

                    // alloc internal receive buffer
                    RcvBuffer[i] = new double[L];
                }
                RcvBufferPin = new GCHandle[RcvBuffer.Length];

                m_Vector_Ext = new double[totL];
            }
        }

        /// <summary>
        /// Internal send buffers
        /// - 1st index: correlates with receiver process, as in <see cref="IParallelization.ProcessesToSendTo"/>.
        /// - 2nd index: determined by the send list, see <see cref="IParallelization.SendCommLists"/>
        /// </summary>
        double[][] SendBuffers;

        GCHandle[] SendBufferPin;

        MPI_Request[] rqst;

        

        /// <summary>
        /// initiates the send/receive - processes and returns immediately;
        /// Every call of this method must be matched by a later call to <see cref="TransceiveFinish"/>;
        /// </summary>
        public void TransceiveStartImReturn() {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);


            var Para = m_master.iParallel;
            


            int[] sndProc = Para.ProcessesToSendTo;


            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);

            Array.Clear(this.rqst, 0, this.rqst.Length);

            unsafe {

                
                
                
                // Sending ...
                // -----------

                // over all processes to which we have to send data to ...
                for (int i = 0; i < sndProc.Length; i++) {

                    // destination processor and comm list
                    int pDest = sndProc[i];
                    int[] commList = Para.SendCommLists[pDest];
                    int Len = commList.Length;

                    // fill send buffer
                    var SendBuffer = SendBuffers[i];

                    int cnt = 0;
                    for (int l = 0; l < Len; l++) {
                        int jCell = commList[l];
                        int N = m_map.GetLength(jCell);
                        int i0 = m_map.LocalUniqueIndex(0, jCell, 0);

                        for (int n = 0; n < N; n++) {
                            SendBuffer[cnt] = this.m_vector[i0 + n];
                            cnt++;
                        }
                    }
                    Debug.Assert(cnt == SendBuffers[i].Length);

                    // MPI send
                    SendBufferPin[i] = GCHandle.Alloc(SendBuffers[i], GCHandleType.Pinned);
                    csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(SendBuffers[i], 0),
                        SendBuffers[i].Length, csMPI.Raw._DATATYPE.DOUBLE, pDest,
                        4442 + MyRank,
                        csMPI.Raw._COMM.WORLD,
                        out rqst[i]);
                }

            }
        }

        double[][] RcvBuffer = null;

        /// <summary>
        /// GC pinning of <see cref="RcvBuffer"/>, index correlates with 1st index of <see cref="RcvBuffer"/>
        /// </summary>
        GCHandle[] RcvBufferPin;



        /// <summary>
        /// Blocks until the send/receive - processes started by <see cref="TransceiveStartImReturn"/> are complete and returns;
        /// The received data is **accumulated** in <see cref="Vector_Ext"/>.
        /// </summary>
        /// <param name="beta">
        /// Pre-scaling of <see cref="Vector_Ext"/>
        /// </param>
        public void TransceiveFinish(double beta) {
            
            var Para = m_master.iParallel;
            
            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);


            var rvcProc = Para.ProcessesToReceiveFrom;
            var sndProc = Para.ProcessesToSendTo;

            unsafe {
                try {

                    // Receiving ...
                    // -------------

                    // over all processes from which we receive data...
                    for (int i = 0; i < rvcProc.Length; i++) {

                        // Source processor and insert index and no of elements to receive ...
                        int pOrigin = rvcProc[i];
                        int Len = Para.RcvCommListsNoOfItems[pOrigin];

                        IntPtr insertAddr;
                        RcvBufferPin[i] = GCHandle.Alloc(RcvBuffer[i]);
                        insertAddr = Marshal.UnsafeAddrOfPinnedArrayElement(RcvBuffer[i], 0);
                        
                        // MPI receive
                        csMPI.Raw.Irecv(insertAddr,
                            RcvBuffer[i].Length, csMPI.Raw._DATATYPE.DOUBLE, pOrigin,
                            4442 + pOrigin,
                            csMPI.Raw._COMM.WORLD,
                            out rqst[i + sndProc.Length]);
                    }

                    // Wait for comm to finish
                    // -----------------------

                    //Array.Clear(staTussies, 0, staTussies.Length);
                    //csMPI.Raw.Waitall(rqst.Length, rqst, staTussies);
#if DEBUG
                    bool[] Markers = new bool[rqst.Length];
#endif

                    for (int iRqs = 0; iRqs < rqst.Length; iRqs++) {
                        csMPI.Raw.Waitany(rqst.Length, rqst, out int index, out MPI_Status stat);
#if DEBUG
                        Debug.Assert(Markers[iRqs] == false);
                        Markers[iRqs] = true;
#endif

                        if(index >= sndProc.Length) {
                            // ++++++++++++++++++++++++++++++++++++++++++
                            // data received: accumulate in output buffer
                            // ++++++++++++++++++++++++++++++++++++++++++
                            int pOrigin = rvcProc[index - sndProc.Length];
                            double[] buf = RcvBuffer[index - sndProc.Length];

                            int accOffset = this.m_map.LocalLength;
                            int J0 = Para.RcvCommListsInsertIndex[pOrigin];
                            int JE = Para.RcvCommListsNoOfItems[pOrigin] + J0;
                            int cnt = 0;
                            for(int jCell = J0; jCell < JE; jCell++) {
                                Debug.Assert(jCell >= m_map.LocalNoOfBlocks);
                                int N = m_map.GetLength(jCell);
                                int i0 = m_map.LocalUniqueIndex(0, jCell, 0);

                                for(int n = 0; n < N; n++) {
                                    int idxAcc = i0 + n - accOffset;
                                    m_Vector_Ext[idxAcc] = m_Vector_Ext[idxAcc] * beta + buf[cnt];
                                    cnt++;
                                }
                            }
                            Debug.Assert(cnt == buf.Length);
                            
                        } else {
                            // sending finished - noop
                        }
                    }

                } finally {
                    // release GC handles
                    // ==================


                    for (int i = 0; i < SendBufferPin.Length; i++)
                        SendBufferPin[i].Free();

                    for (int i = 0; i < RcvBufferPin.Length; i++)
                        RcvBufferPin[i].Free();
                }
            }

        }
    }



    /// <summary>
    /// Parallelization/MPI data exchange for a vector, in *the inverse direction*, i.e. data stored in external cells 
    /// is accumulated in locally owned cells of other processors.
    /// </summary>
    /// <typeparam name="T">vector/array</typeparam>
    public class MPIexchangeInverse<T>
        where T : IList<double> {

        IGridData m_master;
        T m_vector;
        MultigridMapping m_map;
       
        /// <summary>
        /// Locally stored data, handed over by the constructor. 
        /// </summary>
        public T Vector {
            get {
                return m_vector;
            }
        }

        double[] m_Vector_Ext;

        /// <summary>
        /// Data associated with external cells.
        /// </summary>
        public double[] Vector_Ext {
            get {
                return m_Vector_Ext;
            }
        }


        /// <summary>
        /// ctor.
        /// </summary>
        public MPIexchangeInverse(MultigridMapping map, T vector) {

            // misc init
            // =========

            IGridData master = map.AggGrid;
            int J = master.iLogicalCells.Count;
            if (vector.Count != map.LocalLength)
                throw new ArgumentException("wrong length of input vector.");

            
            m_vector = vector;
            m_master = master;
            m_map = map;

            var Para = m_master.iParallel;
            var rvcProc = Para.ProcessesToSendTo; // yes, this is intentional!
            var sndProc = Para.ProcessesToReceiveFrom;  // yes, this is intentional!
            rqst = new MPI_Request[sndProc.Length + rvcProc.Length];
            
            // allocate send buffers
            // =====================
            {
                SendBuffers = new double[sndProc.Length][];
                int totL = 0;
                for (int i = 0; i < SendBuffers.Length; i++) {
                    int p = sndProc[i];

                    // compute length of send list
                    int L = 0;
                    int J0 = Para.RcvCommListsInsertIndex[p];
                    int JE = Para.RcvCommListsNoOfItems[p] + J0;
                    for(int jCell = J0; jCell < JE; jCell++) {
                        Debug.Assert(jCell >= map.LocalNoOfBlocks);
                        L += map.GetLength(jCell);
                    }
                    totL += L;

                    

                    // alloc send buffer
                    SendBuffers[i] = new double[L];
                }
                SendBufferPin = new GCHandle[sndProc.Length];
                m_Vector_Ext = new double[totL];
            }

            // allocate receive buffers
            // ========================

            {
                RcvBuffer = new double[rvcProc.Length][];
                for (int i = 0; i < RcvBuffer.Length; i++) {
                    int p = rvcProc[i];

                    // compute length of receive list
                    int L = 0;
                    foreach (int jCell in Para.SendCommLists[p]) {
                        Debug.Assert(map.IsLocalBlock(jCell + map.FirstBlock));
                        Debug.Assert(map.GetLength(jCell) == map.GetBlockLen(jCell + map.FirstBlock));
                        L += map.GetLength(jCell);
                    }

                    // alloc internal receive buffer
                    RcvBuffer[i] = new double[L];
                }
                RcvBufferPin = new GCHandle[RcvBuffer.Length];
            }
        }

        /// <summary>
        /// Internal send buffers
        /// - 1st index: correlates with receiver process, as in <see cref="IParallelization.ProcessesToReceiveFrom"/>.
        /// - 2nd index: determined by the sequence of external cells
        /// </summary>
        double[][] SendBuffers;

        GCHandle[] SendBufferPin;

        MPI_Request[] rqst;

        

        /// <summary>
        /// initiates the send/receive - processes and returns immediately;
        /// Every call of this method must be matched by a later call to <see cref="TransceiveFinish"/>;
        /// </summary>
        public void TransceiveStartImReturn() {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);


            var Para = m_master.iParallel;
            


            int[] sndProc = Para.ProcessesToReceiveFrom; // yes, intentional


            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);

            Array.Clear(this.rqst, 0, this.rqst.Length);

            unsafe {
                
                
                // Sending ...
                // -----------

                // over all processes to which we have to send data to ...
                for (int i = 0; i < sndProc.Length; i++) {

                    // destination processor and external cell range
                    int pDest = sndProc[i];
                    int J0 = Para.RcvCommListsInsertIndex[pDest];
                    int JE = Para.RcvCommListsNoOfItems[pDest] + J0;
                    

                    // fill send buffer
                    var SendBuffer = SendBuffers[i];
                    int extOffset = m_map.LocalLength;

                    int cnt = 0;
                    for (int jCell = J0; jCell < JE; jCell++) {
                        Debug.Assert(jCell >= m_master.iLogicalCells.NoOfLocalUpdatedCells);

                        int N = m_map.GetLength(jCell);
                        int i0 = m_map.LocalUniqueIndex(0, jCell, 0);

                        for (int n = 0; n < N; n++) {
                            SendBuffer[cnt] = this.m_Vector_Ext[i0 + n - extOffset];
                            cnt++;
                        }
                    }
                    Debug.Assert(cnt == SendBuffers[i].Length);

                    // MPI send
                    SendBufferPin[i] = GCHandle.Alloc(SendBuffers[i], GCHandleType.Pinned);
                    csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(SendBuffers[i], 0),
                        SendBuffers[i].Length, csMPI.Raw._DATATYPE.DOUBLE, pDest,
                        4442 + MyRank,
                        csMPI.Raw._COMM.WORLD,
                        out rqst[i]);
                }

            }
        }

        double[][] RcvBuffer = null;

        /// <summary>
        /// GC pinning of <see cref="RcvBuffer"/>, index correlates with 1st index of <see cref="RcvBuffer"/>
        /// </summary>
        GCHandle[] RcvBufferPin;



        /// <summary>
        /// Blocks until the send/receive - processes started by <see cref="TransceiveStartImReturn"/> are complete and returns;
        /// The received data is **accumulated** in <see cref="Vector"/>.
        /// </summary>
        /// <param name="beta">
        /// Pre-scaling of <see cref="Vector_Ext"/>
        /// </param>
        public void TransceiveFinish(double beta) {
            
            var Para = m_master.iParallel;
            
            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);


            var rvcProc = Para.ProcessesToSendTo; // yes, send and receive roles are exchanged!
            var sndProc = Para.ProcessesToReceiveFrom;

            unsafe {
                try {

                    // Receiving ...
                    // -------------

                    // over all processes from which we receive data...
                    for (int i = 0; i < rvcProc.Length; i++) {

                        // Source processor and insert index and no of elements to receive ...
                        int pOrigin = rvcProc[i];
                        int Len = Para.RcvCommListsNoOfItems[pOrigin];

                        RcvBufferPin[i] = GCHandle.Alloc(RcvBuffer[i]);
                        IntPtr insertAddr = Marshal.UnsafeAddrOfPinnedArrayElement(RcvBuffer[i], 0);
                        
                        // MPI receive
                        csMPI.Raw.Irecv(insertAddr,
                            RcvBuffer[i].Length, csMPI.Raw._DATATYPE.DOUBLE, pOrigin,
                            4442 + pOrigin,
                            csMPI.Raw._COMM.WORLD,
                            out rqst[i + sndProc.Length]);
                    }

                    // Wait for comm to finish
                    // -----------------------

#if DEBUG
                    bool[] Markers = new bool[rqst.Length];
#endif

                    for (int iRqs = 0; iRqs < rqst.Length; iRqs++) {
                        csMPI.Raw.Waitany(rqst.Length, rqst, out int index, out MPI_Status stat);
#if DEBUG
                        Debug.Assert(Markers[iRqs] == false);
                        Markers[iRqs] = true;
#endif

                        if(index >= sndProc.Length) {
                            // ++++++++++++++++++++++++++++++++++++++++++
                            // data received: accumulate in output buffer
                            // ++++++++++++++++++++++++++++++++++++++++++
                            int pOrigin = rvcProc[index - sndProc.Length];
                            double[] buf = RcvBuffer[index - sndProc.Length];

                            int[] InsertIndices = Para.SendCommLists[pOrigin];
                            int L = InsertIndices.Length;


                            int cnt = 0;


                            for (int l = 0; l < L; l++) {
                                int jCell = InsertIndices[l];
                                Debug.Assert(jCell < m_map.LocalNoOfBlocks);

                                int N = m_map.GetLength(jCell);
                                int i0 = m_map.LocalUniqueIndex(0, jCell, 0);

                                for(int n = 0; n < N; n++) {
                                    int idxAcc = i0 + n;
                                    m_vector[idxAcc] = m_vector[idxAcc] * beta + buf[cnt];
                                    cnt++;
                                }
                            }
                            Debug.Assert(cnt == buf.Length);
                            
                        } else {
                            // sending finished - noop
                        }
                    }

                } finally {
                    // release GC handles
                    // ==================


                    for (int i = 0; i < SendBufferPin.Length; i++)
                        SendBufferPin[i].Free();

                    for (int i = 0; i < RcvBufferPin.Length; i++)
                        RcvBufferPin[i].Free();
                }
            }

        }
    }
}
