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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Reflection;
using System.Runtime.InteropServices;
using BoSSS.Foundation.Grid;
using BoSSS.Platform;
using ilPSP.Utils;
using MPI.Wrappers;
using ilPSP;

namespace BoSSS.Foundation.Comm {

    /// <summary>
    /// transceiver for an arbitrary vector
    /// </summary>
    /// <typeparam name="T">vector/array</typeparam>
    /// <typeparam name="V">item type</typeparam>
    public class VectorTransceiver<T, V>
        where T : IList<V>
        where V : struct {

        IGridData m_master;
        T m_vector;
        int m_ItemsPerCell;

        /// <summary>
        /// throws an exception if any member of <paramref name="t"/>
        /// is not suited for transport over the network (e.g. pointers);
        /// </summary>
        [Conditional("DEBUG")]
        private void CheckItemTypeRecursive(Type t) {
            if (t.IsPrimitive)
                return;
            if (!t.IsValueType)
                throw new NotSupportedException("Type T must be composed only of primitive types.");

            FieldInfo[] members = t.GetFields(BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.Instance);

            foreach (FieldInfo member in members) {
                CheckItemTypeRecursive(member.FieldType);
            }
        }


        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="master"></param>
        /// <param name="vector">a vector of length J*<paramref name="ItemsPerCell"/>, where J is the number of local cells (including ghost)</param>
        /// <param name="ItemsPerCell">the number of items that should be transmitted/received per cell</param>
        public VectorTransceiver(IGridData master, T vector, int ItemsPerCell) {
            CheckItemTypeRecursive(typeof(V));

            int J = master.iLogicalCells.NoOfCells;
            if (vector.Count != ItemsPerCell * J)
                throw new ArgumentException("wrong length of input vector.");

            m_ItemsPerCell = ItemsPerCell;
            m_vector = vector;
            m_master = master;

            //sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);

            //// set comm. paths
            //sms.SetCommPathsAndCommit(master.Parallel.ProcessesToSendTo);

            var Para = m_master.iParallel;
            var sndProc = Para.ProcessesToSendTo;
            var rvcProc = Para.ProcessesToReceiveFrom;
            SendBuffers = new V[sndProc.Length][];
            for (int i = 0; i < SendBuffers.Length; i++) {
                int p = sndProc[i];
                SendBuffers[i] = new V[Para.SendCommLists[p].Length * m_ItemsPerCell];
            }
            SendBufferPin = new GCHandle[sndProc.Length];
            rqst = new MPI_Request[sndProc.Length + rvcProc.Length];
            staTussies = new MPI_Status[rqst.Length];
            if (!(typeof(T).IsArray)) {
                RcvBuffer = new V[rvcProc.Length][];
                for (int i = 0; i < RcvBuffer.Length; i++) {
                    int p = rvcProc[i];
                    RcvBuffer[i] = new V[Para.RcvCommListsNoOfItems[p] * m_ItemsPerCell];
                }
                RcvBufferPin = new GCHandle[RcvBuffer.Length];
            }
        }

        ///// <summary>
        ///// it would be much faster to use something else ...
        ///// </summary>
        //SerialisationMessenger sms;

        V[][] SendBuffers;

        GCHandle[] SendBufferPin;

        MPI_Request[] rqst;

        MPI_Status[] staTussies;



        /// <summary>
        /// initiates the send/receive - processes and returns immediately;
        /// Every call of this method must be matched by a later call to <see cref="TransceiveFinish"/>;
        /// </summary>
        public void TransceiveStartImReturn() {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);


            var Para = m_master.iParallel;
            int mul = m_ItemsPerCell;

            /*
            for (int i = 0; i < Para.ProcessesToSendTo.Length; i++) {
                int proc = Para.ProcessesToSendTo[i];
                int[] CellIndexList = Para.m_SendCommLists[proc];
                V[] SendBuf = new V[CellIndexList.Length];
                
                for (int ii = 0; ii < SendBuf.Length; ii++)
                    for( int iii = 0; iii < mul; iii++)
                        SendBuf[ii*mul + iii] = m_vector[CellIndexList[ii]*mul + iii];

                sms.Transmitt(proc, SendBuf);
            }
             */


            int[] sndProc = Para.ProcessesToSendTo;


            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);

            Array.Clear(this.rqst, 0, this.rqst.Length);

            unsafe {

                int cnt = 0;
                int N = m_ItemsPerCell;
                int sz = Marshal.SizeOf(typeof(V));

                // Sending ...
                // -----------

                // over all processes to which we have to send data to ...
                for (int i = 0; i < sndProc.Length; i++) {

                    // destination processor and comm list
                    int pDest = sndProc[i];
                    var commList = Para.SendCommLists[pDest];
                    int Len = commList.Length;

                    // build&fill send buffer
                    var SendBuffer = SendBuffers[i];
                    if (N != 1) {
                        for (int l = 0; l < Len; l++)
                            for (int n = 0; n < N; n++)
                                SendBuffer[l * N + n] = this.m_vector[commList[l] * N + n];
                    } else {
                        for (int l = 0; l < Len; l++)
                            SendBuffer[l] = this.m_vector[commList[l]];
                    }

                    // MPI send
                    SendBufferPin[i] = GCHandle.Alloc(SendBuffers[i], GCHandleType.Pinned);
                    cnt++;
                    csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(SendBuffers[i], 0),
                        Len * sz * N, csMPI.Raw._DATATYPE.BYTE, pDest,
                        4442 + MyRank,
                        csMPI.Raw._COMM.WORLD,
                        out rqst[i]);
                }

            }
        }

        V[][] RcvBuffer = null;
        GCHandle[] RcvBufferPin;



        /// <summary>
        /// blocks until the send/receive - processes started by <see cref="TransceiveStartImReturn"/>
        /// are complete and returns;
        /// </summary>
        public void TransceiveFinish() {
            //ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            //var Para = m_master.Parallel;
            //var mul = this.m_ItemsPerCell;

            //int rcvProc;
            //V[] data;
            //while (sms.GetNext(out rcvProc, out data)) {
            //    int j_insert = Para.m_RcvCommListsInsertIndex[rcvProc];

            //    int L = data.Length;
            //    for( int i = 0; i < L; i++) {
            //        int iDest = i + j_insert;
            //        for( int iii = 0; iii < mul; iii++)
            //            m_vector[iDest*mul + iii] = data[i*mul + iii];
            //    }
            //}


            var Para = m_master.iParallel;
            int mul = m_ItemsPerCell;

            /*
            for (int i = 0; i < Para.ProcessesToSendTo.Length; i++) {
                int proc = Para.ProcessesToSendTo[i];
                int[] CellIndexList = Para.m_SendCommLists[proc];
                V[] SendBuf = new V[CellIndexList.Length];
                
                for (int ii = 0; ii < SendBuf.Length; ii++)
                    for( int iii = 0; iii < mul; iii++)
                        SendBuf[ii*mul + iii] = m_vector[CellIndexList[ii]*mul + iii];

                sms.Transmitt(proc, SendBuf);
            }
             */

            int MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);


            var rvcProc = Para.ProcessesToReceiveFrom;
            var sndProc = Para.ProcessesToSendTo;

            unsafe {
                V[] _vector = m_vector as V[];

                GCHandle _vector_Handle = default(GCHandle);


                try {
                    int N = m_ItemsPerCell;

                    byte* pCA = default(byte*);
                    if (_vector != null) {
                        _vector_Handle = GCHandle.Alloc(_vector, GCHandleType.Pinned);
                        pCA = (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(_vector, 0);
                    }



                    int sz = Marshal.SizeOf(typeof(V));

                    // Receiving ...
                    // -------------

                    // over all processes from which we receive data...
                    for (int i = 0; i < rvcProc.Length; i++) {

                        // Source processor and insert index and no of elements to receive ...
                        int pOrigin = rvcProc[i];
                        int iInsert = Para.RcvCommListsInsertIndex[pOrigin];
                        int Len = Para.RcvCommListsNoOfItems[pOrigin];

                        IntPtr insertAddr;
                        if (_vector != null)
                            insertAddr = (IntPtr)(pCA + iInsert * N * sz);
                        else {
                            RcvBufferPin[i] = GCHandle.Alloc(RcvBuffer[i]);
                            insertAddr = Marshal.UnsafeAddrOfPinnedArrayElement(RcvBuffer[i], 0);
                        }

                        // MPI receive
                        csMPI.Raw.Irecv(insertAddr,
                            Len * N * sz, csMPI.Raw._DATATYPE.BYTE, pOrigin,
                            4442 + pOrigin,
                            csMPI.Raw._COMM.WORLD,
                            out rqst[i + sndProc.Length]);
                    }

                    // Wait for comm to finish
                    // -----------------------

                    Array.Clear(staTussies, 0, staTussies.Length);
                    csMPI.Raw.Waitall(rqst.Length, rqst, staTussies);
                } finally {
                    // release GC handles
                    // ==================


                    for (int i = 0; i < SendBufferPin.Length; i++)
                        SendBufferPin[i].Free();


                    if (_vector != null) {
                        _vector_Handle.Free();
                    } else {
                        for (int i = 0; i < RcvBufferPin.Length; i++)
                            RcvBufferPin[i].Free();
                    }
                }
            }

        }
    }

    /// <summary>
    /// extension methods for MPI communication.
    /// </summary>
    public static class VectorTransceiver_Ext {

        /// <summary>
        /// exchange of an <see cref="MultidimensionalArray"/>
        /// </summary>
        /// <param name="master"></param>
        /// <param name="vector">
        /// The array must be continuous and have zero offset;
        /// first dimension must match local number of cells (including external).
        /// </param>
        static public void MPIExchange(this MultidimensionalArray vector, IGridData master) {

            int[] i0 = new int[vector.Dimension];
            int offset = vector.Index(i0);

            if (!vector.IsContinious || !(offset == 0) || !(vector.Storage.Length == vector.Length))
                throw new NotSupportedException();

            if (vector.GetLength(0) != master.iLogicalCells.NoOfCells)
                throw new ArgumentException("fist dimension must match number of local cells (including external)");

            double[] Stor = vector.Storage;
            int J = master.iLogicalCells.NoOfCells;
            int L = Stor.Length;
            Debug.Assert(L % J == 0);

            var Trx = new VectorTransceiver<double[], double>(master, Stor, L / J);
            Trx.TransceiveStartImReturn();
            Trx.TransceiveFinish();
        }

        /// <summary>
        /// most elegant way of MPI vector exchange....
        /// </summary>
        static public void MPIExchange(this double[] vector, IGridData master) {
            MPIExchange<double[], double>(vector, master);
        }

        /// <summary>
        /// most elegant way of MPI vector exchange....
        /// </summary>
        /// <typeparam name="Tlst"></typeparam>
        /// <typeparam name="Titm"></typeparam>
        /// <param name="vector"></param>
        /// <param name="master"></param>
        static public void MPIExchange<Tlst, Titm>(this Tlst vector, IGridData master)
            where Tlst : IList<Titm>
            where Titm : struct {
            int J = master.iLogicalCells.NoOfCells;
            int L = vector.Count;
            if (L % J != 0)
                throw new ArgumentException("Length of vector must be equal or a multiple of number of cells (including external).");

            var Trx = new VectorTransceiver<Tlst, Titm>(master, vector, L / J);
            Trx.TransceiveStartImReturn();
            Trx.TransceiveFinish();
        }


        /// <summary>
        /// MPI update of a bit-array
        /// </summary>
        /// <param name="b"></param>
        /// <param name="GridData"></param>
        static public void MPIExchange(this BitArray b, IGridData GridData) {
            if (b.Length != GridData.iLogicalCells.NoOfCells) {
                throw new ArgumentException("length must be equal to number of cells.", "b");
            }


            if (GridData.CellPartitioning.MpiSize > 1) {

                // external 
                {
                    int rank, size;
                    csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
                    csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);

                    // setup messenger
                    SerialisationMessenger sms = new SerialisationMessenger(csMPI.Raw._COMM.WORLD);
                    sms.SetCommPathsAndCommit(GridData.iParallel.ProcessesToSendTo);

                    // send data
                    for (int p = 0; p < size; p++) {
                        int[] sendlist = GridData.iParallel.SendCommLists[p];
                        if (sendlist == null)
                            continue;

                        int L = sendlist.Length;
                        System.Collections.BitArray packet_for_p = new System.Collections.BitArray(L, false);
                        for (int l = 0; l < L; l++) {
                            packet_for_p[l] = b[sendlist[l]];
                        }

                        sms.Transmitt(p, packet_for_p);
                    }

                    // receive data
                    System.Collections.BitArray rcv_dat;
                    int rcv_rank;
                    while (sms.GetNext(out rcv_rank, out rcv_dat)) {
                        int insertAt = GridData.iParallel.RcvCommListsInsertIndex[rcv_rank];
                        if (GridData.iParallel.RcvCommListsNoOfItems[rcv_rank] != rcv_dat.Count)
                            throw new ApplicationException("internal error.");

                        int C = rcv_dat.Count;
                        for (int i = 0; i < C; i++) {
                            b[insertAt + i] = rcv_dat[i];
                        }
                    }

                    // dispose
                    sms.Dispose();
                }
            }
        }
    }
}
