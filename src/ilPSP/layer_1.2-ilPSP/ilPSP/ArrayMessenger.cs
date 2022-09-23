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
using System.Collections;
using System.IO;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;
using System.Runtime.InteropServices;
using System.Diagnostics;

using MPI.Wrappers;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;
using System.Drawing;
using System.Reflection;

namespace ilPSP.Utils {
    
    
    /// <summary>
    /// This class provides an easy way of arrays of value types of arbitrary length
    /// objects between MPI processes. Internally, it uses nonblocking MPI point-to-point 
    /// routines.
    /// </summary>
    /// <remarks>
    /// This object is used in the following way:
    /// <list type="bullet">
    ///   <item>
    ///   First, a process needs to specify to which processes it wants to send objects to;<br/>
    ///   This is done by multiple calls to <see cref="SetCommPath"/>;
    ///   After all target processors are specified, the so-called setup phase is finished by
    ///   a single call to <see cref="CommitCommPaths"/>;
    ///   </item>
    ///   <item>
    ///   Now, the messenger is ready to send objects to other processes, by multiple calls
    ///   to <see cref="Transmit"/>; Of course, these calls must match the <see cref="SetCommPath"/>-calls
    ///   in the previous step;
    ///   </item>
    ///   <item>
    ///   After the transmission of objects has been initialized, this class is ready to wait for incomming objects;
    ///   This is done by multiple calls to <see cref="GetNext"/>; Usually, <see cref="GetNext"/> is called in an
    ///   loop that terminates when <see cref="GetNext"/> returns null;
    ///   </item>
    /// </list>
    /// The last two steps can be repeated;
    /// </remarks>
    public class ArrayMessenger<T> : IDisposable 
        where T : struct    
    {

        static int TagCnt = 4567231;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_MPI_comm">the MPI communicator for this messenger (see <see cref="MPI_comm"/>);</param>
        public ArrayMessenger(MPI_Comm _MPI_comm) {
            MPICollectiveWatchDog.Watch(_MPI_comm);


            m_MPI_comm = _MPI_comm;
            csMPI.Raw.Comm_Rank(_MPI_comm, out m_Rank);
            csMPI.Raw.Comm_Size(_MPI_comm, out m_Size);

            m_MyCommPaths = new int[m_Size];
            m_AllCommPaths = new int[m_Size, m_Size];


            m_Requests = new MPI_Request[m_Size * 4];
            m_TransmittCalled = new BitArray(m_Size, false);

            unsafe {
                int myTag = TagCnt;
                TagCnt += 17;

                if (m_Rank == 0)
                    m_MyTagOffset = myTag;

                csMPI.Raw.Bcast((IntPtr)(&myTag), 1, csMPI.Raw._DATATYPE.INT, 0, m_MPI_comm);

                if (m_Rank > 0)  
                    m_MyTagOffset = myTag;
            }
        }

        /// <summary>
        /// offset to MPI send - tags, so that MPI can distinct between more tan one instance of
        /// this object at one time
        /// </summary>
        int m_MyTagOffset;


        /// <summary>
        /// the MPI communicator in which this messenger is acting;
        /// </summary>
        public MPI_Comm MPI_comm {
            get { return m_MPI_comm; }
        }

        /// <summary>
        /// see <see cref="MPI_comm"/>;
        /// </summary>
        MPI_Comm m_MPI_comm;

        /// <summary>
        /// see <see cref="Rank"/>;
        /// </summary>
        int m_Rank;

        /// <summary>
        /// rank of the current process within  the MPI communicator <see cref="MPI_comm"/>;
        /// </summary>
        int Rank {
            get { return m_Rank; }
        }

        /// <summary>
        /// see <see cref="Size"/>;
        /// </summary>
        int m_Size;

        /// <summary>
        /// size (i.e. number of processors) in the MPI communicator <see cref="MPI_comm"/>;
        /// </summary>
        public int Size {
            get { return m_Size; }
        }



        /// <summary>
        /// setup procedure;
        /// Informs the messenger that objects should be transmitted to process <paramref name="TargetProcRank"/>;
        ///  1st function (multiple calls) in the setup process; If all communication paths
        /// are set, <see cref="CommitCommPaths"/> is the next method to call;
        /// </summary>
        /// <param name="TargetProcRank"></param>
        public void SetCommPath(int TargetProcRank, int BufferLength) {
            if (m_CommPathsCommited)
                throw new ApplicationException("setup phase is already finished.");
            if (TargetProcRank == m_Rank)
                throw new ArgumentException("sending data to process itself is not possible.");
            m_MyCommPaths[TargetProcRank] = BufferLength;
        }

        /*
        /// <summary>
        /// performs <see cref="SetCommPath"/> for each element in <paramref name="TargetProcRanks"/>
        /// and calls <see cref="CommitCommPaths"/>;
        /// </summary>
        /// <param name="TargetProcRanks"></param>
        public void SetCommPathsAndCommit(IEnumerable<int> TargetProcRanks) {
            foreach (int p in TargetProcRanks)
                this.SetCommPath(p);
            this.CommitCommPaths();
        }
        */

        /// <summary>
        /// an entry unequal to 0 at index p indicates that objects should be transmitted from 
        /// this process to process p.
        /// </summary>
        int[] m_MyCommPaths;


        /// <summary>
        /// A value unequal to 0 at entry [i,j] indicates that processor j (target process) will receive
        /// an object from processor i (source process);
        /// 1st index: source process
        /// 2nd index: target process
        /// </summary>
        private int[,] m_AllCommPaths;

        
        /// <summary>
        /// setup procedure;
        /// 2nd function (exactly one call in the whole lifecycle of this object is permitted and 
        /// required) in the setup process;
        /// In normal use, called after <see cref="SetCommPath"/>;
        /// </summary>
        public void CommitCommPaths() {

            // -------------------------------------------
            // send communication paths to other processes 
            // -------------------------------------------
            if (m_CommPathsCommited)
                throw new ApplicationException("setup phase is already finished.");

            GCHandle myCommPathsPin = GCHandle.Alloc(m_MyCommPaths, GCHandleType.Pinned);
            GCHandle allCommPathsPin = GCHandle.Alloc(m_AllCommPaths, GCHandleType.Pinned);

            IntPtr pSndBuf = Marshal.UnsafeAddrOfPinnedArrayElement(m_MyCommPaths, 0);
            IntPtr pRcvBuf = Marshal.UnsafeAddrOfPinnedArrayElement(m_AllCommPaths, 0);

            //if (m_Rank == 0)
            //    Debugger.Break();

            csMPI.Raw.Allgather(pSndBuf, m_MyCommPaths.Length, csMPI.Raw._DATATYPE.INT,
                             pRcvBuf, m_MyCommPaths.Length, csMPI.Raw._DATATYPE.INT,
                             this.m_MPI_comm);

            myCommPathsPin.Free();
            allCommPathsPin.Free();
            m_CommPathsCommited = true;

            // -------------------
            // create send buffers
            // -------------------

            for (int p = 0; p < m_Size; p++)
                if (m_MyCommPaths[p] != 0) {
                    m_SendBuffers.Add(p, null);
                }


            // -------------------------------------
            // create key values for receive buffers
            // -------------------------------------

            for (int p = 0; p < m_Size; p++)
                if (m_AllCommPaths[p, m_Rank] != 0)
                    m_ReceiveBuffers.Add(p, null);
        }



        /// <summary>
        /// - keys: processor rank "p";
        /// - values: memory stream for the object that has to be send to process "p";
        /// </summary>
        SortedDictionary<int, T[]> m_SendBuffers = new SortedDictionary<int, T[]>();

        /// <summary>
        /// - keys: processor rank "p";
        /// - values: garbage collector pin handle for the buffer of memory stream at key "p" in <see cref="m_SendBuffers"/>;
        /// </summary>
        SortedDictionary<int, GCHandle> m_SendBuffersPin = new SortedDictionary<int, GCHandle>();


        /// <summary>
        /// - keys: processor rank "p";
        /// - values: memory block for the object that is received from process "p";
        /// </summary>
        SortedDictionary<int, T[]> m_ReceiveBuffers = new SortedDictionary<int,T[]>();
        
        /// <summary>
        /// - keys: processor rank "p";
        /// - values: garbage collector pin handle for the buffer memory at key "p" in <see cref="m_ReceiveBuffers"/>;
        /// </summary>
        SortedDictionary<int, GCHandle> m_ReceiveBuffersPin = new SortedDictionary<int, GCHandle>();

        
        /// <summary>
        /// true indicates that the setup is finished
        /// </summary>
        bool m_CommPathsCommited;

        /// <summary>
        /// see <see cref="TransmissionInProgress"/>;
        /// is set to true by <see cref="InitTransmission"/>
        /// </summary>
        bool m_TransmissionInProgress = false;

        /// <summary>
        /// true indicates that transmission is in progress, and the object 
        /// cannot be disposed (calling <see cref="Dispose"/> will result in an exception);
        /// </summary>
        public bool TransmissionInProgress { get { return m_TransmissionInProgress; } }


        /// <summary>
        /// when the transmission is in progress (<see cref="m_TransmissionInProgress"/> is true),
        /// a true entry at index p indicates that <see cref="Transmit"/> for process p has allready been called.
        /// The only purpose of this array is to detect double <see cref="Transmit"/>-calls to the same target process.
        /// </summary>
        BitArray m_TransmittCalled;


        /// <summary>
        /// 1st phase of the transmission process; must be called after <see cref="CommitCommPaths"/>.
        /// </summary>
        /// <param name="TargetProc"></param>
        /// <param name="data">
        /// array, which will be send to process <paramref name="TargetProc"/>
        /// by using serialization.
        /// </param>
        public void Transmit(int TargetProc, T[] data) {

            // ------------------------------
            // check for correctness of usage
            // ------------------------------

            if (data == null)
                throw new ArgumentNullException("cannot sent 'null' - objects.");
            if (!m_CommPathsCommited) 
                throw new ApplicationException("communication paths have to be committed first.");
            if (!m_SendBuffers.ContainsKey(TargetProc))
                throw new ArgumentException("no communication path was set for specified processor.");
            if (m_MyCommPaths[TargetProc] != data.Length)
                throw new ArgumentException("length of `data` array does not match the size specified before.");
            if ( m_TransmittCalled[TargetProc] == true)
                throw new ApplicationException("Transmit was already called for specified target process.");
            m_TransmittCalled[TargetProc] = true;

            // -------
            // startup
            // -------

            if (!m_TransmissionInProgress)
                InitTransmission();

            // ------------
            // send objects
            // ------------


            // send buffer content
            if (m_SendBuffersPin.ContainsKey(TargetProc))
                m_SendBuffersPin[TargetProc] = GCHandle.Alloc(data, GCHandleType.Pinned);
            else
                m_SendBuffersPin.Add(TargetProc, GCHandle.Alloc(data, GCHandleType.Pinned));
            
            
            if (data.Length > 0) {
                csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(data, 0),
                              data.Length*sizeof_T, csMPI.Raw._DATATYPE.BYTE,
                              TargetProc, TagBufferContent + m_MyTagOffset,
                              m_MPI_comm,
                              out m_Requests[TargetProc]);
            }

            //if (DiagnosisFile != null) {
            //    File.WriteAllBytes(DiagnosisFile + "-smsSend-" + m_Rank + "-" + TargetProc + ".bin", Buffer);
            //}
        }

        /*

        /// <summary>
        /// entry at index p is reserved for the size of the object that is send to process p in bytes;
        /// </summary>
        int[] m_SendObjectSizes;


        /// <summary>
        /// pin handle for <see cref="m_SendObjectSizes"/>;
        /// </summary>
        GCHandle m_SendObjectSizesPin;

        /// <summary>
        /// entry at index p is reseved for the size of object received from process p in bytes;
        /// </summary>
        int[] m_ReceiveObjectSizes;


        /// <summary>
        /// pin handle for <see cref="m_ReceiveObjectSizes"/>;
        /// </summary>
        GCHandle m_ReceiveObjectSizesPin;


        /// <summary>
        /// MPI tag indicating that the size of a buffer is transmitted
        /// </summary>
        const int TagBufferSize = 123;

        */

        /// <summary>
        /// MPI tag indicating that the content of a buffer is transmitted
        /// </summary>
        const int TagBufferContent = 234;

        /// <summary>
        /// MPI request handles for nonblocking send and receives.
        /// Array of length of 2*<see cref="Size"/>:
        /// - 1st half: (index 0 to <see cref="Size"/>-1): send handles for buffer content
        /// - 2nd half: (index <see cref="Size"/> to 2*<see cref="Size"/>-1): receive handles for buffer content;
        /// </summary>
        MPI_Request[] m_Requests;

        int sizeof_T {
            get {
                int sz = Marshal.SizeOf<T>();
                if (typeof(T) == typeof(int))
                    if (sz != sizeof(int))
                        throw new ApplicationException("wrong size for int");
                if (typeof(T) == typeof(long))
                    if (sz != sizeof(long))
                        throw new ApplicationException("wrong size for long");

                return sz;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        private void InitTransmission() {

            if (m_TransmissionInProgress)
                throw new ApplicationException("internal error: transmission already in progress.");
            m_TransmissionInProgress = true;

            //m_SendObjectSizesPin = GCHandle.Alloc(m_SendObjectSizes, GCHandleType.Pinned);
            ReceiveCount = 0;
            // -------------------
            // clear request array
            // -------------------

            for (int i = 0; i < m_Requests.Length; i++) 
                m_Requests[i] = csMPI.Raw.MiscConstants.MPI_REQUEST_NULL;

            // ---------------------------
            // initiate receive of buffers
            // ---------------------------

            for (int p = 0; p < m_Size; p++) {
                if (m_AllCommPaths[p, m_Rank] != 0) {

                    int proc = p;
                    int Len = m_AllCommPaths[p, m_Rank];

                    // reallocate memory, if necessary
                    T[] buffer = m_ReceiveBuffers[proc];
                    if (buffer == null || buffer.Length != Len) {
                        buffer = new T[Len];
                        m_ReceiveBuffers[proc] = buffer;
                    }

                    // pin buffer
                    GCHandle pin = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                    if (m_ReceiveBuffersPin.ContainsKey(proc))
                        m_ReceiveBuffersPin[proc] = pin;
                    else
                        m_ReceiveBuffersPin.Add(proc, pin);

                    // start nonblocking receive
                    csMPI.Raw.Irecv(pin.AddrOfPinnedObject(),
                                     Len * sizeof_T, csMPI.Raw._DATATYPE.BYTE,
                                     proc, TagBufferContent + m_MyTagOffset,
                                     m_MPI_comm,
                                     out m_Requests[m_Size + proc]);
                }

            }
        }

        int ReceiveCount = 0;

        /// <summary>
        /// 2nd and last phase of the transmission process; must be called after all 
        /// calls to <see cref="Transmit"/> have been done;
        /// Must be called multiple times until <paramref name="o"/> is null
        /// after return, which guarantees that all send/receive's are finished
        /// (otherwise is very likely that MPI is in an undefined state).
        /// After this, the object is ready for another send/receive cycle.
        /// </summary>
        /// <param name="TargetProc">
        /// process rank which has send <paramref name="o"/> to this process.
        /// </param>
        /// <param name="o">
        /// on exit, the received array or null, if all objects were received;
        /// and the communication is finished.
        /// </param>
        /// <returns>
        /// true, if <paramref name="o"/> contains a received object;<br/>
        /// false if the communication is finished.
        /// </returns>
        public bool GetNext(out int TargetProc, out T[] o) {

            int size = m_Size;

            // --------------------------
            // check correctness of usage
            // --------------------------

            if (!m_CommPathsCommited)
                throw new ApplicationException("communication paths must be committed first.");

            if (!m_TransmissionInProgress) {
                // call InitTransmission for processes which don't send anything at all

                bool nosend = true;
                for (int p = 0; p < m_MyCommPaths.Length; p++)
                    if (m_MyCommPaths[p] != 0)
                        nosend = false;
                if (nosend)
                    InitTransmission();
            }


            if (!m_TransmissionInProgress)
                throw new ApplicationException("communication must be started.");
            foreach (int p in m_SendBuffers.Keys) {
                if (m_TransmittCalled[p] != true)
                    throw new ApplicationException("before >GetNext< can be called, >Transmit< must be called for every set communication path.");
            }


            // ---------------
            // processing loop
            // ---------------

            while (true) {

                int index;
                MPI_Status status;
                csMPI.Raw.Waitany(m_Requests.Length, m_Requests, out index, out status);
                ReceiveCount++;
                Debug.Assert(index < 0 || m_Requests[index] == csMPI.Raw.MiscConstants.MPI_REQUEST_NULL);
                

                if (index >= 0 && index < size) {
                    
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // buffer content send completed - release gc handle
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    m_SendBuffersPin[index].Free();
                    m_SendBuffersPin.Remove(index);

                
                    
                } else if (index >= size  && index < size * 2) {
                    // +++++++++++++++++++++++++++++++++++++++++
                    // object received - de-serialize and return
                    // +++++++++++++++++++++++++++++++++++++++++

                    int proc = index - size;

                    // free pin
                    m_ReceiveBuffersPin[proc].Free();
                    m_ReceiveBuffersPin.Remove(proc);

                    //if (DiagnosisFile != null) {
                    //    File.WriteAllBytes(DiagnosisFile + "-smsReceive-" + proc + "-" + m_Rank + ".bin", m_ReceiveBuffers[proc]);
                    //}

                    // deserialize
                    o = m_ReceiveBuffers[proc];
                    TargetProc = proc;
                    return true;
                    
                } else if (index == csMPI.Raw.MiscConstants.UNDEFINED) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // we are finished with all - no more objects to receive
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                    m_TransmissionInProgress = false;
                    m_TransmittCalled.SetAll(false);

                    Debug.Assert(m_SendBuffersPin.Count == 0);
                    Debug.Assert(m_ReceiveBuffersPin.Count == 0);


                    o = default(T[]);
                    TargetProc = Int32.MinValue;
                    return false;
                } else {
                   
                    throw new ApplicationException("internal error: index of request out of range.");
                }
            }

            
        }

        /// <summary>
        /// terminates all ongoing communication
        /// </summary>
        public void Dispose() {
            if (m_TransmissionInProgress != false) {
                //throw new ApplicationException("unable to dispose now: communication is not finished yet - application in undefined state.");
                for(int i = 0; i < m_Requests.Length; i++) {
                    if (m_Requests[i] != csMPI.Raw.MiscConstants.MPI_REQUEST_NULL)
                        csMPI.Raw.Cancel(ref m_Requests[i]);
                }
            }

            foreach (var h in m_SendBuffersPin) {
                if (h.Value.IsAllocated)
                    h.Value.Free();
            }
            foreach (var h in m_ReceiveBuffersPin) {
                if (h.Value.IsAllocated)
                    h.Value.Free();
            }
            m_SendBuffersPin = null;
            m_ReceiveBuffersPin = null;
            
            /*
            if (m_SendObjectSizesPin.IsAllocated)
                m_SendObjectSizesPin.Free();

            if (m_ReceiveObjectSizesPin.IsAllocated)
                m_ReceiveObjectSizesPin.Free();
            */

            m_ReceiveBuffers = null;
            m_SendBuffers = null;
        }


        /// <summary>
        /// Exchanges serialize-able data objects in between all MPI processes.
        /// </summary>
        /// <param name="objects_to_send">
        /// Data to send.
        /// - keys: MPI rank of the process to which <em>o</em> should be send to<br/>
        /// - values: some object <em>o</em>
        /// </param>
        /// <param name="comm"></param>
        /// <returns>
        /// Received data.
        /// - keys: MPI rank of the process from which <em>q</em> has been received.<br/>
        /// - values: some object <em>q</em>
        /// </returns>

        public static IDictionary<int, T[]> ExchangeData(IDictionary<int, T[]> objects_to_send, MPI_Comm comm) {
            using (var ams = new ArrayMessenger<T>(comm)) {
                

                foreach(var kv in objects_to_send) {
                    ams.SetCommPath(kv.Key, kv.Value.Length);
                }
                ams.CommitCommPaths();

                foreach (var kv in objects_to_send) {
                    ams.Transmit(kv.Key, kv.Value);
                }

                var R = new Dictionary<int, T[]>();
                T[] obj;
                int rcv_rank;
                while (ams.GetNext(out rcv_rank, out obj)) {
                    R.Add(rcv_rank, obj);
                }

               

                return R;
            }
        }

        /// <summary>
        /// Exchanges serialize-able data objects in between all MPI processes.
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="get_SendData">
        /// returns the data to send for a specific rank
        /// </param>
        /// <param name="receiver_ranks">mpi ranks to send data to</param>
        /// <returns>
        /// Received data.
        /// - keys: MPI rank of the process from which <em>q</em> has been received.<br/>
        /// - values: some object <em>q</em>
        /// </returns>
        public static IDictionary<int, T[]> ExchangeData(IEnumerable<int> receiver_ranks, Func<int, T[]> get_SendData, MPI_Comm comm) {
            var dict = new Dictionary<int, T[]>();
            foreach (var rnk in receiver_ranks)
                dict.Add(rnk, get_SendData(rnk));

            return ExchangeData(dict, comm);
        }


        /// <summary>
        /// equal to <see cref="ExchangeData{T}(IDictionary{int,T},MPI_Comm)"/>, acting on the WORLD-communicator
        /// </summary>
        public static IDictionary<int, T[]> ExchangeData(IDictionary<int, T[]> objects_to_send) {
            return ExchangeData(objects_to_send, csMPI.Raw._COMM.WORLD);
        }

    }
}
