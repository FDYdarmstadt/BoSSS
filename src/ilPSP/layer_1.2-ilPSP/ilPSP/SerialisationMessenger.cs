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

namespace ilPSP.Utils {
    
    
    /// <summary>
    /// This class provides an easy way of exchanging arbitrary (serializeable)
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
    public class SerialisationMessenger : IDisposable {

        static int TagCnt = 1234;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_MPI_comm">the MPI communicator for this messenger (see <see cref="MPI_comm"/>);</param>
        public SerialisationMessenger(MPI_Comm _MPI_comm) {
            MPICollectiveWatchDog.Watch(_MPI_comm);


            m_MPI_comm = _MPI_comm;
            csMPI.Raw.Comm_Rank(_MPI_comm, out m_Rank);
            csMPI.Raw.Comm_Size(_MPI_comm, out m_Size);

            m_MyCommPaths = new byte[m_Size];
            m_AllCommPaths = new byte[m_Size, m_Size];
            m_ReceiveObjectSizes = new int[m_Size];
            m_SendObjectSizes = new int[m_Size];
            m_Requests = new MPI_Request[m_Size * 4];
            m_TransmittCalled = new BitArray(m_Size, false);

            unsafe {
                int myTag = TagCnt;
                TagCnt += 13;

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
        /// rank of the actual process within  the MPI communicator <see cref="MPI_comm"/>;
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
        public void SetCommPath(int TargetProcRank) {
            if (m_CommPathsCommited)
                throw new ApplicationException("setup phase is already finished.");
            if (TargetProcRank == m_Rank)
                throw new ArgumentException("sending data to process itself is not possible.");
            m_MyCommPaths[TargetProcRank] = Byte.MaxValue;
        }

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

        /// <summary>
        /// an entry unequal to 0 at index p indicates that objects should be transmitted from 
        /// this process to process p.
        /// </summary>
        byte[] m_MyCommPaths;


        /// <summary>
        /// A value unequal to 0 at entry [i,j] indicates that processor j (target process) will receive
        /// an object from processor i (source process);
        /// 1st index: source process
        /// 2nd index: target process
        /// </summary>
        private byte[,] m_AllCommPaths;

        
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

            csMPI.Raw.Allgather(pSndBuf, m_MyCommPaths.Length, csMPI.Raw._DATATYPE.BYTE,
                             pRcvBuf, m_MyCommPaths.Length, csMPI.Raw._DATATYPE.BYTE,
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
        /// keys: processor rank "p";
        /// values: memory stream for the object that has to be send to process "p";
        /// </summary>
        SortedDictionary<int, byte[]> m_SendBuffers = new SortedDictionary<int, byte[]>();

        /// <summary>
        /// keys: processor rank "p";
        /// values: garbage collector pin handle for the buffer of memory stream at key "p" in <see cref="m_SendBuffers"/>;
        /// </summary>
        SortedDictionary<int, GCHandle> m_SendBuffersPin = new SortedDictionary<int, GCHandle>();


        /// <summary>
        /// keys: processor rank "p";
        /// values: memory block for the object that is received from process "p";
        /// </summary>
        SortedDictionary<int, byte[]> m_ReceiveBuffers = new SortedDictionary<int,byte[]>();
        
        /// <summary>
        /// keys: processor rank "p";
        /// values: garbage collector pin handle for the buffer memory at key "p" in <see cref="m_ReceiveBuffers"/>;
        /// </summary>
        SortedDictionary<int, GCHandle> m_ReceiveBuffersPin = new SortedDictionary<int, GCHandle>();

        
        /// <summary>
        /// true indicates that the setup is finished
        /// </summary>
        bool m_CommPathsCommited;

        /// <summary>
        /// formatter used for all serialization/de-serialization (legacy)
        /// </summary>
        BinaryFormatter m_Formatter = new BinaryFormatter();

        /// <summary>
        /// formatter used for all serialization/de-serialization
        /// </summary>
        JsonSerializer jsonFormatter = new JsonSerializer() {
            NullValueHandling = NullValueHandling.Include,
            TypeNameHandling = TypeNameHandling.All,
            ConstructorHandling = ConstructorHandling.AllowNonPublicDefaultConstructor,
            ReferenceLoopHandling = ReferenceLoopHandling.Serialize,
            //TypeNameAssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Full
            TypeNameAssemblyFormatHandling = TypeNameAssemblyFormatHandling.Full
        };
        
        bool m_UseJson = false;

        /// <summary>
        /// Whether to use a JSON formatter or the binary formatter
        /// </summary>
        public bool UseJson {
            get {
                return m_UseJson;
            }
            set {
                if (m_CommPathsCommited || m_TransmissionInProgress)
                    throw new NotSupportedException("Cannot be changed after setup phase.");
                m_UseJson = value;
            }

        }

        [Serializable]
        class JsonContainer {
            public object PayLoad;
        }

        public static bool writeSize = false;

        byte[] SerializeObject(object o) {
            using (var ms = new MemoryStream()) {
                Stopwatch stw = new Stopwatch();
                stw.Start();
                if (m_UseJson) {
                    // see: https://stackoverflow.com/questions/25741895/error-serialising-simple-string-to-bson-using-newtonsoft-json-net
                    //var containerArray = Array.CreateInstance(o.GetType(), 1);
                    //containerArray.SetValue(o, 0);
                    var containerObj = new JsonContainer() { PayLoad = o };

                    using (var w = new BsonDataWriter(ms)) {
                        jsonFormatter.Serialize(w, containerObj);
                    }
                } else {
                    m_Formatter.Serialize(ms, o);
                }
                stw.Stop(); 

                byte[] ret = ms.GetBuffer();
                //if (ret.Length > ms.Position) {
                //    Array.Resize(ref ret, (int)ms.Position);
                //}

                if (writeSize) {
                    Console.Error.WriteLine($"r{this.m_Rank}: serialized to {ret.Length/(1024.0*1024.0)} MEGAbyte! took {stw.Elapsed.TotalMilliseconds} msec(UesJSON = {m_UseJson})");
                }

                return ret;
            }

        }

        object DeserializeObject(byte[] data, Type t) {
            using (var ms = new MemoryStream(data)) {
                if (m_UseJson) {
                    //Type ArrayType = Array.CreateInstance(t, 0).GetType();

                    using (var w = new BsonDataReader(ms)) {
                        var containerObj = (JsonContainer) jsonFormatter.Deserialize(w, typeof(JsonContainer));
                        return containerObj.PayLoad;
                    }
                } else {
                    return m_Formatter.Deserialize(ms);
                }
            }
        }



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
        /// 1st phase of the transmission process; must be called after 
        /// <see cref="CommitCommPaths"/>.
        /// </summary>
        /// <param name="TargetProc"></param>
        /// <param name="graph">
        /// a graph of objects, which will be send to process <paramref name="TargetProc"/>
        /// by using serialization.
        /// </param>
        public void Transmit(int TargetProc, object graph) {

            // ------------------------------
            // check for correctness of usage
            // ------------------------------

            if (graph == null)
                throw new ArgumentNullException("cannot sent 'null' - objects.");
            if (!m_CommPathsCommited) 
                throw new ApplicationException("communication paths have to be committed first.");
            if (!m_SendBuffers.ContainsKey(TargetProc))
                throw new ArgumentException("no communication path was set for specified processor.");
            if( m_TransmittCalled[TargetProc] == true)
                throw new ApplicationException("Transmit was already called for specified target process.");
            m_TransmittCalled[TargetProc] = true;


            // -------
            // startup
            // -------

            if (!m_TransmissionInProgress)
                InitTransmission();

            // ---------------
            // serialize object
            // ---------------

            byte[] Buffer = SerializeObject(graph);
            m_SendBuffers[TargetProc] = Buffer;
#if DEBUG
            DeserializeObject(Buffer, graph.GetType());
#endif

            // ------------
            // send objects
            // ------------

            // send size of object 
            m_SendObjectSizes[TargetProc] = m_SendBuffers[TargetProc].Length;
            csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(m_SendObjectSizes, TargetProc),
                          4, csMPI.Raw._DATATYPE.BYTE,
                          TargetProc, TagBufferSize + m_MyTagOffset,
                          m_MPI_comm,
                          out m_Requests[TargetProc]);

            // send buffer content
            if (m_SendBuffersPin.ContainsKey(TargetProc))
                m_SendBuffersPin[TargetProc] = GCHandle.Alloc(Buffer, GCHandleType.Pinned);
            else
                m_SendBuffersPin.Add(TargetProc, GCHandle.Alloc(Buffer, GCHandleType.Pinned));
            if (Buffer.Length > 0) {
                csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(Buffer, 0),
                              Buffer.Length, csMPI.Raw._DATATYPE.BYTE,
                              TargetProc, TagBufferContent + m_MyTagOffset,
                              m_MPI_comm,
                              out m_Requests[m_Size + TargetProc]);
            }

            //if (DiagnosisFile != null) {
            //    File.WriteAllBytes(DiagnosisFile + "-smsSend-" + m_Rank + "-" + TargetProc + ".bin", Buffer);
            //}
        }

        /*
        /// <summary>
        /// Dumping of send/receive blocks for debugging purpose
        /// </summary>
        static public string DiagnosisFile = null;
        */

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

        /// <summary>
        /// MPI tag indicating that the content of a buffer is transmitted
        /// </summary>
        const int TagBufferContent = 234;

        /// <summary>
        /// MPI request handles for nonblocking send and receives.
        /// Array of length of 4*<see cref="m_Size"/>:
        /// <list type="bullet">
        ///   <item>
        ///   1st quarter (index 0 to <see cref="m_Size"/>-1): send handles for buffer size transmission;
        ///   </item>
        ///   <item>
        ///   2nd quarter (index <see cref="m_Size"/> to 2*<see cref="m_Size"/>-1): send handles for buffer content;
        ///   </item>
        ///   <item>
        ///   3rd quarter (index 2*<see cref="m_Size"/> to 3*<see cref="m_Size"/>-1): receive handles for buffer size;
        ///   </item>
        ///   <item>
        ///   4th quarter (index 3*<see cref="m_Size"/> to 4*<see cref="m_Size"/>-1): receive handles for buffer content;
        ///   </item>
        /// </list>
        /// </summary>
        MPI_Request[] m_Requests;



        /// <summary>
        /// 
        /// </summary>
        private void InitTransmission() {

            if (m_TransmissionInProgress)
                throw new ApplicationException("internal error: transmission already in progress.");
            m_TransmissionInProgress = true;

            m_SendObjectSizesPin = GCHandle.Alloc(m_SendObjectSizes, GCHandleType.Pinned);

            // -------------------
            // clear request array
            // -------------------

            for (int i = 0; i < m_Requests.Length; i++) 
                m_Requests[i] = csMPI.Raw.MiscConstants.MPI_REQUEST_NULL;

            // --------------------------------
            // initiate receive of buffer sizes
            // --------------------------------

            if( m_ReceiveObjectSizesPin.IsAllocated)
                throw new ApplicationException("internal error: gc handle already exists.");

            m_ReceiveObjectSizesPin = GCHandle.Alloc(m_ReceiveObjectSizes,GCHandleType.Pinned);

            for( int p = 0; p < m_Size; p++) {
                if (m_AllCommPaths[p, m_Rank] != 0) {
                    csMPI.Raw.Irecv(Marshal.UnsafeAddrOfPinnedArrayElement(m_ReceiveObjectSizes,p),
                                 sizeof(int), csMPI.Raw._DATATYPE.BYTE,
                                 p, TagBufferSize + m_MyTagOffset,
                                 m_MPI_comm,
                                 out m_Requests[2*m_Size + p]);
                }
            }


            
        }

        int ReceiveCount = 0;

        /// <summary>
        /// 2nd and last phase of the transmission process; must be called after all 
        /// calls to <see cref="Transmit"/> have been done;
        /// Must be called multiple times until <paramref name="o"/> is null
        /// after return, which guarantees that all send/receive's are finished
        /// (otherwise is very likly that MPI is in an undefined state).
        /// After this, the object is ready for another send/receive cycle.
        /// </summary>
        /// <param name="TargetProc">
        /// process rank which has send <paramref name="o"/> to this process.
        /// </param>
        /// <param name="o">
        /// on exit, the received object or null, if all objects were received;
        /// and the communication is finished.
        /// </param>
        /// <typeparam name="T">
        /// type of received object
        /// </typeparam>
        /// <returns>
        /// true, if <paramref name="o"/> contains a received object;<br/>
        /// false if the communication is finished.
        /// </returns>
        public bool GetNext<T>(out int TargetProc, out T o) {

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
                    // ++++++++++++++++++++++++++++++++++++++++++++
                    // some send has been completed - nothing to do
                    // ++++++++++++++++++++++++++++++++++++++++++++

                } else if (index >= size && index < size*2) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++
                    // buffer content send completed - release gc handle
                    // +++++++++++++++++++++++++++++++++++++++++++++++++

                    m_SendBuffersPin[index - size].Free();
                    m_SendBuffersPin.Remove(index - size);

                } else if (index >= size * 2 && index < size * 3) {
                    // ++++++++++++++++++++++++++++++++++++++++++++
                    // object size received - start content receive
                    // ++++++++++++++++++++++++++++++++++++++++++++

                    int proc = index - 2*size;
                    
                    // reallocate memory, if necessary
                    byte[] buffer = m_ReceiveBuffers[proc];
                    if (buffer == null) {
                        buffer = new byte[m_ReceiveObjectSizes[proc]];
                        m_ReceiveBuffers[proc] = buffer;
                    } else {
                        if (buffer.Length < m_ReceiveObjectSizes[proc]) {
                            buffer = new byte[m_ReceiveObjectSizes[proc]];
                            m_ReceiveBuffers[proc] = buffer; 
                        }
                    }

                    // pin buffer
                    GCHandle pin = GCHandle.Alloc(buffer, GCHandleType.Pinned);
                    if (m_ReceiveBuffersPin.ContainsKey(proc))
                        m_ReceiveBuffersPin[proc] = pin;
                    else
                        m_ReceiveBuffersPin.Add(proc, pin);

                    // start nonblocking receive
                    if (m_ReceiveObjectSizes[proc] > 0) {
                        csMPI.Raw.Irecv(pin.AddrOfPinnedObject(),
                                     m_ReceiveObjectSizes[proc], csMPI.Raw._DATATYPE.BYTE,
                                     proc, TagBufferContent + m_MyTagOffset,
                                     m_MPI_comm,
                                     out m_Requests[size * 3 + proc]);
                    }
                } else if (index >= size * 3 && index < size * 4) {
                    // +++++++++++++++++++++++++++++++++++++++++
                    // object received - de-serialize and return
                    // +++++++++++++++++++++++++++++++++++++++++

                    int proc = index - size * 3;

                    // free pin
                    m_ReceiveBuffersPin[proc].Free();
                    m_ReceiveBuffersPin.Remove(proc);

                    //if (DiagnosisFile != null) {
                    //    File.WriteAllBytes(DiagnosisFile + "-smsReceive-" + proc + "-" + m_Rank + ".bin", m_ReceiveBuffers[proc]);
                    //}

                    // deserialize
                    o = (T)DeserializeObject(m_ReceiveBuffers[proc], typeof(T));
                    TargetProc = proc;
                    return true;
                    
                } else if (index == csMPI.Raw.MiscConstants.UNDEFINED) {
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++
                    // we are finished with all - no more objects to receive
                    // +++++++++++++++++++++++++++++++++++++++++++++++++++++

                    m_TransmissionInProgress = false;
                    m_TransmittCalled.SetAll(false);

                    m_SendObjectSizesPin.Free();
                    m_ReceiveObjectSizesPin.Free();
                    Debug.Assert(m_SendBuffersPin.Count == 0);
                    Debug.Assert(m_ReceiveBuffersPin.Count == 0);


                    o = default(T);
                    TargetProc = Int32.MinValue;
                    return false;
                } else {
                    //if(PoorManDebugger != null)
                    //{
                    //    PoorManDebugger.WriteLine("are you kidding me?");
                    //    PoorManDebugger.Flush();
                    //}
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
            
            if (m_SendObjectSizesPin.IsAllocated)
                m_SendObjectSizesPin.Free();

            if (m_ReceiveObjectSizesPin.IsAllocated)
                m_ReceiveObjectSizesPin.Free();

            m_ReceiveBuffers = null;
            m_SendBuffers = null;
        }


        /// <summary>
        /// Exchanges serialize-able data objects in between all MPI processes.
        /// </summary>
        /// <param name="objects_to_send">
        /// Data to send.<br/>
        /// keys: MPI rank of the process to which <em>o</em> should be send to<br/>
        /// values: some object <em>o</em>
        /// </param>
        /// <param name="comm"></param>
        /// <returns>
        /// Received data.<br/>
        /// keys: MPI rank of the process from which <em>q</em> has been received.<br/>
        /// values: some object <em>q</em>
        /// </returns>
        /// <param name="__UseJSON"></param>
        public static IDictionary<int, T> ExchangeData<T>(IDictionary<int, T> objects_to_send, MPI_Comm comm, bool __UseJSON = false) {
            using (var sms = new SerialisationMessenger(comm)) {
                sms.UseJson = __UseJSON;
                //if (PoorManDebugger != null) {
                //    PoorManDebugger.WriteLine("tag offset is " + sms.m_MyTagOffset);
                //    PoorManDebugger.Flush();
                //}
                sms.SetCommPathsAndCommit(objects_to_send.Keys);

                foreach (var kv in objects_to_send) {
                    sms.Transmit(kv.Key, kv.Value);
                }

                var R = new Dictionary<int, T>();
                T obj;
                int rcv_rank;
                while (sms.GetNext(out rcv_rank, out obj)) {
                    R.Add(rcv_rank, obj);
                }

                //if(PoorManDebugger != null)
                //{
                //    PoorManDebugger.WriteLine(" leaving ExchangeData");
                //    PoorManDebugger.Flush();
                //}

                return R;
            }
        }

        /*
        public static StreamWriter PoorManDebugger; 
        */

        /// <summary>
        /// equal to <see cref="ExchangeData{T}(IDictionary{int,T},MPI_Comm)"/>, acting on the WORLD-communicator
        /// </summary>
        public static IDictionary<int, T> ExchangeData<T>(IDictionary<int, T> objects_to_send) {
            return ExchangeData<T>(objects_to_send, csMPI.Raw._COMM.WORLD);
        }

    }
}
