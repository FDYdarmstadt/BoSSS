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
using System.Reflection;
using System.Runtime.InteropServices;
using MPI.Wrappers;

namespace BoSSS.Foundation.Comm {
    
    
    /// <summary>
    /// Performs the exchange of general messages (i.e. an array of <typeparamref name="ItemType"/>'s)
    /// from each process to various other processes;
    /// </summary>
    /// <typeparam name="ItemType">must be a "full value-type" i.e. a primitive type or 
    /// a struct which contains no references at all.
    /// Other types produce some exception during construction;
    /// </typeparam>
    public class Many2ManyMessenger<ItemType> : IDisposable
        where ItemType : struct {
        
        static void Memcpy(IntPtr dst, IntPtr src, int NoOfBytes) {
            //Debugger.Break();

            unsafe {
                byte* pDstRst = (byte*)dst;
                byte* pSrcRst = (byte*)src;
                int i = 0;
                for (; i < NoOfBytes; i += sizeof(byte)) {
                    *pDstRst = *pSrcRst;
                    pDstRst++;
                    pSrcRst++;
                }
                
            }
        }


        /// <summary>
        /// calls <see cref="Dispose"/>
        /// </summary>
        ~Many2ManyMessenger() {
            Dispose();
        }

        MPI.Wrappers.MPI_Comm m_Comm;

        int size;
        int MyRank;

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="comm"></param>
        //// <typeparam name="t">the primitive type of the message; must be a value-type;</typeparam>
        public Many2ManyMessenger(MPI.Wrappers.MPI_Comm comm) {
            m_Comm = comm;
            MPI.Wrappers.csMPI.Raw.Comm_Size(m_Comm, out size);
            MPI.Wrappers.csMPI.Raw.Comm_Rank(m_Comm, out MyRank);
            
            m_MyCommPaths = new int[size];
            m_AllCommPaths = new int[size, size];

            m_ReceiveBuffers = new Buffer[size];
            m_SendBuffers = new Buffer[size];

            m_Requests = new MPI_Request[size * 2];
            m_ArrayOfStatuses = new MPI_Status[m_Requests.Length];
            

            // check item type
            m_ItemType = typeof(ItemType);
            CheckItemTypeRecursive(m_ItemType);
            m_ItemSize = Marshal.SizeOf(m_ItemType);
        }

        /// <summary>
        /// throws an exception if any member of <paramref name="t"/>
        /// is not suited for transport over the network (e.g. pointers);
        /// </summary>
        private void CheckItemTypeRecursive(Type t) {
            if (t.IsPrimitive)
                return;
            if (!t.IsValueType)
                throw new NotSupportedException("Type T must be composed only of primitive types.");

            FieldInfo[] members = t.GetFields(BindingFlags.NonPublic|BindingFlags.Public|BindingFlags.Instance);

            foreach (FieldInfo member in members) {
                CheckItemTypeRecursive(member.FieldType);
            }
        }



        /// <summary>
        /// size of one data item in bytes
        /// </summary>
        private int m_ItemSize;

        /// <summary>
        /// type of item
        /// </summary>
        Type m_ItemType;


        /// <summary>
        /// setup procedure;
        /// Informs the messenger that <paramref name="NoOfItems"/> items should be transmitted
        ///  1st function (multiple calls) in the setup process; If all communication paths
        /// are set, <see cref="CommitCommPaths"/> is the next method to call;
        /// </summary>
        /// <param name="TargetProcRank"></param>
        /// <param name="NoOfItems"></param>
        public void SetCommPath(int TargetProcRank, int NoOfItems) {
            if (m_CommPathsCommited)
                throw new ApplicationException("setup phase is already finished.");
            if (TargetProcRank == MyRank)
                throw new ArgumentException("sending data to process itself is not possible.");
            m_MyCommPaths[TargetProcRank] = NoOfItems;
        }


        /// <summary>
        /// A value x at entry m indicates that processor m will receive x items
        /// from this processor;
        /// </summary>
        private int[] m_MyCommPaths;

        /// <summary>
        /// A value x at entry [i,j] indicates that processor j (target process) will receive
        /// x items from processor i (source process);
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
            if(m_CommPathsCommited)
                throw new ApplicationException("setup phase is already finished.");

            GCHandle myCommPathsPin = GCHandle.Alloc(GCHandleType.Pinned);
            GCHandle allCommPathsPin = GCHandle.Alloc(GCHandleType.Pinned);

            IntPtr pSndBuf = Marshal.UnsafeAddrOfPinnedArrayElement(m_MyCommPaths, 0);
            IntPtr pRcvBuf = Marshal.UnsafeAddrOfPinnedArrayElement(m_AllCommPaths, 0);

            csMPI.Raw.Allgather(pSndBuf, m_MyCommPaths.Length * 4, csMPI.Raw._DATATYPE.BYTE,
                             pRcvBuf, m_MyCommPaths.Length * 4, csMPI.Raw._DATATYPE.BYTE,
                             m_Comm);

            myCommPathsPin.Free();
            allCommPathsPin.Free();
            m_CommPathsCommited = true;

            //StringWriter stw = new StringWriter();
            //stw.Write("exch P" + m_Master.MyRank + ": ");
            //for (int j = 0; j < m_Master.Size; j++) {
            //    for (int i = 0; i < m_Master.Size; i++) {
            //        stw.Write(m_AllCommPaths[j, i]);
            //        stw.Write(", ");
            //    }
            //    stw.Write("; ");
            //}
            //Console.WriteLine(stw.ToString());

            // create send buffers...
            for (int p = 0; p < size; p++) {
                if (m_MyCommPaths[p] > 0)
                    m_SendBuffers[p] = new Buffer(this, false, p);
            }

            // create receive buffers...
            for (int p = 0; p < size; p++) {
                if (m_AllCommPaths[p,MyRank] > 0)
                    m_ReceiveBuffers[p] = new Buffer(this, true, p);
            }

        }


        /// <summary>
        /// true indicates that the setup is finished
        /// </summary>
        bool m_CommPathsCommited;


        ///// <summary>
        ///// pointer to unmanaged buffers for sending data items;
        ///// this buffers, if allocated, are released by <see cref="Dispose"/>;
        ///// </summary>
        //IntPtr[] m_SendBuffers;

        ///// <summary>
        ///// size in bytes for the buffers in <see cref="m_SendBufferSize"/>
        ///// </summary>
        //int[] m_SendBufferSize;
        
        ///// <summary>
        ///// pointer to unmanaged buffers for receiving data items;
        ///// this buffers, if allocated, are released by <see cref="Dispose"/>;
        ///// </summary>
        //IntPtr[] m_ReceiveBuffers;

        ///// <summary>
        ///// size in bytes for the buffers in <see cref="m_ReceiveBuffers"/>
        ///// </summary>
        //int[] m_ReceiveBufferSize;


        /// <summary>
        /// 
        /// </summary>
        Buffer[] m_SendBuffers;


        /// <summary>
        /// returns the send buffer to the <paramref name="TargetProcessor"/>'th processor;
        /// if no communication path has been set (<see cref="SetCommPath"/>) for the specified processor,
        /// the buffer is null;
        /// The buffers are not available before <see cref="CommitCommPaths"/> has been called;
        /// </summary>
        /// <param name="TargetProcessor"></param>
        /// <returns></returns>
        public Buffer SendBuffers(int TargetProcessor) {
            if (!m_CommPathsCommited)
                throw new ApplicationException("communication paths have to be commited before buffers can be accessed.");
            return m_SendBuffers[TargetProcessor];
        }

        /// <summary>
        /// 
        /// </summary>
        Buffer[] m_ReceiveBuffers;

        /// <summary>
        /// returns the receive buffer for data from the <paramref name="SourceProcessor"/>'th MPI-process;
        /// if no communication path has been set (ANOTHER process has to do that by calling <see cref="SetCommPath"/>) 
        /// for the specified processor to this processor,
        /// the buffer is null;
        /// The buffers are not available before <see cref="CommitCommPaths"/> has been called;
        /// </summary>
        public Buffer ReceiveBuffers(int SourceProcessor) {
            if (!m_CommPathsCommited)
                throw new ApplicationException("communication paths have to be commited before buffers can be accessed.");
            return m_ReceiveBuffers[SourceProcessor];
        }

        /// <summary>
        /// MPI - request - handles for the unblocking receive and send operations;
        /// the first half of the array contains the requests for receive operations,
        /// the second half requests for send operations.
        /// </summary>
        MPI_Request[] m_Requests;

        /// <summary>
        /// MPI status array
        /// </summary>
        MPI_Status[] m_ArrayOfStatuses;

        /// <summary>
        /// transmission procedure, 1st phase;
        /// makes everything ready to transmit data, which is done
        /// by calling to <see cref="TransmittData"/>
        /// </summary>
        /// <param name="NoOfItemsMultiplyer">
        /// The number of items for each communication path which were specifyed by 
        /// calling <see cref="SetCommPath"/> is multiplyed with this number to get
        /// the actual number of items to sent/receive.
        /// this number must be the same in all processes, otherwise the behavior is undefined
        /// - the implementation doesn't checks for that;
        /// </param>
        public void StartTransmission(int NoOfItemsMultiplyer) {

            if (!m_CommPathsCommited)
                throw new ApplicationException("setup phase must be finished first.");

            m_NoOfItemsMultiplyer = NoOfItemsMultiplyer;
            
            int NoOfProcs = size;
            
            // (re)allocate receive buffers, lock them and start receive 
            // ---------------------------------------------------------
            for (int i = 0; i < NoOfProcs; i++) {
                m_Requests[i] = csMPI.Raw.MiscConstants.MPI_REQUEST_NULL;
                if (m_ReceiveBuffers[i] != null) {


                    m_ReceiveBuffers[i].CheckAllocation();
                    m_ReceiveBuffers[i].SetUnLocked(false);

                    csMPI.Raw.Irecv(m_ReceiveBuffers[i].BufferMemory, m_ReceiveBuffers[i].m_MySize*m_ItemSize,
                                 csMPI.Raw._DATATYPE.BYTE, i, 0, m_Comm, out m_Requests[i]);
                }
            }


            // (re)allocate send buffers, unlock them
            // --------------------------------------
            for (int i = 0; i < NoOfProcs; i++) {
                m_Requests[i + size] = csMPI.Raw.MiscConstants.MPI_REQUEST_NULL;
                if (m_SendBuffers[i] != null) {
                    m_SendBuffers[i].CheckAllocation();
                    m_SendBuffers[i].SetUnLocked(true);
                }

            }


            // finish
            // ------

            m_TransmittDataReady = true;


        }

        /// <summary>
        /// The number of items for each communication path which were specified by 
        /// calling <see cref="SetCommPath"/> is multiplied with this number to get
        /// the actual number of items to sent/receive.
        /// this number must be the same in all processes, otherwise the behavior is undefined
        /// - the implementation doesn't checks for that;
        /// </summary>
        int m_NoOfItemsMultiplyer;

        /// <summary>
        /// true indicates that <see cref="StartTransmission"/> was called and
        /// calls to <see cref="TransmittData"/> are legal;
        /// </summary>
        bool m_TransmittDataReady;


        /// <summary>
        /// 
        /// </summary>
        public sealed class Buffer : IList<ItemType> {

            /// <summary>
            /// constructor;
            /// Cannot be called before communication paths of <paramref name="owner"/> are committed
            /// (<see cref="Many2ManyMessenger{ItemType}.CommitCommPaths"/>);
            /// </summary>
            /// <param name="owner"></param>
            /// <param name="ReceiveOrTransmit">
            /// true: a receive buffer
            /// false: a transmit (or send) buffer
            /// </param>
            /// <param name="ProcRank"></param>
            internal Buffer(Many2ManyMessenger<ItemType> owner, bool ReceiveOrTransmit, int ProcRank) {
                if (owner.m_ItemType != typeof(ItemType))
                    throw new ArgumentException("Item type varies from messenger item type");
                m_Owner = owner;
                m_ProcRank = ProcRank;
                m_ReceiveOrTransmit = ReceiveOrTransmit;

                

                if( !m_Owner.m_CommPathsCommited)
                    throw new ApplicationException("communication paths must be commited first.");
            }

            /// <summary>
            /// owner
            /// </summary>
            Many2ManyMessenger<ItemType> m_Owner;


            /// <summary>
            /// unmanaged memory
            /// </summary>
            IntPtr m_Buffer = IntPtr.Zero;

            /// <summary>
            /// 
            /// </summary>
            internal IntPtr BufferMemory { get { return m_Buffer; } }

            /// <summary>
            /// size of <see cref="m_Buffer"/> in bytes;
            /// maybe more memory than actually needed;
            /// </summary>
            int m_NumberOfAllocatedBytes = 0;


            /// <summary>
            /// allocates/reallocates memory if needed;
            /// sets <see cref="m_MySize"/> to a correct value value;
            /// </summary>
            internal void CheckAllocation() {
                if (m_disposed)
                    throw new ApplicationException("buffer has been disposed.");

                if (m_ReceiveOrTransmit) {
                    // receive buffer
                    m_MySize = m_Owner.m_AllCommPaths[m_ProcRank, m_Owner.MyRank];
                } else {
                    // send buffer
                    m_MySize = m_Owner.m_MyCommPaths[m_ProcRank];
                }

                if (m_MySize <= 0)
                    throw new ApplicationException("internal error: buffer size smaller or equal to 0");

                m_MySize *= m_Owner.m_NoOfItemsMultiplyer;
                int MinBufferSz = m_MySize*m_Owner.m_ItemSize;


                if (m_NumberOfAllocatedBytes < MinBufferSz) {
                    if (m_Buffer == IntPtr.Zero) {
                        m_Buffer = Marshal.AllocHGlobal((IntPtr)MinBufferSz);
                    } else {
                        m_Buffer = Marshal.ReAllocHGlobal(m_Buffer, (IntPtr)MinBufferSz);
                    }
                    m_NumberOfAllocatedBytes = MinBufferSz;
                }

            }



            /// <summary>
            /// 
            /// </summary>
            int m_ProcRank;

            /// <summary>
            /// true: receive buffer;
            /// false: send buffer
            /// </summary>
            bool m_ReceiveOrTransmit;


            /// <summary>
            /// <see cref="UnLocked"/>
            /// </summary>
            bool m_UnLocked = false;

            /// <summary>
            /// sets <see cref="m_UnLocked"/>;
            /// </summary>
            /// <param name="newState"></param>
            internal void SetUnLocked(bool newState) {
                m_UnLocked = newState;
            }


            /// <summary>
            /// if true, the user is allowed to access data;
            /// if false, the buffer memory is looked for MPI communication;
            /// </summary>
            public bool UnLocked { get { return m_UnLocked; } }



            /// <summary>
            /// checks if the user is allowed to access the buffer and throws an exception if he is not
            /// </summary>
            /// <returns></returns>
            void CheckReady() {
                //if (m_Owner.m_disposed == true)
                //    throw new ApplicationException("owner disposed.");

                //if (m_ReceiveOrTransmit) {

                //    if (MySize() > 0)
                //        if (m_Buffer == IntPtr.Zero)
                //            throw new ApplicationException("no internal buffer is allocated-at least one StartTransmission/Finishblocking pair must be called.");
                //    if (m_Owner.m_TransmittDataReady)
                //        throw new ApplicationException("Transmission in progress-receive buffer content is not valid until FinishBlocking has benn called.");
                //} else {
                //    if (!m_Owner.m_TransmittDataReady)
                //        throw new ApplicationException("StartTransmission must be called first.");
                //    if (m_Owner.m_MyCommPaths[m_ProcRank] < 0)
                //        throw new ApplicationException("communication is allready in progress.");
                //}

                if (!m_UnLocked)
                    throw new ApplicationException("Buffer is looked for MPI communication.");
                if (m_disposed)
                    throw new ApplicationException("Buffer has been disposed");
            }


            /// <summary>
            /// size of buffer (in number of items)
            /// </summary>
            /// <returns></returns>
            internal int m_MySize = Int32.MinValue; 
            
            //() {
            //    int sz = m_Owner.m_NoOfItemsMultiplyer;

            //    if (m_ReceiveOrTransmit)
            //        sz *= m_Owner.m_AllCommPaths[m_ProcRank, m_Owner.m_Master.MyRank];
            //    else
            //        sz *= m_Owner.m_MyCommPaths[m_ProcRank];
            //    return sz;
            //}


            /// <summary>
            /// 
            /// </summary>
            /// <param name="i"></param>
            /// <returns></returns>
            void CheckIndex(int i) {
                if (i < 0) throw new IndexOutOfRangeException();
                if (i >= m_MySize) throw new IndexOutOfRangeException();
            }



            /// <summary>
            /// not supported
            /// </summary>
            /// <param name="item"></param>
            /// <returns></returns>
            public int IndexOf(ItemType item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported
            /// </summary>
            /// <param name="index"></param>
            /// <param name="item"></param>
            public void Insert(int index, ItemType item) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// not supported
            /// </summary>
            /// <param name="index"></param>
            public void RemoveAt(int index) {
                throw new NotSupportedException();
            }

            /// <summary>
            /// helper memory
            /// </summary>
            ItemType[] m_itemMemory = new ItemType[1];

            /// <summary>
            /// set/get one item;
            /// causes some overhead - use <see cref="CopyTo(ItemType[],int)"/> and <see cref="CopyFrom(ItemType[],int)"/> to avoid that;
            /// </summary>
            /// <param name="index">zero-based index</param>
            /// <returns>da item!</returns>
            public ItemType this[int index] {
                get {
                    CheckReady();
                    CheckIndex(index);

                   
                    GCHandle OutputPin = GCHandle.Alloc(m_itemMemory, GCHandleType.Pinned);

                    unsafe {
                        IntPtr pItem = Marshal.UnsafeAddrOfPinnedArrayElement(m_itemMemory, 0);
                        int itmsz = m_Owner.m_ItemSize;
                        byte* pSrc = (byte*)m_Buffer + (index*itmsz);

                        Memcpy(pItem, (IntPtr)pSrc, itmsz);
                    }
                    OutputPin.Free();
                    return m_itemMemory[0];

                }
                set {
                    CheckReady();
                    CheckIndex(index);

                    m_itemMemory[0] = value;
                    GCHandle OutputPin = GCHandle.Alloc(m_itemMemory, GCHandleType.Pinned);

                    unsafe {
                        IntPtr pItem = Marshal.UnsafeAddrOfPinnedArrayElement(m_itemMemory, 0);
                        int itmsz = m_Owner.m_ItemSize;
                        byte* pDst = (byte*)m_Buffer + index*itmsz;

                        Memcpy((IntPtr)pDst, pItem, itmsz);
                    }
                    OutputPin.Free();
                }
            }


            /// <summary>
            /// not supported
            /// </summary>
            /// <param name="item"></param>
            public void Add(ItemType item) {
                throw new NotSupportedException("buffer is fixed size");
            }

            /// <summary>
            /// sets all entries 0;
            /// </summary>
            public void Clear() {
                CheckReady();
                
                unsafe {
                    byte* pB = (byte*) m_Buffer;
                    for (int i = 0; i < m_NumberOfAllocatedBytes; i++) {
                        *pB = (byte)0;
                        pB++;
                    }
                }
            }

            /// <summary>
            /// not supported
            /// </summary>
            /// <param name="item"></param>
            /// <returns></returns>
            public bool Contains(ItemType item) {
                throw new NotSupportedException();
            }




            /// <summary>
            /// copys the content of the unmanged buffer into a managed array
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndex">offset into <paramref name="array"/></param>
            public void CopyTo(ItemType[] array, int arrayIndex) {
                CheckReady();
                int cnt = m_MySize;
                if ((array.Length - arrayIndex) < cnt)
                    throw new ArgumentException("target array to short");

                int BytesTocopy = cnt*m_Owner.m_ItemSize;
                GCHandle outPin = GCHandle.Alloc(array, GCHandleType.Pinned);
                unsafe {
                    IntPtr pDst = Marshal.UnsafeAddrOfPinnedArrayElement(array, arrayIndex);
                    Memcpy(pDst, m_Buffer, BytesTocopy);
                }
                outPin.Free();
            }
            
            
            /// <summary>
            /// copys <paramref name="len"/> items form this buffer, starting at index <paramref name="srcOffset"/>
            /// into <paramref name="array"/>, starting at index <paramref name="arrayIndex"/>
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndex">offset into <paramref name="array"/></param>
            /// <param name="len">number of items to copy</param>
            /// <param name="srcOffset">offset into this array</param>
            internal void CopyTo(ItemType[] array, int arrayIndex, int len, int srcOffset) {
                CheckReady();
                if ((array.Length - arrayIndex) < len)
                    throw new ArgumentException("target array to short");

                int BytesTocopy = len*m_Owner.m_ItemSize;
                GCHandle outPin = GCHandle.Alloc(array, GCHandleType.Pinned);
                unsafe {
                    IntPtr pDst = Marshal.UnsafeAddrOfPinnedArrayElement(array, arrayIndex);
                    Memcpy(pDst, (IntPtr)((byte*)m_Buffer + srcOffset*m_Owner.m_ItemSize), BytesTocopy);
                }
                outPin.Free();
            }
            
            
            /// <summary>
            /// copys from an (managed) array into the internal unmanaged buffer
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndex"></param>
            public void CopyFrom(ItemType[] array, int arrayIndex) {
                CheckReady();
                int cnt = m_MySize;
                //Console.WriteLine("P" + m_Owner.m_Master.MyRank + ", al: " + array.Length + " mysize: " + m_MySize + " rcvtrm: " + m_ReceiveOrTransmit);
                if ((array.Length - arrayIndex) < m_MySize)
                    throw new ArgumentException("target array to short");

                int BytesTocopy = cnt*m_Owner.m_ItemSize;
                GCHandle outPin = GCHandle.Alloc(array, GCHandleType.Pinned);
                unsafe {
                    IntPtr pSrc = Marshal.UnsafeAddrOfPinnedArrayElement(array, arrayIndex);
                    Memcpy(m_Buffer, pSrc, BytesTocopy);
                }
                outPin.Free();
            }

            /// <summary>
            /// copys from an (managed) List into the internal unmanaged buffer
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndex"></param>
            public void CopyFromList(IList<ItemType> array, int arrayIndex) {
                CheckReady();
                if ((array.Count - arrayIndex) < m_MySize)
                    throw new ArgumentException("target array to short");

                int itmSz = m_Owner.m_ItemSize;
                ItemType[] mem = new ItemType[1];

                GCHandle memPin = GCHandle.Alloc(mem, GCHandleType.Pinned);
                unsafe {
                    byte* pMem =  (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(mem, 0);
                    byte* pDst = (byte*)m_Buffer;

                    for (int i = 0; i < m_MySize; i++) {
                        mem[0] = array[i+arrayIndex];

                        // i think that for small items a loop is faster than a
                        // memcpy which is called by P/Invoke
                        byte* pSrc = pMem;
                        for (int j = 0; j < itmSz; j++) {
                            *pDst = *pSrc;
                            pDst++;
                            pSrc++;
                        }
                        //csMPI.Raw.Memcpy(m_Buffer, pMem, itmSz); pDst += itmSz;
                    }

                }
                memPin.Free();
            }
            
            
            /// <summary>
            /// copys from an (managed) array into the internal unmanaged buffer;
            /// For each entry i in <paramref name="arrayIndices"/>, this method
            /// copies the entries i*<paramref name="l"/> to i*<paramref name="l"/>+<paramref name="l"/>-1 of <paramref name="array"/>
            /// into the send buffer in subsequent order, starting at index <paramref name="targetoffset"/>.
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndices">indices of the entries in <paramref name="array"/> that should be copied</param>
            /// <param name="l">multiplyer for the indices in <paramref name="arrayIndices"/></param>
            /// <param name="targetoffset"></param>
            internal void CopyFrom(ItemType[] array, int[] arrayIndices, int l, int targetoffset) {
                CheckReady();
                if (arrayIndices.Length*l > (m_MySize-targetoffset))
                    throw new ArgumentException("array bounds exceeded.");

                int itmChunkSz = m_Owner.m_ItemSize*l;
                ItemType[] mem = new ItemType[l];

                GCHandle memPin = GCHandle.Alloc(mem, GCHandleType.Pinned);
                try {
                    unsafe {
                        IntPtr pMem = Marshal.UnsafeAddrOfPinnedArrayElement(mem, 0);
                        byte* pDst = (byte*)m_Buffer + targetoffset*m_Owner.m_ItemSize;

                        int l0 = arrayIndices.Length;
                        for (int i = 0; i < l0; i++) {

                            Array.Copy(array, arrayIndices[i]*l, mem, 0, l);
                            Memcpy((IntPtr)pDst, pMem, itmChunkSz); pDst += itmChunkSz;
                        }
                    }
                } catch (Exception e) {
                    throw e;
                } finally {
                    memPin.Free();
                }
            }

            /// <summary>
            /// opys the content of the unmanged buffer into a managed list
            /// </summary>
            /// <param name="array"></param>
            /// <param name="arrayIndex"></param>
            public void CopyToList(IList<ItemType> array, int arrayIndex) {
                CheckReady();
                int cnt = m_MySize;
                //Console.WriteLine("P" + m_Owner.m_Master.MyRank + ", al: " + array.Length + " mysize: " + m_MySize + " rcvtrm: " + m_ReceiveOrTransmit);
                if ((array.Count - arrayIndex) < m_MySize)
                    throw new ArgumentException("target array to short");

                int itmSz = m_Owner.m_ItemSize;
                ItemType[] mem = new ItemType[1];

                GCHandle memPin = GCHandle.Alloc(mem, GCHandleType.Pinned);
                unsafe {
                    byte* pMem =  (byte*)Marshal.UnsafeAddrOfPinnedArrayElement(mem, 0);
                    byte* pSrc = (byte*)m_Buffer;

                    for (int i = 0; i < m_MySize; i++) {

                        // i think that for small items a loop is faster than a
                        // memcpy which is called by P/Invoke
                        byte* pDst = pMem;
                        for (int j = 0; j < itmSz; j++) {
                            *pDst = *pSrc;
                            pDst++;
                            pSrc++;
                        }
                        //csMPI.Raw.Memcpy(m_Buffer, pMem, itmSz); pDst += itmSz;

                        array[i+arrayIndex] = mem[0];
                    }

                }
                memPin.Free();
            }


            /// <summary>
            /// Number of items in buffer
            /// </summary>
            public int Count {
                get {
                    CheckReady();
                    return m_MySize;
                }
            }

            /// <summary>
            /// always false
            /// </summary>
            public bool IsReadOnly { get { return false; } }

            /// <summary>
            /// not supported;
            /// </summary>
            /// <param name="item"></param>
            /// <returns></returns>
            public bool Remove(ItemType item) {
                throw new NotSupportedException("buffer is fixed size");
            }

            /// <summary>
            /// not implemented
            /// </summary>
            /// <returns></returns>
            public IEnumerator<ItemType> GetEnumerator() {
                throw new Exception("The method or operation is not implemented.");
            }

            /// <summary>
            /// not implemented
            /// </summary>
            /// <returns></returns>
            System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator() {
                throw new Exception("The method or operation is not implemented.");
            }

            /// <summary>
            /// frees unmanaged memory
            /// </summary>
            internal void Dispose() {
                m_UnLocked = false;
                if (m_Buffer != IntPtr.Zero) {
                    Marshal.FreeHGlobal(m_Buffer);
                }
                m_Buffer = IntPtr.Zero;
                m_disposed = true;
            }

            /// <summary>
            /// true after <see cref="Dispose"/> has been called
            /// </summary>
            bool m_disposed = false;
        }




        /// <summary>
        /// transmission procedure, 2nd phase, multiple calls;
        /// calls must exactly match the comm paths which were set up by
        /// calling <see cref="SetCommPath"/>;
        /// After all data items have ben set, a call to <see cref="FinishBlocking"/>
        /// completes the data transmission;
        /// This method ONLY initiates the transmission procedure, the data to transmit must be set
        /// by manipulating the send buffers;
        /// </summary>
        /// <param name="TargetProcRank">index of target process</param>
        public void TransmittData(int TargetProcRank) {
            if (!m_TransmittDataReady)
                throw new ApplicationException("StartTransmission must be called prior to TransmittData");
            if (m_SendBuffers[TargetProcRank] == null)
                throw new ArgumentException("wrong target process - no communication path has been set for that target");
            int sz = size;
            if (m_Requests[sz + TargetProcRank] != csMPI.Raw.MiscConstants.MPI_REQUEST_NULL)
                throw new ApplicationException("TransmitData has already been called for specified target processor.");


            // start nonblocking send
            int NoOfBytes = m_ItemSize*m_SendBuffers[TargetProcRank].m_MySize;
            csMPI.Raw.Issend(m_SendBuffers[TargetProcRank].BufferMemory, NoOfBytes,
                                        csMPI.Raw._DATATYPE.BYTE, TargetProcRank, 0, m_Comm, out m_Requests[sz + TargetProcRank]);
            m_SendBuffers[TargetProcRank].SetUnLocked(false);

            // use -1 to set a flag for the com path
            m_MyCommPaths[TargetProcRank] *= -1;
        }


        /// <summary>
        /// transmission procedure, 3rd (and last) phase;
        /// receives data and blocks until the transmission is finished;
        /// after this method returns, another transmission can be started with
        /// <see cref="StartTransmission"/>;
        /// after this method retunes, all buffers are unlocked;
        /// </summary>
        public void FinishBlocking() {
            m_TransmittDataReady = false;

            for (int i = 0; i < size; i++) {
                if (m_MyCommPaths[i] > 0)
                    throw new ApplicationException("TransmittData for communication path to processor " + i + " was not started.");
            }

            csMPI.Raw.Waitall(m_Requests.Length, m_Requests, m_ArrayOfStatuses);

            for (int i = 0; i < size; i++) 
                m_MyCommPaths[i] *= -1;

            for (int p = 0; p < size; p++) {
                if (m_SendBuffers[p] != null)
                    m_SendBuffers[p].SetUnLocked(true);
                if (m_ReceiveBuffers[p] != null)
                    m_ReceiveBuffers[p].SetUnLocked(true);
            }
        }


        /// <summary>
        /// frees the internal non-managed buffers
        /// </summary>
        public void Dispose() {
            if (m_SendBuffers != null) {
                foreach (Buffer b in m_SendBuffers) {
                    if (b != null)
                        b.Dispose();
                }
                foreach (Buffer b in m_ReceiveBuffers) {
                    if (b != null)
                        b.Dispose();
                }

                this.m_CommPathsCommited = false;
                this.m_ItemSize = Int32.MinValue;
                this.m_ItemType = null;
                this.m_ReceiveBuffers = null;
                this.m_SendBuffers = null;
            }
        }
    }
}
