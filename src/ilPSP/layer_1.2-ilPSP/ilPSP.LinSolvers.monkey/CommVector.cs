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
using MPI.Wrappers;

namespace ilPSP.LinSolvers.monkey {
    
    /// <summary>
    /// 
    /// </summary>
    /// <param name="procRank">rank of MPI-process from which the data was received</param>
    /// <param name="values">pointer (to a double-array) of received data</param>
    internal delegate void OnExternalBlockReceived(int procRank, IntPtr values);
    
    partial class VectorBase {

        /// <summary>
        /// creates ma communication vector object for the given matrix <paramref name="M"/>, i.e.
        /// an vector that can be multiplied with <paramref name="M"/>
        /// </summary>
        abstract public CommVector CreateCommVector(MatrixBase M);

        /// <summary>
        /// a class which contains memory for the external parts 
        /// (those parts which this processor receives from other processors)
        /// of some 
        /// distributed vector;<br/>
        /// Second, this class contains the code to initialize the communication;
        /// </summary>
        /// <remarks>
        /// To multiply a vector with a sparse matrix (in an MPI-parallel context), it is necessary to
        /// exchange parts of the vector between MPI-processes. 
        /// These parts are determined by the popoulation pattern of a sparse matrix,
        /// see <see cref="MatrixBase._SpmvCommPattern"/>.
        /// Therefore, a <see cref="CommVector"/> always depends on a matrix <see cref="Mtx"/>
        /// </remarks>
        abstract public class CommVector : IDisposable {
            
            internal CommVector(MatrixBase M, VectorBase owner) {
                if (!owner.Part.Equals(M.ColPartition)) {
                    throw new ArgumentException("column partition of matrix must be equal to partion of vector", "M");
                }
                m_Mtx = M;
                m_Owner = owner;
                
                // allocate receive buffers
                int i = 0;
                RecvBuffersLock = new GCHandle[m_Mtx._SpmvCommPattern.NoOfReceivedEntries.Keys.Count];
                foreach (int proc in m_Mtx._SpmvCommPattern.NoOfReceivedEntries.Keys) {
                    double[] RecvBuffer = new double[m_Mtx._SpmvCommPattern.NoOfReceivedEntries[proc]];
                    RecvBuffers.Add(proc, RecvBuffer);
                    RecvBuffersLock[i] = GCHandle.Alloc(RecvBuffer, GCHandleType.Pinned);
                    i++;
                }
            }

            /// <summary>
            /// send buffers (starting addresses - can be either unmanaged memory, or pinned managed memory)
            /// </summary>
            protected SortedList<int, IntPtr> SendBuffers = new SortedList<int, IntPtr>();

            /// <summary>
            /// send buffers length (number of items)
            /// </summary>
            protected IDictionary<int, int> SendBuffersLengths = new Dictionary<int, int>();
            
            /// <summary>
            /// receive buffers
            /// </summary>
            internal SortedList<int, double[]> RecvBuffers = new SortedList<int, double[]>();

            VectorBase m_Owner;

            /// <summary>
            /// the vector which is exchanged via this communication object.
            /// </summary>
            public VectorBase Owner { get { return m_Owner; } }

            /// <summary>
            /// <see cref="Mtx"/>;
            /// </summary>
            MatrixBase m_Mtx;

            /// <summary>
            /// the matrix that defines the MPI communication pattern for
            /// this object.
            /// </summary>
            public MatrixBase Mtx { get { return m_Mtx; } }
            
            int Phase = -1;

            /// <summary>
            /// Phase 1 (of 4) of a complete transmission: filling the send buffers.
            /// </summary>
            /// <remarks>
            /// Calling sequence for complete transmission:
            /// <list type="number">
            ///   <item><b>Phase 1/4:</b>Call <see cref="FillSendBuffer"/></item>
            ///   <item><b>Phase 2/4:</b>Call <see cref="StartTransmissionImReturn"/></item>
            ///   <item><b>Phase 3/4:</b>Call <see cref="FillSendBuffer"/></item>
            ///   <item><b>Phase 4/4:</b>Call <see cref="WaitCommFinish"/></item>
            /// </list>
            /// </remarks>
            virtual public void FillSendBuffer() {
                if (Phase >= 0)
                    throw new ApplicationException("transmission in progress");
                Phase++;
            }

            /// <summary>
            /// During transmission, the garbage collector handles (for pinning) for the receive buffers.
            /// index: into the keys-list of <see cref="RecvBuffers"/>
            /// </summary>
            GCHandle[] RecvBuffersLock;

            /// <summary>
            /// MPI requests of send operations (count is number of elements in <see cref="SendBuffers"/>),
            /// followed by
            /// MPI requests of receive operations (count is number of elements in <see cref="RecvBuffers"/>);
            /// </summary>
            MPI_Request[] MPIRequests;

            /// <summary>
            /// Phase 2 (of 4) of a complete transmission: starting the MPI - transmission.
            /// </summary>
            /// <remarks>
            ///  this method returns immediately, i.e. it does not block until the communication is finished
            /// </remarks>
            virtual public void StartTransmissionImReturn() {
                if (Phase != 0)
                    throw new ApplicationException("FillSendBuffer() must be called first!");
                Phase++;

                int MyRank;
                csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out MyRank);

                // init Requests - Array
                {
                    if(MPIRequests == null)
                        MPIRequests = new MPI_Request[SendBuffers.Count + RecvBuffers.Count];

                    for (int i = MPIRequests.Length - 1; i >= 0; i--) {
                        MPIRequests[i] = csMPI.Raw.MiscConstants.MPI_REQUEST_NULL;
                    }
                }

                // Garbage - collector lock & transmission init
                IList<int> procs = SendBuffers.Keys;
                for (int i = 0; i < procs.Count; i++) {
                    csMPI.Raw.Issend(SendBuffers[procs[i]], // Marshal.UnsafeAddrOfPinnedArrayElement(SendBuffer, 0),
                        SendBuffersLengths[procs[i]],    // SendBuffer.Length,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        procs[i],
                        666 + MyRank * 2,
                        csMPI.Raw._COMM.WORLD,
                        out MPIRequests[i]);
                    //Console.WriteLine("Proc #" + MyRank + ": sending " + SendBuffer.Length + " items to proc. #" + procs[i]);
                }
            }

            /// <summary>
            /// Phase 3 (of 4) of a complete transmission: init receive buffers
            /// </summary>
            public void InitReceiveImReturn() {
                if (Phase != 1)
                    throw new ApplicationException("StartTransmissionImReturn() must be called first!");
                Phase++;
                
                // Garbage - collector lock & transmission init
                IList<int> procs = RecvBuffers.Keys;
                int off = SendBuffers.Count;
                for (int i = 0; i < procs.Count; i++) {
                    double[] RecvBuffer = RecvBuffers[procs[i]];
                    csMPI.Raw.Irecv(Marshal.UnsafeAddrOfPinnedArrayElement(RecvBuffer, 0),
                        RecvBuffer.Length,
                        csMPI.Raw._DATATYPE.DOUBLE,
                        procs[i],
                        666 + procs[i] * 2,
                        csMPI.Raw._COMM.WORLD,
                        out MPIRequests[i+off]);
                }

            }

            /// <summary>
            /// Phase 4 (of 4) of a complete transmission: Blocks until communication is complete
            /// </summary>
            /// <param name="OnReceive">
            /// optional callback function, which is called when a send buffer has been received.
            /// </param>
            internal void WaitCommFinish(OnExternalBlockReceived OnReceive) {
                if (Phase != 2)
                    throw new ApplicationException("InitReceiveImReturn() must be called first!");
                Phase = -1;
                
                // wait for incoming transmission
                int index;
                int off = SendBuffers.Count;
                MPI_Status status;
                IList<int> procs = RecvBuffers.Keys;
                while (true) {
                    csMPI.Raw.Waitany(MPIRequests.Length, MPIRequests, out index, out status);

                    if (index == csMPI.Raw.MiscConstants.UNDEFINED) {
                        return;
                    } else if (index >= 0 && index < off) {
                        // some send has been completed
                    } else {
                        int ind = index - off;
                        int proc = procs[ind];
                        OnReceive(proc, Marshal.UnsafeAddrOfPinnedArrayElement(RecvBuffers[proc],0));
                    }
                } 

                //Waitall - version: a bit less efficient
                //if( statussies == null)
                //    statussies =  new MPI_Status[MPIRequests.Length];
                //csMPI.Raw.Waitall(MPIRequests.Length, MPIRequests, statussies);
                //for (int l = SendBuffersLock.Length - 1; l >= 0; l--)
                //    SendBuffersLock[l].Free();
                //for (int l = RecvBuffersLock.Length - 1; l >= 0; l--)
                //    RecvBuffersLock[l].Free();
                //for (int l = procs.Count - 1; l >= 0; l--) {
                //    int proc = procs[l];
                //    OnReceive(proc, RecvBuffers[proc]);
                //}
            }

            

            //MPI_Status[] statussies;

            #region IDisposable Members

            /// <summary>
            /// release the buffers, if necessary
            /// </summary>
            virtual public void Dispose() {
                if (RecvBuffersLock != null) {
                    foreach (GCHandle rbh in RecvBuffersLock)
                        rbh.Free();
                    RecvBuffersLock = null;
                }
            }

            #endregion
        }
    }
}
