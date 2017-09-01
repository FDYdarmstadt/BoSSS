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
using System.Runtime.InteropServices;
using System.Linq;

using MPI.Wrappers;
using BoSSS.Foundation.Grid;

namespace BoSSS.Foundation.Comm {
    
    /// <summary>
    /// performs transmission of DG field coordinates (in the border cells of this process) to 
    /// other processes (Transmitter)
    /// and receives field coordinates from other processes and putts the received values
    /// into the corresponding memory for external cells.
    /// </summary>
    public class Transceiver {

        /// <summary>
        /// checks if the call to <see cref="TransceiveStartImReturn"/>
        /// is matched by a call to <see cref="TransceiveFinish"/>
        /// </summary>
        ~Transceiver() {
            if (started)
                throw new ApplicationException("communication was not finished.");
        }

        int m_MyRank;

        int m_Size;

        IParallelization m_parallel;


        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="TRXfields">
        /// fields which should be exchanged
        /// </param>
        public Transceiver(params DGField[] TRXfields)
            : this((ICollection<DGField>)TRXfields) { }


        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="TRXfields">
        /// fields which should be exchanged
        /// </param>
        public Transceiver(ICollection<DGField> TRXfields) {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out m_MyRank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out m_Size);


            var grd = TRXfields.First().GridDat;
            foreach (var x in TRXfields) {
                if (!object.ReferenceEquals(x.GridDat, grd))
                    throw new ArgumentException("all fields must be assigned to the same grid.");
            }

            this.m_parallel = grd.iParallel;


            if (m_Size == 1) return; // nothing to communicate

            // alloc some arrays
            // =================

            m_ReceiveBuffers = new double[m_Size][];
            m_SendBuffers = new double[m_Size][];
            m_ReceiveBufferPin = new GCHandle[m_parallel.ProcessesToReceiveFrom.Length];
            m_SendBufferPin = new GCHandle[m_parallel.ProcessesToSendTo.Length];
            m_Requests = new MPI_Request[m_parallel.ProcessesToSendTo.Length + m_parallel.ProcessesToReceiveFrom.Length];

            // collect fields
            // ==============

            //m_NumberOfItemsPerCell = 0;
            foreach (DGField f in TRXfields) {
                if (m_TRXFields.Contains(f))
                    continue;
                    //throw new ArgumentException("each field can occur only once.");

                m_TRXFields.Add(f);
                //m_NumberOfItemsPerCell += f.NoOfCoordinatesPerCell;
            }

            //// compute offset indices
            //// ======================

            //m_Offset = new int[master.Size][];
            
            //for( int p = 0; p < master.Size; p++) {
            //    int[] commList = master.m_SendCommLists[p];

            //    if (commList != null) {
            //        m_Offset[p] = new int[m_TRXFields.Count];

            //        int offcur = 0;
            //        for (int i = 0; i < m_TRXFields.Count; i++) {
            //            m_Offset[p][i] = offcur;
            //            offcur += m_TRXFields[i].NoOfCoordinatesPerCell*commList.Length;
            //        }
            //    }
            //}
        }



        int FillSendBuffer(int p) {
            int i0 = 0;

            int TotLen = 0;
            foreach(DGField fld in m_TRXFields) TotLen += fld.GetMPISendBufferSize(p);
            if (m_SendBuffers[p] == null)
                m_SendBuffers[p] = new double[TotLen];
            if (m_SendBuffers[p].Length < TotLen)
                m_SendBuffers[p] = new double[TotLen];
            
            for (int f = 0; f < m_TRXFields.Count; f++) {
                DGField fld = m_TRXFields[f];
                int Len = fld.GetMPISendBufferSize(p);
                fld.FillMPISendBuffer(p, m_SendBuffers[p], i0);
                i0 += Len;
            }

            return TotLen;
        }


        int GetReceiveBufferLen(int p) {
            int TotLen = 0;
            foreach (DGField fld in m_TRXFields) TotLen += fld.GetMPIRecvBufferSize(p);
            
            if (m_ReceiveBuffers[p] == null)
                m_ReceiveBuffers[p] = new double[TotLen];
            if (m_ReceiveBuffers[p].Length < TotLen)
                m_ReceiveBuffers[p] = new double[TotLen];

            return TotLen;
        }

        int WriteReceivedData(int p) {

            int i0 = 0;
            for (int f = 0; f < m_TRXFields.Count; f++) {
                DGField fld = m_TRXFields[f];
                int Len = fld.CopyFromMPIrecvBuffer(p, m_ReceiveBuffers[p],i0);
                i0 += Len;
            }

            return i0;
        }



        /// <summary>
        /// 
        /// </summary>
        List<DGField> m_TRXFields = new List<DGField>();


        double[][] m_SendBuffers;

        double[][] m_ReceiveBuffers;

        GCHandle[] m_SendBufferPin;
        
        GCHandle[] m_ReceiveBufferPin;

        MPI_Request[] m_Requests;

        bool started = false;

        /// <summary>
        /// initiates the send/receive - processes and returns immediately;
        /// Every call of this method must be matched by a later call to <see cref="TransceiveFinish"/>;
        /// </summary>
        public void TransceiveStartImReturn() {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            if (started)
                throw new ApplicationException("transceive must be finished first.");
            started = true;

            if (m_Size == 1) return; // nothing to communicate

            // start sending
            // =============
            int[] sndProc;
            {
                sndProc = m_parallel.ProcessesToSendTo;

                for (int i = 0; i < sndProc.Length; i++) {
                    int proc = sndProc[i];

                    int BufLen = FillSendBuffer(proc);

                    m_SendBufferPin[i] = GCHandle.Alloc(m_SendBuffers[proc],GCHandleType.Pinned);
                    csMPI.Raw.Issend(Marshal.UnsafeAddrOfPinnedArrayElement(m_SendBuffers[proc], 0),
                        BufLen, csMPI.Raw._DATATYPE.DOUBLE, proc,
                        42 + m_MyRank,
                        csMPI.Raw._COMM.WORLD,
                        out m_Requests[i]);
                }
            }

            // init receiving
            // ==============
            {
                int[] rcvProc = m_parallel.ProcessesToReceiveFrom;

                for (int i = 0; i < rcvProc.Length; i++) {
                    int proc = rcvProc[i];

                    int BufLen = GetReceiveBufferLen(proc);

                    m_ReceiveBufferPin[i] = GCHandle.Alloc(m_ReceiveBuffers[proc], GCHandleType.Pinned);
                    csMPI.Raw.Irecv(Marshal.UnsafeAddrOfPinnedArrayElement(m_ReceiveBuffers[proc], 0),
                        BufLen, csMPI.Raw._DATATYPE.DOUBLE, proc,
                        42 + proc,
                        csMPI.Raw._COMM.WORLD,
                        out m_Requests[i + sndProc.Length]);
                }
            }
        }

        /// <summary>
        /// blocks until the send/receive - processes started by <see cref="TransceiveStartImReturn"/>
        /// are complete and returns;
        /// </summary>
        public void TransceiveFinish() {
            ilPSP.MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);
            if (!started)
                throw new ApplicationException("transceive must be started first.");
            started = false;
           
            if (m_Size == 1) return; // nothing to communicate

            int[] sndPrc = m_parallel.ProcessesToSendTo;
            int[] rcvPrc = m_parallel.ProcessesToReceiveFrom;

            int index;
            MPI_Status status;
            for (int i = 0; i < (sndPrc.Length + rcvPrc.Length); i++) {
                csMPI.Raw.Waitany(m_Requests.Length, m_Requests, out index, out status);
                if (index < sndPrc.Length) {
                    m_SendBufferPin[index].Free();
                } else {
                    index -= sndPrc.Length;
                    m_ReceiveBufferPin[index].Free();
                    WriteReceivedData(rcvPrc[index]);
                }
            }
        }
    }
}
