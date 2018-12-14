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

using ilPSP;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;

namespace BoSSS.Foundation.SpecFEM {

    /// <summary>
    /// MPI data exchange for Spectral Elements
    /// </summary>
    public class Transceiver : IDisposable {

        public Transceiver(BoSSS.Foundation.SpecFEM.SpecFemBasis basis) {
            m_Basis = basis;

            List<Tuple<int, IntPtr>> _SendBuffers = new List<Tuple<int, IntPtr>>();
            List<Tuple<int, IntPtr>> _ReceiveBuffers = new List<Tuple<int, IntPtr>>();

            int Size = basis.GridDat.MpiSize;
            for (int rnk = 0; rnk < Size; rnk++) {
                if (basis.MPI_SendLists[rnk] != null) {
                    int L = basis.MPI_SendLists[rnk].Length;
                    _SendBuffers.Add(new Tuple<int, IntPtr>(rnk, Marshal.AllocHGlobal(sizeof(double) * L)));
                }
                if (basis.MPI_InsertLists[rnk] != null) {
                    int L = basis.MPI_InsertLists[rnk].Length;
                    _ReceiveBuffers.Add(new Tuple<int, IntPtr>(rnk, Marshal.AllocHGlobal(sizeof(double) * L)));
                }

            }
            this.SendBuffers = _SendBuffers.ToArray();
            this.ReceiveBuffers = _ReceiveBuffers.ToArray();
        }
        

        ~Transceiver() {
            this.Dispose();
        }


        SpecFemBasis m_Basis;

        Tuple<int, IntPtr>[] SendBuffers;
        Tuple<int, IntPtr>[] ReceiveBuffers;
        public void AccumulateGather(MultidimensionalArray A) {
            if (A.Dimension != 1)
                throw new ArgumentException();
            int K = m_Basis.NoOfLocalNodes;
            if (A.GetLength(0) != K)
                throw new ArgumentException();

            // for accumulate-gather, the role of send and insert lists is reversed!
            int[][] SendLists = m_Basis.MPI_InsertLists;
            int[][] InsertLists = m_Basis.MPI_SendLists;
            var _SendBuffers = this.ReceiveBuffers;
            var _ReceiveBuffers = this.SendBuffers;


            int NoOf_PsendTo = _SendBuffers.Length;    // Number of processes to send to
            int NoOf_PrvcFrm = _ReceiveBuffers.Length; // Number of processes to receive from
            Debug.Assert(NoOf_PsendTo == SendLists.Where(sl => sl != null).Count());
            Debug.Assert(NoOf_PrvcFrm == InsertLists.Where(sl => sl != null).Count());


            MPI_Request[] req = new MPI_Request[NoOf_PsendTo + NoOf_PrvcFrm];

            // set up non-blocking receive
            // ===========================

            for (int i = 0; i < NoOf_PrvcFrm; i++) {
                Tuple<int, IntPtr> kv = _ReceiveBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                csMPI.Raw.Irecv(Buffer, InsertLists[Rank].Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[NoOf_PsendTo + i]);
            }


            // initiate sending
            // ================

            for (int i = 0; i < NoOf_PsendTo; i++) {
                Tuple<int, IntPtr> kv = _SendBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                int[] SndList = SendLists[Rank];

                unsafe {
                    int L = SndList.Length;
                    double* p0 = (double*)Buffer;

                    for (int l = 0; l < L; l++) {
                        *p0 = A[SndList[l]];
                        p0++;
                    }
                }

                csMPI.Raw.Issend(Buffer, SndList.Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[i]);
            }


            // wait for operations to complete
            // ===============================

            int ReqCount = req.Length;
            while (ReqCount > 0) {
                int index;
                MPI_Status status;
                csMPI.Raw.Waitany(NoOf_PsendTo + NoOf_PrvcFrm, req, out index, out status);
                ReqCount--;

                if (index < NoOf_PsendTo)
                    // send finished
                    continue;

                // else: receive finished
                Tuple<int, IntPtr> kv = _ReceiveBuffers[index - NoOf_PsendTo];
                int[] InsList = InsertLists[kv.Item1];

                unsafe {
                    int L = InsList.Length;
                    double* p0 = (double*)kv.Item2;

                    for (int l = 0; l < L; l++) {
                        A[InsList[l]] += *p0;
                        p0++;
                    }
                }
            }
        }


        public void MaxGather(MultidimensionalArray A) {
            if (A.Dimension != 1)
                throw new ArgumentException();
            int K = m_Basis.NoOfLocalNodes;
            if (A.GetLength(0) != K)
                throw new ArgumentException();

            // for accumulate-gather, the role of send and insert lists is reversed!
            int[][] SendLists = m_Basis.MPI_InsertLists;
            int[][] InsertLists = m_Basis.MPI_SendLists;
            var _SendBuffers = this.ReceiveBuffers;
            var _ReceiveBuffers = this.SendBuffers;


            int NoOf_PsendTo = _SendBuffers.Length;    // Number of processes to send to
            int NoOf_PrvcFrm = _ReceiveBuffers.Length; // Number of processes to receive from
            Debug.Assert(NoOf_PsendTo == SendLists.Where(sl => sl != null).Count());
            Debug.Assert(NoOf_PrvcFrm == InsertLists.Where(sl => sl != null).Count());


            MPI_Request[] req = new MPI_Request[NoOf_PsendTo + NoOf_PrvcFrm];

            // set up non-blocking receive
            // ===========================

            for (int i = 0; i < NoOf_PrvcFrm; i++) {
                Tuple<int, IntPtr> kv = _ReceiveBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                csMPI.Raw.Irecv(Buffer, InsertLists[Rank].Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[NoOf_PsendTo + i]);
            }


            // initiate sending
            // ================

            for (int i = 0; i < NoOf_PsendTo; i++) {
                Tuple<int, IntPtr> kv = _SendBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                int[] SndList = SendLists[Rank];

                unsafe {
                    int L = SndList.Length;
                    double* p0 = (double*)Buffer;

                    for (int l = 0; l < L; l++) {
                        *p0 = A[SndList[l]];
                        p0++;
                    }
                }

                csMPI.Raw.Issend(Buffer, SndList.Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[i]);
            }


            // wait for operations to complete
            // ===============================

            int ReqCount = req.Length;
            while (ReqCount > 0) {
                int index;
                MPI_Status status;
                csMPI.Raw.Waitany(NoOf_PsendTo + NoOf_PrvcFrm, req, out index, out status);
                ReqCount--;

                if (index < NoOf_PsendTo)
                    // send finished
                    continue;

                // else: receive finished
                Tuple<int, IntPtr> kv = _ReceiveBuffers[index - NoOf_PsendTo];
                int[] InsList = InsertLists[kv.Item1];

                unsafe {
                    int L = InsList.Length;
                    double* p0 = (double*)kv.Item2;

                    for (int l = 0; l < L; l++) {
                        A[InsList[l]] = Math.Max(A[InsList[l]], *p0);
                        p0++;
                    }
                }
            }
        }

        public void Scatter(MultidimensionalArray A) {
            if (A.Dimension != 1)
                throw new ArgumentException();
            int K = m_Basis.NoOfLocalNodes;
            if (A.GetLength(0) != K)
                throw new ArgumentException();

            int[][] SendLists = m_Basis.MPI_SendLists;
            int[][] InsertLists = m_Basis.MPI_InsertLists;


            int NoOf_PsendTo = this.SendBuffers.Length;    // Number of processes to send to
            int NoOf_PrvcFrm = this.ReceiveBuffers.Length; // Number of processes to receive from

            Debug.Assert(NoOf_PsendTo == SendLists.Where(sl => sl != null).Count());
            Debug.Assert(NoOf_PrvcFrm == InsertLists.Where(sl => sl != null).Count());

            MPI_Request[] req = new MPI_Request[NoOf_PsendTo + NoOf_PrvcFrm];

            // set up non-blocking receive
            // ===========================

            for (int i = 0; i < NoOf_PrvcFrm; i++) {
                Tuple<int, IntPtr> kv = this.ReceiveBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                csMPI.Raw.Irecv(Buffer, InsertLists[Rank].Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[NoOf_PsendTo + i]);
            }


            // initiate sending
            // ================

            for (int i = 0; i < NoOf_PsendTo; i++) {
                Tuple<int, IntPtr> kv = this.SendBuffers[i];
                int Rank = kv.Item1;
                IntPtr Buffer = kv.Item2;

                int[] SndList = SendLists[Rank];

                unsafe {
                    int L = SndList.Length;
                    double* p0 = (double*)Buffer;

                    for (int l = 0; l < L; l++) {
                        *p0 = A[SndList[l]];
                        p0++;
                    }
                }

                csMPI.Raw.Issend(Buffer, SndList.Length, csMPI.Raw._DATATYPE.DOUBLE, Rank, 2341, csMPI.Raw._COMM.WORLD, out req[i]);
            }


            // wait for operations to complete
            // ===============================

            int ReqCount = req.Length;
            while (ReqCount > 0) {
                int index;
                MPI_Status status;
                csMPI.Raw.Waitany(NoOf_PsendTo + NoOf_PrvcFrm, req, out index, out status);
                ReqCount--;

                if (index < NoOf_PsendTo)
                    // send finished
                    continue;

                // else: receive finished
                Tuple<int, IntPtr> kv = this.ReceiveBuffers[index - NoOf_PsendTo];
                int[] InsList = InsertLists[kv.Item1];

                unsafe {
                    int L = InsList.Length;
                    double* p0 = (double*)kv.Item2;

                    for (int l = 0; l < L; l++) {
                        A[InsList[l]] = *p0;
                        p0++;
                    }
                }
            }
        }





        public void Dispose() {
            if (SendBuffers != null) {
                foreach (var kv in SendBuffers) {
                    Marshal.FreeHGlobal(kv.Item2);
                }
                SendBuffers = null;
            }
            if (ReceiveBuffers != null) {
                foreach (var kv in ReceiveBuffers) {
                    Marshal.FreeHGlobal(kv.Item2);
                }
                ReceiveBuffers = null;
            }
        }
    }
}
