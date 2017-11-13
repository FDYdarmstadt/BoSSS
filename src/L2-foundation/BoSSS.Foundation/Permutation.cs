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
using System.Diagnostics;
using System.Runtime.InteropServices;
using ilPSP;
using ilPSP.Utils;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Comm {

    /// <summary>
    /// An MPI-parallel data structure to store permutations.
    /// </summary>
    public sealed class Permutation : ICloneable {

        /// <summary>
        /// creates an empty Permutation (filled with invalid, negative entries)
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="localLength">number of items that are stored in this process</param>
        public Permutation(int localLength, MPI.Wrappers.MPI_Comm comm)
            : this(new long[localLength], comm) {
            for (int i = 0; i < m_Values.Length; i++)
                m_Values[i] = -1;
        }

        /// <summary>
        /// creates a permutation with user-defined entries;
        /// The "correctness" (each entry occurs exactly once over all processes,
        /// and the lowest entry is 0 an no number is skipped) is NOT tested;
        /// </summary>
        /// <param name="comm"></param>
        /// <param name="values">permutation values that are stored on this process</param>
        public Permutation(long[] values, MPI.Wrappers.MPI_Comm comm) {
            ilPSP.MPICollectiveWatchDog.Watch(comm);

            m_Values = values;

            m_Comm = comm;
            m_Partition = new Partitioning(values.Length, comm);

#if DEBUG
            long max = m_Partition.TotalLength;
            int I = values.Length;
            for (int i = 0; i < I; i++) {
                if (values[i] < 0 || values[i] >= max)
                    throw new ArgumentOutOfRangeException();
            }
#endif

            //m_LocalLengths = new int[master.Size];
            //int[] ll = { values.Length };

            //unsafe {
            //    fixed (void* pSndBuf = ll, pRcvBuf = m_LocalLengths) {
            //        csMPI.Raw.Allgather((IntPtr) pSndBuf, 4, MPI_Datatype.BYTE,
            //                         (IntPtr) pRcvBuf, 4, MPI_Datatype.BYTE,
            //                         csMPI.Raw.MPI_COMM_WORLD);
            //    }
            //}

            //int size = m_Master.Size;
            //m_i0Offset = new long[size+1];
            //m_i0Offset[0] = 0;
            //for (int i = 1; i <= size; i++)
            //    m_i0Offset[i] = m_i0Offset[i-1] + m_LocalLengths[i-1];
            //m_TotalLength = m_i0Offset[size];
        }

        /// <summary>
        /// <see cref="Partitioning"/>
        /// </summary>
        private Partitioning m_Partition;

        /// <summary>
        /// partition over MPI processes for this permutation
        /// </summary>
        public Partitioning Partitioning {
            get {
                return m_Partition;
            }
        }

        /// <summary>
        /// constructor used by <see cref="Invert"/>, creates a NON-shallow copy;
        /// </summary>
        /// <param name="o"></param>
        private Permutation(Permutation o) {

            m_Comm = o.m_Comm;
            m_Partition = o.m_Partition;
            //m_TotalLength = o.m_TotalLength;
            //m_i0Offset = (long[]) o.m_i0Offset.Clone();
            //m_LocalLengths = (int[]) o.m_LocalLengths.Clone();
        }

        /*

        /// <summary>
        /// <see cref="StorageGuid"/>
        /// </summary>
        Guid m_StorageGuid = Guid.Empty;


        /// <summary>
        /// Identifier of this permutation within the storage system
        /// (<see cref="BoSSS.Foundation.IO.DatabaseDriver"/>). If this object is just created normally
        /// (i.e. by a new-operation) this Guid is empty.
        /// This value is initialized if this object is stored or loaded.
        /// </summary>
        public Guid StorageGuid {
            get {
                return m_StorageGuid;
            }
        }


        /// <summary>
        /// sets <see cref="StorageGuid"/> to <paramref name="g"/>;
        /// </summary>
        /// <param name="g"></param>
        internal void SetStorageGuid(Guid g) {
            m_StorageGuid = g;
        }
        */

        /// <summary>
        /// MPI communicator on which this permutation 'lives' on
        /// </summary>
        MPI.Wrappers.MPI_Comm m_Comm;

        /// <summary>
        /// total length of the permutation over all processes
        /// </summary>
        public long TotalLength {
            get {
                return m_Partition.TotalLength;
            }
        }

        /// <summary>
        /// index offset for the first entry of the permutation stored by
        /// this process
        /// </summary>
        public long i0Offset {
            get {
                return m_Partition.i0;// m_i0Offset[m_Master.MyRank]; 
            }
        }

        /// <summary>
        /// returns the number of permutation entries which are stored by
        /// this processor;
        /// </summary>
        public int LocalLength {
            get {
                return m_Partition.LocalLength;
            }
        }

        /// <summary>
        /// <see cref="Values"/>
        /// </summary>
        long[] m_Values;

        /// <summary>
        /// permutation values in the actual process; user must ensure that
        /// this values represent a valid permutation among all processes,
        /// otherwise the behavior of the operations is undefined;
        /// </summary>
        public long[] Values {
            get {
                return m_Values;
            }
        }

        /// <summary>
        /// Helper structure to transmit permutation entries over the 
        /// network
        /// </summary>
        [StructLayout(LayoutKind.Sequential)]
        private struct PermutationEntry {

            /// <summary>
            /// Argument 
            /// </summary>
            public long Index;

            /// <summary>
            /// value of the permutation for argument <see cref="Index"/>
            /// </summary>
            public long PermVal;
        }

        /// <summary>
        /// performs a parallel inversion of this permutation
        /// </summary>
        /// <returns></returns>
        public Permutation Invert() {
            using(new FuncTrace()) {
                Permutation inv = new Permutation(this);
                int localLength = m_Partition.LocalLength; //m_LocalLengths[m_Master.MyRank];
                int size;
                MPI.Wrappers.csMPI.Raw.Comm_Size(m_Comm, out size);


                inv.m_Values = new long[localLength];
                long myi0Offset = m_Partition.i0; //m_i0Offset[m_Master.MyRank];

                Many2ManyMessenger<PermutationEntry> m2m = new Many2ManyMessenger<PermutationEntry>(m_Comm);
                int[] NoOfItemsToSent = new int[size];


                // calc i0's for each process
                // ==========================

                long[] i0 = new long[size + 1];
                for(int i = 1; i <= size; i++)
                    i0[i] = i0[i - 1] + m_Partition.GetLocalLength(i - 1);// m_LocalLengths[i - 1];



                // do local inversion and collect items 
                // ====================================
                List<PermutationEntry>[] itemsToSend = new List<PermutationEntry>[size];
                for(int p = 0; p < size; p++)
                    itemsToSend[p] = new List<PermutationEntry>();

                for(int i = 0; i < localLength; i++) {

                    // decide wether inversion of entry i is local or has to be transmitted ...
                    if(m_Values[i] >= myi0Offset && m_Values[i] < (myi0Offset + localLength)) {
                        // local 
                        inv.m_Values[this.m_Values[i] - myi0Offset] = i + myi0Offset;

                    } else {
                        // transmit

                        // find target processor
                        int targProc = m_Partition.FindProcess(m_Values[i]);

                        NoOfItemsToSent[targProc]++;

                        // collect item
                        PermutationEntry pe;
                        pe.Index = myi0Offset + i;
                        pe.PermVal = m_Values[i];
                        itemsToSend[targProc].Add(pe);
                    }
                }


                // setup messenger
                // ===============

                for(int i = 0; i < size; i++) {
                    if(NoOfItemsToSent[i] > 0)
                        m2m.SetCommPath(i, NoOfItemsToSent[i]);
                }
                m2m.CommitCommPaths();

                // transmit data
                // =============

                //csMPI.Raw.Barrier(csMPI.Raw.MPI_COMM_WORLD);
                //Console.WriteLine("one");
                //csMPI.Raw.Barrier(csMPI.Raw.MPI_COMM_WORLD);
                //Console.WriteLine("two");
                //csMPI.Raw.Barrier(csMPI.Raw.MPI_COMM_WORLD);
                //Console.WriteLine("three");

                m2m.StartTransmission(1);

                for(int p = 0; p < size; p++) {
                    if(NoOfItemsToSent[p] <= 0)
                        continue;

                    //Many2ManyMessenger.Buffer<PermEntry> sndbuf = new Many2ManyMessenger.Buffer<PermEntry>(m2m, false, p);
                    PermutationEntry[] items = itemsToSend[p].ToArray();
                    m2m.SendBuffers(p).CopyFrom(items, 0);

                    m2m.TransmittData(p);
                }

                m2m.FinishBlocking();

                // invert received items
                // =====================


                for(int p = 0; p < size; p++) {
                    Many2ManyMessenger<PermutationEntry>.Buffer rcvbuf = m2m.ReceiveBuffers(p);
                    if(rcvbuf == null)
                        continue; // no data from process no. p

                    int cnt = rcvbuf.Count;

                    PermutationEntry[] items = new PermutationEntry[cnt];
                    rcvbuf.CopyTo(items, 0);


                    for(int i = 0; i < cnt; i++) {
                        inv.m_Values[items[i].PermVal - myi0Offset] = items[i].Index;
                    }
                }


                // finalize
                // ========
                m2m.Dispose();
                return inv;
            }
        }

        /// <summary>
        /// Permutation composition.
        /// </summary>
        /// <param name="left">left operand</param>
        /// <param name="right">right operand</param>
        /// <returns>
        /// The product of <paramref name="left"/>- and <paramref name="right"/>-permutation;
        /// the partitioning (<see cref="Permutation.Partitioning"/> is equal to the <paramref name="right"/> operand.
        /// </returns>
        public static Permutation operator *(Permutation left, Permutation right) {

            if (left.m_Comm != right.m_Comm) {
                throw new ArgumentException("left and right operand must share the same communication master.");
            }

            if (left.TotalLength != right.TotalLength) {
                throw new ArgumentException("left and right operand must have the same total length.");
            }

            Permutation result = new Permutation(right.m_Partition.LocalLength, right.m_Comm);
            left.EvaluatePermutation(right.m_Values, result.m_Values);

            return result;
        }

        /// <summary>
        /// parallel evaluator; For arguments which are not stored on this process,
        /// other processes are asked.
        /// </summary>
        /// <param name="Keys">input; indices at which the permutation should be evaluated;</param>
        /// <param name="Result">
        /// output; the i-th entry is the value of the permutation for index 
        /// <paramref name="Keys"/>[i];
        /// </param>
        public void EvaluatePermutation<T1, T2>(T1 Keys, T2 Result)
            where T1 : IList<long>
            where T2 : IList<long> //
        {
            using(new FuncTrace()) {
                if(Keys.Count != Result.Count) {
                    throw new ArgumentException("Keys and Result array must have the same number of elements");
                }

                int cnt = Keys.Count;
                int myrank;
                int size;
                MPI.Wrappers.csMPI.Raw.Comm_Size(this.m_Comm, out size);
                MPI.Wrappers.csMPI.Raw.Comm_Rank(this.m_Comm, out myrank);
                ilPSP.MPICollectiveWatchDog.Watch(this.m_Comm);


                // create objects for questions to other processors
                // ================================================

                List<int>[] Questions = new List<int>[size];
                List<int>[] iQuestions = new List<int>[size];
                for(int p = 0; p < size; p++) {
                    Questions[p] = new List<int>();
                    iQuestions[p] = new List<int>();
                }

                // local evaluation of the permutation mapping
                // ===========================================
                for(int i = 0; i < cnt; i++) {
                    long idx = Keys[i];
                    if(idx < 0 || idx >= TotalLength)
                        throw new ArgumentException("argument out of range; (argument values start at 0)");

                    long idxLocal = idx - i0Offset;
                    if(idxLocal >= 0 && idxLocal < LocalLength) {
                        // local evaluation

                        Result[i] = m_Values[idxLocal]; // that's all, folks
                    } else {
                        // another processor has to be asked

                        int targProcess = m_Partition.FindProcess(idx);
                        int idxLoc = (int)(idx - m_Partition.GetI0Offest(targProcess)); // m_i0Offset[targProcess]);
                        Questions[targProcess].Add(idxLoc);
                        iQuestions[targProcess].Add(i);
                    }
                }

                // --------------------
                // Ask other processors
                // --------------------

                // setup question messenger
                // ========================

                Many2ManyMessenger<int> m2mAsk = new Many2ManyMessenger<int>(m_Comm);
                for(int p = 0; p < size; p++) {
                    if(Questions[p].Count > 0) {
                        m2mAsk.SetCommPath(p, Questions[p].Count);
                    }
                }

                m2mAsk.CommitCommPaths();

                // transmit data
                // =============

                m2mAsk.StartTransmission(1);

                for(int p = 0; p < size; p++) {
                    if(Questions[p].Count > 0) {
                        m2mAsk.SendBuffers(p).CopyFrom(Questions[p].ToArray(), 0);
                        Questions[p] = null;
                        m2mAsk.TransmittData(p);
                    }
                }

                // wait until communication is finished
                // ====================================

                m2mAsk.FinishBlocking();

                // --------------------------
                // answer to other processors
                // --------------------------


                // setup answer messenger
                // ======================

                Many2ManyMessenger<long> m2mAnswer = new Many2ManyMessenger<long>(m_Comm);
                for(int p = 0; p < size; p++) {
                    if(m2mAsk.ReceiveBuffers(p) != null) {
                        m2mAnswer.SetCommPath(p, m2mAsk.ReceiveBuffers(p).Count);
                    }
                }
                m2mAnswer.CommitCommPaths();

                // transmit data
                // =============

                m2mAnswer.StartTransmission(1);

                for(int p = 0; p < size; p++) {
                    Many2ManyMessenger<int>.Buffer q = m2mAsk.ReceiveBuffers(p);
                    Many2ManyMessenger<long>.Buffer r = m2mAnswer.SendBuffers(p);

                    if(q != null) {
                        int _cnt = q.Count;

                        for(int i = 0; i < _cnt; i++) {
                            r[i] = m_Values[q[i]];
                        }

                        m2mAnswer.TransmittData(p);
                    }

                }
                m2mAsk.Dispose();


                // wait until communication is finished
                // ====================================

                m2mAnswer.FinishBlocking();

                // ------------------------------------------
                // evaluate the answers from other processors
                // ------------------------------------------

                for(int p = 0; p < size; p++) {
                    Many2ManyMessenger<long>.Buffer r = m2mAnswer.ReceiveBuffers(p);

                    if(r != null) {
                        // test
                        if(r.Count != iQuestions[p].Count)
                            throw new ApplicationException("internal error: mismatch between number of questions/answers");

                        int c = iQuestions[p].Count;
                        List<int> iQeus = iQuestions[p];

                        for(int i = 0; i < c; i++) {
                            Result[iQeus[i]] = r[i];
                        }


                    }
                }
                m2mAnswer.Dispose();
            }
        }

        [Serializable]
        class ApplyToVector_Helper<I> {
            public List<long> TargetIndices = new List<long>();
            public List<I> Items = new List<I>();
        }

        /// <summary>
        /// Resorts a vector according to this permutation, i.e. the
        /// <em>j</em>-th item of the input vector is copied to the
        /// <see cref="Values"/>[j]-th entry of the output vector.
        /// </summary>
        /// <param name="input">
        /// Input vector, length must be equal to the length of this permutation, unchanged on exit.
        /// </param>
        /// <param name="output">
        /// On exit, <paramref name="output"/>[<see cref="Values"/>[j]] = <paramref name="input"/>[j]
        /// </param>
        public void ApplyToVector<I>(IList<I> input, IList<I> output) {
            ApplyToVector(input, output, this.Partitioning);
        }

        /// <summary>
        /// Resorts a vector according to this permutation, i.e. the
        /// <em>j</em>-th item of the input vector is copied to the
        /// <see cref="Values"/>[j]-th entry of the output vector.
        /// </summary>
        /// <param name="input">
        /// Input vector, length must be equal to the length of this permutation, unchanged on exit.
        /// </param>
        /// <param name="output">
        /// On exit, <paramref name="output"/>[<see cref="Values"/>[j]] = <paramref name="input"/>[j]
        /// </param>
        public void ApplyToVector<I>(IList<I> input, IList<I> output, IPartitioning outputPartitioning) {
            using(new FuncTrace()) {
                if(input.Count != this.LocalLength)
                    throw new ArgumentException("wrong size of input vector.");
                if(output.Count != outputPartitioning.LocalLength)
                    throw new ArgumentException("wrong size of output vector.");

                long[] TargetInd = this.Values;
                // keys: processors which should receive data from this processor
                Dictionary<int, ApplyToVector_Helper<I>> sendData =
                    new Dictionary<int, ApplyToVector_Helper<I>>();

                int out_myI0 = outputPartitioning.i0;
                int out_nextI0 = out_myI0 + outputPartitioning.LocalLength;
                int J = this.Partitioning.LocalLength;

                for(int j = 0; j < J; j++) {
                    if(out_myI0 <= TargetInd[j] && TargetInd[j] < out_nextI0) {
                        // target index located on this processor
                        output[(int)(TargetInd[j] - out_myI0)] = input[j];
                    } else {
                        // target index located on other processor

                        int TargProc = outputPartitioning.FindProcess(TargetInd[j]);

                        ApplyToVector_Helper<I> sendData_TargProc = null;
                        if(!sendData.TryGetValue(TargProc, out sendData_TargProc)) {
                            sendData_TargProc = new ApplyToVector_Helper<I>();
                            sendData.Add(TargProc, sendData_TargProc);
                        }

                        sendData_TargProc.Items.Add(input[j]);
                        sendData_TargProc.TargetIndices.Add(TargetInd[j]);
                    }
                }

                var rcvData = SerialisationMessenger.ExchangeData(
                    sendData, MPI.Wrappers.csMPI.Raw._COMM.WORLD);

                foreach(var rcvPkt in rcvData.Values) {
                    int K = rcvPkt.Items.Count;
                    Debug.Assert(rcvPkt.Items.Count == rcvPkt.TargetIndices.Count);

                    for(int k = 0; k < K; k++) {
                        int locIdx = (int)(rcvPkt.TargetIndices[k]) - out_myI0;
                        output[locIdx] = rcvPkt.Items[k];
                    }
                }
            }
        }

        /// <summary>
        /// equality of permutations;
        /// </summary>
        public override bool Equals(object obj) {
            Permutation othr = obj as Permutation;
            if (othr == null)
                return false;

            if (!this.m_Partition.Equals(othr.m_Partition))
                return false;

            for (int i = this.m_Values.Length - 1; i >= 0; i--) {
                if (this.Values[i] != othr.m_Values[i])
                    return false;
            }


            return true;
        }

        /// <summary>
        /// some hash code
        /// </summary>
        public override int GetHashCode() {
            long code = 0;
            for (int i = Math.Max(20, this.m_Values.Length - 1); i >= 0; i--) {
                code += this.m_Values[i];
            }

            return (int)code;
        }

        /// <summary>
        /// Cloning.
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var R = new Permutation(this.LocalLength, this.m_Comm);
            R.m_Values = (long[]) (this.m_Values.Clone());
            return R;
        }
    }
}
