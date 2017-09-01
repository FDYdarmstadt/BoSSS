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
using System.Linq;
using System.Diagnostics;
using MPI.Wrappers;
using ilPSP.Utils;

namespace ilPSP {

    /// <summary>
    /// This light weighted class describes the partitioning of a vector or a
    /// matrix over a range of MPI processes.
    /// </summary>
    public class Partitioning : ICloneable, IEquatable<IPartitioning>, IPartitioning {

        MPI.Wrappers.MPI_Comm m_comm;

        /// <summary>
        /// ctor; uses the MPI world communicator; 
        /// </summary>
        /// <param name="localsize">
        /// Number of entries that should be stored in this MPI process.
        /// </param>
        /// <param name="_BlockSize">
        /// See <see cref="BlockSize"/>.
        /// </param>
        public Partitioning(int localsize)
            : this(localsize, MPI.Wrappers.csMPI.Raw._COMM.WORLD) {
        }

        /// <summary>
        /// used for cloning
        /// </summary>
        private Partitioning() {
        }

        /// <summary>
        /// ctor.
        /// </summary>
        /// <param name="localsize">
        /// number of entries that should be stored in this MPI process
        /// </param>
        /// <param name="c">
        /// MPI communicator, <see cref="MPI_Comm"/>.
        /// </param>
        /// <param name="_BlockSize">
        /// See <see cref="BlockSize"/>.
        /// </param>
        public Partitioning(int localsize, MPI.Wrappers.MPI_Comm c) {
            m_comm = c;
            MPICollectiveWatchDog.Watch(c);
            csMPI.Raw.Comm_Rank(c, out m_rank);
            csMPI.Raw.Comm_Size(c, out m_size);

            m_LocalLengths = new int[m_size];
            int[] ll = { localsize };

            unsafe {
                fixed (void* pSndBuf = ll, pRcvBuf = m_LocalLengths) {
                    csMPI.Raw.Allgather((IntPtr)pSndBuf, 4, csMPI.Raw._DATATYPE.BYTE,
                                     (IntPtr)pRcvBuf, 4, csMPI.Raw._DATATYPE.BYTE,
                                     m_comm);
                }
            }

            m_i0Offset = new int[m_size + 1];
            m_i0Offset[0] = 0;
            for (int i = 1; i <= m_size; i++)
                m_i0Offset[i] = m_i0Offset[i - 1] + m_LocalLengths[i - 1];
        }

        /// <summary>
        /// rank of actual process in MPI communicator
        /// </summary>
        int m_rank;

        /// <summary>
        /// rank of actual process in MPI communicator
        /// </summary>
        public int MpiRank {
            get {
                return m_rank;
            }
        }

        /// <summary>
        /// no. of MPI processes in communicator
        /// </summary>
        int m_size;

        /// <summary>
        /// no. of MPI processes in communicator
        /// </summary>
        public int MpiSize {
            get {
                return m_size;
            }
        }

        /// <summary>
        /// true, if some index <paramref name="i"/> is within the local range of this MPI Process,
        /// i.e <paramref name="i"/> is greater or equal to <see cref="i0"/> and smaller than <see cref="i0"/>+<see cref="LocalLength"/>;
        /// </summary>
        public bool IsInLocalRange(int i) {
            return (i >= i0 && i < m_i0Offset[m_rank + 1]);
        }

        /// <summary>
        /// true, if some index <paramref name="i"/> is within the local range of this MPI Process,
        /// i.e <paramref name="i"/> is greater or equal to <see cref="i0"/> and smaller than <see cref="i0"/>+<see cref="LocalLength"/>;
        /// </summary>
        public bool IsInLocalRange(long i) {
            return (i >= i0 && i < m_i0Offset[m_rank + 1]);
        }

       

        /// <summary>
        /// subtracts <see cref="i0"/> from <paramref name="iGlob"/>;
        /// </summary>
        public int TransformIndexToLocal(int iGlob) {
            return (iGlob - (int)m_i0Offset[m_rank]);
        }

        
        /// <summary>
        /// Equality of partitioning.
        /// </summary>
        public virtual bool Equals(IPartitioning o) {
            return this.EqualsPartition(o);
        }

       

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override int GetHashCode() {
            return (int)this.TotalLength;
        }

        /// <summary>
        /// MPI Communicator on which this object lives on
        /// </summary>
        public MPI_Comm MPI_Comm {
            get {
                //throw new NotImplementedException("will come soon; currently, it lives on MPI_COMM_WORLD;");
                return m_comm;
            }
        }

        /// <summary>
        /// the first global index that is stored on the actual MPI process
        /// </summary>
        public int i0 {
            get {
                return (int)m_i0Offset[m_rank];
            }
        }

        /// <summary>
        /// the first global index that is stored on the NEXT MPI process
        /// </summary>
        public int iE {
            get {
                return i0 + LocalLength;
            }
        }


        /// <summary>
        /// Returns the number of entries which are stored on
        /// this processor.
        /// </summary>
        public int LocalLength {
            get {
                return m_LocalLengths[m_rank];
            }
        }

        /// <summary>
        /// index offset for some processor
        /// </summary>
        /// <param name="proc">
        /// process index; can be equal to number of processes (i.e. highest process
        /// rank +1); In this case, the <see cref="TotalLength"/> is returned.
        /// </param>
        /// <returns>index of the first permutation entry stored by processor <paramref name="proc"/></returns>
        public int GetI0Offest(int proc) {
            return m_i0Offset[proc];
        }

        /// <summary>
        /// returns the number of entries which are stored by
        /// processor <paramref name="proc"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns></returns>
        public int GetLocalLength(int proc) {
            return m_LocalLengths[proc];
        }

        
        /// <summary>
        /// Total length of the partition over all processes, i.e. sum of <see cref="LocalLength"/> over all MPI processors.
        /// </summary>
        public int TotalLength {
            get {
                Debug.Assert(m_i0Offset.Length == m_LocalLengths.Length + 1);
                return m_i0Offset[m_LocalLengths.Length];
            }
        }

        /// <summary>
        /// <see cref="IPartitioning.IsMutable"/>
        /// </summary>
        virtual public bool IsMutable {
            get {
                return false;
            }
        }

        /// <summary>
        /// offsets for each process; size is equal to number of processors plus 1;
        /// first entry is always 0, last entry is equal to <see cref="m_TotalLength"/>;
        /// </summary>
        int[] m_i0Offset;

        /// <summary>
        /// number of elements stored in each process
        /// </summary>
        int[] m_LocalLengths;

        /// <summary>
        /// returns the process rank which stores the <paramref name="index"/>-th entry;
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public int FindProcess(long index) {
            return FindProcess((int)index);
        }

        /// <summary>
        /// returns the process rank which stores the <paramref name="index"/>-th entry;
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public int FindProcess(int index) {
            Debug.Assert(index >= 0);
            Debug.Assert(index < this.TotalLength);

            // Schuss ins schwarze
            if (index >= this.i0 && index < this.iE)
                return m_rank;

            var p = Array.BinarySearch<int>(m_i0Offset, index);
            if (p < 0) {
                p = (~p) - 1;
            }

            // Play it safe in case of potentially empty partitions
            while (m_i0Offset[p] == m_i0Offset[p + 1]) {
                p++;
            }
            Debug.Assert(p != this.m_rank);
            Debug.Assert(p >= 0);
            Debug.Assert(p < this.m_size);
            Debug.Assert(m_i0Offset[p] <= index);
            Debug.Assert(index < m_i0Offset[p + 1]);

            return p;

            //int targProc = 0;
            //while (index >= m_i0Offset[targProc]) targProc++;
            //targProc--;
            //return targProc;
        }

        /// <summary>
        /// creates a non-shallow copy of this object.
        /// </summary>
        /// <returns></returns>
        public Partitioning CloneAs() {
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            Partitioning ret = new Partitioning();
            ret.m_i0Offset = (int[])this.m_i0Offset.Clone();
            ret.m_LocalLengths = (int[])this.m_LocalLengths.Clone();
            ret.m_rank = this.m_rank;
            ret.m_size = this.m_size;
            ret.m_comm = this.m_comm;
            return ret;
        }

        /// <summary>
        /// creates a non-shallow copy of this object.
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            return CloneAs();
        }

        /// <summary>
        /// equality, see <see cref="Equals(Partitioning)"/>
        /// </summary>
        public override bool Equals(object other) {
            return Equals(other as IPartitioning);
        }

        /// <summary>
        /// <see cref="IPartitioning.GetImmutablePartition"/>
        /// </summary>
        virtual public IPartitioning GetImmutablePartition() {
            if (!this.IsMutable)
                return this;
            else
                return new Partitioning(this.LocalLength, this.MPI_Comm);
        }
    }
}
