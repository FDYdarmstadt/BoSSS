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
using System.Diagnostics;
using MPI.Wrappers;

namespace ilPSP.LinSolvers.monkey.CPU {

    /// <summary>
    /// Reference implementation of <see cref="VectorBase"/>, that works completely on
    /// main processor/main memory.
    /// </summary>
    public class RefVector : VectorBase {

        /// <summary>
        /// constructor that allocates internal memory
        /// </summary>
        /// <param name="P"></param>
        public RefVector(IPartitioning P)
            : base(P) {
            m_Storage = new double[P.LocalLength];
        }

        /// <summary>
        /// constructor which uses memory that is allocated elsewhere
        /// </summary>
        /// <param name="P"></param>
        /// <param name="content">
        /// used to initialize <see cref="Storage"/>
        /// </param>
        public RefVector(IPartitioning P, double[] content)
            : base(P) {
            if (P.LocalLength > content.Length)
                throw new ArgumentException("vector content must match local length of partition", "content");
            m_Storage = content;
        }

        double[] m_Storage;

        /// <summary>
        /// internal storage of the vector;
        /// </summary>
        /// <remarks>
        /// Attention: the length of this array may be grater than the local lenght of this vector;
        /// Only entries in the range 0 to local length - 1 are used. the local length
        /// can be determided by the partition <see cref="VectorBase.Part"/> of this vector.
        /// </remarks>
        public double[] Storage {
            get { return m_Storage; }
        }

        /// <summary>
        /// <see cref="VectorBase.SetValues"/>;
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vals"></param>
        /// <param name="arrayIndex"></param>
        /// <param name="insertAt"></param>
        /// <param name="Length"></param>
        public override void SetValues<T>(T vals, int arrayIndex, int insertAt, int Length) {
            if (typeof(T) == typeof(double[])) {
                // optimized version
                Array.Copy(vals as double[], arrayIndex, m_Storage, insertAt, Length);
            } else {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    m_Storage[i + insertAt] = vals[i + arrayIndex];
            }
        }

        /// <summary>
        /// see <see cref="VectorBase.GetValues"/>;
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="vals"></param>
        /// <param name="arrayIndex"></param>
        /// <param name="readAt"></param>
        /// <param name="Length"></param>
        public override void GetValues<T>(T vals, int arrayIndex, int readAt, int Length) {
            if (typeof(T) == typeof(double[])) {
                // optimized version
                Array.Copy(m_Storage, readAt, vals as double[], arrayIndex, Length);
            } else {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    vals[i + arrayIndex] = m_Storage[i + readAt];
            }
        }

        /// <summary>
        /// gets/sets one entry of the vector
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public override double this[int index] {
            get {
                if (m_IsLocked)
                    throw new ApplicationException("object is locked.");
                return m_Storage[index];
            }
            set {
                if (m_IsLocked)
                    throw new ApplicationException("object is locked.");
                m_Storage[index] = value;
            }
        }

        /// <summary>
        /// sets all entries to 0.0;
        /// </summary>
        public override void Clear() {
            Array.Clear(m_Storage, 0, m_Storage.Length);
        }

        /// <summary>
        /// copies those local indices (indicated by <paramref name="ComList"/>), which must be 
        /// send to another processor, to a send buffer <paramref name="Buffer"/>;
        /// </summary>
        /// <param name="ComList">
        /// an ascending list of entry indices that must be copied to <see cref="Buffer"/>
        /// </param>
        /// <param name="Buffer">
        /// output;
        /// </param>
        /// <remarks>
        /// this method works only in locked mode <see cref="LockAbleObject.Lock"/>;
        /// </remarks>
        private void FillSendBuffer(int[] ComList, IntPtr Buffer) {
            if (m_IsLocked == false)
                throw new ApplicationException("works only in locked mode");
            int L = ComList.Length;
            unsafe {
                double* pBuffer = (double*)Buffer;
                for (int i = 0; i < L; i++)
                    pBuffer[i] = m_Storage[ComList[i]];
            }
        }

        /// <summary>
        /// multiplies all entries by number <paramref name="alpha"/>
        /// </summary>
        public override void Scale(double alpha) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");

            //BLAS.dscal(this.m_Part.LocalLength, alpha, m_Storage, 1);
            int N = this.m_Part.LocalLength;
            unsafe {
                double* pa = this.m_StorageAddr;

                for (int i = 0; i < N; i++)
                    pa[i] *= alpha;
            }

        }

        /// <summary>
        /// this = this + <paramref name="alpha"/>*<paramref name="other"/>;
        /// </summary>
        /// <param name="alpha"></param>
        /// <param name="other"></param>
        public override void Acc(double alpha, VectorBase other) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!other.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            int N = this.m_Part.LocalLength;
            if (other.Part.LocalLength != N)
                throw new ArgumentException("mismatch in vector size.");

            //BLAS.daxpy(N, alpha, (other as RefVector).Storage, 1, this.m_Storage, 1);
            unsafe {
                double* p_tis = this.m_StorageAddr, p_oda = (other as RefVector).m_StorageAddr;
                for (int i = 0; i < N; i++) {
                    p_tis[i] += alpha * p_oda[i];
                }
            }

        }

        //public override void daxpy(VectorBase other, double beta) {
        //    if (!m_IsLocked)
        //        throw new ApplicationException("works only in locked mode.");
        //    if (!other.IsLocked)
        //        throw new ArgumentException("other object must be locked.", "other");
        //    int N = this.m_Part.LocalLength;
        //    if (other.Part.LocalLength != N)
        //        throw new ArgumentException("mismatch in vector size.");

        //    unsafe {
        //        fixed (double* _a = &this.m_Storage[0],
        //                _b = &(other as RefVector).m_Storage[0]) {

        //            //BLAS.daxpy(N, alpha, (other as RefVector).Storage, beta, this.m_Storage, 1);
        //            for (int k = 0; k < N; k++) {
        //                _a[k] *= beta;
        //                _a[k] += _b[k];
        //            }
        //        }
        //    }
        //}


        /// <summary>
        /// the square of the two-norm of this vector;
        /// </summary>
        /// <returns></returns>
        public override double TwoNormSquare() {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            int N = this.m_Part.LocalLength;

            //double NormLoc = BLAS.dnrm2(N, m_Storage, 1);
            //double NormLocPow2 = NormLoc * NormLoc;
            double NormLocPow2 = 0;
            unsafe {
                double* stor = this.m_StorageAddr;
                for (int n = 0; n < N; n++) {
                    double r = stor[n];
                    NormLocPow2 += r * r;
                }
            }


            double NormglobPow2 = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&NormLocPow2), (IntPtr)(&NormglobPow2), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }
            return NormglobPow2;
        }

        /// <summary>
        /// computes the inner product of this timesanother vector;
        /// </summary>
        public override double InnerProd(VectorBase other) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!other.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            int N = this.m_Part.LocalLength;
            if (other.Part.LocalLength != N)
                throw new ArgumentException("mismatch in vector size.");

            //double ddotLoc = BLAS.ddot(N, this.m_Storage, 1, (other as RefVector).m_Storage, 1);

            double ddotLoc = 0;
            unsafe {
                //fixed (double* _a = &m_Storage[0], _b = &(other as RefVector).m_Storage[0]) {
                double* _a = this.m_StorageAddr, _b = (other as RefVector).m_StorageAddr;
                for (int n = 0; n < N; n++) {
                    ddotLoc += _a[n] * _b[n];
                    //    }
                }
            }

            double ddotGlob = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&ddotLoc), (IntPtr) (&ddotGlob), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }
            return ddotGlob;
        }

        /// <summary>
        /// copies the content of <paramref name="other"/> to this vector.
        /// </summary>
        public override void CopyFrom(VectorBase other) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!other.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            int N = this.m_Part.LocalLength;
            if (other.Part.LocalLength != N)
                throw new ArgumentException("mismatch in vector size.");

            Array.Copy((other as RefVector).m_Storage, this.m_Storage, N);
        }

        /// <summary>
        /// exchanges the content of <paramref name="other"/> and this.
        /// </summary>
        public override void Swap(VectorBase other) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!other.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            int N = this.m_Part.LocalLength;
            if (other.Part.LocalLength != N)
                throw new ArgumentException("mismatch in vector size.");

            //BLAS.dswap(N, (other as RefVector).m_Storage, 1, m_Storage, 1);
            unsafe {
                //fixed (double* pa = &(other as RefVector).m_Storage[0], pb = &this.m_Storage[0]) {
                double* pa = (other as RefVector).m_StorageAddr, pb = this.m_StorageAddr;
                double tmp;
                for (int i = 0; i < N; i++) {
                    tmp = pa[i];
                    pa[i] = pb[i];
                    pb[i] = tmp;
                }

            }
        }

        /// <summary>
        /// see <see cref="VectorBase.CopyPartlyFrom"/>
        /// </summary>
        public override void CopyPartlyFrom(VectorBase _src, int[] IdxThis, int PerThis, int[] IdxSrc, int PerSrc) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!_src.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            RefVector src = _src as RefVector;
            if (src == null)
                throw new ArgumentException("type mismatch");

            if (PerThis <= 0)
                throw new ArgumentOutOfRangeException("PerThis", "must be greater than 0.");
            if (PerSrc <= 0)
                throw new ArgumentOutOfRangeException("PerSrc", "must be greater than 0.");


            int Nthis = this.m_Part.LocalLength;
            int Nsrc = src.m_Part.LocalLength;

            if (IdxThis.Length != IdxSrc.Length)
                throw new ArgumentException("length of index arrays must be equal.");
            int L = IdxThis.Length;

            if (Nthis % PerThis != 0)
                throw new ArgumentException("");
            if (Nsrc % PerSrc != 0)
                throw new ArgumentException("");
            
            int NPer = Nthis / PerThis;
            if (NPer != Nsrc / PerSrc)
                throw new ArgumentException("");


            for (int n = 0; n < NPer; n++) {
                for (int l = 0; l < L; l++) {
                    this.m_Storage[n * PerThis + IdxThis[l]] = src.m_Storage[n * PerSrc + IdxSrc[l]];
                }
            }
        }


        GCHandle m_StoragePin;

        unsafe double* m_StorageAddr;

        internal unsafe double* StorageAddr {
            get { return m_StorageAddr; }
        }

        /// <summary>
        /// locks the allocated storrage for the garbage collector
        /// </summary>
        public override void Lock() {
            base.Lock();
            m_StoragePin = GCHandle.Alloc(m_Storage, GCHandleType.Pinned);
            unsafe {
                m_StorageAddr = (double*)Marshal.UnsafeAddrOfPinnedArrayElement(m_Storage, 0);
            }
        }

        /// <summary>
        /// realesed the garbage collector pind handled
        /// </summary>
        public override void Unlock() {
            base.Unlock();
            m_StoragePin.Free();
            unsafe {
                m_StorageAddr = (double*)IntPtr.Zero;
            }
        }
        /// <summary>
        /// For each <em>j</em>, <br/>
        /// this[j] = this[j]*<paramref name="other"/>[j]
        /// </summary>
        /// <param name="other"></param>
        public override void MultiplyElementWise(VectorBase other) {
            if (!m_IsLocked)
                throw new ApplicationException("works only in locked mode.");
            if (!other.IsLocked)
                throw new ArgumentException("other object must be locked.", "other");
            int N = this.m_Part.LocalLength;
            if (other.Part.LocalLength != N)
                throw new ArgumentException("mismatch in vector size.");

            RefVector o = (other as RefVector);

            for (int i = 0; i < N; i++)
                this.m_Storage[i] *= o.m_Storage[i];

        }

        /// <summary>
        /// see <see cref="RefVector.CreateCommVector"/>
        /// </summary>
        public override CommVector CreateCommVector(MatrixBase M) {
            return new RefCommVector(M, this);
        }

        //abstract protected void FillSendBuffer(int[] ComList, double[] Buffer);


        class RefCommVector : CommVector {
            internal RefCommVector(MatrixBase M, RefVector v)
                : base(M, v) {

                m_Owner = v;

                // allocate send buffers
                foreach (int proc in Mtx._SpmvCommPattern.ComLists.Keys) {
                    int len = Mtx._SpmvCommPattern.ComLists[proc].Length;
                    IntPtr pBuffer = Marshal.AllocHGlobal(len * sizeof(double));
                    SendBuffers.Add(proc, pBuffer);
                    SendBuffersLengths.Add(proc, len);
                }
            }

            RefVector m_Owner;

            public override void FillSendBuffer() {
                base.FillSendBuffer();

                // fill the send buffers
                foreach (int targetProc in SendBuffers.Keys) {
                    int[] CommList = Mtx._SpmvCommPattern.ComLists[targetProc];
                    IntPtr SendBuffer = SendBuffers[targetProc];
                    m_Owner.FillSendBuffer(CommList, SendBuffer);

                    //if (Enviroment.MPIEnv.MPI_Rank == 0) {
                    //    //Debugger.Break();
                    //    Console.WriteLine("RefVec.Fill");
                    //    unsafe {
                    //        double* ptr = (double*)SendBuffer;
                    //        System.Console.Error.WriteLine("size:" + CommList.Length);
                    //        for (int i = 0; i < CommList.Length && counter < 1000; i++) {
                    //            System.Console.Error.WriteLine(ptr[i]);
                    //            counter++;
                    //        }
                    //    }
                    //    Console.WriteLine("RefVec.Fill.End");
                    //}
                }


            }

            public override void Dispose() {
                base.Dispose();
                if (base.SendBuffers != null) {
                    foreach (IntPtr pBuffer in base.SendBuffers.Values)
                        Marshal.FreeHGlobal(pBuffer);
                    base.SendBuffers = null;
                }
            }
        }
    }
}
