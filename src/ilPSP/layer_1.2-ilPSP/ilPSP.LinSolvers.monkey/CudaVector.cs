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

namespace ilPSP.LinSolvers.monkey.CUDA {

    class CudaVector : VectorBase {

        CUfunction cuscale;
        CUfunction cuacc;
        CUfunction cudnrm2;
        CUfunction cuinnerprod;
        CUfunction cumew;

        int blocksize = 256;
        int blockcountfull;
        int blockcounthalf;
        int size;

        double[] h_data;
        IntPtr h_result;

        CUdeviceptr d_data;
        CUdeviceptr d_result;

        /// <summary>
        /// constructor that allocates its own memory
        /// </summary>
        /// <param name="p">vector partition among MPI processes</param>
        /// <param name="env"></param>
        public CudaVector(IPartitioning p, CudaEnviroment env)
            : base(p) {
            h_data = new double[p.LocalLength];
            ConstructorCommon(env);
        }

        /// <summary>
        /// constructor which uses memory that is allocated elsewhere
        /// </summary>
        /// <param name="P"></param>
        /// <param name="content">
        /// used to initialize <see cref="h_data"/>
        /// </param>
        /// <param name="env"></param>
        public CudaVector(IPartitioning P, double[] content, CudaEnviroment env) : base(P) {
            if(P.LocalLength > content.Length)
                throw new ArgumentException("vector content must match local length of partition","content");
            h_data = content;
            ConstructorCommon(env);
        }

        CudaEnviroment m_env;
        
        private void ConstructorCommon(CudaEnviroment env) {
            m_env = env;
            cuscale = m_env.Get_CudaVectorKernelDP_Function("scale");
            cuacc = m_env.Get_CudaVectorKernelDP_Function("acc");
            cudnrm2 = m_env.Get_CudaVectorKernelDP_Function("dnrm2");
            cuinnerprod = m_env.Get_CudaVectorKernelDP_Function("innerprod");
            cumew = m_env.Get_CudaVectorKernelDP_Function("mew");

            size = this.Part.LocalLength;
            blockcountfull = (int)Math.Ceiling((decimal)size / blocksize);
            blockcounthalf = (int)Math.Ceiling((decimal)size / (2 * blocksize));
        }

        public CUdeviceptr GetDevicePointer() {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            return d_data;
        }

        public override void SetValues<T>(T vals, int arrayIndex, int insertAt, int Length) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked.");

            if (typeof(T) == typeof(double[])) {
                // optimized version
                Array.Copy(vals as double[], arrayIndex, h_data, insertAt, Length);
            }
            else {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    h_data[i + insertAt] = vals[i + arrayIndex];
            }
        }

        public override void GetValues<T>(T vals, int arrayIndex, int readAt, int Length) {
            if (this.IsLocked)
                throw new ApplicationException("object is locked.");

            if (typeof(T) == typeof(double[])) {
                // optimized version
                Array.Copy(h_data, readAt, vals as double[], arrayIndex, Length);
            }
            else {
                // version which works whith all kinds of IList<double>
                for (int i = 0; i < Length; i++)
                    vals[i + arrayIndex] = h_data[i + readAt];
            }
        }

        public override double this[int index] {
            get {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked.");

                return h_data[index];
            }
            set {
                if (this.IsLocked)
                    throw new ApplicationException("object is locked.");

                h_data[index] = value;
            }
        }

        public override void Clear() {
            if (!this.IsLocked) {
                Array.Clear(h_data, 0, h_data.Length);
            } else {
                cu.MemsetD8(d_data, 0, (uint)(h_data.Length * sizeof(double)));

            }
        }

        public override void Scale(double alpha) {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            int offset = 0;

            cu.ParamSetp(cuscale, offset, d_data);
            offset += sizeof(long);
            cu.ParamSetd(cuscale, offset, alpha);
            offset += sizeof(double); ;
            cu.ParamSeti(cuscale, offset, size);
            offset += sizeof(uint);

            cu.ParamSetSize(cuscale, (uint)offset);
            cu.FuncSetBlockShape(cuscale, blocksize, 1, 1);

            cu.LaunchGrid(cuscale, blockcountfull, 1);
        }

        public override void Acc(double alpha, VectorBase other) {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            CudaVector _other = other as CudaVector;
            if(_other == null)
                throw new ArgumentException("other must be of type CudaVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            int offset = 0;

            cu.ParamSetp(cuacc, offset, d_data);
            offset += sizeof(long);
            cu.ParamSetp(cuacc, offset, _other.GetDevicePointer());
            offset += sizeof(long);
            cu.ParamSetd(cuacc, offset, alpha);
            offset += sizeof(double);
            cu.ParamSeti(cuacc, offset, size);
            offset += sizeof(uint);

            cu.ParamSetSize(cuacc, (uint)offset);
            cu.FuncSetBlockShape(cuacc, blocksize, 1, 1);

            cu.LaunchGrid(cuacc, blockcountfull, 1);
        }

        public override double TwoNormSquare() {
            if (!this.IsLocked)
                throw new ApplicationException("works only in locked mode");

            int offset = 0;
            double finalResult = 0.0;
            
            cu.ParamSetp(cudnrm2, offset, d_data);
            offset += sizeof(long);
            cu.ParamSetp(cudnrm2, offset, d_result);
            offset += sizeof(long);
            cu.ParamSeti(cudnrm2, offset, size);
            offset += sizeof(uint);

            cu.ParamSetSize(cudnrm2, (uint)offset);
            cu.FuncSetBlockShape(cudnrm2, blocksize, 1, 1);
            cu.FuncSetSharedSize(cudnrm2, (uint)(blocksize * sizeof(double)));

            cu.LaunchGrid(cudnrm2, blockcounthalf, 1);
            cu.CtxSynchronize();

            unsafe {
                double* ptr = (double*)h_result;
                for (int i = 0; i < blockcounthalf; i++) {
                    finalResult += ptr[i];
                }
            }

            double NormglobPow2 = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr) (&finalResult), (IntPtr) (&NormglobPow2), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }
            return NormglobPow2;
        }

        public override double InnerProd(VectorBase other) {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            CudaVector _other = other as CudaVector;
            if (_other == null)
                throw new ArgumentException("other must be of type CudaVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            int offset = 0;
            double finalResult = 0.0;

            cu.ParamSetp(cuinnerprod, offset, d_data);
            offset += sizeof(long);
            cu.ParamSetp(cuinnerprod, offset, _other.GetDevicePointer());
            offset += sizeof(long);
            cu.ParamSetp(cuinnerprod, offset, d_result);
            offset += sizeof(long);
            cu.ParamSeti(cuinnerprod, offset, size);
            offset += sizeof(uint);

            cu.ParamSetSize(cuinnerprod, (uint)offset);
            cu.FuncSetBlockShape(cuinnerprod, blocksize, 1, 1);
            cu.FuncSetSharedSize(cuinnerprod, (uint)(blocksize * sizeof(double)));

            cu.LaunchGrid(cuinnerprod, blockcounthalf, 1);
            cu.CtxSynchronize();

            unsafe {
                double* ptr = (double*)h_result;
                for (int i = 0; i < blockcounthalf; i++) {
                    finalResult += ptr[i];
                }
            }

            double dotProdGlobal = double.NaN;
            unsafe {
                csMPI.Raw.Allreduce((IntPtr)(&finalResult), (IntPtr) (&dotProdGlobal), 1, csMPI.Raw._DATATYPE.DOUBLE, csMPI.Raw._OP.SUM, csMPI.Raw._COMM.WORLD);
            }
            return dotProdGlobal;
        }

        public override void CopyFrom(VectorBase other) {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            CudaVector _other = other as CudaVector;
            if(_other == null)
                throw new ArgumentException("other must be of type CudaVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            cu.MemcpyDtoD(d_data, _other.GetDevicePointer(), (uint)(size * sizeof(double)));
        }

        public override void Swap(VectorBase other) {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            CudaVector _other = other as CudaVector;
            if(_other == null)
                throw new ArgumentException("other must be of type CudaVector.", "other");

            if(_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            CUdeviceptr temp = _other.d_data;
            _other.d_data = this.d_data;
            this.d_data = temp;
        }

        /// <summary>
        /// For each <em>j</em>, <br/>
        /// this[j] = this[j]*<paramref name="other"/>[j]
        /// </summary>
        /// <param name="other"></param>
        public override void MultiplyElementWise(VectorBase other) {
            if (!this.IsLocked || !other.IsLocked)
                throw new ApplicationException("works only in locked mode");

            CudaVector _other = other as CudaVector;
            if (_other == null)
                throw new ArgumentException("other must be of type CudaVector.", "other");

            if (_other.Part.LocalLength != this.Part.LocalLength)
                throw new ArgumentException("mismatch in vector size.");

            int offset = 0;

            cu.ParamSetp(cumew, offset, d_data);
            offset += sizeof(long);
            cu.ParamSetp(cumew, offset, _other.GetDevicePointer());
            offset += sizeof(long);
            cu.ParamSeti(cumew, offset, size);
            offset += sizeof(uint);

            cu.ParamSetSize(cumew, (uint)offset);
            cu.FuncSetBlockShape(cumew, blocksize, 1, 1);

            cu.LaunchGrid(cumew, blockcountfull, 1);
        }

        public override void Lock() {
            base.Lock();

            cu.MemAlloc(out d_data, (uint)size * sizeof(double));
            cu.MemcpyHtoD(d_data, h_data, (uint)size * sizeof(double));

            cu.MemHostAlloc(out h_result, (uint)blockcounthalf * sizeof(double), CUmem_host_alloc.CU_MEMHOSTALLOC_DEVICEMAP);
            cu.MemHostGetDevicePointer(out d_result, h_result, 0);
        }

        public override void Unlock() {
            base.Unlock();

            cu.MemcpyDtoH(h_data, d_data, (uint)size * sizeof(double));
            cu.MemFree(d_data);

            cu.MemFreeHost(h_result);
        }


        public override CommVector CreateCommVector(MatrixBase M) {
            CudaMatrix cM = (CudaMatrix)M;
            return new CudaCommVector(M, this, cM.extStream);
        }


        class CudaCommVector : CommVector {

            //CUmodule module;
            CUfunction cufill;
            CUstream stream;

            int size;
            int blocksize = 128;
            int blockcount;

            /// <summary>
            /// starting address for the send buffer - page locked memory
            /// </summary>
            IntPtr h_SendBuffer;

            /// <summary>
            /// Device memory pointer for <see cref="h_SendBuffer"/>
            /// </summary>
            CUdeviceptr d_SendBuffer;

            /// <summary>
            /// for the SpMv y=beta*y + alpha*M*x, the indices of those elements of x which must 
            /// be send to other processors.
            /// </summary>
            int[] h_IndicesToSend;

            /// <summary>
            /// Device memory pointer for <see cref="h_IndicesToSend"/>
            /// </summary>
            CUdeviceptr d_IndicesToSend;

            /// <summary>
            /// CudaVector that this CommVector belongs to
            /// </summary>
            CudaVector owner;

            internal CudaCommVector(MatrixBase M, CudaVector v, CUstream stream)
                : base(M, v) {


                this.owner = v;
                this.stream = stream;
                cufill = owner.m_env.Get_CudaVectorKernelDP_Function("fillSendBuffer");

                IDictionary<int, int[]> comLists = M._SpmvCommPattern.ComLists;
                //int[] procranks = new int[comLists.Count]; // put all proccessor ranks in one list to have a unique ordering
                
                int totLen = 0;
                foreach(int procRnk in comLists.Keys) {
                    int l = comLists[procRnk].Length;
                    base.SendBuffersLengths[procRnk] = l;
                    totLen += l;
                }

                size = totLen;
                blockcount = (int)Math.Ceiling((decimal)size / blocksize);
                if (size > 0) {
                    // alloc
                    h_IndicesToSend = new int[size];
                    cu.MemAlloc(out d_IndicesToSend, (uint)size * sizeof(int));

                    cu.MemHostAlloc(out h_SendBuffer, sizeof(double) * (uint)size, CUmem_host_alloc.CU_MEMHOSTALLOC_DEVICEMAP);
                    cu.MemHostGetDevicePointer(out d_SendBuffer, h_SendBuffer, 0);

                    // concat lists: 
                    int i0 = 0;
                    unsafe {
                        double* P0 = (double*)h_SendBuffer;

                        foreach (int procRnk in comLists.Keys) {
                            base.SendBuffers[procRnk] = (IntPtr)P0;  // startaddres for sending to process 'procRnk'

                            int l = base.SendBuffersLengths[procRnk];
                            P0 += l;
                            Array.Copy(comLists[procRnk], 0, h_IndicesToSend, i0, l); // concat comm list
                            i0 += l;
                        }
                    }

                    cu.MemcpyHtoD(d_IndicesToSend, h_IndicesToSend, (uint)size * sizeof(int));
                }
            }

            public override void FillSendBuffer() {
 	            if(!owner.IsLocked)
                    throw new ApplicationException("works only in locked mode");

                base.FillSendBuffer();

                if (size > 0) {
                    int offset = 0;

                    cu.ParamSetp(cufill, offset, d_SendBuffer);
                    offset += sizeof(long);
                    cu.ParamSetp(cufill, offset, d_IndicesToSend);
                    offset += sizeof(long);
                    cu.ParamSetp(cufill, offset, owner.GetDevicePointer());
                    offset += sizeof(long);
                    cu.ParamSeti(cufill, offset, size);
                    offset += sizeof(uint);

                    cu.ParamSetSize(cufill, (uint)offset);
                    cu.FuncSetBlockShape(cufill, blocksize, 1, 1);
                    //{
                    //    int major, minor;
                    //    cu.DeviceComputeCapability(out major, out minor, this.m_Cu m_CUDAdev);
                    //    if (major >= 2)
                    //        cu.FuncSetCacheConfig(cufill, CUfunc_cache.CU_FUNC_CACHE_PREFER_L1);
                    //}
                    cu.LaunchGridAsync(cufill, blockcount, 1, stream);
                    cu.StreamSynchronize(stream);
                }
            }

            public override void Dispose() {
                base.Dispose();

                if (size > 0) {
                    cu.MemFree(d_IndicesToSend);
                    if (h_SendBuffer != IntPtr.Zero) {
                        cu.MemFreeHost(h_SendBuffer);
                        h_SendBuffer = IntPtr.Zero;
                    }
                }
            }
        }
        
        /// <summary>
        /// see <see cref="VectorBase.CopyPartlyFrom"/>
        /// </summary>
        public override void CopyPartlyFrom(VectorBase _src, int[] IdxThis, int PerThis, int[] IdxSrc, int PerSrc) {
            throw new NotImplementedException();
        }
    }
}
