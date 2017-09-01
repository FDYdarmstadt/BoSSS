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
using MPI.Wrappers;
using ilPSP.Utils;
using System.Runtime.InteropServices;

namespace ilPSP.LinSolvers.monkey.CUDA {
    
    /// <summary>
    /// base class for all CUDA matrices
    /// </summary>
    abstract public class CudaMatrix : MatrixBase, IDisposable {

        /// <summary>
        /// CUDA environment
        /// </summary>
        protected CudaEnviroment m_CudaEnv;
        
        ///// <summary>
        ///// Module handle for the PTX file containing the CUDA kernels
        ///// </summary>
        //protected CUmodule module;

        /// <summary>
        /// CUDA function for accumulating internal and external part
        /// </summary>
        private CUfunction cuaccext;

        /// <summary>
        /// CUDA function for SpMV (sparse matrix - vector product)
        /// </summary>
        protected CUfunction sparseMultiply;

        /// <summary>
        /// used by <see cref="SpMV_External_RecvCallBack"/>;
        /// key: MPI processor rank; <br/>
        /// content: indices into <see cref="h_ElementsToAcc"/>;
        /// </summary>
        SortedDictionary<int, int[]> acc_ind = new SortedDictionary<int, int[]>();
        
        /// <summary>
        /// Number of external elements needed by this part
        /// </summary>
        private int extSize;

        /// <summary>
        /// <see cref="h_ElementsToAcc"/>
        /// </summary>
        private int[] h_IndicesToAccumulate;

        /// <summary>
        /// Device memory pointer for <see cref="h_IndicesToAccumulate"/>
        /// </summary>
        private CUdeviceptr d_IndicesToAccumulate;

        /// <summary>
        /// This are the elements of the 
        /// page-locked memory, allocated by 'cuMemAllocHost';<br/>
        /// Let be y = beta*y + alpha*M*x the SpMV, with scalars alpha and beta, vectors x and y and matrix M;<br/>
        /// Then, the 
        /// </summary>
        private IntPtr h_ElementsToAcc;

        /// <summary>
        /// Device memory pointer for <see cref="h_ElementsToAcc"/>
        /// </summary>
        private CUdeviceptr d_ElementsToAcc;

        /// <summary>
        /// Device memory pointer for memory of result vector
        /// </summary>
        private CUdeviceptr d_acc;

        private double m_alpha = double.NaN;

        private bool disposed = false;

        /// <summary>
        /// CUDA stream for asynchron processing of internal matrix on GPU
        /// </summary>
        internal CUstream intStream;

        /// <summary>
        /// and external MPI on CPU
        /// </summary>
        internal CUstream extStream;
        
        private void LMAA() {
            
            List<int> ins_ind = new List<int>();
            foreach (int proc in ExtMatrix.Keys) {
                External ext = ExtMatrix[proc];
                ins_ind.AddRange(ext.rowInd);
            }

            // remove duplicate entries
            ins_ind.Sort();
            if (ins_ind.Count > 1) {
                for (int i = 1; i < ins_ind.Count; i++) {
                    while (i < ins_ind.Count && ins_ind[i] == ins_ind[i - 1])
                        ins_ind.RemoveAt(i);
                }
            }

            //build new lists
            foreach (int proc in ExtMatrix.Keys) {
                External ext = ExtMatrix[proc];
                int[] accIndProc = new int[ext.rowInd.Length];
                for (int i = 0; i < accIndProc.Length; i++) {
                    accIndProc[i] = ins_ind.BinarySearch(ext.rowInd[i]);
                    if (accIndProc[i] < 0)
                        // should never happen
                        throw new ApplicationException("error in alg");
                }
                acc_ind.Add(proc, accIndProc);
            }

            h_IndicesToAccumulate = ins_ind.ToArray();
            extSize = h_IndicesToAccumulate.Length;
        }


        /// <summary>
        /// ctor
        /// </summary>
        public CudaMatrix(MsrMatrix M, string funcName, CudaEnviroment CudaEnv)
            : base(M) {

            m_CudaEnv = CudaEnv;
            base.PackMatrix(M);
            
            cu.StreamCreate(out intStream, 0);
            cu.StreamCreate(out extStream, 0);
            disposed = false;

            sparseMultiply = CudaEnv.Get_CudaMatrixKernelDP_Function(funcName);
            cuaccext = CudaEnv.Get_CudaMatrixKernelDP_Function("accumulateExternal");

            //int numreg;
            //cu.FuncGetAttribute(out numreg, CUfunction_attribute.CU_FUNC_ATTRIBUTE_NUM_REGS, sparseMultiply);
            //int version;
            //cu.FuncGetAttribute(out version, CUfunction_attribute.CU_FUNC_ATTRIBUTE_BINARY_VERSION, sparseMultiply);
            //System.Console.WriteLine("Number of registers: " + numreg + ", version: " + version);

            LMAA();

            if (extSize > 0) {
                // allocate page-locked mem
                cu.MemHostAlloc(out h_ElementsToAcc, sizeof(double) * (uint)extSize, CUmem_host_alloc.CU_MEMHOSTALLOC_DEVICEMAP);
                //test_ext = new double[totLen];
                cu.MemHostGetDevicePointer(out d_ElementsToAcc, h_ElementsToAcc, 0);

                cu.MemAlloc(out d_IndicesToAccumulate, (uint)extSize * sizeof(int));

                // Copy indices for combining external and internal part to GPU as they don't change over execution
                cu.MemcpyHtoD(d_IndicesToAccumulate, h_IndicesToAccumulate, (uint)extSize * sizeof(int));
            }
        }

        /// <summary>
        /// dtor
        /// </summary>
        ~CudaMatrix() {
            //Dispose();
        }

        /// <summary>
        /// disp
        /// </summary>
        public override void Dispose() {
            base.Dispose();

            if (disposed)
                return;

            cu.MemFreeHost(h_ElementsToAcc);
            h_ElementsToAcc = IntPtr.Zero;
            cu.MemFree(d_IndicesToAccumulate);

            cu.StreamDestroy(intStream);
            intStream = default(CUstream);
            cu.StreamDestroy(extStream);
            extStream = default(CUstream);

            disposed = true;
        }

        /// <summary>
        /// See <see cref="MatrixBase.CreateVec{T}(T, IPartitioning, out bool)"/>.
        /// </summary>
        protected override VectorBase CreateVec<T>(T a, IPartitioning len, out bool CopyIsShallow) {
            if (a.Count != len.LocalLength)
                throw new ArgumentException("count of a must be at least 'len'!", "len");
            if (len.MPI_Comm != this.ColPartition.MPI_Comm)
                throw new ArgumentException();


            IPartitioning part_a = len;
            double[] _a_stor = ArrayTools.List2Array<double>(a, 0, len.LocalLength);
            CopyIsShallow = false;
            return new CudaVector(part_a, _a_stor, this.m_CudaEnv);
        }

        internal override void SpMV_Local_Start(double alpha, VectorBase a, double beta, VectorBase acc) {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");

            CudaVector _a = a as CudaVector;
            CudaVector _acc = acc as CudaVector;

            if(_a == null)
                throw new ArgumentException("a must be of type CudaVector.", "a");
            if (_acc == null)
                throw new ArgumentException("acc must be of type CudaVector.", "acc");

            CallDriver(intStream, alpha, _a, beta, _acc);
        }

        abstract internal void CallDriver(CUstream stream, double alpha, CudaVector a, double beta, CudaVector acc);

        /// <summary>
        /// no operation
        /// </summary>
        internal override void SpMV_Local_Middle(double alpha, VectorBase a, double beta, VectorBase acc) {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");
        }
        
        /// <summary>
        /// blocks until the computation of the internal part is finished.
        /// </summary>
        internal override void SpMV_Local_End(double alpha, VectorBase a, double beta, VectorBase acc) {
            if (!m_IsLocked)
                throw new ApplicationException("object must be locked.");
            
            cu.StreamSynchronize(intStream);
        }

        internal override void SpMV_External_Begin(double alpha, double beta, VectorBase acc) {
            m_alpha = alpha;
            CudaVector _acc = (CudaVector)acc;
            d_acc = _acc.GetDevicePointer();

            unsafe {
                double* _acc_stor = (double*)h_ElementsToAcc;
                for (int i = (int)extSize - 1; i >= 0; i--) {
                    *_acc_stor = 0;
                    _acc_stor++;
                }
            }
        }

        /// <summary>
        /// multiplies the external vector (received from MPI-process <paramref name="procRank"/>)
        /// with the corresponding part of the matrix.
        /// This is always done on the CPU.
        /// The result is stored in the buffer <see cref="h_ElementsToAcc"/>
        /// </summary>
        internal override void SpMV_External_RecvCallBack(int procRank, IntPtr values) {
            External externalPart = ExtMatrix[procRank];
            
            int m_MyRank;
            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out m_MyRank);

            int[] accIndProc = acc_ind[procRank];

            unsafe {
                double* _acc_stor = (double*)h_ElementsToAcc;
                double* _values = (double*)values;
                fixed (double* _Val = &externalPart.Val[0]) {
                    fixed (int* _rowSt = &externalPart.RowStart[0],
                                _colInd = &externalPart.ColInd[0]) {

                        int K = externalPart.rowInd.Length;

                        int* Col = _colInd;
                        double* MtxEntry = _Val;

                        for (int k = 0; k < K; k++) {

                            int RowStart = _rowSt[k];
                            int RowEnd = _rowSt[k + 1];

                            double rowacc = 0;
                            for (int i = RowStart; i < RowEnd; i++) {

                                rowacc += (*MtxEntry) * _values[*Col];
                                Col++;
                                MtxEntry++;
                            }

                            _acc_stor[accIndProc[k]] += rowacc;
                        }
                    }
                }
            }
        }

        internal override void SpMV_External_Finalize() {
            if (extSize > 0) {
                int blocksize = 128;
                int blockcount = (int)Math.Ceiling((decimal)extSize / blocksize);
                int offset = 0;

                cu.ParamSetp(cuaccext, offset, d_acc);
                offset += sizeof(long);
                cu.ParamSetp(cuaccext, offset, d_IndicesToAccumulate);
                offset += sizeof(long);
                cu.ParamSetp(cuaccext, offset, d_ElementsToAcc);
                offset += sizeof(long);
                cu.ParamSetd(cuaccext, offset, m_alpha);
                offset += sizeof(double);
                cu.ParamSeti(cuaccext, offset, extSize);
                offset += sizeof(uint);

                cu.ParamSetSize(cuaccext, (uint)offset);
                cu.FuncSetBlockShape(cuaccext, blocksize, 1, 1);
                //{
                //    int major, minor;
                //    cu.DeviceComputeCapability(out major, out minor, this.m_CudaEnv.CUDAdev);
                //    if (major >= 2)
                //        cu.FuncSetCacheConfig(cuaccext, CUfunc_cache.CU_FUNC_CACHE_PREFER_L1);
                //}
                cu.LaunchGridAsync(cuaccext, blockcount, 1, extStream);
                cu.StreamSynchronize(extStream);
            }
        }
    }
}
