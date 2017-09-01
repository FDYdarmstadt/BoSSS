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
using System.Linq;
using ilPSP.Tracing;

namespace ilPSP.LinSolvers.monkey.CUDA
{
    /// <summary>
    /// Version of ELLPACK with shared memory usage
    /// </summary>
    public class CudaELLPACKcacheMatrix : CudaMatrix
    {
        private CUdeviceptr d_val;
        private CUdeviceptr d_colIdx;
        private CUdeviceptr d_xSubStart;
        private CUdeviceptr d_blockSubVector;

        private ManualCacheELLPACK m_internalData;

        private int blocksize = 256;
        private int maxBlockSubVectorLength = -1;
        private int size;
        private int colCount;
        private int valStride;
        private int colStride;
        private int blockcount;


        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="M">Sparse matrix in MSR format</param>
        /// <param name="CudaEnv"></param>
        public CudaELLPACKcacheMatrix(MsrMatrix M, CudaEnviroment CudaEnv)
            : base(M, "mcellMultiply", CudaEnv)
        {
            using (new FuncTrace()) {
                m_internalData = (ManualCacheELLPACK)m_LocalMtx;
                size = m_internalData.NoOfRows;
                colCount = m_internalData.NoOfPackedCols;
                valStride = m_internalData.MtxEntries.ColStride;
                colStride = m_internalData.ColIndBlock.ColStride;
                blockcount = (int)Math.Ceiling((Decimal)size / blocksize);
            }
        }

        private void CreateDataStructures()
        {
            int[] h_xSubStart = new int[m_internalData.BlockSubVectorLength.Length + 1];
            int len = 0;
            for (int i = 0; i < m_internalData.BlockSubVectorLength.Length; i++)
            {
                h_xSubStart[i] = len;
                len += m_internalData.BlockSubVectorLength[i];

                if (m_internalData.BlockSubVectorLength[i] > maxBlockSubVectorLength)
                {
                    maxBlockSubVectorLength = m_internalData.BlockSubVectorLength[i];
                }
            }
            h_xSubStart[m_internalData.BlockSubVectorLength.Length] = len;

            int[] h_blockSubVector = new int[len];
            int idx = 0;
            for (int i = 0; i < m_internalData.BlockSubVectorLength.Length; i++)
            {
                for (int j = 0; j < m_internalData.BlockSubVectorLength[i]; j++)
                {
                    h_blockSubVector[idx] = m_internalData.BlockSubVector[i, j];
                    idx++;
                }
            }

            cu.MemAlloc(out d_xSubStart, (uint)(h_xSubStart.Length + 1) * sizeof(int));
            cu.MemAlloc(out d_blockSubVector, (uint)h_blockSubVector.Length * sizeof(int));
            cu.MemcpyHtoD(d_xSubStart, h_xSubStart, (uint)(h_xSubStart.Length + 1) * sizeof(int));
            cu.MemcpyHtoD(d_blockSubVector, h_blockSubVector, (uint)h_blockSubVector.Length * sizeof(int));
        }

        /// <summary>
        /// see <see cref="MatrixBase.AssembleFinalFormat"/>
        /// </summary>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp) {
            using (new FuncTrace()) {
                ManualCacheELLPACK ret = null;
                bool TooBig = true;
                int maxSharesSize;
                cu.DeviceGetAttribute(out maxSharesSize, CUdevice_attribute.CU_DEVICE_ATTRIBUTE_MAX_SHARED_MEMORY_PER_BLOCK, m_CudaEnv.CUDAdev);

                while (TooBig) {
                    ret = new ManualCacheELLPACK(tmp, 1, blocksize, 32);
                    int ReqSharedSize = ret.BlockSubVectorLength.Max() * sizeof(double);

                    if (ReqSharedSize > maxSharesSize) {
                        blocksize /= 2;
                        if (blocksize < 16)
                            throw new ApplicationException("unable to create ELLPACKcache - matrix; not enough shared memory available;");
                    } else {
                        TooBig = false;
                    }

                }
                return ret;
            }
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            CreateDataStructures();

            cu.MemAlloc(out d_val, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double));
            cu.MemAlloc(out d_colIdx, (uint)m_internalData.ColIndBlock.Values.Length * sizeof(ushort));
            cu.MemcpyHtoD(d_val, m_internalData.MtxEntries.Values, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double));
            cu.MemcpyHtoD(d_colIdx, m_internalData.ColIndBlock.Values, (uint)m_internalData.ColIndBlock.Values.Length * sizeof(ushort));
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cu.MemFree(d_val);
            cu.MemFree(d_colIdx);
            cu.MemFree(d_blockSubVector);
            cu.MemFree(d_xSubStart);
        }

        internal override void CallDriver(CUstream stream, double alpha, CudaVector a, double beta, CudaVector acc)
        {
            CUdeviceptr d_x = a.GetDevicePointer();
            CUdeviceptr d_result = acc.GetDevicePointer();

            int offset = 0;

            cu.ParamSetp(sparseMultiply, offset, d_val);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_colIdx);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_xSubStart);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_blockSubVector);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_x);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_result);
            offset += sizeof(long);
            cu.ParamSetd(sparseMultiply, offset, alpha);
            offset += sizeof(double);
            cu.ParamSetd(sparseMultiply, offset, beta);
            offset += sizeof(double);
            cu.ParamSeti(sparseMultiply, offset, size);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, colCount);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, valStride);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, colStride);
            offset += sizeof(uint);

            cu.ParamSetSize(sparseMultiply, (uint)offset);
            cu.FuncSetBlockShape(sparseMultiply, blocksize, 1, 1);
            cu.FuncSetSharedSize(sparseMultiply, (uint)maxBlockSubVectorLength * sizeof(double));

            cu.LaunchGridAsync(sparseMultiply, blockcount, 1, stream);
        }
    }
}
