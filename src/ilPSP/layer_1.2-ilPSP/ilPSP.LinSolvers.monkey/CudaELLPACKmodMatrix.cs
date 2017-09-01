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

namespace ilPSP.LinSolvers.monkey.CUDA
{
    /// <summary>
    /// Implements matrix-vector multiplication for matrix in ELLPACKmod format.
    /// </summary>
    public class CudaELLPACKmodMatrix : CudaMatrix
    {
        private CUdeviceptr d_val;
        private CUdeviceptr d_colIdx;

        private ELLPACKmod m_internalData;

        private int blocksize = 256;
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
        public CudaELLPACKmodMatrix(MsrMatrix M, CudaEnviroment CudaEnv)
            : base(M, "ellMultiply", CudaEnv)
        {
            m_internalData = (ELLPACKmod)m_LocalMtx;
            size = m_internalData.NoOfRows;
            colCount = m_internalData.NoOfPackedCols;
            valStride = m_internalData.MtxEntries.ColStride;
            colStride = m_internalData.ColInd.ColStride;
            blockcount = (int)Math.Ceiling((Decimal)size / blocksize);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="tmp"></param>
        /// <returns></returns>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp)
        {
            return new ELLPACKmod(tmp, 1, blocksize);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            cu.MemAlloc(out d_val, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double));
            cu.MemAlloc(out d_colIdx, (uint)m_internalData.ColInd.Values.Length * sizeof(int));
            cu.MemcpyHtoD(d_val, m_internalData.MtxEntries.Values, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double));
            cu.MemcpyHtoD(d_colIdx, m_internalData.ColInd.Values, (uint)m_internalData.ColInd.Values.Length * sizeof(int));
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cu.MemFree(d_val);
            cu.MemFree(d_colIdx);
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

            cu.LaunchGridAsync(sparseMultiply, blockcount, 1, stream);
        }
    }
}
