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
using System.Linq;
using System.Text;

namespace ilPSP.LinSolvers.monkey.CL
{
    class clELLPACKcacheMatrix : clMatrix
    {
        private cl_mem d_val;
        private cl_mem d_colIdx;
        private cl_mem d_xSubStart;
        private cl_mem d_blockSubVector;

        private ManualCacheELLPACK m_internalData;

        private int maxBlockSubVectorLength = -1;
        private int colCount;
        private int valStride;
        private int colStride;
        private int size;

        public clELLPACKcacheMatrix(MsrMatrix M, clDevice device)
            : base(M, device, "mcellMultiply")
        {
            m_internalData = (ManualCacheELLPACK)m_LocalMtx;

            size = m_internalData.NoOfRows;
            colCount = m_internalData.NoOfPackedCols;
            valStride = m_internalData.MtxEntries.ColStride;
            colStride = m_internalData.ColIndBlock.ColStride;

            localsize = 256;
            globalsize = size;
            int m = size % localsize;
            if (m > 0)
            {
                globalsize += localsize - m;
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

            d_xSubStart = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)(h_xSubStart.Length + 1) * sizeof(int), h_xSubStart);
            d_blockSubVector = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)h_blockSubVector.Length * sizeof(int), h_blockSubVector);
        }

        /// <summary>
        /// returns a <see cref="MatrixBase.CSR"/>-object.
        /// </summary>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp)
        {
            return new MatrixBase.ManualCacheELLPACK(tmp, 1, localsize, 32);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            CreateDataStructures();

            d_val = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double), m_internalData.MtxEntries.Values);
            d_colIdx = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.ColIndBlock.Values.Length * sizeof(ushort), m_internalData.ColIndBlock.Values);
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cl.ReleaseMemObject(d_val);
            cl.ReleaseMemObject(d_colIdx);
            cl.ReleaseMemObject(d_blockSubVector);
            cl.ReleaseMemObject(d_xSubStart);
        }

        internal override void SetArguments(double alpha, clVector a, double beta, clVector acc)
        {
            cl_mem d_x = a.GetDevicePointer();
            cl_mem d_result = acc.GetDevicePointer();

            cl.SetKernelArg(clmultiply, 0, d_val);
            cl.SetKernelArg(clmultiply, 1, d_colIdx);
            cl.SetKernelArg(clmultiply, 2, d_xSubStart);
            cl.SetKernelArg(clmultiply, 3, d_blockSubVector);
            cl.SetKernelArg(clmultiply, 4, d_x);
            cl.SetKernelArg(clmultiply, 5, d_result);
            cl.SetKernelArgLocalSize(clmultiply, 6, (uint)maxBlockSubVectorLength * sizeof(double));
            cl.SetKernelArg(clmultiply, 7, alpha);
            cl.SetKernelArg(clmultiply, 8, beta);
            cl.SetKernelArg(clmultiply, 9, size);
            cl.SetKernelArg(clmultiply, 10, colCount);
            cl.SetKernelArg(clmultiply, 11, valStride);
            cl.SetKernelArg(clmultiply, 12, colStride);
        }
    }
}
