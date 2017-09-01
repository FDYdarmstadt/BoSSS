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
    class clELLPACKmodMatrix : clMatrix
    {
        private cl_mem d_val;
        private cl_mem d_colIdx;

        private ELLPACKmod m_internalData;

        private int colCount;
        private int valStride;
        private int colStride;
        private int size;

        public clELLPACKmodMatrix(MsrMatrix M, clDevice device)
            : base(M, device, "ellMultiply")
        {
            m_internalData = (ELLPACKmod)m_LocalMtx;

            size = m_internalData.NoOfRows;
            colCount = m_internalData.NoOfPackedCols;
            valStride = m_internalData.MtxEntries.ColStride;
            colStride = m_internalData.ColInd.ColStride;

            localsize = 256;
            globalsize = size;
            int m = size % localsize;
            if (m > 0)
            {
                globalsize += localsize - m;
            }
        }

        /// <summary>
        /// returns a <see cref="MatrixBase.CSR"/>-object.
        /// </summary>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp)
        {
            return new MatrixBase.ELLPACKmod(tmp, 1, localsize);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            d_val = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.MtxEntries.Values.Length * sizeof(double), m_internalData.MtxEntries.Values);
            d_colIdx = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.ColInd.Values.Length * sizeof(int), m_internalData.ColInd.Values);
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cl.ReleaseMemObject(d_val);
            cl.ReleaseMemObject(d_colIdx);
        }

        internal override void SetArguments(double alpha, clVector a, double beta, clVector acc)
        {
            cl_mem d_x = a.GetDevicePointer();
            cl_mem d_result = acc.GetDevicePointer();

            cl.SetKernelArg(clmultiply, 0, d_val);
            cl.SetKernelArg(clmultiply, 1, d_colIdx);
            cl.SetKernelArg(clmultiply, 2, d_x);
            cl.SetKernelArg(clmultiply, 3, d_result);
            cl.SetKernelArg(clmultiply, 4, alpha);
            cl.SetKernelArg(clmultiply, 5, beta);
            cl.SetKernelArg(clmultiply, 6, size);
            cl.SetKernelArg(clmultiply, 7, colCount);
            cl.SetKernelArg(clmultiply, 8, valStride);
            cl.SetKernelArg(clmultiply, 9, colStride);
        }
    }
}
