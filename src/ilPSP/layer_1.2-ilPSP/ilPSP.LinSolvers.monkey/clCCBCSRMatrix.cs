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
    class clCCBCSRMatrix : clMatrix
    {
        private cl_mem d_cellData;
        private cl_mem d_cellColIdx;

        private CCBCSR m_internalData;

        private int cellsize;
        private int cellrowsperblock;
        private int cellsperrow;
        private int stride;
        private int size;

        public clCCBCSRMatrix(MsrMatrix M, clDevice device)
            : base(M, device, "ccbcsrMultiply")
        {
            m_internalData = (CCBCSR)m_LocalMtx;

            size = base.RowPartitioning.LocalLength;
            cellsize = m_internalData.CellSize;
            // Number of cells per block, choose so that it is around 256 threads per block
            cellrowsperblock = (int)Math.Ceiling(128.0 / cellsize);
            cellsperrow = m_internalData.NoOfCellsPerRow;
            stride = m_internalData.CellStride;

            // Number of threads per block
            localsize = cellsize * cellrowsperblock;
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
            return new MatrixBase.CCBCSR(tmp, 1, base.CellSize);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            d_cellData = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.Val.Length * sizeof(double), m_internalData.Val);
            d_cellColIdx = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)m_internalData.CellColumn.Length * sizeof(double), m_internalData.CellColumn);
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cl.ReleaseMemObject(d_cellData);
            cl.ReleaseMemObject(d_cellColIdx);
        }

        internal override void SetArguments(double alpha, clVector a, double beta, clVector acc)
        {
            cl_mem d_x = a.GetDevicePointer();
            cl_mem d_result = acc.GetDevicePointer();

            cl.SetKernelArg(clmultiply, 0, d_cellData);
            cl.SetKernelArg(clmultiply, 1, d_x);
            cl.SetKernelArg(clmultiply, 2, d_cellColIdx);
            cl.SetKernelArg(clmultiply, 3, d_result);
            cl.SetKernelArgLocalSize(clmultiply, 4, (uint)(localsize * sizeof(double)));
            cl.SetKernelArgLocalSize(clmultiply, 5, (uint)(cellrowsperblock * sizeof(int)));
            cl.SetKernelArgLocalSize(clmultiply, 6, (uint)(cellrowsperblock * sizeof(int)));
            cl.SetKernelArg(clmultiply, 7, alpha);
            cl.SetKernelArg(clmultiply, 8, beta);
            cl.SetKernelArg(clmultiply, 9, cellsize);
            cl.SetKernelArg(clmultiply, 10, cellrowsperblock);
            cl.SetKernelArg(clmultiply, 11, cellsperrow);
            cl.SetKernelArg(clmultiply, 12, stride);
            cl.SetKernelArg(clmultiply, 13, size);
        }
    }
}
