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
    class clCSRMatrix : clMatrix
    {
        private cl_mem d_val;
        private cl_mem d_rowStart;
        private cl_mem d_colIdx;

        private int size;

        public clCSRMatrix(MsrMatrix M, clDevice device)
            : base(M, device, "csrMultiply")
        {
            size = base.RowPartitioning.LocalLength;
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
            return new MatrixBase.CSR(tmp);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            MatrixBase.CSR LocalMatrix = (MatrixBase.CSR)base.m_LocalMtx;

            d_val = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)LocalMatrix.Val.Length * sizeof(double), LocalMatrix.Val);
            d_rowStart = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)LocalMatrix.RowStart.Length * sizeof(double), LocalMatrix.RowStart);
            d_colIdx = cl.CreateBuffer(device.env.context, cl_mem_flags.CL_MEM_COPY_HOST_PTR | cl_mem_flags.CL_MEM_READ_ONLY, (uint)LocalMatrix.ColInd.Length * sizeof(double), LocalMatrix.ColInd);
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cl.ReleaseMemObject(d_val);
            cl.ReleaseMemObject(d_rowStart);
            cl.ReleaseMemObject(d_colIdx);
        }

        internal override void SetArguments(double alpha, clVector a, double beta, clVector acc)
        {
            cl_mem d_x = a.GetDevicePointer();
            cl_mem d_result = acc.GetDevicePointer();

            cl.SetKernelArg(clmultiply, 0, d_val);
            cl.SetKernelArg(clmultiply, 1, d_colIdx);
            cl.SetKernelArg(clmultiply, 2, d_rowStart);
            cl.SetKernelArg(clmultiply, 3, d_result);
            cl.SetKernelArg(clmultiply, 4, d_x);
            cl.SetKernelArgLocalSize(clmultiply, 5, (uint)(localsize + 1) * sizeof(int));
            cl.SetKernelArg(clmultiply, 6, alpha);
            cl.SetKernelArg(clmultiply, 7, beta);
            cl.SetKernelArg(clmultiply, 8, size);
        }
    }
}
