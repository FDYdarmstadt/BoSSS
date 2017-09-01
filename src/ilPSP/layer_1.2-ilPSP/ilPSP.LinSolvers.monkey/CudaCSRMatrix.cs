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

namespace ilPSP.LinSolvers.monkey.CUDA {
    
    /// <summary>
    /// CUDA matrix in standart CSR format (not very performant, only for reference purposes)
    /// </summary>
    public class CudaCSRMatrix : CudaMatrix
    {
        private CUdeviceptr d_val;
        private CUdeviceptr d_rowStart;
        private CUdeviceptr d_colIdx;

        private int blocksize = 256;
        private int size;
        private int blockcount;

        /// <summary>
        /// Constructor
        /// </summary>
        public CudaCSRMatrix(MsrMatrix M, CudaEnviroment CudaEnv)
            : base(M, "sparseMultiply", CudaEnv) {
            size = base.RowPartitioning.LocalLength;
            blockcount = (int)Math.Ceiling((decimal)size / blocksize);
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

            cu.MemAlloc(out d_val, (uint)LocalMatrix.Val.Length * sizeof(double));
            cu.MemAlloc(out d_rowStart, (uint)LocalMatrix.RowStart.Length * sizeof(int));
            cu.MemAlloc(out d_colIdx, (uint)LocalMatrix.ColInd.Length * sizeof(int));
            cu.MemcpyHtoD(d_val, LocalMatrix.Val, (uint)LocalMatrix.Val.Length * sizeof(double));
            cu.MemcpyHtoD(d_rowStart, LocalMatrix.RowStart, (uint)LocalMatrix.RowStart.Length * sizeof(int));
            cu.MemcpyHtoD(d_colIdx, LocalMatrix.ColInd, (uint)LocalMatrix.ColInd.Length * sizeof(int));

        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cu.MemFree(d_val);
            cu.MemFree(d_rowStart);
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
            cu.ParamSetp(sparseMultiply, offset, d_rowStart);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_result);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_x);
            offset += sizeof(long);
            cu.ParamSetd(sparseMultiply, offset, alpha);
            offset += sizeof(double);
            cu.ParamSetd(sparseMultiply, offset, beta);
            offset += sizeof(double);
            cu.ParamSeti(sparseMultiply, offset, size);
            offset += sizeof(uint);
            cu.ParamSetSize(sparseMultiply, (uint)offset);
            cu.FuncSetBlockShape(sparseMultiply, blocksize, 1, 1);
            cu.FuncSetSharedSize(sparseMultiply, (uint)((blocksize + 1) * sizeof(int)));

            cu.LaunchGridAsync(sparseMultiply, blockcount, 1, stream);
        }
    }
}
