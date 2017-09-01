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
    /// CUDA driver for <see cref="MatrixBase.CCBCSR"/> matrix format
    /// </summary>
    public class CudaCCBCSRMatrix : CudaMatrix
    {
        private CUdeviceptr d_cellData;
        private CUdeviceptr d_cellColIdx;

        private CCBCSR m_internalData;

        private int rowcount;
        private int cellsize;
        private int cellrowsperblock;
        private int cellsperrow;
        private int stride;
        private int blocksize;
        private int blockcount;

        /// <summary>
        /// 
        /// </summary>
        internal CudaCCBCSRMatrix(MsrMatrix M, CudaEnviroment CudaEnv)
            : base(M, "blockMultiply2", CudaEnv)
        {
            m_internalData = (CCBCSR)m_LocalMtx;

            rowcount = base.RowPartitioning.LocalLength;
            cellsize = m_internalData.CellSize;
            // Number of cells per block, choose so that it is around 128 threads per block
            cellrowsperblock = (int)Math.Ceiling(128.0 / cellsize);
            cellsperrow = m_internalData.NoOfCellsPerRow;
            stride = m_internalData.CellStride;
            // Number of threads per block
            blocksize = cellsize * cellrowsperblock;
            // Number of blocks
            blockcount = (int)Math.Ceiling((Decimal)rowcount / blocksize);
        }

        /// <summary>
        /// Allocate memory and copy matrix data on GPU
        /// </summary>
        public override void Lock()
        {
            base.Lock();

            cu.MemAlloc(out d_cellData, (uint)m_internalData.Val.Length * sizeof(double));
            cu.MemAlloc(out d_cellColIdx, (uint)m_internalData.CellColumn.Length * sizeof(int));
            cu.MemcpyHtoD(d_cellData, m_internalData.Val, (uint)m_internalData.Val.Length * sizeof(double));
            cu.MemcpyHtoD(d_cellColIdx, m_internalData.CellColumn, (uint)m_internalData.CellColumn.Length * sizeof(int));
        }

        /// <summary>
        /// Release memory resources on GPU
        /// </summary>
        public override void Unlock()
        {
            base.Unlock();

            cu.MemFree(d_cellData);
            cu.MemFree(d_cellColIdx);
        }

        /// <summary>
        /// 
        /// </summary>
        protected override MatrixBase.FormatBase AssembleFinalFormat(MatrixBase.TempCSR tmp)
        {
            return new CCBCSR(tmp, 1, base.CellSize);
        }

        internal override void CallDriver(CUstream stream, double alpha, CudaVector a, double beta, CudaVector acc)
        {
            CUdeviceptr d_x = a.GetDevicePointer();
            CUdeviceptr d_result = acc.GetDevicePointer();

            int offset = 0;

            cu.ParamSetp(sparseMultiply, offset, d_cellData);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_x);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_cellColIdx);
            offset += sizeof(long);
            cu.ParamSetp(sparseMultiply, offset, d_result);
            offset += sizeof(long);
            cu.ParamSetd(sparseMultiply, offset, alpha);
            offset += sizeof(double);
            cu.ParamSetd(sparseMultiply, offset, beta);
            offset += sizeof(double);
            cu.ParamSeti(sparseMultiply, offset, cellsize);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, cellrowsperblock);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, cellsperrow);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, stride);
            offset += sizeof(uint);
            cu.ParamSeti(sparseMultiply, offset, rowcount);
            offset += sizeof(uint);

            cu.ParamSetSize(sparseMultiply, (uint)offset);
            cu.FuncSetBlockShape(sparseMultiply, blocksize, 1, 1);
            cu.FuncSetSharedSize(sparseMultiply, (uint)(blocksize * sizeof(double) + 2 * cellrowsperblock * sizeof(int)));

            cu.LaunchGridAsync(sparseMultiply, blockcount, 1, stream);
        }
    }
}
