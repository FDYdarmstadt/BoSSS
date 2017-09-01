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
using ilPSP.LinSolvers;
using BoSSS.Platform;

namespace NSE_SIMPLE {

    /// <summary>
    /// Inverse of approximation of some matrix <see cref="MatrixAssemblyApprox"/>.
    /// </summary>
    public class MatrixAssemblyApproxInv : SIMPLEMatrixAssembly {

        SIMPLEMatrixAssembly m_Approx;
        SIMPLEControl m_SolverConf;
        int m_LocalNoOfCells;

        /// <summary>
        /// Ctor.
        /// </summary>
        public MatrixAssemblyApproxInv(SIMPLEControl SolverConf, int LocalNoOfCells, SIMPLEMatrixAssembly Approx, int MaxUseMatrix = 1)
            : base(Approx.IsConstant, false, MaxUsePerIterMatrix: MaxUseMatrix) {

            m_Approx = Approx;
            m_SolverConf = SolverConf;
            m_LocalNoOfCells = LocalNoOfCells;
            base.Initialize();
        }

        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Approx = m_Approx.AssemblyMatrix;
            MsrMatrix ApproxInv = new MsrMatrix(Approx.RowPartitioning, Approx.ColPartition);

            switch (m_SolverConf.PredictorApproximation) {
                case PredictorApproximations.Identity:
                case PredictorApproximations.Identity_IP1:
                case PredictorApproximations.Diagonal:
                    int i0 = Approx.RowPartitioning.i0;
                    int LocalLength = Approx.RowPartitioning.LocalLength;

                    for (int row = 0; row < LocalLength; row++) {
                        double Approx_ii = Approx[row + i0, row + i0];
                        ApproxInv[row + i0, row + i0] = 1.0 / Approx_ii;
                    }
                    break;
                case PredictorApproximations.BlockDiagonal:
                    BlockDiagonalMatrix ApproxBlock = new BlockDiagonalMatrix(Approx, Approx.RowPartitioning.LocalLength / m_LocalNoOfCells, Approx.ColPartition.LocalLength / m_LocalNoOfCells);
                    BlockDiagonalMatrix ApproxBlockInv = ApproxBlock.Invert();
                    ApproxInv.Acc(1.0, ApproxBlockInv);
                    break;
                default:
                    throw new ArgumentException();
            }

            return ApproxInv;
        }

        protected override double[] ComputeAffine() {
            double[] ZeroAffine = new double[m_Approx.LocalLength];
            return ZeroAffine;
        }
    }
}
