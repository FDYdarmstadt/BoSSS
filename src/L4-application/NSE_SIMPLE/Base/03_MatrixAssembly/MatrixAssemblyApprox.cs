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
    /// Approximation of a matrix, e.g. Identity, Diagonal, BlockDiagonal.
    /// Used e.g. to approximate the Predictor in Corrector step (the inverse of this approximation <see cref="MatrixAssemblyApproxInv"/>)
    /// and for under-relaxation.
    /// </summary>
    public class MatrixAssemblyApprox : SIMPLEMatrixAssembly {

        SolverConfiguration m_SolverConf;
        SIMPLEMatrixAssembly m_Src;
        BDFScheme m_BDF;
        int m_NoOfCells;
        BlockDiagonalMatrix m_Rho;

        /// <summary>
        /// [VariableDensity] Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="Src"></param>
        /// <param name="BDF"></param>
        /// <param name="Rho"></param>
        /// <param name="MaxUseMatrix"></param>
        public MatrixAssemblyApprox(SolverConfiguration SolverConf, int NoOfCells, SIMPLEMatrixAssembly Src, BDFScheme BDF, BlockDiagonalMatrix Rho, int MaxUseMatrix = 1)
            : base(SolverConf.Control.PredictorApproximationIsConstant, false, MaxUsePerIterMatrix: MaxUseMatrix) {

            m_SolverConf = SolverConf;
            m_Src = Src;
            m_BDF = BDF;
            m_NoOfCells = NoOfCells;

            if (Rho != null) {
                // VariableDensity
                m_Rho = Rho;
                if (m_Src.i0 != m_Rho.RowPartitioning.i0)
                    throw new NotImplementedException("Different row partitions, which should be equal.");
            } else {
                // Incompressible
                m_Rho = new BlockDiagonalMatrix(m_Src.LocalLength, 1);
                m_Rho.AccEyeSp();
            }

            base.Initialize();
        }

        /// <summary>
        /// [Incompressible] Ctor.
        /// </summary>
        public MatrixAssemblyApprox(SolverConfiguration SolverConf, int NoOfCells, SIMPLEMatrixAssembly Src, BDFScheme BDF, int MaxUseMatrix = 1)
            : this(SolverConf, NoOfCells, Src, BDF, null, MaxUseMatrix) {
        }

        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Src = m_Src.AssemblyMatrix;
            MsrMatrix Approx = new MsrMatrix(Src.RowPartitioning, Src.ColPartition);
            
            int i0 = Src.RowPartitioning.i0;
            int LocalLength = Src.RowPartitioning.LocalLength;

            double BDFfactor;
            switch (m_SolverConf.Control.Algorithm) {
                case SolutionAlgorithms.Steady_SIMPLE:
                    BDFfactor = 0.0;
                    break;
                case SolutionAlgorithms.Unsteady_SIMPLE:
                    int BDFOrder = m_SolverConf.BDFOrder;
                    double dt = m_SolverConf.dt;
                    BDFfactor = m_BDF.beta[BDFOrder - 1][0] / (m_BDF.gamma[BDFOrder - 1] * dt);
                    break;
                default:
                    throw new ArgumentException();
            }

            switch (m_SolverConf.Control.PredictorApproximation) {
                case PredictorApproximations.Identity:
                case PredictorApproximations.Identity_IP1:
                    Approx.AccEyeSp();
                    Approx.AccEyeSp(BDFfactor);
                    break;
                case PredictorApproximations.Diagonal:
                    for (int row = 0; row < LocalLength; row++) {
                        double Src_ii = Src[row + i0, row + i0];
                        double Rho_ii = m_Rho[row + i0, row + i0];
                        Approx[row + i0, row + i0] = BDFfactor * Rho_ii + Src_ii;
                    }
                    break;
                case PredictorApproximations.BlockDiagonal:
                    BlockDiagonalMatrix SrcBlock = new BlockDiagonalMatrix(Src, Src.RowPartitioning.LocalLength / m_NoOfCells, Src.ColPartition.LocalLength / m_NoOfCells );
                    Approx.Acc(1.0, SrcBlock);
                    Approx.Acc(BDFfactor, m_Rho);
                    break;
                default:
                    throw new ArgumentException();
            }

            return Approx;
        }

        protected override double[] ComputeAffine() {
            double[] ZeroAffine = new double[m_Src.LocalLength];
            return ZeroAffine;
        }
    }
}
