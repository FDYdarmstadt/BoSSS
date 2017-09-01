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

namespace NSE_SIMPLE {

    /// <summary>
    /// Corrector for variable density flows.
    /// </summary>
    public class MatrixAssemblyVariableDensityCorrector : SIMPLEMatrixAssembly {

        SIMPLEOperator[] m_VelocityDivergence;
        SIMPLEMatrixAssembly[] m_PredictorApproxInv;
        SIMPLEOperator[] m_PressureGradient;

        /// <summary>
        /// Ctor for multiphase.
        /// </summary>
        /// <param name="VelocityDivergence"></param>
        /// <param name="PredictorApproxInv"></param>
        /// <param name="PressureGradient"></param>
        /// <param name="UpdateCycleAppPred"></param>
        public MatrixAssemblyVariableDensityCorrector(SIMPLEOperator[] VelocityDivergence, SIMPLEMatrixAssembly[] PredictorApproxInv, SIMPLEOperator[] PressureGradient,
            int UpdateCycleAppPred)
            : base(PredictorApproxInv[0].IsConstant, false, MaxUsePerIterMatrix: UpdateCycleAppPred) {

            if (VelocityDivergence.Length != PressureGradient.Length)
                throw new ArgumentException("Mismatch in dimensions");

            m_VelocityDivergence = VelocityDivergence;
            m_PredictorApproxInv = PredictorApproxInv;
            m_PressureGradient = PressureGradient;

            base.Initialize();
        }

        /// <summary>
        /// Ctor for Low-Mach.
        /// </summary>
        /// <param name="VelocityDivergence"></param>
        /// <param name="PredictorApproxInv"></param>
        /// <param name="PressureGradient"></param>
        public MatrixAssemblyVariableDensityCorrector(SIMPLEOperator[] VelocityDivergence, SIMPLEMatrixAssembly[] PredictorApproxInv, SIMPLEOperator[] PressureGradient)
            : base(false, false) {

            if (VelocityDivergence.Length != PressureGradient.Length)
                throw new ArgumentException("Mismatch in dimensions");

            m_VelocityDivergence = VelocityDivergence;
            m_PredictorApproxInv = PredictorApproxInv;
            m_PressureGradient = PressureGradient;

            base.Initialize();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix CorrectorMatrix = new MsrMatrix(m_VelocityDivergence[0].LocalLength);

            for (int comp = 0; comp < m_VelocityDivergence.Length; comp++) {
                MsrMatrix prod1 = MsrMatrix.Multiply(m_PredictorApproxInv[comp].AssemblyMatrix, m_PressureGradient[comp].OperatorMatrix);
                MsrMatrix prod2 = MsrMatrix.Multiply(m_VelocityDivergence[comp].OperatorMatrix, prod1);
                CorrectorMatrix.Acc(1.0, prod2);
            }

            //CorrectorMatrix.AssumeSymmetric = false;

            return CorrectorMatrix;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override double[] ComputeAffine() {
            //Corrector has got no affine part
            return null;
        }

    }
}
