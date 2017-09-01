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
    /// Matrix assembly for corrector.
    /// </summary>
    public class MatrixAssemblyIncompressibleCorrector : SIMPLEMatrixAssembly {

        SIMPLEOperator[] m_DivergenceConti;
        SIMPLEMatrixAssembly m_PredictorApproxInv;
        SIMPLEOperator[] m_PressureGradient;

        SIMPLEOperator m_PressureStabilization = null;        

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="DivergenceConti"></param>
        /// <param name="PredictorApproxInv"></param>
        /// <param name="PressureGradient"></param>        
        /// <param name="UpdateCycleAppPred"></param>
        /// <param name="PressureStabilization">
        /// Can be null,
        /// i.e. no pressure stabilization is used.
        /// </param>        
        public MatrixAssemblyIncompressibleCorrector(SIMPLEOperator[] DivergenceConti, SIMPLEMatrixAssembly PredictorApproxInv, SIMPLEOperator[] PressureGradient,
            int UpdateCycleAppPred, SIMPLEOperator PressureStabilization)
            : base(PredictorApproxInv.IsConstant, false) {

            if (DivergenceConti.Length != PressureGradient.Length)
                throw new ArgumentException("Mismatch in dimensions");

            m_DivergenceConti = DivergenceConti;
            m_PredictorApproxInv = PredictorApproxInv;
            m_PressureGradient = PressureGradient;

            m_PressureStabilization = PressureStabilization;            

            base.Initialize();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix CorrectorMatrix = new MsrMatrix(m_DivergenceConti[0].LocalLength);

            for (int comp = 0; comp < m_DivergenceConti.Length; comp++) {
                MsrMatrix prod1 = MsrMatrix.Multiply(m_PredictorApproxInv.AssemblyMatrix, m_PressureGradient[comp].OperatorMatrix);
                MsrMatrix prod2 = MsrMatrix.Multiply(m_DivergenceConti[comp].OperatorMatrix, prod1);
                CorrectorMatrix.Acc(1.0, prod2);
            }

            if (m_PressureStabilization != null)
                CorrectorMatrix.Acc(-1.0, m_PressureStabilization.OperatorMatrix);            
            
            //CorrectorMatrix.AssumeSymmetric = true;

            return CorrectorMatrix;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override double[] ComputeAffine() {
            //Corrector has got no affine part
            double[] ZeroAffine = new double[m_DivergenceConti[0].LocalLength];
            return ZeroAffine;
        }
    }
}