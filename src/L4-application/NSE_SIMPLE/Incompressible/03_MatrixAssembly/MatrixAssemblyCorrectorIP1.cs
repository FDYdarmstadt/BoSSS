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
using BoSSS.Foundation;

namespace NSE_SIMPLE {

    /// <summary>
    /// Matrix assembly for corrector of incompressible flows using SIP discretization.
    /// </summary>
    public class MatrixAssemblyCorrectorIP1 : SIMPLEMatrixAssembly {

        SIMPLEOperator m_IPOperator;        

        SIMPLEOperator m_PressureStabilization = null;
        SolverConfiguration m_SolverConf;
        BDFScheme m_BDF;

        /// <summary>
        /// Ctor.
        /// </summary>        
        /// <param name="IPOperator"></param>                
        /// <param name="PressureStabilization"></param>
        /// <param name="SolverConf"></param>
        /// <param name="BDF"></param>
        public MatrixAssemblyCorrectorIP1(SIMPLEOperator IPOperator, SIMPLEOperator PressureStabilization, SolverConfiguration SolverConf, BDFScheme BDF)
            : base((SolverConf.Control.Algorithm == SolutionAlgorithms.Steady_SIMPLE), false) {

            m_IPOperator = IPOperator;            

            m_PressureStabilization = PressureStabilization;
            m_SolverConf = SolverConf;
            m_BDF = BDF;

            base.Initialize();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix CorrectorMatrix = new MsrMatrix(m_IPOperator.OperatorMatrix);

            switch (m_SolverConf.Control.Algorithm) {
                case SolutionAlgorithms.Steady_SIMPLE:
                    break;
                case SolutionAlgorithms.Unsteady_SIMPLE:
                    // gamma * dt / (beta_0 + gamma * dt)
                    double UnsteadyFactor = m_BDF.gamma[m_SolverConf.BDFOrder - 1] * m_SolverConf.dt / (m_BDF.beta[m_SolverConf.BDFOrder - 1][0] + m_BDF.gamma[m_SolverConf.BDFOrder - 1] * m_SolverConf.dt);
                    CorrectorMatrix.Scale(UnsteadyFactor);
                    break;
                default:
                    throw new NotImplementedException();
            }

            if (m_PressureStabilization != null) {
                CorrectorMatrix.Acc(-1.0, m_PressureStabilization.OperatorMatrix);
            }                

            CorrectorMatrix.AssumeSymmetric = true;
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
