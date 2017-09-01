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
    /// Matrix assembly for Predictor (i.e. convection and viscous operator of momentum equation)
    /// of incompressible flows.
    /// </summary>
    public class MatrixAssemblyIncompressiblePredictor : SIMPLEMatrixAssembly {

        SIMPLEOperator m_Convection;
        SIMPLEOperator m_Visc;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="Convection">
        /// Convection operator for momentum equation.
        /// </param>
        /// <param name="Visc">
        /// Viscous (Laplace) operator  for momentum equation.
        /// </param>        
        /// <param name="MaxUsePerIterMatrix"></param>
        /// <param name="MaxUsePerIterAffine"></param>
        public MatrixAssemblyIncompressiblePredictor(SIMPLEOperator Convection, SIMPLEOperator Visc, int MaxUsePerIterMatrix, int MaxUsePerIterAffine)
            : base(false, Convection.OnlyAffine, MaxUsePerIterMatrix, MaxUsePerIterAffine) {
            m_Convection = Convection;
            m_Visc = Visc;

            base.Initialize();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Matrix = new MsrMatrix(m_Convection.OperatorMatrix);
            Matrix.Acc(m_Visc.OperatorMatrix, 1.0);
            return Matrix;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override double[] ComputeAffine() {
            double[] Affine = new double[m_Convection.LocalLength];
            double[] AffineConvection = m_Convection.OperatorAffine;
            double[] AffineVisc = m_Visc.OperatorAffine;
            for (int i = 0; i < AffineConvection.Length; i++) {
                Affine[i] = AffineConvection[i] + AffineVisc[i];
            }
            return Affine;
        }
    }
}