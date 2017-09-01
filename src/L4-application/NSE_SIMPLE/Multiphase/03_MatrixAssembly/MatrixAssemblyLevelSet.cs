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
using BoSSS.Solution.NSECommon;

namespace NSE_SIMPLE {

    /// <summary>
    /// Matrix assembly for level-set equation.
    /// </summary>
    public class MatrixAssemblyLevelSet : SIMPLEMatrixAssembly {

        SIMPLEOperator m_Convection;
        
        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="Convection">
        /// Convection operator for scalar.
        /// </param>
        public MatrixAssemblyLevelSet(SIMPLEOperator Convection) : base(false, false, MaxUsePerIterMatrix: 2) {
            m_Convection = Convection;

            base.Initialize();
        }

        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Matrix = m_Convection.OperatorMatrix;
            return Matrix;
        }

        protected override double[] ComputeAffine() {
            double[] Affine = m_Convection.OperatorAffine;
            return Affine;
        }
    }
}