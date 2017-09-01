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
    /// Matrix assembly for temperature equation in Low-Mach flows
    /// (i.e. convection and heat conduction).
    /// </summary>
    public class MatrixAssemblyTemperature : SIMPLEMatrixAssembly {

        SIMPLEOperator Convection;
        SIMPLEOperator HeatConduction;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Convection"></param>
        /// <param name="HeatConduction"></param>
        public MatrixAssemblyTemperature(SIMPLEOperator Convection, SIMPLEOperator HeatConduction)
            : base(false, false, MaxUsePerIterMatrix: 2) {

            this.Convection = Convection;
            this.HeatConduction = HeatConduction;

            base.Initialize();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Matrix = new MsrMatrix(Convection.OperatorMatrix);
            Matrix.Acc(HeatConduction.OperatorMatrix, 1.0);
            return Matrix;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override double[] ComputeAffine() {
            double[] Affine = new double[Convection.LocalLength];
            double[] AffineConvection = Convection.OperatorAffine;
            double[] AffineHeatConduction = HeatConduction.OperatorAffine;
            for (int i = 0; i < Affine.Length; i++) {
                Affine[i] = AffineConvection[i] + AffineHeatConduction[i];
            }
            return Affine;
        }
    }
}
