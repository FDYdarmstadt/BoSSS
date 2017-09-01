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
    /// [VariableDensity] Matrix assembly for Predictor,
    /// i.e. convection, viscosity term1
    /// and implicit parts of viscosity term2 and term3.
    /// (Term 3 is zero for multiphase flows)
    /// </summary>
    public class MatrixAssemblyVariableDensityPredictor : SIMPLEMatrixAssembly {

        int Component;
        SIMPLEOperator[] Convection;
        SIMPLEOperator[] ViscTerm1;
        SIMPLEMatrixAssembly[,] ViscTerm2Term3;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Component">Spatial component of velocity vector.</param>
        /// <param name="Convection"></param>
        /// <param name="ViscTerm1"></param>
        /// <param name="ViscTerm2Term3"></param>
        public MatrixAssemblyVariableDensityPredictor(int Component, SIMPLEOperator[] Convection, SIMPLEOperator[] ViscTerm1, SIMPLEMatrixAssembly[,] ViscTerm2Term3)
            : base(false, false, 2, 1) {

                this.Component = Component;
                this.Convection = Convection;
                this.ViscTerm1 = ViscTerm1;
                this.ViscTerm2Term3 = ViscTerm2Term3;

                base.Initialize();
        }

        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Matrix = new MsrMatrix(Convection[0].OperatorMatrix);
            Matrix.Acc(ViscTerm1[0].OperatorMatrix, 1.0);
            Matrix.Acc(ViscTerm2Term3[Component, Component].AssemblyMatrix, 1.0);
            return Matrix;
        }

        protected override double[] ComputeAffine() {
            double[] Affine = new double[Convection[Component].LocalLength];
            double[] AffineConvection = Convection[Component].OperatorAffine;
            double[] AffineViscTerm1 = ViscTerm1[Component].OperatorAffine;
            double[] AffineViscTerm2Term3 = ViscTerm2Term3[Component, Component].AssemblyAffine;
            for (int i = 0; i < Affine.Length; i++) {
                Affine[i] = AffineConvection[i] + AffineViscTerm1[i] + AffineViscTerm2Term3[i];
            }
            return Affine;
        }
    }
}
