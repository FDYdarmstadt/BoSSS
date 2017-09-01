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
    /// Splitting of Swip2 and Swip3 term for decoupled solver.
    /// </summary>
    public class MatrixAssemblyViscSplit : SIMPLEMatrixAssembly {

        SIMPLEOperator Swip2;
        SIMPLEOperator Swip3;
        int[] RowIndicesSource;
        int[] ColumnIndicesSource;
        int[] ColumnIndicesTarget;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="Swip2"></param>
        /// <param name="Swip3">Can be null for multiphase flows.</param>
        /// <param name="Component">Component of velocity vector.</param>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        public MatrixAssemblyViscSplit(SIMPLEOperator Swip2, SIMPLEOperator Swip3, int Component,
            UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping)
            : base(false, false) {

            this.Swip2 = Swip2;
            this.Swip3 = Swip3;

            RowIndicesSource = new int[Swip2.LocalLength];
            int i0 = Swip2.RowPartition.i0;
            for (int i = 0; i < Swip2.LocalLength; i++) {
                RowIndicesSource[i] = i + i0;
            }

            ColumnIndicesSource = VelocityVectorMapping.GetSubvectorIndices(true, Component);
            ColumnIndicesTarget = VelocityMapping.GetSubvectorIndices(true, 0);

            base.Initialize();
        }

        protected override MsrMatrix ComputeMatrix() {
            MsrMatrix Matrix = new MsrMatrix(Swip2.RowPartition);

            Swip2.OperatorMatrix.WriteSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(Matrix,
                RowIndicesSource,
                null,
                ColumnIndicesSource,
                ColumnIndicesTarget);

            if (Swip3 != null) {
                Swip3.OperatorMatrix.AccSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(1.0,
                    Matrix,
                    RowIndicesSource,
                    null,
                    ColumnIndicesSource,
                    ColumnIndicesTarget);
            }

            return Matrix;
        }

        protected override double[] ComputeAffine() {            
            if (Swip3 != null) {
                double[] Affine = new double[Swip2.LocalLength];
                double[] AffineSwip2 = Swip2.OperatorAffine;
                double[] AffineSwip3 = Swip3.OperatorAffine;
                for (int i = 0; i < Affine.Length; i++) {
                    Affine[i] = AffineSwip2[i] + AffineSwip3[i];
                }
                return Affine;
            } else {
                return Swip2.OperatorAffine;
            }
        }
    }
}
