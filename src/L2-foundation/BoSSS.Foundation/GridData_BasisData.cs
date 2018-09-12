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
using System.Diagnostics;
using System.Linq;
using BoSSS.Foundation.Caching;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;

namespace BoSSS.Foundation.Grid.Classic {

    partial class GridData {

        /// <summary>
        /// see <see cref="BasisData"/>
        /// </summary>
        public BasisData ChefBasis {
            private set;
            get;
        }

        /// <summary>
        /// common, heavy-weighted data for all <see cref="Basis"/>-objects.
        /// </summary>
        class _BasisData : BasisData {

            GridData m_Owner; 

            /// <summary>
            /// ctor
            /// </summary>
            internal _BasisData(GridData o) : base(o) {
                m_Owner = o;
            }

           
            /// <summary>
            /// Orthonormalization for curved elements
            /// </summary>
            protected override MultidimensionalArray Compute_OrthonormalizationTrafo(int j0, int Len, int Degree) {

                // init
                // ====

                int iKref = this.m_Owner.iGeomCells.GetRefElementIndex(j0);
                PolynomialList Polys = this.GetOrthonormalPolynomials(Degree)[iKref];
                int N = Polys.Count;
                CellType cellType = this.m_Owner.iGeomCells.GetCellType(j0);

#if DEBUG
                // checking
                for(int j = 1; j < Len; j++) {
                    int jCell = j + j0;

                    if(this.m_Owner.iGeomCells.GetCellType(jCell) != cellType)
                        throw new NotSupportedException("All cells in chunk must have same type.");
                }
#endif
                // storage for result
                MultidimensionalArray NonlinOrtho = MultidimensionalArray.Create(Len, N, N);


                if(m_Owner.iGeomCells.IsCellAffineLinear(j0)) {
                    // affine-linear branch
                    // ++++++++++++++++++++

                    MultidimensionalArray scl = this.Scaling;

                    for(int j = 0; j < Len; j++) {
                        int jCell = j + j0;

                        double scl_j = scl[jCell];

                        for(int n = 0; n < N; n++) {
                            NonlinOrtho[j, n, n] = scl_j;
                        }
                    }
                    
                } else {
                    // nonlinear cells branch
                    // ++++++++++++++++++++++

                    // init
                    // ====

                    var Kref = m_Owner.iGeomCells.GetRefElement(j0);
                    int deg; // polynomial degree of integrand: Degree of Jacobi determinat + 2* degree of basis polynomials in ref.-space.
                    {
                        int D = m_Owner.SpatialDimension;
                        deg = Kref.GetInterpolationDegree(cellType);
                        if(deg > 1)
                            deg -= 1;
                        deg *= D;
                        deg += 2 * Degree;
                    }
                    var qr = Kref.GetQuadratureRule((int)deg);
                    int K = qr.NoOfNodes;

                    // evaluate basis polys in ref space
                    // =================================

                    MultidimensionalArray BasisValues = this.EvaluateBasis(qr.Nodes, Degree);
                    Debug.Assert(BasisValues.GetLength(0) == K);
                    if(BasisValues.GetLength(0) > N) {
                        BasisValues = BasisValues.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, N - 1 });
                    }

                    // compute \f$ A_{k n m} = \phi_{k n} \phi_{k m} w_{k} \f$
                    // (cell-INdependentd part of mass matrix computation)
                    // ==================================================================

                    MultidimensionalArray A = MultidimensionalArray.Create(K, N, N);
                    A.Multiply(1.0, qr.Weights, BasisValues, BasisValues, 0.0, "knm", "k", "kn", "km");

                    // Determine mass matrix  \f$ M \f$ in all cells by quadrature
                    // \f$ M_{n m} = \sum_{k} J_k A_{k n m} \f$
                    // =====================================================

                    MultidimensionalArray J = m_Owner.JacobianDeterminat.GetValue_Cell(qr.Nodes, j0, Len);

                    // store mass-matrix in 'NonlinOrtho' to save mem alloc
                    NonlinOrtho.Multiply(1.0, A, J, 0.0, "jmn", "knm", "jk");


                    // Compute change-of-basis for all cells
                    // (invese of cholesky)
                    // =====================================


                    for(int j = 0; j < Len; j++) {
                        MultidimensionalArray Mj = NonlinOrtho.ExtractSubArrayShallow(j, -1, -1); // mass-matrix of basis on refrence element in cell j+j0


#if DEBUG
                        MultidimensionalArray MjClone = Mj.CloneAs();

#endif

                        //Mj.InvertSymmetrical();
                        Mj.SymmetricLDLInversion(Mj, null);

                        // clear lower triangular part
                        // Debug.Assert(Mj.NoOfCols == N);
                        for(int n = 1; n < N; n++) {
                            for(int m = 0; m < n; m++) {
                                Mj[n, m] = 0.0;
                            }
                        }
#if DEBUG

                        MultidimensionalArray B = NonlinOrtho.ExtractSubArrayShallow(j, -1, -1);
                        MultidimensionalArray Bt = B.Transpose();
                        double MjNorm = MjClone.InfNorm();
                        double Bnorm = B.InfNorm();


                        MultidimensionalArray check = IMatrixExtensions.GEMM(Bt, MjClone, B);
                        check.AccEye(-1.0);
                        double checkNorm = check.InfNorm();

                        double RelErr = checkNorm / Math.Max(MjNorm, Bnorm);

                        Debug.Assert(RelErr < 1.0e-5, "Fatal error in numerical orthonomalization on nonlinear cell.");
#endif

                    }
                }

                // return
                // =========

                return NonlinOrtho;
            }


           
        }
    }
}
