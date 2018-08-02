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
        public class BasisData {

            /// <summary>
            /// reference to owner object
            /// </summary>
            private IGridData m_Owner;

            /// <summary>
            /// ctor
            /// </summary>
            internal BasisData(GridData o) {
                m_Owner = o;

                MultidimensionalArray JacDet = m_Owner.iGeomCells.JacobiDet;

                int J = m_Owner.iGeomCells.Count;
                var sc = MultidimensionalArray.Create(J);
                Scaling = sc;
                for (int j = 0; j < J; j++) {
                    if (m_Owner.iGeomCells.IsCellAffineLinear(j))
                        sc[j] = 1.0 / Math.Sqrt(JacDet[j]);
                    else
                        sc[j] = double.NaN;
                }

                this.OrthonormalizationTrafo = new CacheLogicImpl_CP(this.m_Owner, this.Compute_OrthonormalizationTrafo);

                // caches for reference eval
                // -------------------------

                this.BasisValues = new CacheLogicImpl_NsP(this.EvaluateBasis);
                this.BasisGradientValues = new CacheLogicImpl_NsP(this.EvaluateBasisGradient);
                this.BasisHessianValues = new CacheLogicImpl_NsP(this.EvaluateBasisHessian);

                this.EdgeEval = new BasisEdgeValuesCacheLogic(this.m_Owner, false);
                this.EdgeGradientEval = new BasisEdgeValuesCacheLogic(this.m_Owner, true);

                // caches for cell/physical eval
                // -----------------------------

                this.CellBasisValues = new CacheLogicImpl_CNsP(this.m_Owner, this.EvaluateCellBasis,
                    delegate(int j0, int Len, int Degree, NodeSet Ns) {
                        int Ndim = GetNdim(j0, Len, Degree, Ns);
                        return MultidimensionalArray.Create(Len, Ns.NoOfNodes, Ndim);
                    });
                this.CellBasisGradientValues = new CacheLogicImpl_CNsP(this.m_Owner, this.EvaluateCellBasisGradient,
                    delegate(int j0, int Len, int Degree, NodeSet Ns) {
                        int Ndim = GetNdim(j0, Len, Degree, Ns);
                        int D = this.m_Owner.SpatialDimension;
                        return MultidimensionalArray.Create(Len, Ns.NoOfNodes, Ndim, D);
                    });
                this.CellBasisHessianValues = new CacheLogicImpl_CNsP(this.m_Owner, this.EvaluateCellBasisHessian,
                    delegate(int j0, int Len, int Degree, NodeSet Ns) {
                        int Ndim = GetNdim(j0, Len, Degree, Ns);
                        int D = this.m_Owner.SpatialDimension;
                        return MultidimensionalArray.Create(Len, Ns.NoOfNodes, Ndim, D, D);
                    });

                

            }

            /// <summary>
            /// Number of degrees-of-freedom for polynomial degree <paramref name="Degree"/>.
            /// </summary>
            /// <param name="j0"></param>
            /// <param name="Len"></param>
            /// <param name="Degree"></param>
            /// <param name="Ns"></param>
            /// <returns></returns>
            private int GetNdim(int j0, int Len, int Degree, NodeSet Ns) {
                int Ndim;
                NodeCoordinateSystem NsSys = Ns.GetNodeCoordinateSystem(this.m_Owner);
                switch(NsSys) {
                    case NodeCoordinateSystem.CellCoord: {
                        int iKref = this.m_Owner.iGeomEdges.GetRefElementIndex(j0);
                        Ndim = this.GetOrthonormalPolynomials(Degree)[iKref].Count;
#if DEBUG
                        for(int j = 0; j < Len; j++) {
                            int jCell = j + j0;
                            Debug.Assert(iKref == this.m_Owner.iGeomEdges.GetRefElementIndex(jCell));
                        }
#endif
                        break;
                    }

                    case NodeCoordinateSystem.EdgeCoord: {
                        int iKref = this.m_Owner.iGeomCells.GetRefElementIndex(this.m_Owner.iGeomEdges.CellIndices[j0, 0]);
                        Ndim = this.GetOrthonormalPolynomials(Degree)[iKref].Count;

#if DEBUG
                        for(int j = 0; j < Len; j++) {
                            int iEdge = j + j0;


                            Debug.Assert(iKref == this.m_Owner.iGeomCells.GetRefElementIndex(this.m_Owner.iGeomEdges.CellIndices[iEdge, 0]));
                            if(this.m_Owner.iGeomEdges.CellIndices[iEdge, 1] >= 0) {
                                Debug.Assert(iKref == this.m_Owner.iGeomCells.GetRefElementIndex(this.m_Owner.iGeomEdges.CellIndices[iEdge, 1]));
                            }
                        }


#endif
                        break;
                    }

                    default: throw new NotImplementedException("Todo.");
                }
                return Ndim;
            }

            /// <summary>
            /// For a cell <i>K<sub>j</sub></i>, with an affine-linear
            /// transformation to the reference cell, the basis polynomials are
            /// scaled by the factor <see cref="Scaling"/>[<i><sub>j</sub></i>]
            /// in order to preserve orthonormality under the transformation.
            /// </summary>
            /// <remarks>
            /// index: cell index.
            /// </remarks>
            public MultidimensionalArray Scaling {
                private set;
                get;
            }


            SortedDictionary<int, PolynomialList[]> m_PolynomialLists = new SortedDictionary<int, PolynomialList[]>();

            /// <summary>
            /// Returns the orthonormal approximation polynomials \f$ \phi_n \f$ up to a specific degree.
            /// </summary>
            public PolynomialList[] GetOrthonormalPolynomials(int Degree) {
                PolynomialList[] R;
                if(!this.m_PolynomialLists.TryGetValue(Degree, out R)) {
                    var Krefs = this.m_Owner.iGeomCells.RefElements;
                    R = new PolynomialList[Krefs.Length];
                    for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                        R[iKref] = Krefs[iKref].GetOrthonormalPolynomials(Degree);
                    }
                    this.m_PolynomialLists.Add(Degree, R);
                }
                return R;
            }

            SortedDictionary<int, PolynomialList[,]> m_Polynomial1stDerivLists = new SortedDictionary<int, PolynomialList[,]>();

            /// <summary>
            /// Returns the 1st derivatives of the orthonormal approximation polynomials, 
            /// i.e. \f$ \nabla_{\vec{\xi}} \phi_n \f$ up to a specific degree.
            /// </summary>
            /// <returns>
            /// 1st index: reference element;
            /// 2nd index: spatial direction
            /// </returns>
            public PolynomialList[,] GetOrthonormalPolynomials1stDeriv(int Degree) {
                PolynomialList[,] R;
                if(!this.m_Polynomial1stDerivLists.TryGetValue(Degree, out R)) {
                    var Krefs = this.m_Owner.iGeomCells.RefElements;
                    int D = this.m_Owner.SpatialDimension;

                    R = new PolynomialList[Krefs.Length, D];
                    PolynomialList[] Polys = this.GetOrthonormalPolynomials(Degree);

                    int[] DerivExp = new int[D];

                    for(int d = 0; d < D; d++) {
                        DerivExp[d] = 1;
                        for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                            int N = Polys[iKref].Count;
                            Polynomial[] tmp = new Polynomial[N];
                            for(int n = 0; n < N; n++) {
                                tmp[n] = Polys[iKref][n].Derive(DerivExp);
                            }

                            R[iKref, d] = new PolynomialList(tmp);
                        }
                        DerivExp[d] = 0;
                    }
                    this.m_Polynomial1stDerivLists.Add(Degree, R);
                }
                return R;
            }

            SortedDictionary<int, PolynomialList[,,]> m_Polynomial2ndDerivLists = new SortedDictionary<int, PolynomialList[,,]>();


            /// <summary>
            /// Returns the 2nd derivatives of the orthonormal approximation polynomials, 
            /// i.e. \f$ \nabla_{\vec{\xi}} \phi_n \f$ up to a specific degree.
            /// </summary>
            /// <returns>
            /// 1st index: reference element;
            /// 2nd index: spatial direction
            /// 2nd index: spatial direction
            /// </returns>
            public PolynomialList[,,] GetOrthonormalPolynomials2ndDeriv(int Degree) {
                PolynomialList[,,] R;
                if(!this.m_Polynomial2ndDerivLists.TryGetValue(Degree, out R)) {
                    var Krefs = this.m_Owner.iGeomCells.RefElements;
                    int D = this.m_Owner.SpatialDimension;

                    R = new PolynomialList[Krefs.Length, D, D];
                    PolynomialList[] Polys = this.GetOrthonormalPolynomials(Degree);

                    int[] DerivExp = new int[D];

                    for(int d1 = 0; d1 < D; d1++) {
                        DerivExp[d1]++;
                        for(int d2 = 0; d2 < D; d2++) {
                            DerivExp[d2]++;
                            
                            for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                                int N = Polys[iKref].Count;
                                Polynomial[] tmp = new Polynomial[N];
                                for(int n = 0; n < N; n++) {
                                    tmp[n] = Polys[iKref][n].Derive(DerivExp);
                                }

                                R[iKref, d1, d2] = new PolynomialList(tmp);
                            }
                            DerivExp[d2]--;
                        }
                        DerivExp[d1]--;
                    }

#if DEBUG
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = d1 + 1; d2 < D; d2++) {
                            for(int iKref = 0; iKref < Krefs.Length; iKref++) {
                                int N = Polys[iKref].Count;

                                for(int n =0; n < N; n++) {
                                    Polynomial P_d1d2 = R[iKref, d1, d2][n];
                                    Polynomial P_d2d1 = R[iKref, d2, d1][n];

                                    Debug.Assert(P_d1d2.Equals(P_d2d1),"Hessian seems unsymmetric.");
                                }
                            }  
                        }
                    }
#endif


                    this.m_Polynomial2ndDerivLists.Add(Degree, R);
                }
                return R;
            }


            /// <summary>
            /// Change of basis to provide orthonormality in physical coordinates.
            /// The orthonormalization transformation depends on cell and polynomial degree of the basis,
            /// but is independent from the used node set.
            /// </summary>
            public CacheLogicImpl_CP OrthonormalizationTrafo {
                get;
                private set;
            }

            MultidimensionalArray Compute_OrthonormalizationTrafo(int j0, int Len, int Degree) {

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



            /// <summary>
            /// Cached evaluation of basis polynomials \f$ \phi_n \f$ in reference space.
            /// </summary>
            public CacheLogic_NsP BasisValues {
                get;
                private set;
            }
                        
            /// <summary>
            /// Used by the <see cref="BasisValues"/>-cache.
            /// </summary>
            MultidimensionalArray EvaluateBasis(NodeSet NS, int MinDegree) {
                int iKref = NS.GetVolumeRefElementIndex(this.m_Owner);
                PolynomialList Polys = this.GetOrthonormalPolynomials(MinDegree)[iKref];
                MultidimensionalArray R = MultidimensionalArray.Create(NS.NoOfNodes, Polys.Count);
                Polys.Evaluate(NS, R);
                return R;
            }

            /// <summary>
            /// Cached evaluation of basis polynomials gradients \f$ \nabla_{\vec{xi}} \phi_n \f$ in reference space.
            /// </summary>
            public CacheLogic_NsP BasisGradientValues {
                get;
                private set;
            }

            /// <summary>
            /// Used by the <see cref="BasisGradientValues"/>-cache.
            /// </summary>
            MultidimensionalArray EvaluateBasisGradient(NodeSet NS, int MinDegree) {
                int iKref = NS.GetVolumeRefElementIndex(this.m_Owner);
                int D = this.m_Owner.SpatialDimension;
                PolynomialList[,] Polys = this.GetOrthonormalPolynomials1stDeriv(MinDegree);
                int N = Polys[iKref, 0].Count;

                MultidimensionalArray R = MultidimensionalArray.Create(NS.NoOfNodes, N, D);
                for(int d = 0; d < D; d++) {
                    Debug.Assert(Polys[iKref, d].Count == N);
                    //Debug.Assert(Polys[iKref, d].MaxAbsoluteDegree == MinDegree);
                        
                    Polys[iKref, d].Evaluate(NS, R.ExtractSubArrayShallow(-1, -1, d));
                }
                return R;
            }

            /// <summary>
            /// Cached evaluation of basis polynomials Hessian \f$ \partial_{\vec{xi}}^2 \phi_n \f$ in reference space.
            /// </summary>
            public CacheLogic_NsP BasisHessianValues {
                get;
                private set;
            }


            /// <summary>
            /// Used by the <see cref="BasisHessianValues"/>-cache.
            /// </summary>
            MultidimensionalArray EvaluateBasisHessian(NodeSet NS, int MinDegree) {
                int iKref = NS.GetVolumeRefElementIndex(this.m_Owner);
                int D = this.m_Owner.SpatialDimension;
                PolynomialList[,,] Polys = this.GetOrthonormalPolynomials2ndDeriv(MinDegree);
                int N = Polys[0, 0, 0].Count;

                MultidimensionalArray R = MultidimensionalArray.Create(NS.NoOfNodes, N, D, D);
                for(int d1 = 0; d1 < D; d1++) {
                    for(int d2 = 0; d2 < D; d2++) {
                        Debug.Assert(Polys[iKref, d1, d2].Count == N);
                        //Debug.Assert(Polys[iKref, d1, d2].MaxAbsoluteDegree == MinDegree);
                        Polys[iKref, d1, d2].Evaluate(NS, R.ExtractSubArrayShallow(-1, -1, d1, d2));
                    }
                }
                return R;
            }

            void EvaluateCellBasis(NodeSet NodeSet, int j0, int Len, int degree, MultidimensionalArray output) {
                
                MultidimensionalArray BasisValRef = this.BasisValues.GetValues(NodeSet, degree); // basis values on reference element
                Debug.Assert(BasisValRef.Dimension == 2);
                Debug.Assert(BasisValRef.GetLength(0) == NodeSet.NoOfNodes);

                Debug.Assert(output.Dimension == 3);
                Debug.Assert(output.GetLength(0) == Len);
                Debug.Assert(output.GetLength(1) == NodeSet.NoOfNodes);

                int N = output.GetLength(2);
                Debug.Assert(BasisValRef.GetLength(1) >= N);
                if(BasisValRef.GetLength(1) > N)
                    BasisValRef = BasisValRef.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { NodeSet.NoOfNodes - 1, N - 1 });

                bool AffineLinear = m_Owner.iGeomCells.IsCellAffineLinear(j0);
#if DEBUG
                for(int ii = 1; ii < Len; ii++) {
                    if(m_Owner.iGeomCells.IsCellAffineLinear(j0 + ii) != AffineLinear)
                        throw new NotSupportedException("Switching cell type in one chunk is not supported.");
                }
#endif


                if(AffineLinear) {

                    for(int i = 0; i < Len; i++) {
                        int jCell = i + j0;
                        var res_j = output.ExtractSubArrayShallow(i, -1, -1);

                        double scale = this.Scaling[jCell];

                        Debug.Assert(res_j.Dimension == 2 && BasisValRef.Dimension == 2);
                        Debug.Assert(res_j.GetLength(0) == BasisValRef.GetLength(0));
                        Debug.Assert(res_j.GetLength(1) == BasisValRef.GetLength(1));

                        res_j.Clear();
                        res_j.Acc(scale, BasisValRef);
                    }

                } else {

                    MultidimensionalArray Trafo = this.OrthonormalizationTrafo.GetValue_Cell(j0, Len, degree);
                    if(Trafo.GetLength(1) != N)
                        Trafo = Trafo.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, N - 1, N - 1 });

                    output.Multiply(1.0, BasisValRef, Trafo, 0.0, "jkn", "km", "jmn");
                }
            }

            /// <summary>
            /// Values of polynomials in cells,  \f$ \phi_{j n} \f$
            /// </summary>
            public Caching.CacheLogicImpl_CNsP CellBasisValues {
                get;
                private set;
            }

            /// <summary>
            /// Values of polynomials gradients in cells, \f$ \nabla_{\vec{x}} \phi_{j n} \f$
            /// </summary>
            public Caching.CacheLogicImpl_CNsP CellBasisGradientValues {
                get;
                private set;
            }

            void EvaluateCellBasisGradient(NodeSet NodeSet, int j0, int Len, int degree, MultidimensionalArray output) {
                MultidimensionalArray resRef = this.BasisGradientValues.GetValues(NodeSet, degree);
                Debug.Assert(resRef.GetLength(0) == NodeSet.NoOfNodes);
                Debug.Assert(resRef.Dimension == 3);
                
                int N = output.GetLength(2);
                int D = this.m_Owner.SpatialDimension;
                //MultidimensionalArray output = MultidimensionalArray.Create(Len, NodeSet.NoOfNodes, N);
                Debug.Assert(output.Dimension == 4);
                Debug.Assert(output.GetLength(0) == Len);
                Debug.Assert(output.GetLength(1) == NodeSet.NoOfNodes);
                Debug.Assert(output.GetLength(3) == D);

                Debug.Assert(resRef.GetLength(1) >= N);
                if(resRef.GetLength(1) > N)
                    resRef = resRef.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { NodeSet.NoOfNodes - 1, N - 1, D - 1 });


                bool AffineLinear = m_Owner.iGeomCells.IsCellAffineLinear(j0);
#if DEBUG
                for(int ii = 1; ii < Len; ii++) {
                    if(m_Owner.iGeomCells.IsCellAffineLinear(j0 + ii) != AffineLinear)
                        throw new NotSupportedException("Switching cell type in one chunk is not supported.");
                }
#endif


                

                if(AffineLinear) {

                    var R = this.m_Owner.iGeomCells.InverseTransformation;
                    var scales = this.Scaling;

                    for(int i = 0; i < Len; i++) {
                        int jCell = j0 + i;
                        double scale = scales[jCell];

                        var output_i = output.ExtractSubArrayShallow(i, -1, -1, -1);
                        var R_i = R.ExtractSubArrayShallow(j0 + i, -1, -1);
                        output_i.Multiply(scale, R_i, resRef, 0.0, "kmd", "ld", "kml");
                    }

                } else {
                    var weights = this.OrthonormalizationTrafo.GetValue_Cell(j0, Len, degree);
                    if(weights.GetLength(1) != N)
                        weights = weights.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, N - 1, N - 1 });

                    MultidimensionalArray JacInverse = this.m_Owner.InverseJacobian.GetValue_Cell(NodeSet, j0, Len);
                    output.Multiply(1.0, JacInverse, weights, resRef, 0.0, "jknd", "jked", "jmn", "kme");
                }
            }

            /// <summary>
            /// Values of polynomials gradients in cells, \f$ \partial^2_{\vec{x}} \phi_{j n} \f$
            /// </summary>
            public Caching.CacheLogicImpl_CNsP CellBasisHessianValues {
                get;
                private set;
            }


            void EvaluateCellBasisHessian(NodeSet NodeSet, int j0, int Len, int degree, MultidimensionalArray Hess) {
                MultidimensionalArray gbv = this.BasisHessianValues.GetValues(NodeSet, degree);
                
                int N = Hess.GetLength(2);
                int D = this.m_Owner.SpatialDimension;
                int K = NodeSet.NoOfNodes;
                //MultidimensionalArray output = MultidimensionalArray.Create(Len, NodeSet.NoOfNodes, N);
                
                Debug.Assert(Hess.Dimension == 5);
                Debug.Assert(Hess.GetLength(0) == Len);
                Debug.Assert(Hess.GetLength(1) == K);
                Debug.Assert(gbv.GetLength(1) >= N);
                Debug.Assert(Hess.GetLength(3) == D);
                Debug.Assert(Hess.GetLength(4) == D);

                Debug.Assert(gbv.Dimension == 4);
                Debug.Assert(gbv.GetLength(0) == K);
                Debug.Assert(gbv.GetLength(1) >= N);
                Debug.Assert(gbv.GetLength(2) == D);
                Debug.Assert(gbv.GetLength(3) == D);
                if(gbv.GetLength(1) > N)
                    gbv = gbv.ExtractSubArrayShallow(new int[] { 0, 0, 0, 0 }, new int[] { NodeSet.NoOfNodes - 1, N - 1, D - 1, D - 1 });
                
                bool AffineLinear = m_Owner.iGeomCells.IsCellAffineLinear(j0);
#if DEBUG
                for(int ii = 1; ii < Len; ii++) {
                    if(m_Owner.iGeomCells.IsCellAffineLinear(j0 + ii) != AffineLinear)
                        throw new NotSupportedException("Switching cell type in one chunk is not supported.");
                }
#endif

                               
                if(AffineLinear) {
                    MultidimensionalArray scales = this.Scaling;
                    var Tinv = this.m_Owner.iGeomCells.InverseTransformation;
                

                    for(int j = 0; j < Len; j++) { // loop over cells
                        double sc = scales[j0 + j];
                        
                        MultidimensionalArray Tinv_j = Tinv.ExtractSubArrayShallow(j + j0, -1, -1);

                        // loop over nodes ...
                        for(int k = 0; k < K; k++) {


                            // transform...
                            for(int n = 0; n < N; n++) { // loop over basis functions

                                for(int d1 = 0; d1 < D; d1++) {  // 1st derivation index 
                                    for(int d2 = 0; d2 < D; d2++) { // 2nd derivation index

                                        double r = 0.0;
                                        for(int k1 = 0; k1 < D; k1++) {
                                            double rr = 0;
                                            for(int k2 = 0; k2 < D; k2++) {
                                                rr += gbv[k, n, k1, k2] * Tinv_j[k2, d1];
                                            }
                                            r += rr * Tinv_j[k1, d2];
                                        }

                                        Hess[j, k, n, d1, d2] = r * sc;
                                    }
                                }
                            }
                        }
                    }

                } else {
                    throw new NotImplementedException("todo");
                }
            }

            /// <summary>
            /// Evaluation of the reference basis values \f$ \phi_n \f$ for edges.
            /// </summary>
            public Caching.BasisEdgeValuesCacheLogic EdgeEval {
                get;
                private set;
            }

            /// <summary>
            /// Evaluation of the reference basis gradient values \f$ \nabla_{\vec{xi}} \phi_n \f$ for edges.
            /// </summary>
            public Caching.BasisEdgeValuesCacheLogic EdgeGradientEval {
                get;
                private set;
            }
        }
    }
}
