using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.XDG {
    partial class LevelSetTracker {


        /// <summary>
        /// Caches for normals, curvature, etc, for each level-set.
        /// </summary>
        public class LevelSetData {
            LevSetValueCache[] m_LevelSetValc;

            LevelSetGradientCache[] m_LevelSetGradientsCache;

            LevelSetNormalsCache[] m_LevelSetNormalsCache;

            LevelSetReferenceGradientCache[] m_LevelSetReferenceGradientsCache;

            LevelSetReferenceNormalsCache[] m_LevelSetReferenceNormalsCache;

            LevelSetReferenceCurvatureCache[] m_LevelSetReferenceCurvatureCache;

            LevelSetReferenceHessianCache[] m_LevelSetReferenceHessianCache;



            /// <summary>
            /// Cached evaluation of the level set fields (values of the level-set field).
            /// </summary>
            /// <param name="levSetInd"></param>
            /// <param name="NS"></param>
            /// <param name="j0">local index of first cell to evaluate</param>
            /// <param name="Len">number of cells to evaluate</param>
            /// <returns>
            /// values of the level set field at the nodes of the specified node set.
            /// </returns>
            public MultidimensionalArray GetLevSetValues(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetValc[levSetInd].GetValue_Cell(NS, j0, Len);
            }

            /// <summary>
            /// 
            /// </summary>
            class LevSetValueCache : Caching.CacheLogic_CNs {

                internal LevSetValueCache(ILevelSet levSet, GridData gridData)
                    : base(gridData) {
                    m_LevSet = levSet;
                }

                ILevelSet m_LevSet;

                /// <summary>
                /// 
                /// </summary>
                protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                    m_LevSet.Evaluate(j0, Len, N, output);
                }

                /// <summary>
                /// 
                /// </summary>
                protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                    return MultidimensionalArray.Create(Len, NS.NoOfNodes);
                }
            }



            /// <summary>
            /// Caches the gradient of a level set
            /// </summary>
            private class LevelSetGradientCache : Caching.CacheLogic_CNs {

                /// <summary>
                /// The level set in question
                /// </summary>
                private ILevelSet levelSet;

                /// <summary>
                /// Constructs a value cache for the evaluation of the gradient of
                /// the given level set.
                /// </summary>
                /// <param name="levelSet">
                /// The level set in question
                /// </param>
                /// <param name="nsc">
                /// The node set controller holding the evaluation nodes
                /// </param>
                internal LevelSetGradientCache(ILevelSet levelSet, GridData gridData)
                    : base(gridData) //
                {
                    this.levelSet = levelSet;
                }

                /// <summary>
                /// <see cref="ILevelSet.EvaluateGradient"/>
                /// </summary>
                /// <param name="NodeSetIndex">
                /// <see cref="ILevelSet.EvaluateGradient"/>
                /// </param>
                /// <param name="j0">
                /// <see cref="ILevelSet.EvaluateGradient"/>
                /// </param>
                /// <param name="Len">
                /// <see cref="ILevelSet.EvaluateGradient"/>
                /// </param>
                /// <param name="output">
                /// <see cref="ILevelSet.EvaluateGradient"/>
                /// </param>
                protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                    levelSet.EvaluateGradient(j0, Len, N, output);
                }

                /// <summary>
                /// <see cref="NodeSetController.ByCellValueCache{MultidimensionalArray}.Allocate"/>
                /// </summary>
                protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet N) {
                    return MultidimensionalArray.Create(Len, N.NoOfNodes, N.SpatialDimension);
                }


            }

            /// <summary>
            /// Calculates the gradients of the level set in every affected cell
            /// and every point contained the node set identified by
            /// <paramref name="NodeSetIndex"/>. The layout of the resulting array
            /// is equivalent to <see cref="ILevelSet.EvaluateGradient"/>.
            /// </summary>
            /// <param name="levSetInd">The index of the level set</param>
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            /// <returns></returns>
            public MultidimensionalArray GetLevelSetGradients(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetGradientsCache[levSetInd].GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetReferenceGradientCache : Caching.CacheLogic_CNs {

                private ILevelSet levelSet;

                private LevelSetTracker owner;

                internal LevelSetReferenceGradientCache(ILevelSet levelSet, LevelSetTracker owner)
                    : base(owner.GridDat) //
                {
                    this.levelSet = levelSet;
                    this.owner = owner;
                }

                protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                    LevelSet levelSetField = levelSet as LevelSet;
                    if(levelSetField == null) {
                        ComputeValuesNonField(levelSet, N, j0, Len, output);
                    } else {
                        ComputeValuesField(levelSetField, N, j0, Len, output);
                    }
                }

                protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                    int D = base.GridData.SpatialDimension;
                    return MultidimensionalArray.Create(Len, NS.NoOfNodes, D);
                }

                private void ComputeValuesField(LevelSet levelSet, NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    int D = levelSet.Basis.GridDat.SpatialDimension;
                    MultidimensionalArray grad = levelSet.Basis.EvaluateGradient(NS);
                    int noOfNodes = grad.GetLength(0);

                    unsafe {
                        fixed (double* pGrad = &grad.Storage[0], pRes = &output.Storage[0]) {
                            double* pResCur = pRes;
                            double* pGradCur = pGrad;
                            for(int i = 0; i < Len; i++) {
                                for(int j = 0; j < noOfNodes; j++) {
                                    //double norm = 0.0;
                                    for(int d = 0; d < D; d++) {
                                        for(int k = 0; k < levelSet.Basis.MinimalLength; k++) {
                                            *(pResCur + d) += levelSet.Coordinates[i + j0, k] * *(pGradCur + grad.Index(j, k, d));
                                        }
                                    }

                                    pResCur += D;
                                }
                            }
                        }
                    }

                    //// Reference implementation
                    //for (int i = 0; i < Len; i++) {
                    //    for (int j = 0; j < noOfNodes; j++) {
                    //        for (int d = 0; d < D; d++) {
                    //            for (int k = 0; k < levelSetField.Basis.MinimalLength; k++) {
                    //                output[i, j, d] += levelSetField.Coordinates[i + j0, k] * grad[j, k, d];
                    //            }
                    //        }
                    //    }
                    //}
                }

                private void ComputeValuesNonField(ILevelSet levelSet, NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    int noOfNodes = NS.NoOfNodes;
                    int D = NS.SpatialDimension;
                    var R = owner.GridDat.Cells.Transformation;
                    var JacDet = owner.GridDat.ChefBasis.Scaling;
                    var Cells = owner.GridDat.Cells;

                    MultidimensionalArray physGradient = owner.GetLevelSetGradients(0, NS, j0, Len);
                    for(int i = 0; i < Len; i++) {
                        int jCell = j0 + i;
                        if(Cells.IsCellAffineLinear(jCell)) {
                            double det = JacDet[j0 + i];

                            for(int j = 0; j < noOfNodes; j++) {
                                for(int d = 0; d < D; d++) {
                                    double r = 0.0;
                                    for(int dd = 0; dd < D; dd++) {
                                        r += R[i + j0, dd, d] * physGradient[i, j, dd];
                                    }

                                    output[i, j, d] = r / det;
                                }
                            }
                        } else {
                            throw new NotImplementedException("todo: nonlinear cell");
                        }
                    }
                }



            }


            /// <summary>
            /// Calculates the gradients in the reference coordinate system of the level set in every affected cell
            /// and every point contained the node set identified by
            /// <paramref name="NodeSetIndex"/>. The layout of the resulting array
            /// is equivalent to <see cref="ILevelSet.EvaluateGradient"/>.
            /// </summary>
            /// <param name="levSetInd">The index of the level set</param>
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            /// <returns></returns>
            public MultidimensionalArray GetLevelSetReferenceGradients(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceGradientsCache[levSetInd].GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetNormalsCache : Caching.CacheLogic_CNs {

                private LevelSetTracker owner;

                private int levelSetIndex;

                public LevelSetNormalsCache(LevelSetTracker owner, int levelSetIndex)
                    : base(owner.GridDat) //
                {
                    this.owner = owner;
                    this.levelSetIndex = levelSetIndex;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    MultidimensionalArray gradient =
                        owner.m_LevelSetGradientsCache[levelSetIndex].GetValue_Cell(NS, j0, Len);
                    //gradient.Storage.CopyTo(output.Storage, 0);
                    Debug.Assert(gradient.Dimension == 3 && output.Dimension == 3);
                    Debug.Assert(gradient.GetLength(1) == output.GetLength(1));
                    Debug.Assert(gradient.GetLength(2) == output.GetLength(2));
                    Debug.Assert(output.GetLength(0) == Len);
                    output.Set(gradient.ExtractSubArrayShallow(new int[] { 0, 0, 0 }, new int[] { Len - 1, gradient.GetLength(1) - 1, gradient.GetLength(2) - 1 }));

                    int N = gradient.GetLength(1);
                    int D = gradient.GetLength(2);

                    for(int i = 0; i < Len; i++) {
                        for(int j = 0; j < gradient.GetLength(1); j++) {
                            double normOfGradient = 0.0;
                            for(int k = 0; k < D; k++) {
                                normOfGradient += gradient[i, j, k] * gradient[i, j, k];
                            }

                            // Avoid NaN. Normal is zero in this case
                            if(normOfGradient == 0.0) {
                                continue;
                            }

                            double OOnormOfGradient = 1.0 / Math.Sqrt(normOfGradient);

                            for(int k = 0; k < D; k++) {
                                output[i, j, k] *= OOnormOfGradient;
                            }
                        }
                    }
                }

                protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                    return MultidimensionalArray.Create(Len, NS.NoOfNodes, NS.SpatialDimension);
                }
            }

            /// <summary>
            /// Variant of <see cref="GetLevelSetGradients"/> which normalizes the
            /// gradients before returning them
            /// </summary>
            ///  <param name="levSetInd">The index of the level set</param>
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            ///  /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            public MultidimensionalArray GetLevelSetNormals(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetNormalsCache[levSetInd].GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetReferenceNormalsCache : Caching.CacheLogic_CNs {

                private LevelSetTracker owner;

                private int levelSetIndex;

                public LevelSetReferenceNormalsCache(LevelSetTracker owner, int levelSetIndex)
                    : base(owner.GridDat) //
                {
                    this.owner = owner;
                    this.levelSetIndex = levelSetIndex;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    MultidimensionalArray gradient =
                        owner.m_LevelSetReferenceGradientsCache[levelSetIndex].GetValue_Cell(NS, j0, Len);
                    gradient.Storage.CopyTo(output.Storage, 0);

                    for(int i = 0; i < gradient.GetLength(0); i++) {
                        for(int j = 0; j < gradient.GetLength(1); j++) {
                            double normOfGradient = 0.0;
                            for(int d = 0; d < gradient.GetLength(2); d++) {
                                normOfGradient += gradient[i, j, d] * gradient[i, j, d];
                            }

                            // Avoid NaN. Normal is zero in this case
                            if(normOfGradient == 0.0) {
                                continue;
                            }

                            normOfGradient = Math.Sqrt(normOfGradient);

                            for(int d = 0; d < gradient.GetLength(2); d++) {
                                output[i, j, d] = gradient[i, j, d] / normOfGradient;
                            }
                        }
                    }
                }

                protected override MultidimensionalArray Allocate(int i0, int Len, NodeSet NS) {
                    return MultidimensionalArray.Create(Len, NS.NoOfNodes, NS.SpatialDimension);
                }
            }

            /// <summary>
            /// Variant of <see cref="GetLevelSetReferenceGradients(int, NodeSet, int, int)"/> which normalizes the
            /// gradients before returning them.
            /// </summary>
            ///  <param name="levSetInd">The index of the level set</param>
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            ///  /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            public MultidimensionalArray GetLevelSetReferenceNormals(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceNormalsCache[levSetInd].GetValue_Cell(NS, j0, Len);
            }

            public MultidimensionalArray GetLevelSetNormalReferenceToPhysicalMetrics(int levSetInd, NodeSet NS, int j0, int Len) {
                MultidimensionalArray physGradients = GetLevelSetGradients(levSetInd, NS, j0, Len);
                MultidimensionalArray refGradients = GetLevelSetReferenceGradients(levSetInd, NS, j0, Len);
                int noOfNodes = physGradients.GetLength(1);
                int D = GridDat.Grid.SpatialDimension;
                var OneOverSqrt_AbsJacobiDet = GridDat.ChefBasis.Scaling;
                MultidimensionalArray result = MultidimensionalArray.Create(Len, noOfNodes);

                for(int i = 0; i < Len; i++) {
                    int jCell = j0 + i;

                    if(GridDat.Cells.IsCellAffineLinear(jCell)) {
                        double sc = OneOverSqrt_AbsJacobiDet[jCell];
                        for(int j = 0; j < noOfNodes; j++) {
                            double normPhys = 0.0;
                            double normRef = 0.0;

                            for(int d = 0; d < D; d++) {
                                normPhys += physGradients[i, j, d] * physGradients[i, j, d];
                                normRef += refGradients[i, j, d] * refGradients[i, j, d];
                            }

                            result[i, j] = Math.Sqrt(normRef / normPhys) * sc;
                        }
                    } else {

                        throw new NotImplementedException("nonlinear cell: todo");
                    }
                }


                //{
                //    double erracc = 0;

                //    var NormalsRef = this.GetLevelSetReferenceNormals(levSetInd, NodeSetIndex, j0, Len);
                //    var NoramlsPhys = this.GetLevelSetNormals(levSetInd, NodeSetIndex, j0, Len);

                //    MultidimensionalArray Jacobi = MultidimensionalArray.Create(D, D);
                //    double[] v1 = new double[D];


                //    for (int i = 0; i < Len; i++) {
                //        int jCell = i + j0;
                //        for (int d1 = 0; d1 < D; d1++) {
                //            for (int d2 = 0; d2 < D; d2++) {
                //                Jacobi[d1, d2] = Ctx.GridDat.Transformation[jCell, d1 + D*d2];
                //            }
                //        }

                //        for (int iNode = 0; iNode < noOfNodes; iNode++) {

                //            for (int d1 = 0; d1 < D; d1++) {
                //                double acc = 0;
                //                for (int d2 = 0; d2 < D; d2++) {
                //                    acc += Jacobi[d1, d2]*NormalsRef[i, iNode, d2];
                //                }
                //                v1[d1] = acc;
                //            }

                //            double metrix = 0;
                //            for (int d = 0; d < D; d++)
                //                metrix += v1[d]*NoramlsPhys[i, iNode, d];

                //            double anderemetrix = result[i, iNode];

                //            erracc += (metrix - anderemetrix).Pow2();
                //        }
                //    }

                //    Console.WriteLine("metrix test: " + erracc);
                //}




                return result;
            }



            private class LevelSetReferenceHessianCache : Caching.CacheLogic_CNs {

                private LevelSetTracker owner;

                private int levelSetIndex;

                public LevelSetReferenceHessianCache(LevelSetTracker owner, int levelSetIndex)
                    : base(owner.GridDat) //
                {
                    this.owner = owner;
                    this.levelSetIndex = levelSetIndex;
                }



                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    LevelSet LevSet = (LevelSet)this.owner.LevelSets[this.levelSetIndex];

                    var BasisHessian = LevSet.Basis.Evaluate2ndDeriv(NS);

                    Debug.Assert(output.GetLength(0) == Len);
                    int N = LevSet.Basis.Length;
                    Debug.Assert(BasisHessian.GetLength(1) == N);
                    int D = this.owner.m_gDat.SpatialDimension;
                    Debug.Assert(D == BasisHessian.GetLength(2));
                    Debug.Assert(D == BasisHessian.GetLength(3));
                    Debug.Assert(D == output.GetLength(2));
                    Debug.Assert(D == output.GetLength(3));
                    int K = output.GetLength(1); // No of nodes
                    Debug.Assert(K == BasisHessian.GetLength(0));

                    var Coordinates = ((MultidimensionalArray)LevSet.Coordinates).ExtractSubArrayShallow(new int[] { j0, 0 }, new int[] { j0 + Len - 1, N - 1 });
                    output.Multiply(1.0, BasisHessian, Coordinates, 0.0, "jkdr", "kndr", "jn");
                }

                protected override MultidimensionalArray Allocate(int j0, int Len, NodeSet N) {
                    return MultidimensionalArray.Create(Len, N.NoOfNodes, N.SpatialDimension, N.SpatialDimension);
                }
            }

            /// <summary>
            /// The Hessian of the level set field with respect to reference coordinates.
            /// </summary>
            /// <param name="levSetInd"></param>
            /// <param name="NodeSet"></param>
            /// <param name="j0">first cell to evaluate</param>
            /// <param name="Len">number of cells to evaluate</param>
            /// <returns></returns>
            public MultidimensionalArray GetLevelSetReferenceHessian(int levSetInd, NodeSet nodes, int j0, int Len) {
                return m_LevelSetReferenceHessianCache[levSetInd].GetValue_Cell(nodes, j0, Len);
            }

            private class LevelSetReferenceCurvatureCache : Caching.CacheLogic_CNs {

                private LevelSetTracker owner;

                private int levelSetIndex;

                public LevelSetReferenceCurvatureCache(LevelSetTracker owner, int levelSetIndex)
                    : base(owner.GridDat) //
                {
                    this.owner = owner;
                    this.levelSetIndex = levelSetIndex;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {

                    MultidimensionalArray Phi = owner.GetLevSetValues(this.levelSetIndex, NS, j0, Len);
                    MultidimensionalArray GradPhi = owner.GetLevelSetReferenceGradients(this.levelSetIndex, NS, j0, Len);
                    MultidimensionalArray HessPhi = owner.GetLevelSetReferenceHessian(this.levelSetIndex, NS, j0, Len);

                    MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
                    MultidimensionalArray Laplace = new MultidimensionalArray(2);
                    MultidimensionalArray Q = new MultidimensionalArray(3);

                    int K = output.GetLength(1);
                    int D = GradPhi.GetLength(2);
                    Debug.Assert(D == this.owner.GridDat.SpatialDimension);

                    ooNormGrad.Allocate(Len, K);
                    Laplace.Allocate(Len, K);
                    Q.Allocate(Len, K, D);


                    // compute the monstrous formula
                    // -----------------------------

                    // norm of Gradient:
                    for(int d = 0; d < D; d++) {
                        var GradPhi_d = GradPhi.ExtractSubArrayShallow(-1, -1, d);
                        ooNormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    ooNormGrad.ApplyAll(x => 1.0 / Math.Sqrt(x));

                    // laplacian of phi:
                    for(int d = 0; d < D; d++) {
                        var HessPhi_d_d = HessPhi.ExtractSubArrayShallow(-1, -1, d, d);
                        Laplace.Acc(1.0, HessPhi_d_d);
                    }

                    // result = Laplacian(phi)/|Grad phi|
                    output.Multiply(1.0, Laplace, ooNormGrad, 0.0, "ik", "ik", "ik");


                    // result = Grad(1/|Grad(phi)|)
                    for(int d1 = 0; d1 < D; d1++) {
                        var Qd = Q.ExtractSubArrayShallow(-1, -1, d1);

                        for(int d2 = 0; d2 < D; d2++) {
                            var Grad_d2 = GradPhi.ExtractSubArrayShallow(-1, -1, d2);
                            var Hess_d2_d1 = HessPhi.ExtractSubArrayShallow(-1, -1, d2, d1);

                            Qd.Multiply(-1.0, Grad_d2, Hess_d2_d1, 1.0, "ik", "ik", "ik");
                        }
                    }

                    ooNormGrad.ApplyAll(x => x * x * x);

                    output.Multiply(1.0, GradPhi, Q, ooNormGrad, 1.0, "ik", "ikd", "ikd", "ik");
                }

                protected override MultidimensionalArray Allocate(int j0, int Len, NodeSet N) {
                    return MultidimensionalArray.Create(Len, N.NoOfNodes);
                }
            }

            /// <summary>
            /// The curvature of the level set field with respect to reference coordinates.
            /// </summary>
            /// <param name="levSetInd"></param>
            /// <param name="NodeSet"></param>
            /// <param name="j0">first cell to evaluate</param>
            /// <param name="Len">number of cells to evaluate</param>
            /// <returns></returns>
            public MultidimensionalArray GetLevelSetReferenceCurvature(int levSetInd, NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceCurvatureCache[levSetInd].GetValue_Cell(NS, j0, Len);
            }

        }

    }
}
