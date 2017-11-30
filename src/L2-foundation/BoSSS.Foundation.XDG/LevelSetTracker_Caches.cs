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


        HistoryStack<LevelSetData>[] m_DataHistories;

        /// <summary>
        /// Stack for data of previous level-set states. 
        /// </summary>
        public HistoryStack<LevelSetData>[] DataHistories {
            get {
                return m_DataHistories;
            }
        }


        /// <summary>
        /// Caches for normals, curvature, etc, for each level-set.
        /// </summary>
        public class LevelSetData {
            LevSetValueCache m_LevelSetValc;

            LevelSetGradientCache m_LevelSetGradientsCache;

            LevelSetNormalsCache m_LevelSetNormalsCache;

            LevelSetReferenceGradientCache m_LevelSetReferenceGradientsCache;

            LevelSetReferenceNormalsCache m_LevelSetReferenceNormalsCache;

            LevelSetReferenceCurvatureCache m_LevelSetReferenceCurvatureCache;

            LevelSetReferenceHessianCache m_LevelSetReferenceHessianCache;

            internal void ClearCaches() {
                m_LevelSetValc.Clear();
                m_LevelSetGradientsCache.Clear();
                m_LevelSetReferenceGradientsCache.Clear();
                m_LevelSetNormalsCache.Clear();
                m_LevelSetReferenceNormalsCache.Clear();
            }

            /// <summary>
            /// Link to background grid.
            /// </summary>
            public GridData GridDat {
                get {
                    return m_owner.GridDat;
                }
            }

            /// <summary>
            /// Region object which correlates with actual <see cref="HistoryIndex"/>.
            /// </summary>
            public LevelSetTracker.LevelSetRegions Region {
                get {
                    return m_owner.RegionsHistory[HistoryIndex];
                }
            }

            LevelSetTracker m_owner; // don't expose this public, it would circumvent the idea of encapsulating a certain level-set time level.
            internal int m_StackIdx = 1;
            int m_iLevSet;

            /// <summary>
            /// Index into the level-set history, see e.g. <see cref="LevelSetTracker.LevelSetHistories"/>.
            /// </summary>
            public int HistoryIndex {
                get {
                    return m_StackIdx;
                }
            }

            /// <summary>
            /// Index of the corresponding level-set, see e.g. <see cref="LevelSetTracker.LevelSetHistories"/>.
            /// </summary>
            public int LevelSetIndex {
                get {
                    return m_iLevSet;
                }
            }


            ILevelSet GetLevSet() {
                return m_owner.m_LevelSetHistories[m_iLevSet][m_StackIdx];
            }


            internal LevelSetData(LevelSetTracker __owner, int iLevSet) {
                m_owner = __owner;
                m_iLevSet = iLevSet;

                m_LevelSetValc = new LevSetValueCache(this);
                m_LevelSetGradientsCache = new LevelSetGradientCache(this);
                m_LevelSetNormalsCache = new LevelSetNormalsCache(this);
                m_LevelSetReferenceGradientsCache = new LevelSetReferenceGradientCache(this);
                m_LevelSetReferenceNormalsCache = new LevelSetReferenceNormalsCache(this);
                m_LevelSetReferenceCurvatureCache = new LevelSetReferenceCurvatureCache(this);
                m_LevelSetReferenceHessianCache = new LevelSetReferenceHessianCache(this);

            }



            /// <summary>
            /// Cached evaluation of the level set fields (values of the level-set field).
            /// </summary>
            /// <param name="NS"></param>
            /// <param name="j0">local index of first cell to evaluate</param>
            /// <param name="Len">number of cells to evaluate</param>
            /// <returns>
            /// values of the level set field at the nodes of the specified node set.
            /// </returns>
            public MultidimensionalArray GetLevSetValues(NodeSet NS, int j0, int Len) {
                return m_LevelSetValc.GetValue_Cell(NS, j0, Len);
            }

            /// <summary>
            /// 
            /// </summary>
            class LevSetValueCache : Caching.CacheLogic_CNs {

                internal LevSetValueCache(LevelSetData o)
                    : base(o.m_owner.GridDat) //
                {
                    m_owner = o;    
                }

                LevelSetData m_owner;

                /// <summary>
                /// 
                /// </summary>
                protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                    m_owner.GetLevSet().Evaluate(j0, Len, N, output);
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

                LevelSetData m_owner;

                /// <summary>
                /// Constructs a value cache for the evaluation of the gradient of
                /// the given level set.
                /// </summary>
                internal LevelSetGradientCache(LevelSetData o)
                    : base(o.m_owner.GridDat) //
                {
                    m_owner = o;
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
                    m_owner.GetLevSet().EvaluateGradient(j0, Len, N, output);
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
            public MultidimensionalArray GetLevelSetGradients(NodeSet NS, int j0, int Len) {
                return m_LevelSetGradientsCache.GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetReferenceGradientCache : Caching.CacheLogic_CNs {

                LevelSetData m_owner;

                internal LevelSetReferenceGradientCache(LevelSetData o)
                    : base(o.m_owner.GridDat) //
                {
                    this.m_owner = o;
                }

                protected override void ComputeValues(NodeSet N, int j0, int Len, MultidimensionalArray output) {
                    LevelSet levelSetField = m_owner.GetLevSet() as LevelSet;
                    if(levelSetField == null) {
                        ComputeValuesNonField(m_owner.GetLevSet(), N, j0, Len, output);
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
                    var R = m_owner.m_owner.GridDat.Cells.Transformation;
                    var JacDet = m_owner.m_owner.GridDat.ChefBasis.Scaling;
                    var Cells = m_owner.m_owner.GridDat.Cells;

                    MultidimensionalArray physGradient = m_owner.GetLevelSetGradients(NS, j0, Len);
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
            public MultidimensionalArray GetLevelSetReferenceGradients(NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceGradientsCache.GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetNormalsCache : Caching.CacheLogic_CNs {

                private LevelSetData m_owner;


                public LevelSetNormalsCache(LevelSetData o)
                    : base(o.m_owner.GridDat) //
                {
                    this.m_owner = o;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    MultidimensionalArray gradient =
                        m_owner.m_LevelSetGradientsCache.GetValue_Cell(NS, j0, Len);
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
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            ///  /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            public MultidimensionalArray GetLevelSetNormals(NodeSet NS, int j0, int Len) {
                return m_LevelSetNormalsCache.GetValue_Cell(NS, j0, Len);
            }

            private class LevelSetReferenceNormalsCache : Caching.CacheLogic_CNs {

                private LevelSetData m_owner;


                public LevelSetReferenceNormalsCache(LevelSetData o)
                    : base(o.m_owner.GridDat) //
                {
                    m_owner = o;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    MultidimensionalArray gradient =
                        m_owner.m_LevelSetReferenceGradientsCache.GetValue_Cell(NS, j0, Len);
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
            /// <param name="NS">
            /// Nodes to evaluate at.
            /// </param>
            ///  /// <param name="j0">
            /// The first cell to be evaluated
            /// </param>
            /// <param name="Len">
            /// The number of cell to be evaluated
            /// </param>
            public MultidimensionalArray GetLevelSetReferenceNormals(NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceNormalsCache.GetValue_Cell(NS, j0, Len);
            }

            /// <summary>
            /// Some integral transformation metrics....
            /// </summary>
            public MultidimensionalArray GetLevelSetNormalReferenceToPhysicalMetrics(NodeSet NS, int j0, int Len) {
                MultidimensionalArray physGradients = GetLevelSetGradients(NS, j0, Len);
                MultidimensionalArray refGradients = GetLevelSetReferenceGradients(NS, j0, Len);
                int noOfNodes = physGradients.GetLength(1);
                int D = m_owner.GridDat.Grid.SpatialDimension;
                var OneOverSqrt_AbsJacobiDet = m_owner.GridDat.ChefBasis.Scaling;
                MultidimensionalArray result = MultidimensionalArray.Create(Len, noOfNodes);

                for(int i = 0; i < Len; i++) {
                    int jCell = j0 + i;

                    if(m_owner.GridDat.Cells.IsCellAffineLinear(jCell)) {
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

                private LevelSetData m_owner;

                private int levelSetIndex;

                public LevelSetReferenceHessianCache(LevelSetData owner)
                    : base(owner.m_owner.GridDat) //
                {
                    this.m_owner = owner;
                }



                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {
                    LevelSet LevSet = (LevelSet)(this.m_owner.GetLevSet());

                    var BasisHessian = LevSet.Basis.Evaluate2ndDeriv(NS);

                    Debug.Assert(output.GetLength(0) == Len);
                    int N = LevSet.Basis.Length;
                    Debug.Assert(BasisHessian.GetLength(1) == N);
                    int D = this.m_owner.m_owner.GridDat.SpatialDimension;
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
            /// <param name="NodeSet"></param>
            /// <param name="j0">first cell to evaluate</param>
            /// <param name="Len">number of cells to evaluate</param>
            /// <returns></returns>
            public MultidimensionalArray GetLevelSetReferenceHessian(NodeSet nodes, int j0, int Len) {
                return m_LevelSetReferenceHessianCache.GetValue_Cell(nodes, j0, Len);
            }

            private class LevelSetReferenceCurvatureCache : Caching.CacheLogic_CNs {

                private LevelSetData m_owner;

                private int levelSetIndex;

                public LevelSetReferenceCurvatureCache(LevelSetData owner)
                    : base(owner.m_owner.GridDat) //
                {
                    this.m_owner = owner;
                }

                protected override void ComputeValues(NodeSet NS, int j0, int Len, MultidimensionalArray output) {

                    MultidimensionalArray Phi = m_owner.GetLevSetValues(NS, j0, Len);
                    MultidimensionalArray GradPhi = m_owner.GetLevelSetReferenceGradients(NS, j0, Len);
                    MultidimensionalArray HessPhi = m_owner.GetLevelSetReferenceHessian(NS, j0, Len);

                    MultidimensionalArray ooNormGrad = new MultidimensionalArray(2);
                    MultidimensionalArray Laplace = new MultidimensionalArray(2);
                    MultidimensionalArray Q = new MultidimensionalArray(3);

                    int K = output.GetLength(1);
                    int D = GradPhi.GetLength(2);
                    Debug.Assert(D == this.m_owner.m_owner.GridDat.SpatialDimension);

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
            public MultidimensionalArray GetLevelSetReferenceCurvature(NodeSet NS, int j0, int Len) {
                return m_LevelSetReferenceCurvatureCache.GetValue_Cell(NS, j0, Len);
            }

            public object Clone() {
                throw new NotImplementedException();
            }
        }

    }
}
