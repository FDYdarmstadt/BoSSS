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

using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using System.Linq;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// The cubic reference element, i.e. $` K^{\textrm{cube} } = ( -1,1 )^3 $`.
    /// </summary>
    public partial class Cube : RefElement {
        

        /// <summary>
        /// The encoding to identify all six faces of the cube.
        /// </summary>
        public enum Faces {
            /// <summary>
            /// edge between this cell the neighbour cell with \f$ x \f$--coordinates closer to negative infinity.
            /// </summary>
            Left = 0,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ x \f$--coordinates closer to positive infinity.
            /// </summary>
            Right = 1,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ y \f$--coordinates closer to positive infinity.
            /// </summary>
            Top = 2,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ y \f$--coordinates closer to negative infinity.
            /// </summary>
            Bottom = 3,


            /// <summary>
            /// edge between this cell the neighbour cell with \f$ z \f$--coordinates closer to positive infinity
            /// </summary>
            Front = 4,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ z \f$--coordinates closer to negative infinity
            /// </summary>
            Back = 5

        }



        private static Cube instance = null;
        private static readonly object padlock = new object();
        
        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static Cube Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new Cube();
                    }
                    return instance;
                }
            }
        }




        /// <summary>
        /// standard constructor
        /// </summary>
        private Cube() {
            using(new FuncTrace()) { 

                // ===============
                // define vertices
                // ===============

                var _Vertices = new double[8, 3] { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 },
                                                   { 1, 1, -1 }, { 1, 1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }};
                this.m_Vertices = new NodeSet(this, 8, 3, false);
                this.m_Vertices.InitializeFrom(_Vertices);
                this.m_Vertices.LockForever();

                m_NoOfFaces = 6;

                // ============
                // edge simplex
                // ============

                m_FaceRefElement = Square.Instance;
            
                // ===================================
                // define Quadrature nodes and weights
                // ===================================

                {
                    //int mem = 0;
                    var qrTemp1D = QuadRuleResource.DecodeFromBase64( Resource.LineQuadRules_bin);
                    foreach(var q in qrTemp1D) {
                        m_1Drules.Add(q.Item1, (q.Item2, q.Item3));

                        
                    }

                    
                }

                // ==================================
                // define the orthonormal polynomials
                // ==================================
                DefinePolynomials();
            }
        }

        SortedDictionary<int, (double[,] Nodes, double[] Weights)> m_1Drules = new SortedDictionary<int, (double[,] Nodes, double[] Weights)>();

        SortedDictionary<int, QuadRule> m_3drules = new SortedDictionary<int, QuadRule>();


        /// <summary>
        /// <see cref="RefElement.GetQuadratureRule"/>
        /// </summary>
        /// <remarks>
        /// The 3D-Rules occupy a significant amount of memory (about 200 MB, particularly bad when running with lots of MPI cores:
        /// E.g., 100 Processors, 20 GB of memory just for quadrature rules), so we only create those that we need on the fly.
        /// </remarks>
        public override QuadRule GetQuadratureRule(int DesiredOrder) {
            if (DesiredOrder > HighestKnownOrder) {
                throw new ArgumentOutOfRangeException("no quadrature rule for desired order " + DesiredOrder + " available for simplex " + this.GetType().Name + ".", "DesiredOrder");
            }

            lock (m_3drules) {

                QuadRule realQr;
                if (!m_3drules.TryGetValue(DesiredOrder, out realQr)) {

                    //
                    //
                    //

                    int OrderOfPrecision = m_1Drules.Keys.Where(order => order >= DesiredOrder).Min();

                    if (!m_3drules.TryGetValue(OrderOfPrecision, out realQr)) {
                        var _1Drule = m_1Drules[OrderOfPrecision];

                        int NN = _1Drule.Weights.GetLength(0);
                        int D = this.SpatialDimension;
                        realQr = QuadRule.CreateEmpty(this, NN * NN * NN, D, true);

                        for (int i = 0; i < NN; i++) {
                            for (int j = 0; j < NN; j++) {
                                for (int k = 0; k < NN; k++) {
                                    realQr.Nodes[(i * NN + j) * NN + k, 0] = _1Drule.Nodes[k, 0];
                                    realQr.Nodes[(i * NN + j) * NN + k, 1] = _1Drule.Nodes[j, 0];
                                    realQr.Nodes[(i * NN + j) * NN + k, 2] = _1Drule.Nodes[i, 0];
                                    realQr.Weights[(i * NN + j) * NN + k] = _1Drule.Weights[i] * _1Drule.Weights[j] * _1Drule.Weights[k];
                                }
                            }
                        }

                        realQr.OrderOfPrecision = OrderOfPrecision;
                        realQr.Nodes.LockForever();
                        realQr.Weights.LockForever();
                        m_3drules.Add(OrderOfPrecision, realQr); // the rule of order 'OrderOfPrecision' must also be used for order 'DesiredOrder'
                        if (DesiredOrder != OrderOfPrecision)
                            m_3drules.Add(DesiredOrder, realQr);
                    } else {
                        m_3drules.Add(DesiredOrder, realQr);
                    }
                }
                
                return realQr;
            }

        }

        /// <summary>
        /// <see cref="RefElement.HighestKnownOrder"/>
        /// </summary>
        override public int HighestKnownOrder {
            get {
                return m_1Drules.Keys.Max();
            }
        }
        /// <summary>
        /// transforms some vertices (<paramref name="FaceVertices"/>) from the local 2D-coordinate system of either
        /// the top, bottom, left, right, front or back edge (see <see cref="Faces"/>) to the local 
        /// coordinate system of the cube;
        /// </summary>
        /// <param name="FaceIndex">0, 1, 2, 3, 4 or 5; <see cref="Faces"/></param>
        /// <param name="FaceVertices">input;</param>
        /// <param name="VolumeVertices">output;</param>
        public override void TransformFaceCoordinates(int FaceIndex, MultidimensionalArray FaceVertices, MultidimensionalArray VolumeVertices) {
            if (FaceVertices.Dimension != 2)
                throw new ArgumentException("dimension of EdgeVertices must be 2.", "EdgeVertices");
            if (VolumeVertices.Dimension != 2)
                throw new ArgumentException("dimension of VolumeVertices must be 2.", "VolumeVertices");
            if (VolumeVertices.GetLength(1) != 3)
                throw new ArgumentException("wrong spatial dimension of output", "VolumeVertices");
            if (FaceVertices.GetLength(1) != 2)
                throw new ArgumentException("wrong spatial dimension of input", "EdgeVertices");
            if (FaceVertices.GetLength(0) != VolumeVertices.GetLength(0))
                throw new ArgumentException("mismatch in number of vertices between input and output.", "EdgeVertices,VolumeVertices");
            
            int L = FaceVertices.GetLength(0);

            switch (FaceIndex) {
                case (int)Faces.Front:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = FaceVertices[i, 0];
                        VolumeVertices[i, 1] = FaceVertices[i, 1];
                        VolumeVertices[i, 2] = 1.0;
                    }
                    break;

                case (int)Faces.Back:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = FaceVertices[i, 0];
                        VolumeVertices[i, 1] = FaceVertices[i, 1];
                        VolumeVertices[i, 2] = -1.0;
                    }
                    break;


                case (int)Faces.Left:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = -1.0;
                        VolumeVertices[i, 1] = FaceVertices[i, 1];
                        VolumeVertices[i, 2] = FaceVertices[i, 0];
                    }
                    break;

                case (int)Faces.Right:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = 1.0;
                        VolumeVertices[i, 1] = FaceVertices[i, 1];
                        VolumeVertices[i, 2] = FaceVertices[i, 0];
                    }
                    break;

                case (int)Faces.Top:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = FaceVertices[i, 0];
                        VolumeVertices[i, 1] = 1.0;
                        VolumeVertices[i, 2] = FaceVertices[i, 1];
                    }
                    break;

                case (int)Faces.Bottom:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = FaceVertices[i, 0];
                        VolumeVertices[i, 1] = -1.0;
                        VolumeVertices[i, 2] = FaceVertices[i, 1];
                    }
                    break;


                default:
                    throw new ArgumentException("EdgeIndex out of range");
            }
        }

        /// <summary>
        /// partitions this cube into 8 sub-cubes of equal size;
        /// </summary>
        /// <returns></returns>
        public override AffineTrafo[] GetSubdivision() {
            AffineTrafo[] ret = new AffineTrafo[8];

            for (int i = 0; i < 8; i++) {
                ret[i] = new AffineTrafo(3);
                ret[i].Matrix = MultidimensionalArray.Create(3,3); ret[i].Matrix.AccEye(1.0);
                ret[i].Matrix.Scale(0.5);

                ret[i].Affine = Vertices.ExtractSubArrayShallow(i, -1).To1DArray();
                BLAS.dscal(3, 0.5, ret[i].Affine, 1);
            }

            return ret;
        }

        /// <summary>
        /// tests whether <paramref name="pt"/> is within the convex hull of
        /// vertices or not;
        /// </summary>
        /// <param name="pt"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public override bool IsWithin(double[] pt, double tolerance) {
            if (pt.Length != 3)
                throw new ArgumentException("wrong spatial dimension.", "pt");

            if ((pt[0] < -1.0 - tolerance) || (pt[0] > 1.0 + tolerance)
                || (pt[1] < -1.0 - tolerance) || (pt[1] > 1.0 + tolerance)
                || (pt[2] < -1.0 - tolerance) || (pt[2] > 1.0 + tolerance))
                return false;
            else
                return true;
        }


        /// <summary>
        /// see <see cref="RefElement.GetNodeSet(int,out MultidimensionalArray,out int[])"/>
        /// </summary>
        protected override void GetNodeSet(int px, out MultidimensionalArray Nodes, out int[] Type) {
            if (px < 2)
                throw new ArgumentOutOfRangeException("at least two nodes in each direction are required.");

            Nodes = MultidimensionalArray.Create(px * px * px, 3);
            Type = new int[Nodes.GetLength(0)];

            var Nodes1D = GenericBlas.Linspace(-1, 1, px);
            int cnt = 0;
            for (int i = 0; i < px; i++) {
                int i_edge = (i == 0 || i == px - 1) ? 1 : 0;

                for (int j = 0; j < px; j++) {
                    int j_edge = (j == 0 || j == px - 1) ? 1 : 0;

                    for (int k = 0; k < px; k++) {
                        int k_edge = (k == 0 || k == px - 1) ? 1 : 0;

                        Nodes[cnt, 0] = Nodes1D[i];
                        Nodes[cnt, 1] = Nodes1D[j];
                        Nodes[cnt, 2] = Nodes1D[k];

                        Type[cnt] = i_edge + j_edge + k_edge;

                        cnt++;
                    }
                }
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetInterpolationNodes_NonLin"/>
        /// </summary>
        override protected void GetInterpolationNodes_NonLin(CellType Type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex) {
            switch (Type) {
                case CellType.Cube_8: {
                        base.SelectNodalPolynomials(2, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_20: {
                        //MultidimensionalArray _InterpolationNodes;
                        //Polynomial[] _InterpolationPolynomials;
                        //int[] _NodeType;
                        //int[] _EntityIndex;

                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex, NodeTypeFilter: new int[] { 0, 1 },
                            ModalBasisSelector: delegate(Polynomial p) {

                            for (int i = 0; i < p.Coeff.Length; i++) {
                                int p1 = p.Exponents[i, 0];
                                int p2 = p.Exponents[i, 1];
                                int p3 = p.Exponents[i, 2];

                                if (p1 > 2 || p2 > 2 || p3 > 2)
                                    return false;

                                if (p1 == 2 && p2 == 2)
                                    return false;
                                if (p2 == 2 && p3 == 2)
                                    return false;
                                if (p1 == 2 && p3 == 2)
                                    return false;
                            }

                            return true;
                        });
                        Debug.Assert(NodeType.Length == 20);

                        //// should be 27 Nodes, so we have to drop seven:
                        //// it will be the volume node in the center of the cell and all face nodes

                        //int _K = _InterpolationNodes.GetLength(0);
                        //Debug.Assert(_K == 27);
                        //Debug.Assert(_NodeType[0] == 0 && _NodeType[1] == 1 && _NodeType[2] == 1 && _NodeType[3] == 1 && _NodeType[4] == 1 && _NodeType[5] == 1 && _NodeType[6] == 1, "first 7 nodes should be the volume node and 6 face nodes");
                        //int D = _InterpolationNodes.GetLength(1);
                        //Debug.Assert(D == 3, "spatial dimension is expected to be 3.");


                        //int K = 20;
                        //int offset = _K - K;
                        //InterpolationNodes = MultidimensionalArray.Create(K, D);
                        //InterpolationNodes.Set(_InterpolationNodes.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { _K - 1, D - 1 }));
                        //InterpolationPolynomials = new Polynomial[K];
                        //Array.Copy(_InterpolationPolynomials, offset, InterpolationPolynomials, 0, K);
                        //NodeType = new int[K];
                        //Array.Copy(_NodeType, offset, NodeType, 0, K);

                        //EntityIndex = new int[K];
                        //Array.Copy(_EntityIndex, offset, EntityIndex, 0, K);
                        return;
                    }
                case CellType.Cube_27: {
                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_64: {
                        base.SelectNodalPolynomials(4, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_125: {
                        base.SelectNodalPolynomials(5, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_216: {
                        base.SelectNodalPolynomials(6, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                //case CellType.Cube_343: {
                //        base.SelectNodalPolynomials(7, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_512: {
                //        base.SelectNodalPolynomials(8, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_729: {
                //        base.SelectNodalPolynomials(9, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_1000: {
                //        base.SelectNodalPolynomials(10, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementMapping"/>
        /// </summary>
        /// <param name="Type"></param>
        /// <param name="conv"></param>
        /// <returns></returns>
        public override int[] GetForeignElementMapping(CellType Type, RefElement.ExchangeFormats conv) {
            int[] permutationArray = new int[Vertices.GetLength(0)];
            for (int i = 0; i < permutationArray.Length; i++) {
                permutationArray[i] = i;
            }
            if (conv == ExchangeFormats.Gmsh || conv == ExchangeFormats.CGNS) {
                SwitchNode(ref permutationArray[2], ref permutationArray[3]);
                SwitchNode(ref permutationArray[3], ref permutationArray[4]);
                SwitchNode(ref permutationArray[5], ref permutationArray[6]);
                return permutationArray;
            }
            return permutationArray;
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementType"/>
        /// </summary>
        /// <param name="Type"></param>
        /// <param name="conv"></param>
        /// <param name="ForeignName"></param>
        /// <param name="ForeignTypeConstant"></param>
        public override void GetForeignElementType(CellType Type, RefElement.ExchangeFormats conv, out string ForeignName, out int ForeignTypeConstant) {
            ForeignName = "Hexagon";
            ForeignTypeConstant = 0;
            if (conv == ExchangeFormats.Gmsh) {
                if (Type == CellType.Cube_Linear) {
                    ForeignTypeConstant = 5;
                } else if (Type == CellType.Cube_27) {
                    ForeignTypeConstant = 12;
                } else if (Type == CellType.Cube_20) {
                    ForeignTypeConstant = 27;
                } else if (Type == CellType.Cube_64) {
                    ForeignTypeConstant = 92;
                } else if (Type == CellType.Cube_125) {
                    ForeignTypeConstant = 93;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.CGNS) {
                if (Type == 0) {
                    ForeignTypeConstant = 17;
                } else if (Type == CellType.Cube_20) {
                    ForeignTypeConstant = 18;
                } else if (Type == CellType.Cube_27) {
                    ForeignTypeConstant = 19;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.GambitNeutral) {
                ForeignName = "Brick";
                if (Type == CellType.Cube_Linear) {
                    ForeignTypeConstant = 4;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else {
                throw new NotSupportedException("Wrong foreign convention type");
            }
        }
    }
}
