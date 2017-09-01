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

using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// The point reference element, i.e. \f$ K^{\textrm{point}} = \{ 0 \} \f$.
    /// </summary>
    public class Point : RefElement {


        private static Point instance = null;
        private static readonly object padlock = new object();

        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static Point Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new Point();
                    }
                    return instance;
                }
            }
        }

        /// <summary>
        /// default constructor
        /// </summary>
        private Point() {
            using (new FuncTrace()) {
            // ===============
            // define vertices
            // ===============

            m_Vertices = new NodeSet(this, 1, 1);
            m_Vertices[0, 0] = 0.0;
            m_Vertices.LockForever();


            // ==============================
            // define orthonormal polynomials
            // ==============================

            OrthonormalPolynomials = null;

        }
        }


        /// <summary>
        /// not supported - a point has no edges;
        /// </summary>
        /// <param name="EdgeIndex"></param>
        /// <param name="EdgeVertices"></param>
        /// <param name="VolumeVertices"></param>
        public override void TransformFaceCoordinates(int EdgeIndex, MultidimensionalArray EdgeVertices, MultidimensionalArray VolumeVertices) {
            throw new NotSupportedException("A Point has no edges");
        }



        /// <summary>
        /// quadrature for a point is trivial - it's just evaluation a the point,
        /// so there is only one quadrature rule with one node that is valid for any order;
        /// </summary>
        /// <param name="DesiredOrder">
        /// any positive integer value;
        /// </param>
        /// <returns>
        /// An exact quad. rule for arbitrary functions in zero dimensions !!!
        /// </returns>
        public override QuadRule GetQuadratureRule(int DesiredOrder) {
            if (DesiredOrder < 0)
                throw new ArgumentOutOfRangeException("negative polynomial degree.");

            QuadRule qr = new QuadRule();
            qr.OrderOfPrecision = int.MaxValue;
            qr.Nodes = new NodeSet(this, 1, 1); 
            qr.Nodes[0, 0] = 0.0;
            qr.Nodes.LockForever();
            qr.Weights = MultidimensionalArray.Create(1);
            qr.Weights[0] = 1.0;

            
            return qr;
        }

        /// <summary>
        /// dimension of a point is 0;
        /// </summary>
        public override int SpatialDimension {
            get {
                return 0;
            }
        }

        /// <summary>
        /// the maximum positive integer
        /// </summary>
        public override int HighestKnownOrder {
            get {
                return Int32.MaxValue;
            }
        }

        /// <summary>
        /// the Point cannot be subdivided - the subdivision contains only one element;
        /// </summary>
        /// <returns></returns>
        public override AffineTrafo[] GetSubdivision() {
            AffineTrafo[] ret = new AffineTrafo[] { AffineTrafo.Identity(1) };
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
            if (pt[0] != 1)
                throw new ArgumentException("wrong spatial dimension");

            // Minimum tolerance determined experimentally by Florian
            tolerance = Math.Max(tolerance, 1e-8);
            if (Math.Abs(pt[0]) > tolerance)
                return false;
            else
                return true;
        }

        /// <summary>
        /// here, the same issues as for <see cref="GetQuadratureRule"/> apply;
        /// </summary>
        public override QuadRule GetBruteForceQuadRule(int NoOfSubDiv, int BaseRuleOrder) {
            return GetQuadratureRule(BaseRuleOrder);
        }

        /*
        /// <summary>
        /// See <see cref="RefElement.JacobianDetTransformation"/>.
        /// For the point simplex, always 0.
        /// </summary>
        public override void JacobianDetTransformation(MultidimensionalArray ReferenceVerticesIn, MultidimensionalArray JacobianDetOut, int jOffset, CellType MinorCellType, MultidimensionalArray TransformationParams) {
            int N = ReferenceVerticesIn.GetLength(0);
            int D = SpatialDimension;

            Debug.Assert(ReferenceVerticesIn.Dimension == 2);
            Debug.Assert(JacobianDetOut.Dimension == 2);
            Debug.Assert(ReferenceVerticesIn.GetLength(1) == 1, "wrong spatial dimension of ReferenceVerticesIn");
            Debug.Assert(ReferenceVerticesIn.GetLength(0) == JacobianDetOut.GetLength(1), "mismatch in number of vertices per cell.");
            Debug.Assert(TransformationParams.GetLength(0) == 1, "wrong spatial dimension of TransformationParams");

            if(MinorCellType != 0)
                throw new ArgumentException("unknown cell type.");

            JacobianDetOut.Clear();
        }
         */

        /// <summary>
        /// Operation not meaningful for zero-dimensional object.
        /// </summary>
        protected override void GetNodeSet(int px, out MultidimensionalArray Nodes, out int[] Type) {
            if (px != 0)
                throw new ArgumentOutOfRangeException();

            Nodes = MultidimensionalArray.Create(1, 0);
            Type = new int[1];
        }

        /// <summary>
        /// Not implemented.
        /// </summary>
        protected override void GetInterpolationNodes_NonLin(CellType Type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex) {
            throw new NotImplementedException();
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementType"/>
        /// </summary>
        public override void GetForeignElementType(CellType Type, RefElement.ExchangeFormats conv, out string ForeignName, out int ForeignTypeConstant) {
            ForeignName = "Point";
            ForeignTypeConstant = 0;

            if (Type != CellType.Point)
                throw new ArgumentException();

            if (conv == ExchangeFormats.Gmsh) {
                throw new NotSupportedException("Wrong minor cell type");
            } else if (conv == ExchangeFormats.CGNS) {
                ForeignTypeConstant = 2;
            } else if (conv == ExchangeFormats.GambitNeutral) {
                throw new NotSupportedException("Wrong minor cell type");
            } else {
                throw new NotSupportedException("Wrong foreign convention type");
            }
        }
    }
}
