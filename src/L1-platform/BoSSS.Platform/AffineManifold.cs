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
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// Affine, D-1 dimensional manifold M;
    /// x is in M, if <see cref="Normal"/>*x = <see cref="a"/>.
    /// </summary>
    public class AffineManifold : ICloneable {

        /// <summary>
        /// constructs an empty, degenerate manifold
        /// </summary>
        /// <param name="D"></param>
        public AffineManifold(int D) {
            if (D <= 1)
                throw new ArgumentException("1D or lower is not supported");
            Normal = new Vector(D);
        }

        /// <summary>
        /// constructs an affine manifold from normal vector and offset point
        /// </summary>
        /// <param name="_Normal">normal onto the affine manifold</param>
        /// <param name="_Offset">one point in the plane, can be null</param>
        public AffineManifold(double[] _Normal, double[] _Offset) 
            : this(new Vector(_Normal), new Vector(_Offset != null ? _Offset : new double[_Normal.Length]))
        {
#if DEBUG
            Debug.Assert(_Normal.Length == this.Normal.Dim);
            for(int d = 0; d < _Normal.Length; d++) {
                Debug.Assert(_Normal[d] == this.Normal[d]);
            }
#endif
        }

        /// <summary>
        /// constructs an affine manifold from normal vector and offset point
        /// </summary>
        /// <param name="_Normal">normal onto the affine manifold</param>
        /// <param name="_Offset">one point in the plane, can be null</param>
        public AffineManifold(Vector _Normal, Vector _Offset) {
            if (_Normal.Dim != _Offset.Dim)
                throw new ArgumentException();
            if (_Normal.Dim <= 1)
                throw new ArgumentException("1D or lower is not supported");

            Normal = _Normal;
            if(Normal.AbsSquare() <= 0.0)
                throw new ArgumentException();

            int D = Normal.Dim;
            for (int d = 0; d < D; d++) {
                this.a += Normal[d] * _Offset[d];
            }

        }

        /// <summary>
        /// constructs an affine manifold from a set of points
        /// </summary>
        /// <param name="pts">
        /// points on the manifold;
        /// 1st index: point index.
        /// </param>
        static public AffineManifold FromPoints(params Vector[] pts) {
            int D = pts.Length;
            if (D <= 1)
                throw new ArgumentException("1D or lower is not supported");
            for (int i = 0; i < D; i++)
                if (pts[i].Dim != D)
                    throw new ArgumentException();
            var ret = new AffineManifold(D);

            switch (D) {
                case 2:
                ret.Normal[0] = -(pts[1][1] - pts[0][1]);
                ret.Normal[1] = +(pts[1][0] - pts[0][0]);

                break;
                case 3:
                double a_1 = pts[1][0] - pts[0][0];
                double a_2 = pts[1][1] - pts[0][1];
                double a_3 = pts[1][2] - pts[0][2];
                double b_1 = pts[2][0] - pts[0][0];
                double b_2 = pts[2][1] - pts[0][1];
                double b_3 = pts[2][2] - pts[0][2];

                ret.Normal[0] = a_2 * b_3 - a_3 * b_2;
                ret.Normal[1] = a_3 * b_1 - a_1 * b_3;
                ret.Normal[2] = a_1 * b_2 - a_2 * b_1;

                break;
                default:
                throw new NotSupportedException();

            }
            ret.a = ret.Normal * pts[0];
            if(ret.Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();

            return ret;
        }



        /// <summary>
        /// constructs an affine manifold from a set of points
        /// </summary>
        /// <param name="pts">
        /// points on the manifold;
        /// 1st index: point index.
        /// 2nd index: spatial dimension.
        /// </param>
        static public AffineManifold FromPoints(double[,] pts) {
            var _pts = MultidimensionalArray.Create(pts.GetLength(0), pts.GetLength(1));
            _pts.Set2DArray(pts);
            return FromPoints(_pts);
        }


        /// <summary>
        /// constructs an affine manifold from a set of points
        /// </summary>
        /// <param name="pts">
        /// points on the manifold;
        /// 1st index: point index.
        /// 2nd index: spatial dimension.
        /// </param>
        static public AffineManifold FromPoints(MultidimensionalArray pts) {
            int D = pts.GetLength(1);
            var ret = new AffineManifold(D);

            switch (D) {
                case 2:
                ret.Normal[0] = -(pts[1, 1] - pts[0, 1]);
                ret.Normal[1] = +(pts[1, 0] - pts[0, 0]);

                break;
                case 3:
                double a_1 = pts[1, 0] - pts[0, 0];
                double a_2 = pts[1, 1] - pts[0, 1];
                double a_3 = pts[1, 2] - pts[0, 2];
                double b_1 = pts[2, 0] - pts[0, 0];
                double b_2 = pts[2, 1] - pts[0, 1];
                double b_3 = pts[2, 2] - pts[0, 2];

                ret.Normal[0] = a_2 * b_3 - a_3 * b_2;
                ret.Normal[1] = a_3 * b_1 - a_1 * b_3;
                ret.Normal[2] = a_1 * b_2 - a_2 * b_1;

                break;
                default:
                throw new NotSupportedException();

            }
            double[] Offset = pts.ExtractSubArrayShallow(0, -1).To1DArray();
            ret.a = ret.Normal * (new Vector(Offset));
            if(ret.Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();

            return ret;
        }



        /// <summary>
        /// surface normal
        /// </summary>
        public Vector Normal;

        /// <summary>
        /// affine offset (result of scalar product of <see cref="Normal"/>*x, with x being an arbitrary point in the surface)
        /// </summary>
        public double a;

        /// <summary>
        /// the distance of some point to the affine manifold/plane, in multiples of the norm of <see cref="Normal"/>;
        /// </summary>
        public double PointDistance(params double[] pt) {
            if (pt.Length != this.Normal.Dim)
                throw new ArgumentException();
            if(Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();

            double ret = 0;
            for (int d = this.Normal.Dim - 1; d >= 0; d--) {
                ret += pt[d]*this.Normal[d];
            }
            ret -= this.a;

            return ret;
        }

        /// <summary>
        /// the distance of some point to the affine manifold/plane, in multiples of the norm of <see cref="Normal"/>;
        /// </summary>
        public double PointDistance(Vector pt) {
            if (pt.Dim != this.Normal.Dim)
                throw new ArgumentException();
            if (Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();

            double ret = 0;
            for (int d = this.Normal.Dim - 1; d >= 0; d--) {
                ret += pt[d] * this.Normal[d];
            }
            ret -= this.a;

            return ret;
        }

        /// <summary>
        /// returns the point in this affine manifold (plane) which is closest to <paramref name="pt"/>
        /// </summary>
        public Vector ProjectPoint(Vector pt) {
            if (pt.Dim != this.Normal.Dim)
                throw new ArgumentException();
            if (Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();


            double dist = 0;
            for (int d = this.Normal.Dim - 1; d >= 0; d--) {
                dist += pt[d] * this.Normal[d];
            }
            dist -= this.a;

            Vector cp = pt;
            cp.Acc(Normal, -dist);

            return cp;
        }



        /// <summary>
        /// equality check: tests if two objects represent the same affine manifold/plane
        /// </summary>
        public bool Equals(object _other, double RelTol) {
            var other = _other as AffineManifold;
            if (other == null)
                throw new ArgumentException();
            if (other.Normal.Dim != this.Normal.Dim)
                throw new ArgumentException();

            double l1 = this.Normal.AbsSquare();
            double l2 = other.Normal.AbsSquare();
            double Tol = (l1 + l2) * RelTol;

            double inner_prod =   this.Normal * other.Normal;

            double cos_alpha = inner_prod / Math.Sqrt(l1*l2);
            if (Math.Abs(Math.Abs(cos_alpha) - 1) > 1.0E-12)
                // Normals not parallel -> manifolds can't be equal
                return false;


            //double offset = inner_prod/(this.a*other.a);
            double defect = inner_prod - (other.a / this.a) * l1;
            if (Math.Abs(defect) > Tol)
                // manifolds don't intersect
                return false;

            return true;
        }

        /// <summary>
        /// equality check: tests if two objects represent the same affine manifold/plane
        /// </summary>
        public override bool Equals(object _other) {
            return this.Equals(_other, 1.0e-12);
        }

        /// <summary>
        /// hash
        /// </summary>
        public override int GetHashCode() {
            return (int)((a + (this.Normal.AbsSquare()))*1000.0);
        }

        /// <summary>
        /// non-shallow cloning
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            var ret = new AffineManifold(this.Normal.Dim);
            ret.Normal = this.Normal;
            ret.a = this.a;
            return ret;
        }

        /// <summary>
        /// Intersection of two 1D affine manifolds in 2D space
        /// </summary>
        static public Vector Intersect2D(AffineManifold A, AffineManifold B) {
            if (A.Normal.Dim != 2)
                throw new ArgumentException();
            if (B.Normal.Dim != 2)
                throw new ArgumentException();
            if(A.Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();
            if(B.Normal.AbsSquare() <= 0.0)
                throw new ArithmeticException();

            double n1x = A.Normal[0];
            double n1y = A.Normal[1];
            double a1 = A.a;
            double n2x = B.Normal[0];
            double n2y = B.Normal[1];
            double a2 = B.a;

            return new Vector(
                (-n1y*a2+a1*n2y)/(n1x*n2y-n2x*n1y),
               -(-n1x*a2+n2x*a1)/(n1x*n2y-n2x*n1y)
               );
        }
    }
}

