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
using System.Runtime.InteropServices;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// A 2D Vector2D, or Stage 1 Tensor in 3D Space;
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct Vector3D {

        /// <summary>
        /// initializes to (<paramref name="__x"/>,<paramref name="__y"/><paramref name="__z"/>).
        /// </summary>
        /// <param name="__x">x - component</param>
        /// <param name="__y">y - component</param>
        /// <param name="__z">z - component</param>
        public Vector3D(double __x, double __y, double __z) {
            x = __x;
            y = __y;
            z = __z;
        }

        /// <summary>
        /// Constructs a new vector based on the given coordinates
        /// </summary>
        /// <param name="coordinates">
        /// The coordinates of the vector
        /// </param>
        public Vector3D(params double[] coordinates) {
            if (coordinates.Length != 3) {
                throw new ArgumentException("Invalid dimension");
            }

            x = coordinates[0];
            y = coordinates[1];
            z = coordinates[2];
        }

        /// <summary>
        /// x - component
        /// </summary>
        public double x;
        /// <summary>
        /// y - component
        /// </summary>
        public double y;
        /// <summary>
        /// z - component
        /// </summary>
        public double z;

        /// <summary>
        /// adds <paramref name="v"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        public void Acc(Vector3D v) {
            this.x += v.x;
            this.y += v.y;
            this.z += v.z;
        }

        /// <summary>
        /// adds <paramref name="v"/>*<paramref name="s"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        /// <param name="s">the scaleing of the vector to add</param>
        public void Acc(Vector3D v, double s) {
            this.x += v.x * s;
            this.y += v.y * s;
            this.z += v.z * s;
        }

        /// <summary>
        /// Calculates the cross product between this vector and the given
        /// vector <paramref name="v"/>.
        /// </summary>
        /// <param name="v">A vector</param>
        /// <returns>
        /// \f$ \vec{u} \times \vec{v}\f$ 
        /// </returns>
        public Vector3D CrossProduct(Vector3D v) {
            return new Vector3D() {
                x = this.y * v.z - this.z * v.y,
                y = this.z * v.x - this.x * v.z,
                z = this.x * v.y - this.y * v.x
            };
        }

        /// <summary>
        /// the absolute value (length) of this vector 
        /// </summary>
        /// <returns></returns>
        public double Abs() {
            return Math.Sqrt(x * x + y * y + z * z);
        }

        /// <summary>
        /// the squared magnitude of this vector
        /// </summary>
        /// <returns></returns>
        public double AbsSquare() {
            return x * x + y * y + z * z;
        }

        ///// <summary>
        ///// The angle between this vector and the positive x-Axis, counted counterclockwise
        ///// </summary>
        ///// <returns></returns>
        //public double Angle() {
        //    return Math.Atan2(this.y, this.x);
        //}

        ///// <summary>
        ///// sets (the cartesian) coordinates of this vector from polar coordinates
        ///// </summary>
        ///// <param name="r">Distance to origin.</param>
        ///// <param name="phi">angle to the positive x-Axis, counted counterclockwise</param>
        //public void FromPolar(double r, double phi) {
        //    this.x = Math.Cos(phi)*r;
        //    this.y = Math.Sin(phi)*r;
        //}


        /// <summary>
        /// normalizes this vector
        /// </summary>
        public void Normalize() {
            double l = 1.0 / Abs();
            x *= l;
            y *= l;
            z *= l;
        }

        /// <summary>
        /// subtracts <paramref name="v"/> from this vector;
        /// </summary>
        /// <param name="v">the vector that should be substracted</param>
        public void Sub(Vector3D v) {
            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
        }

        /// <summary>
        /// multiplies this vector with  scalar <paramref name="s"/>;
        /// </summary>
        /// <param name="s"></param>
        public void Scale(double s) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
        }

        /// <summary>
        /// sets the components of this object
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public void Set(double x, double y, double z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        /// <summary>
        /// Copies the components 0 to <paramref name="length"/> - 1 of this
        /// object into <paramref name="destination"/>
        /// </summary>
        /// <param name="destination">An array of length three</param>
        /// <param name="length">The number of elements to be copied</param>
        public void CopyTo(double[] destination, int length) {
            if (destination.Length < length) {
                throw new ArgumentException("The destination array is too small", "destination");
            }

            if (length > 3 || length < 1) {
                throw new ArgumentException("Length can only be 0, 1 or 2", "length");
            }

            for (int i = 0; i < length; i++) {
                destination[i] = this[i];
            }
        }

        /// <summary>
        /// Copies the components 0 to <paramref name="length"/> - 1 of this
        /// object into <paramref name="destination"/> starting with <paramref name="destinationIndex" />
        /// </summary>
        /// <param name="destination">The target array</param>
        /// <param name="length">The number of elements to be copied</param>
        /// <param name="destinationIndex">The index of the target array at which copying should be started</param>
        public void CopyTo(double[] destination, int destinationIndex, int length) {
            if (destination.Length < destinationIndex + length) {
                throw new ArgumentException("The destination array is too small", "destination");
            }

            if (length > 3 || length < 1) {
                throw new ArgumentException("Length can only be 0, 1 or 2", "length");
            }

            if (destinationIndex < 0) {
                throw new ArgumentException("The destination index must be >= 0", "destinationIndex");
            }

            for (int i = 0; i < length; i++) {
                destination[i + destinationIndex] = this[i];
            }
        }


        ///// <summary>
        ///// multiplyes this vector by a matrix from the left hand side,
        ///// i.e does "matrix <paramref name="M"/>" * "column vector";
        ///// </summary>
        ///// <param name="M"></param>
        //public void MulLeft(TensorSt2 M) {
        //    double xold = x;

        //    x = M._00*xold + M._01*y;
        //    y = M._10*xold + M._11*y;
        //}

        ///// <summary>
        ///// multiplyes this vector by a matrix from the right hand side,
        ///// i.e does "row vector" * "matrix <paramref name="M"/>";
        ///// </summary>
        ///// <param name="M"></param>
        //public void MulRight(TensorSt2 M) {
        //    double xold = x;

        //    x = M._00*xold + M._10*y;
        //    y = M._01*xold + M._11*y;
        //}



        /// <summary>
        /// standard vector addition
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector3D operator +(Vector3D left, Vector3D right) {
            Vector3D ret = left;
            ret.Acc(right);
            return ret;
        }

        /// <summary>
        /// standard vector addition
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector3D operator -(Vector3D left, Vector3D right) {
            Vector3D ret = left;
            ret.Sub(right);
            return ret;
        }

        /// <summary>
        /// multiplication by a scalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector3D operator *(Vector3D v, double s) {
            Vector3D ret = v;
            ret.Scale(s);
            return ret;
        }

        /// <summary>
        /// multiplication by a scalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector3D operator *(double s, Vector3D v) {
            return (v * s);
        }

        /// <summary>
        /// Division by a scalar
        /// </summary>
        /// <param name="v">The vector</param>
        /// <param name="s">The scalar</param>
        /// <returns>v / s</returns>
        public static Vector3D operator /(Vector3D v, double s) {
            Vector3D result = v;
            result.Scale(1.0 / s);
            return result;
        }

        /// <summary>
        /// Standard Dot/Vector2D/Inner-Product.
        /// </summary>
        /// <param name="a">1st operand</param>
        /// <param name="b">2nd operand</param>
        /// <returns>a*b*</returns>
        public static double operator *(Vector3D a, Vector3D b) {
            return (a.x * b.x + a.y * b.y + a.z * b.z);
        }

        /// <summary>
        /// Implicit conversion to double array of length 3.
        /// </summary>
        /// <param name="v">The vector to be converted.</param>
        /// <returns>
        /// Any array of length 3 with the entries
        /// [<see cref="x"/>, <see cref="y"/>, <see cref="z"/>]
        /// </returns>
        public static implicit operator double[](Vector3D v) {
            return new double[] { v.x, v.y, v.z };
        }

        /// <summary>
        /// a vector notation: (x|y|z);
        /// </summary>
        /// <returns>(x|y|z)</returns>
        public override string ToString() {
            return ("(" + x + "|" + y + "|" + z + ")");
        }

        /// <summary>
        /// the <paramref name="d"/>-th Standard basis
        /// </summary>
        /// <param name="d">Dimension index</param>
        /// <returns>
        /// (1,0,0) if <paramref name="d"/>=0, (0,1,0) if <paramref name="d"/>=1,
        /// (0,0,1) if <paramref name="d"/>=2;
        /// </returns>
        public static Vector3D StdBasis(int d) {
            if (d == 0)
                return new Vector3D(1, 0, 0);
            else if (d == 1)
                return new Vector3D(0, 1, 0);
            else if (d == 2)
                return new Vector3D(0, 0, 1);
            else
                throw new ArgumentOutOfRangeException("in 2D, a Dimension index must be either 0 or 1");
        }

        /// <summary>
        /// set/get entries
        /// </summary>
        /// <param name="i">either 0 (x-component) or 1 (y-component) or 2 (z-component) </param>
        /// <returns></returns>
        public double this[int i] {
            set {
                switch (i) {
                    case 0:
                        x = value;
                        return;
                    case 1:
                        y = value;
                        return;
                    case 2:
                        z = value;
                        return;
                    default:
                        throw new IndexOutOfRangeException("vector component index must be either 0 or 1 or 2.");
                }
            }
            get {
                switch (i) {
                    case 0:
                        return x;
                    case 1:
                        return y;
                    case 2:
                        return z;
                    default:
                        throw new IndexOutOfRangeException("vector component index must be either 0 or 1 or 2.");
                }
            }
        }

        /// <summary>
        /// Euclidean distance between the points <paramref name="a"/> and <paramref name="b"/>
        /// </summary>
        public static double Dist(Vector3D a, Vector3D b) {
            double distPow2 = 0;
            double d;
            d = (a.x - b.x);
            distPow2 += d * d;
            d = (a.y - b.y);
            distPow2 += d * d;
            d = (a.z - b.z);
            distPow2 += d * d;
            return Math.Sqrt(distPow2);
        }

        /// <summary>
        /// Computes the transformation for a arbitrary vector from the
        /// standard coordinate system into a coordinate system with the axes
        /// \f$  \vec{e}_1 = \vec{n}\f$ ,
        /// \f$  \vec{e}_2 = (-n_2, n_1, 0)^T\f$  and
        /// \f$  \vec{e}_3 = \vec{e}_1 \times \vec{e}_2\f$ ,
        /// where \f$  \vec{n} = (n_1, n_2, n_3)^T\f$ 
        /// is given by <paramref name="edgeNormal"/>.
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// An orthonormal 3x3 matrix that, when applied to a vector,
        /// transforms this vector into the above-mentioned coordinate system.
        /// </returns>
        private MultidimensionalArray GetTransformationToEdgeCoordinates(Vector3D edgeNormal) {
            Vector3D t1 = new Vector3D(-edgeNormal[1], edgeNormal[0], 0.0);
            Vector3D t2 = edgeNormal.CrossProduct(t1);
            MultidimensionalArray trafo = MultidimensionalArray.Create(3, 3);
            for (int i = 0; i < 3; i++) {
                trafo[0, i] = edgeNormal[i];
                trafo[1, i] = t1[i];
                trafo[2, i] = t2[i];
            }

            return trafo;
        }

        /// <summary>
        /// Transforms this vector into a new vector in a coordinate system whose
        /// first axis is aligned with the given <paramref name="edgeNormal"/>.
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// This vector in a coordinate system with the axes
        /// \f$  \vec{e}_1 = \vec{n}\f$ ,
        /// \f$  \vec{e}_2 = (-n_2, n_1, 0)^T\f$  and
        /// \f$  \vec{e}_3 = \vec{e}_1 \times \vec{e}_2\f$ ,
        /// where \f$  \vec{n} = (n_1, n_2, n_3)^T\f$ 
        /// is given by <paramref name="edgeNormal"/>
        /// </returns>
        public Vector3D ToEdgeCoordinates(Vector3D edgeNormal) {
            double[] transformedVector = new double[3];
            GetTransformationToEdgeCoordinates(edgeNormal).gemv(
                1.0, (double[])this, 0.0, transformedVector);
            return new Vector3D(transformedVector);
        }

        /// <summary>
        /// Transforms this vector from a coordinate system as defined by
        /// <see cref="ToEdgeCoordinates"/> into a new vector in the standard
        /// coordinate system
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// This vector in the standard coordinate system
        /// </returns>
        public Vector3D FromEdgeCoordinates(Vector3D edgeNormal) {
            double[] transformedVector = new double[3];
            GetTransformationToEdgeCoordinates(edgeNormal).Transpose().gemv(
                1.0, (double[])this, 0.0, transformedVector);
            return new Vector3D(transformedVector);
        }
    }
}
