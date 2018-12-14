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
using System.Text;
using System.Runtime.InteropServices;
using System.Diagnostics;
using ilPSP.Utils;
using ilPSP;

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// A spatial coordinate or vector, in 1D, 2D, 3D
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct Vector {

        /// <summary>
        /// initializes a <paramref name="D"/>-dimensional vector.
        /// </summary>
        public Vector(int D) {
            if (D < 1 || D > 3) {
                throw new ArgumentException("supports only 1D, 2D or 3D");
            }
            x = 0;
            y = 0;
            z = 0;
            Dim = D;
            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// initializes a 1D vector.
        /// </summary>
        public Vector(double __x) {
            x = __x;
            Dim = 1;
            y = 0;
            z = 0;
            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// initializes a 2D vector.
        /// </summary>
        public Vector(double __x, double __y) {
            x = __x;
            y = __y;
            Dim = 2;
            z = 0;
            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// initializes a 3D vector.
        /// </summary>
        public Vector(double __x, double __y, double __z) {
            x = __x;
            y = __y;
            z = __z;
            Dim = 3;
            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// initializes a vector from an array
        /// </summary>
        /// <param name="X"></param>
        public Vector(double[] X) {
            if (X.Length < 1 || X.Length > 3) {
                throw new ArgumentException();
            }

            this.Dim = X.Length;
            x = X[0];
            if(this.Dim > 1)
                y = X[1];
            else
                y = 0;

            if(this.Dim > 2)
                z = X[2];
            else
                z = 0;

            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// x - component
        /// </summary>
        public double x;
        /// <summary>
        /// y - component (used for dimensions >= 2)
        /// </summary>
        public double y;

        /// <summary>
        /// z - component (used for dimensions >= 3)
        /// </summary>
        public double z;

        /// <summary>
        /// spatial dimension
        /// </summary>
        public int Dim;

        /// <summary>
        /// Dummy entry/reserved; enforces the entire structure to be 256 bit in size.
        /// </summary>
        public int Dummy_256bitAlign;


        /// <summary>
        /// set/get entries
        /// </summary>
        /// <param name="i">either 0 (x-component) or 1 (y-component)</param>
        /// <returns></returns>
        public double this[int i] {
            set {
                if(i < 0 || i >= this.Dim)
                    throw new IndexOutOfRangeException("vector component out of range.");

                if (i == 0)
                    x = value;
                else if (i == 1 && Dim > 1)
                    y = value;
                else if (i == 2 && Dim > 2)
                    z = value;
                else
                    throw new IndexOutOfRangeException("vector component index must be either 0 or 1.");
            }
            get {
                if (i < 0 || i >= this.Dim)
                    throw new IndexOutOfRangeException("vector component out of range.");

                if (i == 0)
                    return x;
                else if(i == 1 && Dim > 1)
                    return y;
                else if(i == 2 && Dim > 2)
                    return z;
                else
                    throw new IndexOutOfRangeException("vector component index must be either 0 or 1.");

            }
        }


        /// <summary>
        /// adds <paramref name="v"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        public void Acc(Vector v) {
            if(v.Dim != this.Dim)
                throw new ArgumentException("Dimension mismatch");

            this.x += v.x;
            this.y += v.y;
            this.z += v.z;
        }

        /// <summary>
        /// adds <paramref name="v"/>*<paramref name="s"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        /// <param name="s">the scaling of the vector to add</param>
        public void Acc(Vector v, double s) {
             if(v.Dim != this.Dim)
                throw new ArgumentException("Dimension mismatch");

            this.x += v.x * s;
            this.y += v.y * s;
            this.z += v.z * s;
        }

        /// <summary>
        /// Calculates the cross product between this vector and the given vector <paramref name="v"/>.
        /// </summary>
        /// <param name="v">A vector</param>
        /// <returns>
        /// \f$ \vec{u} \times \vec{v}\f$ 
        /// Note
        /// - return value will always be a 3D vector (<see cref="Dim"/>==3)
        /// - for 2D arguments, only the <see cref="z"/> component will be unequal 0.
        /// </returns>
        public Vector CrossProduct(Vector v) {
            if(v.Dim != this.Dim)
                throw new ArgumentException("Dimension mismatch");

            var R = new Vector(
                __x: this.y * v.z - this.z * v.y,
                __y: this.z * v.x - this.x * v.z,
                __z: this.x * v.y - this.y * v.x
            );
            Debug.Assert(R.Dim == 3);
            return R;
        }

        /// <summary>
        /// z-component of the cross product between 2D-vectors with zero z-component.
        /// </summary>
        public double CrossProduct2D(Vector b) {
            if (this.Dim != 2)
                throw new NotSupportedException("Only supported for 2D vectors.");
            if (b.Dim != 2)
                throw new ArgumentException("Only supported for 2D vectors.");
            return this[0] * b[1] - this[1] * b[0];
        }


        /// <summary>
        /// the absolute value (length) of this vector  (synonym for <see cref="L2Norm"/>)
        /// </summary>
        public double Abs() {
            return Math.Sqrt(x * x + y * y + z * z);
        }

        /// <summary>
        /// the absolute distance (length) to some other point <paramref name="o"/> 
        /// </summary>
        public double Dist(Vector o) {
            return (this - o).L2Norm();
        }
        
        /// <summary>
        /// the absolute value (length) of this vector (synonym for <see cref="Abs"/>)
        /// </summary>
        public double L2Norm() {
            return Math.Sqrt(x * x + y * y + z * z);
        }
        
        /// <summary>
        /// the absolute value (length) of this vector to the power of two
        /// </summary>
        public double AbsSquare() {
            return (x * x + y * y + z * z);
        }

        /// <summary>
        /// The angle (in radians) between this vector and the positive x-Axis, counted counterclockwise
        /// </summary>
        public double Angle2D() {
            if(this.Dim != 2)
                throw new NotSupportedException();
            return Math.Atan2(this.y, this.x);
        }

        /// <summary>
        /// The angle, in radians, between this vector and vector <paramref name="o"/>
        /// </summary>
        public double AngleTo(Vector o) {
            if (o.Dim != this.Dim)
                throw new ArgumentException("spatial dimension mismatch");
            if (o.AbsSquare() <= 0)
                throw new ArithmeticException("other vector is zero - unable to determine angle");
            if (this.AbsSquare() <= 0)
                throw new ArithmeticException("this vector is zero - unable to determine angle");

            Vector tn = this;
            tn.Normalize();

            Vector on = o;
            on.Normalize();

            double inner = tn * on;

            Debug.Assert(inner <= 1.0 + BLAS.MachineEps.Sqrt());

            // clamp value to range [-1..+1]: e.g. an inner product of 1.00000000002 could cause an NAN in Math.Acos
            inner = Math.Max(-1.0, inner); 
            inner = Math.Min(+1.0, inner);

            double angle =  Math.Acos(inner);
            return angle;
        }

        /// <summary>
        /// sets (the Cartesian) coordinates of this vector from polar coordinates
        /// </summary>
        /// <param name="r">Distance to origin.</param>
        /// <param name="phi">angle to the positive x-Axis, counted counterclockwise</param>
        public void FromPolar(double r, double phi) {
            if(this.Dim < 2)
                throw new NotSupportedException();
            this.x = Math.Cos(phi) * r;
            this.y = Math.Sin(phi) * r;
        }


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
        /// <param name="v">the vector that should be subtracted</param>
        public void Sub(Vector v) {

            if(v.Dim != this.Dim)
                throw new ArgumentException("Dimension mismatch");

            this.x -= v.x;
            this.y -= v.y;
            this.z -= v.z;
        }

        /// <summary>
        /// multiples this vector with  scalar <paramref name="s"/>;
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
        public void Set(double x, double y) {
            if(this.Dim != 2)
                throw new NotSupportedException("spatial dimension mismatch");
            this.x = x;
            this.y = y;
        }

        /// <summary>
        /// sets the components of this object
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public void Set(double x, double y, double z) {
            if(this.Dim != 3)
                throw new NotSupportedException("spatial dimension mismatch");
            this.x = x;
            this.y = y;
            this.z = z;
        }


       

        /// <summary>
        /// standard vector addition
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector operator +(Vector left, Vector right) {
            Vector ret = left;
            ret.Acc(right);
            return ret;
        }

        /// <summary>
        /// standard vector subtraction
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector operator -(Vector left, Vector right) {
            Vector ret = left;
            ret.Sub(right);
            return ret;
        }

        /// <summary>
        /// multiplication by a scalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector operator *(Vector v, double s) {
            Vector ret = v;
            ret.Scale(s);
            return ret;
        }

        /// <summary>
        /// multiplication by a scalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector operator *(double s, Vector v) {
            return (v * s);
        }

        /// <summary>
        /// Division by a scalar
        /// </summary>
        /// <param name="v">The vector</param>
        /// <param name="s">The scalar</param>
        /// <returns>v / s</returns>
        public static Vector operator /(Vector v, double s) {
            Vector result = v;
            result.Scale(1.0 / s);
            return result;
        }


        /// <summary>
        /// standard Dot/Vector/Inner-Product.
        /// </summary>
        /// <param name="a">1st operand</param>
        /// <param name="b">2nd operand</param>
        /// <returns>a*b*</returns>
        public static double operator *(Vector a, Vector b) {
             if(a.Dim != b.Dim)
                throw new ArgumentException("Dimension mismatch");

            return (a.x * b.x + a.y * b.y + a.z * b.z);
        }

        /// <summary>
        /// Implicit conversion an array of doubles of length 2
        /// </summary>
        /// <param name="v">The vector to be converted</param>
        /// <returns>An array of doubles of length 2</returns>
        public static implicit operator double[] (Vector v) {
            switch(v.Dim) {
                case 1:
                return new double[] { v.x };
                case 2:
                return new double[] { v.x, v.y };
                case 3:
                return new double[] { v.x, v.y, v.z };
                default:
                throw new NotSupportedException();
            }
        }

        /// <summary>
        /// a vector notation: (x|y);
        /// </summary>
        /// <returns>(x|y)</returns>
        public override string ToString() {
            switch(this.Dim) {
                case 1:
                return ("(" + x + ")");
                case 2:
                return ("(" + x + "|" + y + ")");
                case 3:
                return ("(" + x + "|" + y + "|" + z + ")");
                default:
                throw new NotSupportedException();
            }
        }

        /// <summary>
        /// the <paramref name="d"/>-th standard basis
        /// </summary>
        /// <param name="D">spatial dimension</param>
        /// <param name="d">direction index</param>
        /// <returns>
        /// a <paramref name="D"/>-dimensional vector with <paramref name="d"/>-th coordinate equal to 1.0.
        /// </returns>
        public static Vector StdBasis(int d, int D) {
            Vector e = default(Vector);
            e.Dim = D;
            e[d] = 1.0;
            return e;
        }

        /// <summary>
        /// the <paramref name="d"/>-th standard basis in 2D
        /// </summary>
        /// <param name="d">direction index</param>
        /// <returns>
        /// a two-dimensional vector with <paramref name="d"/>-th coordinate equal to 1.0.
        /// </returns>
        public static Vector StdBasis2D(int d, int D) {
            return StdBasis(d, 2);
        }

        /// <summary>
        /// the <paramref name="d"/>-th standard basis in 3D
        /// </summary>
        /// <param name="d">direction index</param>
        /// <returns>
        /// a three-dimensional vector with <paramref name="d"/>-th coordinate equal to 1.0.
        /// </returns>
        public static Vector StdBasis3D(int d, int D) {
            return StdBasis(d, 3);
        }

        /// <summary>
        /// Euclidean distance between the points <paramref name="a"/> and <paramref name="b"/>
        /// </summary>
        public static double Dist(Vector a, Vector b) {
             if(a.Dim != b.Dim)
                throw new ArgumentException("Dimension mismatch");

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
    }

    /// <summary>
    /// Extension methods for <see cref="Vector"/>
    /// </summary>
    public static class VectorExtensions {

        /// <summary>
        /// extracts the <paramref name="RowNo"/>-th row from
        /// <paramref name="inp"/>.
        /// </summary>
        /// <param name="inp">
        /// input matrix
        /// </param>
        /// <param name="RowNo">
        /// row which should be extracted
        /// </param>
        /// <returns>
        /// A vector with dimension (<see cref="Vector.Dim"/>) equal to 2nd length of <paramref name="inp"/>, containing the
        /// <paramref name="RowNo"/>-th row of <paramref name="inp"/>
        /// </returns>
        public static Vector GetRowPt(this IMatrix inp, int RowNo) {
            switch(inp.NoOfCols) {
                case 1:
                return new Vector(inp[RowNo, 0]);
                case 2:
                return new Vector(inp[RowNo, 0], inp[RowNo, 1]);
                case 3:
                return new Vector(inp[RowNo, 0], inp[RowNo, 1], inp[RowNo, 2]);
                default:
                throw new ArgumentException("Matrix has " + inp.NoOfCols + " columns, this cannot be a spatial dimension.");
            }
        }


        /// <summary>
        /// sets the <paramref name="RowNo"/>-th row from <paramref name="inp"/> to values provided by <paramref name="row"/>.
        /// </summary>
        /// <param name="inp">
        /// matrix that should be altered
        /// </param>
        /// <param name="RowNo">
        /// row index of the row to set
        /// </param>
        /// <param name="row">
        /// a vector 
        /// </param>
        public static void SetRowPt(this IMatrix inp, int RowNo, Vector row) {
            if (row.Dim != inp.NoOfCols)
                throw new ArgumentException("Dimension mismatch.");

            switch(row.Dim) {
                case 1:
                inp[RowNo, 0] = row.x; return;
                case 2:
                inp[RowNo, 0] = row.x; inp[RowNo, 1] = row.y; return;
                case 3:
                inp[RowNo, 0] = row.x; inp[RowNo, 1] = row.y; inp[RowNo, 2] = row.z; return;
                default:
                throw new NotImplementedException();
            }
        }
    }

}
