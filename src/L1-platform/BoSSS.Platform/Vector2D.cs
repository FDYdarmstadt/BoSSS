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

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// A 2D Vector2D, or Stage 1 Tensor in 2D Space;
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct Vector2D {

        /// <summary>
        /// initializes to (<paramref name="__x"/>,<paramref name="__y"/>).
        /// </summary>
        /// <param name="__x">x - component</param>
        /// <param name="__y">y - component</param>
        public Vector2D(double __x, double __y) {
            x = __x;
            y = __y;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="X"></param>
        public Vector2D(double[] X) {
            if (X.Length != 2) {
                throw new ArgumentException();
            }

            x = X[0];
            y = X[1];
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
        /// set/get entries
        /// </summary>
        /// <param name="i">either 0 (x-component) or 1 (y-component)</param>
        /// <returns></returns>
        public double this[int i] {
            set {
                if (i == 0)
                    x = value;
                else if (i == 1)
                    y = value;
                else
                    throw new IndexOutOfRangeException("vector component index must be either 0 or 1.");
            }
            get {
                if (i == 0)
                    return x;
                else if (i == 1)
                    return y;
                else
                    throw new IndexOutOfRangeException("vector component index must be either 0 or 1.");
            }
        }


        /// <summary>
        /// adds <paramref name="v"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        public void Acc(Vector2D v) {
            this.x += v.x;
            this.y += v.y;
        }

        /// <summary>
        /// adds <paramref name="v"/>*<paramref name="s"/> to this vector;
        /// </summary>
        /// <param name="v">the vector that should be added</param>
        /// <param name="s">the scaleing of the vector to add</param>
        public void Acc(Vector2D v, double s) {
            this.x += v.x * s;
            this.y += v.y * s;
        }

        /// <summary>
        /// the absolut value (length) of this vector 
        /// </summary>
        /// <returns></returns>
        public double Abs() {
            return Math.Sqrt(x * x + y * y);
        }

        /// <summary>
        /// The angle between this vector and the positive x-Axis, counted counterclockwise
        /// </summary>
        /// <returns></returns>
        public double Angle() {
            return Math.Atan2(this.y, this.x);
        }

        /// <summary>
        /// sets (the cartesian) coordinates of this vector from polar coordinates
        /// </summary>
        /// <param name="r">Distance to origin.</param>
        /// <param name="phi">angle to the positive x-Axis, counted counterclockwise</param>
        public void FromPolar(double r, double phi) {
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
        }

        /// <summary>
        /// subtracts <paramref name="v"/> from this vector;
        /// </summary>
        /// <param name="v">the vector that should be substracted</param>
        public void Sub(Vector2D v) {
            this.x -= v.x;
            this.y -= v.y;
        }

        /// <summary>
        /// multiplyes this vector with  scalar <paramref name="s"/>;
        /// </summary>
        /// <param name="s"></param>
        public void Scale(double s) {
            this.x *= s;
            this.y *= s;
        }

        /// <summary>
        /// sets the components of this object
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        public void Set(double x, double y) {
            this.x = x;
            this.y = y;
        }

        /// <summary>
        /// multiplyes this vector by a matrix from the left hand side,
        /// i.e does "matrix <paramref name="M"/>" * "column vector";
        /// </summary>
        /// <param name="M"></param>
        public void MulLeft(TensorSt2 M) {
            double xold = x;

            x = M._00 * xold + M._01 * y;
            y = M._10 * xold + M._11 * y;
        }

        /// <summary>
        /// multiplyes this vector by a matrix from the right hand side,
        /// i.e does "row vector" * "matrix <paramref name="M"/>";
        /// </summary>
        /// <param name="M"></param>
        public void MulRight(TensorSt2 M) {
            double xold = x;

            x = M._00 * xold + M._10 * y;
            y = M._01 * xold + M._11 * y;
        }



        /// <summary>
        /// standart vector addition
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector2D operator +(Vector2D left, Vector2D right) {
            Vector2D ret = left;
            ret.Acc(right);
            return ret;
        }

        /// <summary>
        /// standart vector addition
        /// </summary>
        /// <param name="left">op1</param>
        /// <param name="right">op2</param>
        /// <returns>clear;</returns>
        public static Vector2D operator -(Vector2D left, Vector2D right) {
            Vector2D ret = left;
            ret.Sub(right);
            return ret;
        }

        /// <summary>
        /// multiplication by a acalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector2D operator *(Vector2D v, double s) {
            Vector2D ret = v;
            ret.Scale(s);
            return ret;
        }

        /// <summary>
        /// multiplication by a acalar
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="v">the vector</param>
        /// <returns>clear;</returns>
        public static Vector2D operator *(double s, Vector2D v) {
            return (v * s);
        }


        /// <summary>
        /// Standart Dot/Vector2D/Inner-Product.
        /// </summary>
        /// <param name="a">1st operand</param>
        /// <param name="b">2nd operand</param>
        /// <returns>a*b*</returns>
        public static double operator *(Vector2D a, Vector2D b) {
            return (a.x * b.x + a.y * b.y);
        }

        /// <summary>
        /// Implicit conversion an array of doubles of length 2
        /// </summary>
        /// <param name="v">The vector to be converted</param>
        /// <returns>An array of doubles of length 2</returns>
        public static implicit operator double[] (Vector2D v) {
            return new double[] { v[0], v[1] };
        }

        /// <summary>
        /// a vector notation: (x|y);
        /// </summary>
        /// <returns>(x|y)</returns>
        public override string ToString() {
            return ("(" + x + "|" + y + ")");
        }

        /// <summary>
        /// the <paramref name="d"/>-th Standart basis
        /// </summary>
        /// <param name="d">Dimension index</param>
        /// <returns>(1,0) if <paramref name="d"/>=0, (0,1) if <paramref name="d"/>=1;</returns>
        public static Vector2D StdBasis(int d) {
            if (d == 0)
                return new Vector2D(1, 0);
            else if (d == 1)
                return new Vector2D(0, 1);
            else
                throw new ArgumentOutOfRangeException("in 2D, a Dimension index must be either 0 or 1");
        }

    }



    /// <summary>
    /// A stage 2 - Tensor in 2D - Space, or Matrix
    /// </summary>
    public struct TensorSt2 {
        /// <summary> entry in first row, first column </summary>
        public double _00;
        /// <summary> entry in first row, second column </summary>
        public double _01;
        /// <summary> entry in second row, first column </summary>
        public double _10;
        /// <summary> entry in second row, second column </summary>
        public double _11;

        /// <summary>
        /// matrix times column vector
        /// </summary>
        /// <param name="M">Matrix</param>
        /// <param name="v">vector</param>
        /// <returns><paramref name="M"/>*<paramref name="v"/></returns>
        public static Vector2D operator *(TensorSt2 M, Vector2D v) {
            v.MulLeft(M);
            return v;
        }

        /// <summary>
        /// row vector times matrix
        /// </summary>
        /// <param name="M">Matrix</param>
        /// <param name="v">vector</param>
        /// <returns><paramref name="v"/>*<paramref name="M"/></returns>
        public static Vector2D operator *(Vector2D v, TensorSt2 M) {
            v.MulRight(M);
            return v;
        }
    }

}
