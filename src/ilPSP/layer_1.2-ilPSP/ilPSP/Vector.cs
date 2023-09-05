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
using System.Collections;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using Newtonsoft.Json.Linq;
using Newtonsoft.Json;

namespace ilPSP {

    /// <summary>
    /// A spatial coordinate or vector, in 1D, 2D, 3D
    /// </summary>
    [Serializable] 
    [StructLayout(LayoutKind.Sequential)]
    [DataContract]
    public struct Vector : IList<double> {

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
            if (this.Dim > 1)
                y = X[1];
            else
                y = 0;

            if (this.Dim > 2)
                z = X[2];
            else
                z = 0;

            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// initializes a vector some part of an array
        /// </summary>
        /// <param name="X"></param>
        /// <param name="D">spatial dimension (<see cref="Dim"/>), number of entries to take form <paramref name="X"/></param>
        /// <param name="offset">offset into <paramref name="X"/></param>
        public Vector(double[] X, int offset, int D) {
            if (D < 1 || D > 3) {
                throw new ArgumentException();
            }

            this.Dim = D;
            x = X[offset + 0];
            if (this.Dim > 1)
                y = X[offset + 1];
            else
                y = 0;

            if (this.Dim > 2)
                z = X[offset + 2];
            else
                z = 0;

            Dummy_256bitAlign = 0;
        }

        /// <summary>
        /// x - component
        /// </summary>
        [DataMember]
        public double x;

        /// <summary>
        /// y - component (used for dimensions >= 2)
        /// </summary>
        [DataMember]
        public double y;

        /// <summary>
        /// z - component (used for dimensions >= 3)
        /// </summary>
        [DataMember]
        public double z;

        /// <summary>
        /// spatial dimension
        /// </summary>
        [DataMember]
        public int Dim;

        /// <summary>
        /// Dummy entry/reserved; enforces the entire structure to be 256 bit in size.
        /// </summary>
        public int Dummy_256bitAlign;

        /// <summary>
        /// equal to <see cref="Dim"/>
        /// </summary>
        public int Count {
            get {
                return Dim;
            }
        }

        /// <summary>
        /// always true, according to interface definition: elements can be changed, but no elements can be added/removed
        /// </summary>
        public bool IsReadOnly {
            get {
                return true;
            }
        }


        /// <summary>
        /// set/get entries
        /// </summary>
        /// <param name="i">either 0 (x-component) or 1 (y-component)</param>
        /// <returns></returns>
        public double this[int i] {
            //[MethodImpl(MethodImplOptions.AggressiveInlining)]
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
                
                //unsafe {
                //    fixed(Vector* d = &this) {
                //        ((double*)d)[i] = value;
                //    }
                //}
            }
            //[MethodImpl(MethodImplOptions.AggressiveInlining)]
            get {
                if(i < 0 || i >= this.Dim)
                    throw new IndexOutOfRangeException("vector component out of range.");

                if (i == 0)
                    return x;
                else if(i == 1 && Dim > 1)
                    return y;
                else if(i == 2 && Dim > 2)
                    return z;
                else
                    throw new IndexOutOfRangeException("vector component index must be either 0 or 1.");
                    
                //unsafe {
                //    fixed(Vector* d = &this) {
                //        return ((double*)d)[i];
                //    }
                //}

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
        /// sets all entries to 0.0
        /// </summary>
        public void Clearentries() {
            this.x = 0.0;
            this.y = 0.0;
            this.z = 0.0;
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
                throw new NotSupportedException("Only supported for 2D vectors.");
            return Math.Atan2(this.y, this.x);
        }

        /// <summary>
        /// Rotates a 2D vector in the xy-plane.
        /// </summary>
        public Vector Rotate2D(double angle) {
            if(this.Dim != 2)
                throw new NotSupportedException("Only supported for 2D vectors.");

            double cs = Math.Cos(angle);
            double ss = Math.Sin(angle);

            Vector r = default(Vector);
            r.Dim = 2;
            r.x = cs * this.x - ss * this.y;
            r.y = ss * this.x + cs * this.y;
            return r;
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
            tn.NormalizeInPlace();

            Vector on = o;
            on.NormalizeInPlace();

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
        /// normalizes this vector (same direction, length is 1); overwrites components
        /// </summary>
        public void NormalizeInPlace() {
            double l = 1.0 / Abs();
            x *= l;
            y *= l;
            z *= l;
        }

        /// <summary>
        /// returns a normalized (same direction, length is 1) copy of this vector
        /// </summary>
        public Vector Normalize() {
            Vector r = this;
            r.NormalizeInPlace();
            return r;
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
        /// multiples this vector with scalar <paramref name="s"/>;
        /// </summary>
        /// <param name="s">scaling factor</param>
        public void ScaleInPlace(double s) {
            this.x *= s;
            this.y *= s;
            this.z *= s;
        }
        
        /// <summary>
        /// multiples this vector with scalar <paramref name="s"/>, returns the result;
        /// this vector will be un-changed.
        /// </summary>
        /// <param name="s"></param>
        public Vector Scale(double s) {
            Vector r = this;
            r.ScaleInPlace(s);
            return r;
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
            ret.ScaleInPlace(s);
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
        /// multiplication by a matrix
        /// </summary>
        /// <param name="s">the scalar</param>
        /// <param name="m">the matrix</param>
        /// <returns><paramref name="m"/>*<paramref name="s"/></returns>
        public static Vector operator *(MultidimensionalArray m, Vector s)
        {
            Vector result = new Vector(s.Dim);
            for (int i = 0; i < s.Dim; ++i)
            {
                for(int j = 0; j < s.Dim; ++j)
                {
                    result[i] += m[i, j] * s[j];
                }
            }
            return result;
        }

        /// <summary>
        /// Division by a scalar
        /// </summary>
        /// <param name="v">The vector</param>
        /// <param name="s">The scalar</param>
        /// <returns>v / s</returns>
        public static Vector operator /(Vector v, double s) {
            Vector result = v;
            result.ScaleInPlace(1.0 / s);
            return result;
        }


        /// <summary>
        /// standard Dot/Vector/Inner-Product.
        /// </summary>
        /// <param name="a">1st operand</param>
        /// <param name="b">2nd operand</param>
        /// <returns>a*b</returns>
        public static double operator *(Vector a, Vector b) {
             if(a.Dim != b.Dim)
                throw new ArgumentException("Dimension mismatch");

            return (a.x * b.x + a.y * b.y + a.z * b.z);
        }
        
        /// <summary>
        /// Implicit conversion an array of doubles of length <see cref="Dim"/>
        /// </summary>
        /// <param name="v">The vector to be converted</param>
        /// <returns>An array of doubles of length <see cref="Dim"/></returns>
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
        /// Implicit conversion an array of doubles of length <see cref="Dim"/>
        /// </summary>
        /// <param name="v">The vector to be converted</param>
        /// <returns>An array of doubles of length <see cref="Dim"/></returns>
        public static implicit operator Vector( double[]  v) {
            return new Vector(v);
        }
        

        /// <summary>
        /// a vector notation: (x|y);
        /// </summary>
        /// <returns>(x|y)</returns>
        public override string ToString() {
            switch(this.Dim) {
                case 0:
                return "nüll";
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
        /// empty/default vector of dimension <paramref name="D"/>
        /// </summary>
        public static Vector Empty(int D) {
            if(D < 0 || D > 3)
                throw new ArgumentOutOfRangeException("Invalid spatial dimension for vector: " + D);
            Vector e = default(Vector);
            e.Dim = D;
            return e;
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
        public static Vector StdBasis2D(int d) {
            return StdBasis(d, 2);
        }

        /// <summary>
        /// the <paramref name="d"/>-th standard basis in 3D
        /// </summary>
        /// <param name="d">direction index</param>
        /// <returns>
        /// a three-dimensional vector with <paramref name="d"/>-th coordinate equal to 1.0.
        /// </returns>
        public static Vector StdBasis3D(int d) {
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
                throw new ArgumentException("Length can only be 1, 2 or 3, not " + length, "length");
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
                // throw new ArgumentException("Length can only be 0, 1 or 2", "length");
                throw new ArgumentException("Length can only be 1, 2 or 3, not " + length, "length");
            }

            if (destinationIndex < 0) {
                throw new ArgumentException("The destination index must be >= 0", "destinationIndex");
            }

            for (int i = 0; i < length; i++) {
                destination[i + destinationIndex] = this[i];
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public int IndexOf(double item) {
            for(int d = 0; d < Dim; d++) {
                if (item == this[d])
                    return d;
            }
            return -1;
        }

        /// <summary>
        /// inserting an entry, increasing vector dimension until 3D is reached
        /// </summary>
        public void Insert(int index, double item) {
            if(Dim >= 3)
                throw new NotSupportedException("Vector is already 3D - adding another dimension is not supported.");
            if(index < 0 || index >= Dim)
                throw new ArgumentOutOfRangeException();

            Dim++;
            for(int i = Dim - 1; i > index; i--) {
                this[i] = this[i - 1];
            }
            this[index] = item;
        }

        /// <summary>
        /// removing some entry, reducing vector dimension
        /// </summary>
        public void RemoveAt(int index) {
            if(index < 0 || index >= Dim)
                throw new ArgumentOutOfRangeException();

            for(int i = index; i < Dim - 1; i++) {
                this[i] = this[i + 1];
            }
            this[Dim - 1] = 0;
            Dim--;
        }

        /// <summary>
        /// adding until 3D is reached, dimension of vector changes
        /// </summary>
        public void Add(double item) {
            if(Dim >= 3)
                throw new NotSupportedException("Vector is already 3D - adding another dimension is not supported.");

            Dim++;
            this[Dim - 1] = item;
        }

        /// <summary>
        /// not supported - read-only (<see cref="Dim"/> could be changed)
        /// </summary>
        public void Clear() {
            throw new NotSupportedException("Not supported, since it could be misleading (removing all entries and setting dimension to zero vs. setting all entries to 0.0 and keeping spatial dimension).");
        }

        /// <summary>
        /// %
        /// </summary>
        public bool Contains(double item) {
            return this.IndexOf(item) >= 0;
        }

        /// <summary>
        /// not supported - read-only (<see cref="Dim"/> could be changed)
        /// </summary>
        public bool Remove(double item) {
            throw new NotSupportedException("Read-Only.");
        }

        public IEnumerator<double> GetEnumerator() {
            return new MyEnum(this);
        }

        IEnumerator IEnumerable.GetEnumerator() {
            return GetEnumerator();
        }

        class MyEnum : IEnumerator<double> {
            internal MyEnum(Vector __vec) {
                vec = __vec;
            }

            private Vector vec;
            int pos = -1;

            public double Current {
                get {
                    if (pos < 0)
                        throw new InvalidOperationException();
                    if (pos >= vec.Dim)
                        throw new InvalidOperationException();
                    return vec[pos];
                }
            }

            object IEnumerator.Current {
                get {
                    return Current;
                }
            }

            public void Dispose() {
            }

            public bool MoveNext() {
                pos++;
                return (pos < vec.Dim);
            }

            public void Reset() {
                pos = -1;
            }
        }

        /// <summary>
        /// Initializes this from the second dimension of a 2D array
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="i1">index into 1st dimension of <paramref name="mda"/></param>
        public void SetFrom(MultidimensionalArray mda, int i1) {
#if DEBUG
            if (mda.Dimension != 2)
                throw new ArgumentException("Expecting a 2D-array."); ;
            if (mda.GetLength(1) != Dim)
                throw new ArgumentException("Second dimension mismatch.");
#endif
            x = mda[i1, 0];
            if(Dim > 1) 
                y = mda[i1, 1];
            if(Dim > 2) 
                z = mda[i1, 2];
        }

        /// <summary>
        /// Initializes this from the third dimension of a 3D array
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="i1">index into 1st dimension of <paramref name="mda"/></param>
        /// <param name="i2">index into 2nd dimension of <paramref name="mda"/></param>
        public void SetFrom(MultidimensionalArray mda, int i1, int i2) {
#if DEBUG
            if (mda.Dimension != 3)
                throw new ArgumentException("Expecting a 3D-array."); ;
            if (mda.GetLength(2) != Dim)
                throw new ArgumentException("Second dimension mismatch.");
#endif
            x = mda[i1, i2, 0];
            if(Dim > 1) 
                y = mda[i1, i2, 1];
            if(Dim > 2) 
                z = mda[i1, i2, 2];
        }

        /// <summary>
        /// Initializes this from an array with arbitrary dimensions
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="IndexOffset">the index into <paramref name="mda"/>, where to start reading data</param>
        /// <param name="OriginDim">the dimension from which to take the vector entries</param>
        public void SetFrom(int OriginDim, MultidimensionalArray mda, params int[] IndexOffset) {
#if DEBUG
            if (mda.Dimension != IndexOffset.Length)
                throw new ArgumentException("Expecting an array of dimension equal to length of IndexOffset.");
            if (mda.GetLength(OriginDim) != Dim)
                throw new ArgumentException("Origin dimension has mismatching length."); 
#endif
            x = mda[IndexOffset];
            if (Dim > 1) {
                IndexOffset[OriginDim]++;
                y = mda[IndexOffset];
            }
            if (Dim > 2) {
                IndexOffset[OriginDim]++;
                z = mda[IndexOffset];
            }
        }

        /// <summary>
        /// Initializes this from an array with arbitrary dimensions
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="inc">index increase</param>
        /// <param name="IndexOffset">the index into <paramref name="mda"/>, where to start reading data</param>
        /// <param name="OriginDim">the dimension from which to take the vector entries</param>
        public void SetFrom(int OriginDim, int inc, MultidimensionalArray mda, params int[] IndexOffset) {
#if DEBUG
            if (mda.Dimension != IndexOffset.Length)
                throw new ArgumentException("Expecting an array of dimension equal to length of IndexOffset.");
#endif
            x = mda[IndexOffset];
            if (Dim > 1) {
                IndexOffset[OriginDim] += inc;
                y = mda[IndexOffset];
            }
            if (Dim > 2) {
                IndexOffset[OriginDim] += inc;
                z = mda[IndexOffset];
            }
        }



        /// <summary>
        /// Writes the components of this vector to a 2D array
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="i1">index into 1st dimension of <paramref name="mda"/></param>
        public void WriteTo(MultidimensionalArray mda, int i1) {
#if DEBUG
            if (mda.Dimension != 2)
                throw new ArgumentException("Expecting a 2D-array."); ;
            if (mda.GetLength(1) != Dim)
                throw new ArgumentException("Second dimension mismatch.");
#endif
            mda[i1, 0] = x;
            if(Dim > 1) 
                mda[i1, 1] = y;
            if(Dim > 2) 
                mda[i1, 2] = z;
        }

        /// <summary>
        /// Writes the components of this vector to a 3D array
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="i1">index into 1st dimension of <paramref name="mda"/></param>
        /// <param name="i2">index into 2nd dimension of <paramref name="mda"/></param>
        public void WriteTo(MultidimensionalArray mda, int i1, int i2) {
#if DEBUG
            if (mda.Dimension != 3)
                throw new ArgumentException("Expecting a 3D-array."); ;
            if (mda.GetLength(2) != Dim)
                throw new ArgumentException("Second dimension mismatch.");
#endif
            mda[i1, i2, 0] = x;
            if(Dim > 1) 
                mda[i1, i2, 1] = y;
            if(Dim > 2) 
                mda[i1, i2, 2] = z;
        }

        /// <summary>
        /// Writes the components of this vector to an array with arbitrary dimensions
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="IndexOffset">the index into <paramref name="mda"/>, where to start writing data</param>
        /// <param name="DestDim">the dimension from which to take the vector entries</param>
        public void WriteTo(int DestDim, MultidimensionalArray mda, params int[] IndexOffset) {
#if DEBUG
            if (mda.Dimension != IndexOffset.Length)
                throw new ArgumentException("Expecting an array of dimension equal to length of IndexOffset."); ;
            if (mda.GetLength(DestDim) != Dim)
                throw new ArgumentException("Destination dimension has mismatching length.");
#endif

            mda[IndexOffset] = x;
            if (Dim > 1) {
                IndexOffset[DestDim]++;
                mda[IndexOffset] = y;
            }
            if (Dim > 2) {
                IndexOffset[DestDim]++;
                mda[IndexOffset] = z;
            }
        }

        /// <summary>
        /// Writes the components of this vector to an array with arbitrary dimensions
        /// </summary>
        /// <param name="mda">data origin</param>
        /// <param name="inc">index increase</param>
        /// <param name="IndexOffset">the index into <paramref name="mda"/>, where to start writing data</param>
        /// <param name="DestDim">the dimension from which to take the vector entries</param>
        public void WriteTo(int DestDim, int inc, MultidimensionalArray mda, params int[] IndexOffset) {
#if DEBUG
            if (mda.Dimension != IndexOffset.Length)
                throw new ArgumentException("Expecting an array of dimension equal to length of IndexOffset."); ;
#endif

            mda[IndexOffset] = x;
            if (Dim > 1) {
                IndexOffset[DestDim] += inc;
                mda[IndexOffset] = y;
            }
            if (Dim > 2) {
                IndexOffset[DestDim] += inc;
                mda[IndexOffset] = z;
            }
        }

        /// <summary>
        /// Helps with the Serialization/Deserialization of Vectors in control files
        /// </summary>
        public class VectorConverter : JsonConverter<Vector> {
            public override bool CanWrite => true;
            public override bool CanRead => true;

            public override void WriteJson(JsonWriter writer, Vector value, JsonSerializer serializer) {
                // Create a JObject and write the x and y properties
                var jObject = new JObject();
                jObject["Dim"] = value.Dim;
                jObject["x"] = value.x;
                if (value.Dim > 1)
                    jObject["y"] = value.y;
                if (value.Dim > 2)
                    jObject["z"] = value.z;

                jObject.WriteTo(writer);
            }

            public override Vector ReadJson(JsonReader reader, Type objectType, Vector existingValue, bool hasExistingValue, JsonSerializer serializer) {


                // Read JObject from the reader
                var jObject = JObject.Load(reader);

                // Deserialize x and y from JObject
                var Dim = jObject["Dim"].Value<int>();

                var ret = new Vector(Dim);
                ret.x = jObject["x"].Value<double>();
                if (Dim > 1)
                    ret.y = jObject["y"].Value<double>();
                if (Dim > 2)
                    ret.z = jObject["z"].Value<double>();

                return ret;
            }
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
        /// extracts a vector from a <see cref="MultidimensionalArray"/> <paramref name="inp"/>;
        /// </summary>
        /// <param name="inp">
        /// input array
        /// </param>
        /// <param name="i0">
        /// first index into <paramref name="inp"/>
        /// </param>
        /// <param name="i1">
        /// second index into <paramref name="inp"/>
        /// </param>
        /// <returns>
        /// A vector with dimension (<see cref="Vector.Dim"/>) equal to 3rd length of <paramref name="inp"/>
        /// </returns>
        public static Vector GetRowPt(this MultidimensionalArray inp, int i0, int i1) {
            switch(inp.GetLength(2)) {
                case 1:
                return new Vector(inp[i0, i1, 0]);
                case 2:
                return new Vector(inp[i0, i1, 0], inp[i0, i1, 1]);
                case 3:
                return new Vector(inp[i0, i1, 0], inp[i0, i1, 1], inp[i0, i1, 2]);
                default:
                throw new ArgumentException("2nd length of input is " + inp.GetLength(2) + ", this cannot be a spatial dimension.");
            }
        }


        /// <summary>
        /// Matrix-Vector product
        /// </summary>
        public static Vector MtxVecMul(this IMatrix M, Vector v) {
            if(M.NoOfCols != v.Dim)
                throw new ArgumentException();

            var R = new Vector(M.NoOfRows);
            for(int i = M.NoOfRows - 1; i >= 0; i--) {
                double acc = 0;
                for(int j = 0; j < v.Dim; j++) {
                    acc += M[i, j] * v[j];
                }

                R[i] = acc;
            }

            return R;
        }


        /// <summary>
        /// sets the (<paramref name="i0"/>,<paramref name="i1"/>)-th row from <paramref name="inp"/> to values provided by <paramref name="row"/>.
        /// </summary>
        /// <param name="inp">
        /// matrix that should be altered
        /// </param>
        /// <param name="i0">
        /// first index into <paramref name="inp"/>
        /// </param>
        /// <param name="i1">
        /// second index into <paramref name="inp"/>
        /// </param>
        /// <param name="row">
        /// a vector 
        /// </param>
        public static void SetRowPt(this MultidimensionalArray inp, int i0, int i1, Vector row) {
            if (row.Dim != inp.NoOfCols)
                throw new ArgumentException("Dimension mismatch.");

            switch(row.Dim) {
                case 1:
                inp[i0, i1, 0] = row.x; return;
                case 2:
                inp[i0, i1, 0] = row.x; inp[i0, i1, 1] = row.y; return;
                case 3:
                inp[i0, i1, 0] = row.x; inp[i0, i1, 1] = row.y; inp[i0, i1, 2] = row.z; return;
                default:
                throw new NotImplementedException();
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
        public static Vector GetRowPt(this double[,] inp, int RowNo) {
            switch(inp.GetLength(1)) {
                case 1:
                return new Vector(inp[RowNo, 0]);
                case 2:
                return new Vector(inp[RowNo, 0], inp[RowNo, 1]);
                case 3:
                return new Vector(inp[RowNo, 0], inp[RowNo, 1], inp[RowNo, 2]);
                default:
                throw new ArgumentException("Matrix has " + inp.GetLength(1) + " columns, this cannot be a spatial dimension.");
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
        public static void SetRowPt(this double[,] inp, int RowNo, Vector row) {
            if (row.Dim != inp.GetLength(1))
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
