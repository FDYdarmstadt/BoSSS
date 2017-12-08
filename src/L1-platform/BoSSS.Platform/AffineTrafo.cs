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
using System.Runtime.Serialization;
using ilPSP.Utils;
using ilPSP;
using System.Diagnostics;

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// represents an affine - linear transformation;
    /// </summary>
    /// <remarks> 
    /// The transformation of some point \f$ \vec{x} \f$ into it's 
    /// image \f$ \vec{y} \f$
    /// is defined as
    /// 
    /// <i>y</i> = <see cref="Matrix"/>*\f$ \vec{x} \f$ + <see cref="Affine"/>
    /// </remarks>
    [Serializable]
    [DataContract]
    public class AffineTrafo : ICloneable, IEquatable<AffineTrafo> {

        /// <summary>
        /// initializes a zero transformation
        /// </summary>
        /// <param name="D_cod">
        /// spatial dimension of codomain
        /// </param>
        /// <param name="D_dom">
        /// spatial dimension of domain
        /// </param>
        public AffineTrafo(int D_dom, int D_cod) {
            Matrix = MultidimensionalArray.Create(D_cod, D_dom);
            Affine = new double[D_cod];
        }

        /// <summary>
        /// initializes a zero transformation
        /// </summary>
        /// <param name="D">spatial dimension</param>
        public AffineTrafo(int D) : this(D, D) {
        }

        /// <summary>
        /// initializes a transformation with matrix and affine vector.
        /// </summary>
        public AffineTrafo(IMatrix __Matrix, double[] __Affine) {
            if (__Affine.Length != __Matrix.NoOfRows)
                throw new ArgumentException();
            this.Matrix = MultidimensionalArray.Create(__Matrix.NoOfRows, __Matrix.NoOfCols); 
            this.Matrix.SetMatrix(__Matrix);
            this.Affine = __Affine;
        }

        /// <summary>
        /// empty ctor.
        /// </summary>
        private AffineTrafo() {
        }

        /// <summary>
        /// linear part of the transformation;
        /// </summary>
        [DataMember]
        public MultidimensionalArray Matrix;

        /// <summary>
        /// affine part of the transformation;
        /// </summary>
        [DataMember]
        public double[] Affine;

        /// <summary>
        /// Computes <see cref="Matrix"/>*<paramref name="vtx"/> + <see cref="Affine"/>;
        /// </summary>
        /// <param name="vtx">
        /// the vector/point/vertex that should be transformed;
        /// </param>
        /// <returns></returns>
        /// <remarks>
        /// This method should not be used for performance-critical
        /// issues; One reason for this is that it allocates a small array
        /// for return value.
        /// </remarks>
        public double[] Transform(params double[] vtx) {
            if (vtx.Length != Matrix.NoOfCols)
                throw new ArgumentException("dimension mismatch");

            double[] ret = (double[])Affine.Clone();
            Matrix.gemv(1.0, vtx, 1.0, ret);

            return ret;
        }

        /// <summary>
        /// applies this transformation to a series of vectors
        /// </summary>
        public void Transform(double[] vtxIn, double[] vtxOut) {
            if (vtxIn.GetLength(0) != DomainDimension)
                throw new ArgumentException("spatial dimension mismatch");
            if (vtxOut.GetLength(0) != CodomainDimension)
                throw new ArgumentException("spatial dimension mismatch");

            int M = DomainDimension; // number of columns
            int N = CodomainDimension; // number of rows

            for (int n = 0; n < N; n++) {
                double acc = Affine[n];

                for (int m = 0; m < M; m++)
                    acc += Matrix[n, m] * vtxIn[m];

                vtxOut[n] = acc;
            }

        }

        /// <summary>
        /// Applies this transformation to a series of vectors.
        /// </summary>
        public void Transform(MultidimensionalArray vtxIn, MultidimensionalArray vtxOut) {
            if (vtxIn.GetLength(1) != DomainDimension)
                throw new ArgumentException("spatial dimension mismatch");
            if (vtxOut.GetLength(1) != CodomainDimension)
                throw new ArgumentException("spatial dimension mismatch");
            if (vtxIn.GetLength(0) != vtxOut.GetLength(0))
                throw new ArgumentException("mismatch in number of vectors between in and out;");

            int L = vtxIn.GetLength(0);
            int M = DomainDimension; // number of columns
            int N = CodomainDimension; // number of rows
            for (int l = 0; l < L; l++) {
                for (int n = 0; n < N; n++) {
                    double acc = Affine[n];

                    for (int m = 0; m < M; m++)
                        acc += Matrix[n, m] * vtxIn[l, m];

                    vtxOut[l, n] = acc;
                }
            }
        }
        

        /// <summary>
        /// applies this transformation to a series of vectors
        /// </summary>
        public double[,] Transform(double[,] vtx) {
            if (vtx.GetLength(1) != DomainDimension)
                throw new ArgumentException("spatial dimension mismatch");

            int L = vtx.GetLength(0);
            var ret = new double[L, CodomainDimension];
            double[] x = new double[DomainDimension];

            for (int l = 0; l < L; l++) {
                ArrayTools.GetRow(vtx, l, x);
                var yl = Transform(x);
                ret.SetRow(l, yl);
            }

            return ret;
        }

        /// <summary>
        /// applies this transformation to a series of vectors
        /// </summary>
        public MultidimensionalArray Transform(MultidimensionalArray vtx) {
            if (vtx.Dimension != 2)
                throw new ArgumentException();
            if (vtx.GetLength(1) != DomainDimension)
                throw new ArgumentException("spatial dimension mismatch");

            int L = vtx.GetLength(0);
            var ret = MultidimensionalArray.Create(L, CodomainDimension);
            double[] x = new double[DomainDimension];

            for (int l = 0; l < L; l++) {
                vtx.GetRow(l, x);
                var yl = Transform(x);
                ret.SetSubVector(yl, l, -1);
            }

            return ret;
        }

        /// <summary>
        /// computes the inverse transformation.<br/>
        /// Quite surprising, ey?
        /// </summary>
        /// <returns></returns>
        public AffineTrafo Invert() {
            AffineTrafo inv = new AffineTrafo();
            int D = this.Matrix.NoOfCols;

            inv.Matrix = MultidimensionalArray.Create(D, D);
            this.Matrix.InvertTo(inv.Matrix);

            inv.Affine = new double[D];
            inv.Matrix.gemv(-1.0, this.Affine, 0.0, inv.Affine);

            return inv;
        }

        /// <summary>
        /// dimension of domain
        /// </summary>
        public int DomainDimension {
            get {
                return Matrix.NoOfCols;
            }
        }

        /// <summary>
        /// dimension of codomain
        /// </summary>
        public int CodomainDimension {
            get {
                return Matrix.NoOfRows;
            }
        }

        /// <summary>
        /// composition of two transformations;
        /// </summary>
        /// <param name="left"></param>
        /// <param name="right"></param>
        /// <returns>
        /// returns a transformation, which's application is equivalent to
        /// first applying <paramref name="right"/> and secondly applying <paramref name="left"/>;
        /// </returns>
        /// <remarks>
        /// let be <i>M</i> and <i>o</i> matrix (see <see cref="Matrix"/>) and affine vector (see <see cref="Affine"/>)
        /// of the <paramref name="left"/> Transformation and, respectively <i>N</i> and <i>p</i>
        /// matrix and vector of <paramref name="right"/>.<br/>
        /// Then the matrix of the returned value is <i>M</i>*<i>N</i>, and the vector of the returned 
        /// transformation is <i>M</i>*<i>p</i> + <i>o</i>;
        /// </remarks>
        public static AffineTrafo operator *(AffineTrafo left, AffineTrafo right) {
            if (right.CodomainDimension != left.DomainDimension)
                throw new ArgumentException("mismatch in spatial dimension");

            var ret = new AffineTrafo(right.DomainDimension, left.CodomainDimension);

            ret.Matrix.GEMM(1.0, left.Matrix, right.Matrix, 0.0);
            Array.Copy(left.Affine, ret.Affine, left.Affine.Length);
            left.Matrix.gemv(1.0, right.Affine, 1.0, ret.Affine);

            return ret;
        }

        /// <summary>
        /// creates an identity transformation
        /// </summary>
        /// <param name="D">spatial dimension</param>
        /// <returns></returns>
        public static AffineTrafo Identity(int D) {
            AffineTrafo tr = new AffineTrafo();
            tr.Matrix = MultidimensionalArray.Create(D, D); 
            tr.Matrix.AccEye(1.0);
            tr.Affine = new double[D];
            return tr;
        }

        /// <summary>
        /// computes matrix and affine vector of the affine-linear transformation
        /// that maps the vectors in the <paramref name="preimage"/> to <paramref name="image"/>;
        /// </summary>
        /// <param name="preimage">
        /// The preimage: \f$  (D+1) \f$  vectors of dimension \f$ D \f$;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is \f$ D+1 \f$</item>
        ///   <item>2nd index: spatial dimension, lenght is D</item>
        /// </list>
        /// Tip: use <see cref="ilPSP.Utils.ArrayTools.Transpose{t}(t[,],t[,])"/>
        /// if the sequence of indices is not appropriate;
        /// </param>
        /// <param name="image">
        /// The image: \f$  (D+1) \f$ vectors of dimension \f$ D \f$;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is \f$ D + 1\f$</item>
        ///   <item>2nd index: spatial dimension, length is \f$ D \f$</item>
        /// </list>
        /// </param>
        public static AffineTrafo FromPoints(MultidimensionalArray preimage, MultidimensionalArray image) {
            if (preimage.Dimension != 2)
                throw new ArgumentException("expecting 2D - array.", "preimage");
            if (image.Dimension != 2)
                throw new ArgumentException("expecting 2D - array.", "image");


            return FromPoints(preimage.To2DArray(), image.To2DArray());
        }

        /// <summary>
        /// Returns a 2D transformation.
        /// </summary>
        public static AffineTrafo Some2DRotation(double angle) {
            AffineTrafo Trafo = new AffineTrafo();

            MultidimensionalArray dreh = MultidimensionalArray.Create(2, 2);
            dreh[0, 0] = Math.Cos(angle);
            dreh[0, 1] = -Math.Sin(angle);
            dreh[1, 0] = Math.Sin(angle);
            dreh[1, 1] = Math.Cos(angle);
            Trafo.Affine = new double[2];
            Trafo.Matrix = dreh;

            return Trafo;
        }

        /// <summary>
        /// computes matrix and affine vector of the affine-linear transformation
        /// that maps the vectors in the <paramref name="preimage"/> to <paramref name="image"/>;
        /// </summary>
        /// <param name="preimage">
        /// The preimage: L vectors of dimension D_dom;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is L</item>
        ///   <item>2nd index: spatial dimension, length is D_dom</item>
        /// </list>
        /// Tip: use <see cref="ilPSP.Utils.ArrayTools.Transpose{t}(t[,],t[,])"/>
        /// if the sequence of indices is not appropriate;
        /// </param>
        /// <param name="image">
        /// The image: L vectors of dimension D_cod;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is L</item>
        ///   <item>2nd index: spatial dimension, length is D_cod</item>
        /// </list>
        /// </param>
        /// <remarks>
        /// L*D_cod == D_cod + L*D_dom
        /// </remarks>
        public static AffineTrafo FromPoints(double[,] preimage, double[,] image) {
            int D_dom = preimage.GetLength(1);
            int D_cod = image.GetLength(1);
            int L = preimage.GetLength(0);  // number of points

            int rows = L * D_cod;
            int cols = D_cod + D_cod * D_dom;

            if (image.GetLength(0) != preimage.GetLength(0))
                throw new ArgumentException("number of points in image and preimage must be equal.");
            if (rows < cols)
                throw new ArgumentException("insufficient information - not enough points.");
            if (rows > cols)
                throw new ArgumentException("over-determined information - to many points.");



            MultidimensionalArray M = MultidimensionalArray.Create(rows, cols);
            double[] rhs = new double[M.NoOfCols];
            double[] x = (double[])rhs.Clone();

            // gesucht: affine lineare Transformation: "xi -> x"
            // preimage = xi, image = x;
            // 
            //          x = Tr*xi + o
            //
            // Sei z.B. D=2 (andere D's ergeben sich analog)
            // x1 = (x11,x12)^T, x2 = (x21,x22)^T und x3 seien die Bilder von
            // xi1 = (xi11,xi12), xi2 sowie xi3.
            //
            // gesucht sind die Matrix Tr und der Vektor o = (o1,o2)^T,
            //
            //                   [ Tr11  Tr12 ]
            //              Tr = [            ],
            //                   [ Tr21  Tr22 ]
            //
            // welche mit einem Glsys. der folgenden Form berechnet werden:
            //
            // [ 1 0   xi11 xi12     0    0  ]   [ o1  ]   [ x11 ]
            // [ 0 1    0    0     xi11 xi12 ]   [ o2  ]   [ x12 ]
            // [                             ]   [     ]   [     ]
            // [ 1 0   xi21 xi22     0    0  ] * [ T11 ] = [ x21 ]
            // [ 0 1    0    0     xi21 xi22 ]   [ T12 ]   [ x22 ]
            // [                             ]   [     ]   [     ]
            // [ 1 0   xi31 xi32     0    0  ]   [ T21 ]   [ x31 ]
            // [ 0 1    0    0     xi31 xi32 ]   [ T22 ]   [ x31 ]
            //
            //

            // build rhs
            for (int l = 0; l < L; l++) { // loop over image points
                for (int d = 0; d < D_cod; d++) {
                    rhs[l * D_cod + d] = image[l, d];
                }
            }

            // build matrix
            M.Clear();
            for (int l = 0; l < L; l++) { // loop over image points
                for (int d = 0; d < D_cod; d++) {
                    M[l * D_cod + d, d] = 1.0;

                    for (int ddd = 0; ddd < D_dom; ddd++) {
                        M[l * D_cod + d, D_cod + D_dom * d + ddd] = preimage[l, ddd];
                    }

                }
            }

            // solve Matrix
            M.Solve(x, rhs);

            // save results, return
            AffineTrafo Trafo = new AffineTrafo(D_dom, D_cod);
            for (int d = 0; d < D_cod; d++) {
                Trafo.Affine[d] = x[d];

                for (int dd = 0; dd < D_dom; dd++) {
                    Trafo.Matrix[d, dd] = x[D_cod + d * D_dom + dd];
                }
            }
            return Trafo;
        }

        /// <summary>
        /// creates a non-shallow copy
        /// </summary>
        public AffineTrafo CloneAs() {
            AffineTrafo ret = new AffineTrafo();
            ret.Matrix = this.Matrix.CloneAs();
            ret.Affine = this.Affine.CloneAs();
            return ret;
        }

        /// <summary>
        /// creates a non-shallow copy
        /// </summary>
        public object Clone() {
            return CloneAs();
        }

        /// <summary>
        /// See <see cref="Equals(AffineTrafo)"/>
        /// </summary>
        /// <param name="obj">See <see cref="object.Equals(object)"/></param>
        /// <returns>See <see cref="Equals(AffineTrafo)"/></returns>
        public override bool Equals(object obj) {
            if (obj is AffineTrafo) {
                return Equals((AffineTrafo)obj);
            } else {
                return false;
            }
        }

        /// <summary>
        /// Builds a hash code based on the hashes of <see cref="Affine"/> and
        /// <see cref="Matrix"/>.
        /// </summary>
        /// <returns>A hash of this object</returns>
        public override int GetHashCode() {
            unchecked {
                int hash = 17;
                hash = hash * 23 + Affine.GetHashCode();
                hash = hash * 23 + Matrix.GetHashCode();
                return hash;
            }
        }

        #region IEquatable<AffineTrafo> Members

        /// <summary>
        /// Checks this object for equality to the given transformation based
        /// on <see cref="Affine"/> and <see cref="Matrix"/>.
        /// </summary>
        /// <param name="other">
        /// The object to be compared with.
        /// </param>
        /// <returns>
        /// True, if <see cref="Affine"/> and <see cref="Matrix"/> of the given
        /// transformation are equal to the entities in this object. False,
        /// otherwise.
        /// </returns>
        public bool Equals(AffineTrafo other) {
            return ApproximateEquals(other, 1.0e-12);
        }

        #endregion

        /// <summary>
        /// approximately equal
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        /// <param name="RelTol">relative tolerance for comparison</param>
        public bool ApproximateEquals(AffineTrafo other, double RelTol = 1.0e-9) {
            int D = this.Affine.Length;

            if (other.Affine.Length != D || other.Matrix.NoOfCols != this.Matrix.NoOfCols || other.Matrix.NoOfRows != this.Matrix.NoOfRows)
                throw new ArgumentException("mismatch in spatial dimension.");

            double tolAffine = GenericBlas.L2NormPow2(this.Affine) * RelTol;
            if (GenericBlas.L2DistPow2(this.Affine, other.Affine) > tolAffine)
                return false;

            Debug.Assert(this.Matrix.IsContinious && this.Matrix.Index(0, 0) == 0);
            Debug.Assert(other.Matrix.IsContinious && other.Matrix.Index(0, 0) == 0);
            double tolTrafo = GenericBlas.L2NormPow2(this.Matrix.Storage) * RelTol;
            if(GenericBlas.L2DistPow2(other.Matrix.Storage, this.Matrix.Storage) > tolTrafo)
                return false;

            return true;
        }
    }
}

