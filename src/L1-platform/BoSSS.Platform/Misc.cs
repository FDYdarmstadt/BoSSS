/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;
using System.Collections.Generic;
using System.Text;
using ilPSP.Utils;

namespace BoSSS.Platform.LinAlg {

    /// <summary>
    /// represents an affine - linear transformation;
    /// </summary>
    /// <remarks> 
    /// The transformation of some point <i>x</i> into it's 
    /// image <i>y</i>
    /// is defined as follows: <br/>
    /// <i>y</i> = <see cref="Matrix"/>*<i>x</i> + <see cref="Affine"/>
    /// </remarks>
    public struct AffineTrafo : ICloneable, IEquatable<AffineTrafo> {

        /// <summary>
        /// linear part of the transformation;
        /// </summary>
        public FullMatrix Matrix;

        /// <summary>
        /// affine part of the transformation;
        /// </summary>
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
            Matrix.gemv<double[], double[]>(1.0, vtx, 1.0, ret);

            return ret;
        }

        /// <summary>
        /// Computes: <paramref name="res"/> := <see cref="Matrix"/>*<paramref name="vtx"/> + <see cref="Affine"/>;
        /// </summary>
        /// <param name="vtx">
        /// the vector/point/vertex that should be transformed;
        /// </param>
        /// <param name="res">output</param>
        public void Transform(double[] vtx, double[] res) {
            if (vtx.Length != Matrix.NoOfCols)
                throw new ArgumentException("dimension mismatch");
            if (res.Length != Matrix.NoOfRows)
                throw new ArgumentException("dimension mismatch");

            Array.Copy(this.Affine, res, this.Affine.Length);
            Matrix.gemv<double[], double[]>(1.0, vtx, 1.0, res);
        }

        /// <summary>
        /// Transforms the list of points given in <paramref name="nodes"/> and
        /// stores the result in <paramref name="result"/>.
        /// </summary>
        /// <typeparam name="Tin"></typeparam>
        /// <typeparam name="Tout"></typeparam>
        /// <param name="nodes">
        /// The nodes to be transformed.
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Spatial dimension</item>
        /// </list>
        /// </param>
        /// <param name="result">
        /// On exit: The transformed nodes
        /// <list type="bullet">
        ///     <item>1st index: Node index</item>
        ///     <item>2nd index: Spatial dimension</item>
        /// </list>
        /// </param>
        public void Transform<Tin, Tout>(Tin nodes, Tout result)
            where Tin : IMatrix
            where Tout : IMatrix {
            int D = nodes.NoOfCols;
            int noOfNodes = nodes.NoOfRows;

            if (Matrix.NoOfCols != D) {
                throw new ArgumentException("Invalid spatial dimension", "nodes");
            }

            if (result.NoOfCols != D) {
                throw new ArgumentException("Invalid spatial dimension", "result");
            }

            if (result.NoOfRows != noOfNodes) {
                throw new ArgumentException("Invalid number of vertices", "result");
            }

            for (int i = 0; i < noOfNodes; i++) {
                for (int d = 0; d < D; d++) {
                    result[i, d] = Affine[d];

                    for (int dd = 0; dd < D; dd++) {
                        result[i, d] += Matrix[d, dd] * nodes[i, dd];
                    }
                }
            }
        }

        /// <summary>
        /// computes the inverse transformation.<br/>
        /// Quite surprising, ey?
        /// </summary>
        /// <returns></returns>
        public AffineTrafo Invert() {
            AffineTrafo inv = new AffineTrafo();
            int D = this.Matrix.NoOfCols;

            inv.Matrix = new FullMatrix(D, D);
            this.Matrix.Invert(inv.Matrix);

            inv.Affine = new double[D];
            inv.Matrix.gemv<double[], double[]>(-1.0, this.Affine, 0.0, inv.Affine);

            return inv;
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
            if (left.Affine.Length != right.Affine.Length)
                throw new ArgumentException("dimensional mismatch.");

            //int D = left.Affine.Length;
            AffineTrafo ret = new AffineTrafo();
            ret.Matrix = left.Matrix * right.Matrix;

            ret.Affine = (double[])left.Affine.Clone();
            left.Matrix.gemv<double[], double[]>(1.0, right.Affine, 1.0, ret.Affine);

            return ret;
        }


        /// <summary>
        /// creates an identity transformation
        /// </summary>
        /// <param name="D">spatial dimension</param>
        /// <returns></returns>
        public static AffineTrafo ID(int D) {
            AffineTrafo tr = new AffineTrafo();
            tr.Matrix = FullMatrix.Eye(D);
            tr.Affine = new double[D];
            return tr;
        }

        /// <summary>
        /// returns a 2D transformation
        /// </summary>
        public static AffineTrafo Some2DRotation(double angle) {
            AffineTrafo Trafo = new AffineTrafo();
           

            FullMatrix dreh = new FullMatrix(2, 2);
            dreh[0, 0] = Math.Cos(angle); dreh[0, 1] = -Math.Sin(angle);
            dreh[1, 0] = Math.Sin(angle); dreh[1, 1] = Math.Cos(angle);
            Trafo.Affine = new double[2];
            Trafo.Matrix = dreh;

            return Trafo;
        }


        /// <summary>
        /// computes matrix and affine vector of the affine-linear transformation
        /// that maps the vectors in the <paramref name="preimage"/> to <paramref name="image"/>;
        /// </summary>
        /// <param name="preimage">
        /// The preimage: (D+1) vectors of dimension D;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is D+1</item>
        ///   <item>2nd index: spatial dimension, length is D</item>
        /// </list>
        /// Tip: use <see cref="ilPSP.Utils.ArrayTools.Transpose{t}(t[,],t[,])"/>
        /// if the sequence of indices is not appropriate;
        /// </param>
        /// <param name="image">
        /// The image: (D+1) vectors of dimension D;
        /// <list type="bullet">
        ///   <item>1st index: vector/vertex index, length is D+1</item>
        ///   <item>2nd index: spatial dimension, length is D</item>
        /// </list>
        /// </param>
        /// <returns></returns>
        public static AffineTrafo FromPoints(double[,] preimage, double[,] image) {
            int D = preimage.GetLength(1);
            if (image.GetLength(1) != D || image.GetLength(0) != D + 1)
                throw new ArgumentException("wrong size; expecting " + D + "x" + (D + 1) + " - array;", "image");
            if (preimage.GetLength(0) != D + 1)
                throw new ArgumentException("wrong size; expecting " + D + "x" + (D + 1) + " - array;", "preimage");


            FullMatrix M = new FullMatrix(D * (D + 1), D * (D + 1));
            double[] rhs = new double[M.NoOfCols];
            double[] x = (double[])rhs.Clone();

            // gesucht: affine lineare Transformation: "xi -> x"
            // preimage = xi, image = x;
            // 
            //          x = Tr*xi + o
            //
            // Sei z.B. D=2 (andere D's ergeben sich analog)
            // x1 = (x11,x12)^T, x2 = (x21,x22)^T, x3 seien die Bilder von
            // xi1 = (xi11,xi12), xi2, xi3.
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
            // [                             ]   [  0  ]   [     ]
            // [ 1 0   xi21 xi22     0    0  ] * [ T11 ] = [ x21 ]
            // [ 0 1    0    0     xi21 xi22 ]   [ T12 ]   [ x22 ]
            // [                             ]   [     ]   [     ]
            // [ 1 0   xi31 xi32     0    0  ]   [ T21 ]   [ x31 ]
            // [ 0 1    0    0     xi31 xi32 ]   [ T22 ]   [ x31 ]
            //
            //


            // build rhs
            for (int d = 0; d < (D + 1); d++) {
                for (int dd = 0; dd < D; dd++) {
                    rhs[d * D + dd] = image[d, dd];
                }
            }


            // build matrix
            M.Clear();
            for (int d = 0; d < (D + 1); d++) {
                for (int dd = 0; dd < D; dd++) {
                    M[d * D + dd, dd] = 1.0;

                    for (int ddd = 0; ddd < D; ddd++) {
                        M[d * D + dd, D + D * dd + ddd] = preimage[d, ddd];
                    }

                }
            }

            // solve Matrix
            M.Solve(x, rhs);

            // save results, return
            AffineTrafo Trafo = new AffineTrafo();
            Trafo.Affine = new double[D];
            Trafo.Matrix = new FullMatrix(D, D);
            for (int d = 0; d < D; d++) {
                Trafo.Affine[d] = x[d];

                for (int dd = 0; dd < D; dd++) {
                    Trafo.Matrix[d, dd] = x[D + d * D + dd];
                }
            }
            return Trafo;
        }


        /// <summary>
        /// creates a non-shallow copy
        /// </summary>
        /// <returns></returns>
        public object Clone() {
            AffineTrafo ret;
            ret.Matrix = (FullMatrix)this.Matrix.Clone();
            ret.Affine = (double[])this.Affine.Clone();
            return ret;
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
            if (ArrayTools.Equals(Affine, other.Affine) && Matrix.Equals(other.Matrix)) {
                return true;
            }

            return false;
        }

        #endregion
    }
}
