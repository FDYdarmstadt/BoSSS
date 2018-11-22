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
using System.Diagnostics;
using BoSSS.Platform;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Foundation.Grid {

    /// <summary>
    /// Just from its name, you might guess what this class is good for.
    /// </summary>
    public sealed class Polynomial : IComparable<Polynomial>, ICloneable {

        /// <summary>
        /// 1st index: Coefficient index; 2nd index: spatial dimension index;
        /// See <see cref="Coeff"/> for further information;
        /// </summary>
        public int[,] Exponents = new int[0, 0];

        /// <summary>
        /// Polynomial coefficients;
        /// The polynomial value at point (x,y) is evaluated as the sum 
        /// <see cref="Coeff"/>[i]*x^<see cref="Exponents"/>[i,0]*y^<see cref="Exponents"/>[i,1]
        /// over all i.
        /// </summary>
        public double[] Coeff = new double[0];

        /// <summary>
        /// Obsolete.
        /// </summary>
        public Polynomial() {
        }

        /// <summary>
        /// Obsolete. The polynomial GUID is no longer used.
        /// </summary>
        [Obsolete]
        public Polynomial(Guid __Guid) {
        }


        /// <summary>
        /// adds a coefficient/exponent pair to this polynomial
        /// </summary>
        public void AddCoeff(double _Coeff, int[] _Exponents) {
            if (Coeff == null || Coeff.Length == 0) {
                Coeff = new double[1];
                Exponents = new int[1, _Exponents.Length];
            } else {


                Array.Resize<double>(ref Coeff, Coeff.Length + 1);
                int[,] exp_new = new int[Exponents.GetLength(0) + 1, Exponents.GetLength(1)];
                Array.Copy(Exponents, exp_new, Exponents.Length);
                Exponents = exp_new;
            }

            int n = Coeff.Length - 1;
            Coeff[n] = _Coeff;
            int l = _Exponents.Length;
            if (l != Exponents.GetLength(1))
                throw new ArgumentException("wrong number of exponents");
            for (int i = 0; i < l; i++) {
                Exponents[n, i] = _Exponents[i];
            }
        }

        /// <summary>
        /// returns the coefficient index (index into <see cref="Coeff"/>) for a specific exponent.
        /// </summary>
        public int GetCoefficientIndex(int[] exponent) {
            if (exponent.Length != Exponents.GetLength(1)) {
                throw new ArgumentException(
                    String.Format(
                        "Invalid exponent. Length must {0}, but {1} was given",
                        Exponents.GetLength(1),
                        exponent.Length),
                    "exponent");
            }

            for (int i = 0; i < Coeff.Length; i++) {
                bool equals = true;
                for (int j = 0; j < Exponents.GetLength(1); j++) {
                    if (Exponents[i, j] != exponent[j]) {
                        equals = false;
                        break;
                    }
                }

                if (equals) {
                    return i;
                }
            }

            return -1;
        }

        /// <summary>
        /// polynomial degree
        /// </summary>
        public int AbsoluteDegree {
            get {
                Debug.Assert(this.Coeff.Length == this.Exponents.GetLength(0));
            
                int deg = 0;
                for (int i = 0; i < Exponents.GetLength(0); i++) {
                    int d = 0;
                    for (int j = 0; j < Exponents.GetLength(1); j++) {
                        d += Exponents[i, j];
                    }

                    if (d > deg)
                        deg = d;
                }
                return deg;
            }
        }

        /// <summary>
        /// 2 for polynomials in x,y and 3 for polynomials in x,y,z;
        /// </summary>
        public int SpatialDimension {
            get {
                Debug.Assert(this.Coeff.Length == this.Exponents.GetLength(0));
                return Exponents.GetLength(1);
            }
        }

        /// <summary>
        /// Creates a new polynomial by multiplying the given polynomial by -1
        /// </summary>
        /// <param name="p">
        /// A polynomial to be multiplied by -1
        /// </param>
        /// <returns>
        /// A new polynomial q with q(x) = -<paramref name="p"/>(x)
        /// </returns>
        public static Polynomial operator -(Polynomial p) {
            Polynomial negation = (-1.0) * p;
            return negation;
        }

        /// <summary>
        /// Creates a new polynomial by multiplying the given polynomial by <paramref name="scale"/>
        /// </summary>
        public static Polynomial operator *(double scale, Polynomial p) {
            Polynomial negation = p.CloneAs();
            for (int i = 0; i < negation.Coeff.Length; i++) {
                negation.Coeff[i] *= scale;
            }
            return negation;
        }

        /// <summary>
        /// Creates a new polynomial by multiplying the given polynomial by <paramref name="scale"/>
        /// </summary>
        public static Polynomial operator *(Polynomial p, double scale) {
            return scale * p;
        }

        /// <summary>
        /// Symbolic multiplication of two polynomials
        /// </summary>
        public static Polynomial operator *(Polynomial p, Polynomial q) {
            if (p.SpatialDimension != q.SpatialDimension)
                throw new ArgumentException("Spatial dimension mismatch.");
            Polynomial R = new Polynomial();
            if(p.Coeff.Length > 0 && q.Coeff.Length > 0) {
                for( int i = 0; i < q.Coeff.Length; i++) {
                    Polynomial Ri = p.MultiplyByMonomial(q.Exponents.GetRow(i), q.Coeff[i]);
                    if (i == 0)
                        R = Ri;
                    else
                        R = R + Ri;
                }
            }
            return R;
        }

        /// <summary>
        /// multiplies this polynomial with a monomial expression
        /// </summary>
        public Polynomial MultiplyByMonomial(int[] _Exponents, double alpha) {
            if (_Exponents.Length != this.SpatialDimension)
                throw new ArgumentException("Spatial dimension mismatch.");
            var R = this.CloneAs();
            int D = this.SpatialDimension;

            for( int i = 0; i < R.Coeff.Length; i++) {
                R.Coeff[i] *= alpha;

                for (int d = 0; d < D; d++) {
                    R.Exponents[i, d] *= _Exponents[d];
                }
            }
            
            return R;
        }



        /// <summary>
        /// Calculates the sum of <paramref name="p"/> and <paramref name="q"/>
        /// </summary>
        /// <param name="p">First operand</param>
        /// <param name="q">First operand</param>
        /// <returns>
        /// A new polynomial r with
        /// r(x) = <paramref name="p"/>(x) + <paramref name="q"/>(x)
        /// </returns>
        public static Polynomial operator +(Polynomial p, Polynomial q) {
            Polynomial sumPolynomial = p.CloneAs();
            int D = p.Exponents.GetLength(1);

            for (int i = 0; i < q.Coeff.Length; i++) {
                int[] exponent = new int[D];
                for (int d = 0; d < D; d++) {
                    exponent[d] = q.Exponents[i, d];
                }

                int index = sumPolynomial.GetCoefficientIndex(exponent);
                if (index < 0) {
                    sumPolynomial.AddCoeff(q.Coeff[i], exponent);
                } else {
                    sumPolynomial.Coeff[index] += q.Coeff[i];
                }
            }

            return sumPolynomial;
        }

        /// <summary>
        /// Calculates the difference between <paramref name="p"/> and
        /// <paramref name="q"/>
        /// </summary>
        /// <param name="p">First operand</param>
        /// <param name="q">First operand</param>
        /// <returns>
        /// A new polynomial r with
        /// r(x) = <paramref name="p"/>(x) - <paramref name="q"/>(x)
        /// </returns>
        public static Polynomial operator -(Polynomial p, Polynomial q) {
            return p + (-q);
        }


        /// <summary>
        /// Computes the derivative of a monomial.
        /// </summary>
        /// <param name="MonomExp">Monomial exponents.</param>
        /// <param name="DerivExp">Derivative symbol exponents.</param>
        static Tuple<int,int[]> DeriveMonomial(int[] MonomExp, int[] DerivExp) {
            if(MonomExp.Length != DerivExp.Length)
                throw new ArgumentException();
            int D = MonomExp.Length;
            

            for(int d = 0; d < D; d++) {
                if(MonomExp[d] < 0)
                    throw new ArgumentException("All monomial exponents must be non-negative.");
                if(DerivExp[d] < 0)
                    throw new ArgumentException("All derivative exponents must be non-negative.");
            }

            int ret = 1;
            int[] RetExp = DerivExp.CloneAs();

            for(int d = 0; d < D; d++) {
                int n = DerivExp[d];
                int m = MonomExp[d];

                if(n > m) {
                    ret = 0;
                    break;
                } else {

                    for(int i = m; i >= (m - n + 1); i--) {
                        ret *= i;
                    }

                    RetExp[d] = m - n;
                }

            }

            return new Tuple<int, int[]>(ret, RetExp);
        }


        /// <summary>
        /// Returns a derivative of this polynomial.
        /// </summary>
        /// <param name="Deriv">Derivative symbol exponents.</param>
        public Polynomial Derive(params int[] Deriv) {
            Polynomial R = new Polynomial();

            Debug.Assert(this.Coeff.Length == this.Exponents.GetLength(0));
            for(int i = 0; i < this.Coeff.Length; i++) {
                var MonomDeriv = DeriveMonomial(this.Exponents.GetRow(i), Deriv);
                if(MonomDeriv.Item1 != 0) {
                    R.AddCoeff(this.Coeff[i] * MonomDeriv.Item1, MonomDeriv.Item2);
                }
            }

            if(R.Coeff.Length == 0) {
                // empty polynomial: add dummy coefficient 
                R.AddCoeff(0.0, new int[this.SpatialDimension]);
            }

            R.Collect();
            
            return R;
        }

        private void Collect() {
            Debug.Assert(this.Coeff.Length == this.Exponents.GetLength(0));
            for(int i = 0; i < this.Coeff.Length; i++) {
                int[] Exp_i = this.Exponents.GetRow(i);
                Debug.Assert(Exp_i.Length == this.SpatialDimension);

                for(int j = i + 1; j < this.Coeff.Length; j++) {
                    int[] Exp_j = this.Exponents.GetRow(j);
                    Debug.Assert(Exp_j.Length == this.SpatialDimension);

                    if(ArrayTools.Equals(Exp_i, Exp_j)) {
                        //Console.WriteLine("gtcha");
                        throw new NotImplementedException("todo");
                    } else {
                        //
                    }

                }
            }

        }

        
        /*

        /// <summary>
        /// a vectorized evaluation of this polynomial
        /// </summary>
        /// <param name="result">
        /// Dimension must be 1;
        /// On exit, the k-th entry contains the value of the polynomial 
        /// at the point <paramref name="LocalPoints"/>[k,:];
        /// </param>
        /// <param name="LocalPoints">
        /// Points at which the polynomial should be evaluated;
        /// 1st index: Point index; 2nd index: spatial coordinate index,
        /// 0,1 for 2D and 0,1,2 for 3D;
        /// </param>
        public void Evaluate(MultidimensionalArray result, MultidimensionalArray LocalPoints) {
            Evaluate(result, LocalPoints, GetMonomials(LocalPoints));
        }

         */


        /// <summary>
        /// A vectorized evaluation of this polynomial at specific nodes.
        /// </summary>
        /// <param name="result">
        /// Dimension must be 1;
        /// On exit, the k-th entry contains the value of the polynomial 
        /// at the point <paramref name="LocalPoints"/>[k,:];
        /// </param>
        /// <param name="LocalPoints">
        /// Points at which the polynomial should be evaluated;
        /// 1st index: Point index; 2nd index: spatial coordinate index,
        /// 0,1 for 2D and 0,1,2 for 3D;
        /// </param>
        public void Evaluate(MultidimensionalArray result, NodeSet LocalPoints) {
            if (result.Dimension != 1) {
                throw new ArgumentException("dimension of result must be 1");
            }

            if (result.GetLength(0) != LocalPoints.GetLength(0)) {
                throw new ArgumentException("length of result array must be equal to number of coordinates.");
            }

            result.Clear();

            // Empty polynomial, always zero
            if (this.Coeff.Length == 0) {
                return;
            }

            if (LocalPoints.GetLength(1) != this.SpatialDimension) {
                throw new ArgumentException("wrong spatial dimension;");
            }

            int iE = Coeff.Length;
            int L = result.GetLength(0);
            int D = SpatialDimension;

            MultidimensionalArray monomials = Caching.MonomialCache.Instance.GetMonomials(LocalPoints, this.AbsoluteDegree);

            // Required because the number of $monomials per point may be
            // higher than required by the $Degree of this polynomial
            int maxMonomialExponent = monomials.GetLength(2);

            unsafe {
                fixed (double* pResult = result.Storage, pCoeff = Coeff, pMonomials = monomials.Storage) {
                    fixed (int* pExponents = Exponents) {
                        // loop over local points ...
                        double* pMonomialsCur = pMonomials;
                        for (int j = 0; j < L; j++) {
                            // Beware: $result is not necessarily continuous so
                            // we cannot safely optimize this pointer evaluation
                            double* pResultCur = pResult + result.Index(j);

                            // loop over all exponents (with nonzero coefficients) ...
                            double* pCoeffCur = pCoeff;
                            int* pExponentsCur = pExponents;
                            double acc = 0;
                            for (int i = 0; i < iE; i++) {
                                double monom = 1.0;

                                // loop over space dimensions ...
                                for (int k = 0; k < D; k++) {
                                    monom *= *(pMonomialsCur + k * maxMonomialExponent + *(pExponentsCur++));
                                }

                                acc += monom * *(pCoeffCur++);
                            }
                            *pResultCur = acc;

                            pMonomialsCur += D * maxMonomialExponent;
                        }
                    }
                }
            }

            //// Reference implementation
            //// loop over local points ...
            //for (int j = 0; j < L; j++) {

            //    // loop over all exponents (with nonzero coeficients) ...
            //    for (int i = 0; i < iE; i++) {
            //        double monom = 1.0;

            //        // loop over space dimensions ...
            //        for (int k = 0; k < D; k++) {
            //            monom *= monomials[j, k, Exponents[i, k]];
            //        }

            //        monom *= Coeff[i];
            //        result[j] += monom;
            //    }
            //}
        }

       

        /// <summary>
        /// compares according to polynomial <see cref="AbsoluteDegree"/>;
        /// </summary>
        /// <param name="other"></param>
        /// <returns></returns>
        public int CompareTo(Polynomial other) {
            return (this.AbsoluteDegree - other.AbsoluteDegree);
        }
        
        /// <summary>
        /// type-safe cloning
        /// </summary>
        public Polynomial CloneAs() {
            Polynomial p = new Polynomial();
            p.Coeff = (double[])this.Coeff.Clone();
            p.Exponents = (int[,])this.Exponents.Clone();
            return p;
        }

        #region ICloneable Members

        /// <summary>
        /// cloning
        /// </summary>
        public object Clone() {
            return CloneAs();
        }

        #endregion

        /// <summary>
        /// If all exponents and coefficients are the same, two polynomials are equal.
        /// </summary>
        /// <param name="obj"></param>
        /// <returns></returns>
        public override bool Equals(object obj) {
            Polynomial other = obj as Polynomial;

            if(other == null)
                return false;

            int L = this.Coeff.Length;
            int D = this.SpatialDimension;
            Debug.Assert(this.Coeff.Length == this.Exponents.GetLength(0));
            Debug.Assert(other.Coeff.Length == other.Exponents.GetLength(0));

            if(L != other.Coeff.Length)
                return false;
            if(D != other.SpatialDimension)
                return false;

            if(!ArrayTools.AreEqual(this.Exponents, other.Exponents))
                return false;
            if(!ArrayTools.AreEqual(this.Coeff, other.Coeff))
                return false;

            return true;

        }

        /// <summary>
        /// Hash code based on the exponents.
        /// </summary>
        public override int GetHashCode() {
            int R = 879;
            for(int i = 0; i < this.Exponents.GetLength(0); i++)
                for(int j = 0; j < this.Exponents.GetLength(1); j++)
                    R += this.Exponents[i, j] * (i + 13);
            return R;
        }
    }
}
