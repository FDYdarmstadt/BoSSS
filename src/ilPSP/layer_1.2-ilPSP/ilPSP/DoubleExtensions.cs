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
using System.Linq;
using System.Text;
using ilPSP.Utils;

namespace ilPSP {

    /// <summary>
    /// Extension methods for <see cref="double"/> values.
    /// </summary>
    public static class DoubleExtensions {

        /// <summary>
        /// to string with invariant number formatting;
        /// (see <see cref="System.Globalization.NumberFormatInfo.InvariantInfo"/>)
        /// </summary>
        static public string ToStringDot(this double d) {
            return d.ToString(System.Globalization.NumberFormatInfo.InvariantInfo);
        }

        /// <summary>
        /// to string with invariant number formatting;
        /// (see <see cref="System.Globalization.NumberFormatInfo.InvariantInfo"/>)
        /// </summary>
        static public string ToStringDot(this double d, string format) {
            return d.ToString(format, System.Globalization.NumberFormatInfo.InvariantInfo);
        }

        /// <summary>
        /// <see cref="Math.Abs(double)"/>
        /// </summary>
        /// <param name="d"><see cref="Math.Abs(double)"/></param>
        /// <returns><see cref="Math.Abs(double)"/></returns>
        public static double Abs(this double d) {
            return Math.Abs(d);
        }

        /// <summary>
        /// <see cref="double.IsNaN"/>
        /// </summary>
        /// <param name="d"><see cref="double.IsNaN"/></param>
        /// <returns><see cref="double.IsNaN"/></returns>
        public static bool IsNaN(this double d) {
            return double.IsNaN(d);
        }

        /// <summary>
        /// <see cref="Math.Sign(double)"/>
        /// </summary>
        /// <param name="d"><see cref="Math.Sign(double)"/></param>
        /// <returns><see cref="Math.Sign(double)"/></returns>
        public static int Sign(this double d) {
            return Math.Sign(d);
        }

        /// <summary>
        /// Heaviside
        /// </summary>
        public static int Heaviside(this double d) {
            if (d>0)
                return 1;
            else
                return 0;
        }

        /// <summary>
        /// <see cref="Math.Sqrt(double)"/>
        /// </summary>
        /// <param name="d"><see cref="Math.Sqrt(double)"/></param>
        /// <returns><see cref="Math.Sqrt(double)"/></returns>
        public static double Sqrt(this double d) {
            return Math.Sqrt(d);
        }

        /// <summary>
        /// the square of <paramref name="d"/>
        /// </summary>
        public static double Pow2(this double d) {
            return d * d;
        }

        /// <summary>
        /// logarithm at basis 10 for <paramref name="d"/>
        /// </summary>
        public static double Log10(this double d) {
            return Math.Log10(d);
        }

        /// <summary>
        /// logarithmus naturalis for see <paramref name="d"/>
        /// </summary>
        public static double Log(this double d) {
            return Math.Log(d, Math.E);
        }

        /// <summary>
        /// <see cref="Math.Pow"/>
        /// </summary>
        /// <param name="d"><see cref="Math.Pow"/></param>
        /// <param name="exponent"><see cref="Math.Pow"/></param>
        /// <returns><see cref="Math.Pow"/></returns>
        public static double Pow(this double d, double exponent) {
            return Math.Pow(d, exponent);
        }

        /// <summary>
        /// Checks <paramref name="d1"/> and <paramref name="d2"/> for equality
        /// within the given tolerance <paramref name="tol"/>.
        /// </summary>
        /// <param name="d1">First operand</param>
        /// <param name="d2">Second operand</param>
        /// <param name="tol">Tolerance</param>
        /// <returns>
        /// <see cref="Math.Abs(double)"/>(<paramref name="d1"/>, <paramref name="d2"/>)
        /// &lt;=<paramref name="tol"/>
        /// </returns>
        public static bool RoughlyEquals(this double d1, double d2, double tol) {
            return ((d1 - d2).Abs() <= tol);
        }

        /// <summary>
        /// Checks whether the elements of <paramref name="d1"/> and
        /// <paramref name="d2"/> are equal within the given tolerance
        /// <paramref name="tol"/>.
        /// </summary>
        /// <param name="d1">First sequence</param>
        /// <param name="d2">Second sequence</param>
        /// <param name="tol">Tolerance</param>
        /// <returns>
        /// True, if <paramref name="d1"/> and <paramref name="d2"/> are
        /// not null, have the same length and if the associated elements of
        /// both sequences are roughly equal (see
        /// <see cref="RoughlyEquals(double, double, double)"/>)
        /// </returns>
        public static bool RoughlyEquals(this IEnumerable<double> d1, IEnumerable<double> d2, double tol) {
            if (d1 == null || d2 == null) {
                return false;
            }

            int count = d1.Count();
            if (count != d2.Count()) {
                throw new ArgumentException("Both sequences must have the same length");
            }

            using (var d1Enumerator = d1.GetEnumerator()) {
                using (var d2Enumerator = d2.GetEnumerator()) {
                    while (d1Enumerator.MoveNext()) {
                        d2Enumerator.MoveNext();

                        if (!d1Enumerator.Current.RoughlyEquals(d2Enumerator.Current, tol)) {
                            return false;
                        }
                    }
                }
            }

            return true;
        }

        /// <summary>
        /// Calculates the L2-norm of the elements of
        /// <paramref name="source"/>.
        /// </summary>
        /// <param name="source">
        /// A generic sequence of numbers.
        /// </param>
        /// <returns>
        /// \f$ \sqrt{\sum \limits_{e \in S} e^2}\f$  where S is
        /// <paramref name="source"/>
        /// </returns>
        public static double L2Norm(this IEnumerable<double> source) {
            if (source.IsNullOrEmpty()) {
                throw new ArgumentException("Sequence empty");
            }

            double normSquared = 0.0;
            using (var enumerator = source.GetEnumerator()) {
                while (enumerator.MoveNext()) {
                    normSquared += enumerator.Current * enumerator.Current;
                }
            }
            return normSquared.Sqrt();
        }

        /// <summary>
        /// normalizes some vector 
        /// </summary>
        public static T Normalize<T>(this T vec) where T : IList<double> {
            double ooN = 1.0/GenericBlas.L2Norm(vec);
            int N = vec.Count;
            for (int i = 0; i < N; i++) vec[i] *= ooN;
            return vec;
        }


        /// <summary>
        /// Calculates the distance between the given sequences in the L2-norm.
        /// </summary>
        /// <param name="source">
        /// A generic sequence of numbers.
        /// </param>
        /// <param name="other">
        /// Another generic sequence of numbers. Must have the same length as
        /// <paramref name="source"/>.
        /// </param>
        /// <returns>
        /// \f$ \sqrt{\sum \limits (s_i - o_i)^2}\f$ 
        /// </returns>
        public static double L2Distance(this IEnumerable<double> source, IEnumerable<double> other) {
            if (source == null) {
                throw new ArgumentException("Sequence must not be null", "source");
            }

            if (other == null) {
                throw new ArgumentException("Sequence must not be null", "other");
            }

            if (source.Count() != other.Count()) {
                throw new ArgumentException("Both sequences must have the same length");
            }

            double normSquared = 0.0;
            using (var sourceEnumerator = source.GetEnumerator()) {
                using (var otherEnumerator = other.GetEnumerator()) {

                    while (sourceEnumerator.MoveNext()) {
                        otherEnumerator.MoveNext();

                        double diff = (sourceEnumerator.Current - otherEnumerator.Current);
                        normSquared += diff * diff;
                    }
                }
            }

            return normSquared.Sqrt();
        }


       

        

    }
}
