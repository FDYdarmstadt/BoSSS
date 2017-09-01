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

namespace CutCellQuadrature {

    /// <summary>
    /// A polynomial regularization of the delta distribution.
    /// </summary>
    /// <remarks>
    /// For the definition and the meaning of the parameters, see Tornberg2000.
    /// </remarks>
    public class DeltaPolynomial : RegularizationPolynomoial {

        /// <summary>
        /// All the polynomials.
        /// </summary>
        private static DeltaPolynomial[,] polynomials = new DeltaPolynomial[,] {
            { // One vanishing moment
                new DeltaPolynomial(
                    (xi) => -3.0 / 4.0 * xi * xi + 3.0 / 4.0,
                    1,
                    0),
                new DeltaPolynomial(
                    (xi) => 15.0 / 16.0 + (-15.0 / 8.0 + 15.0 / 16.0 * xi * xi) * xi * xi,
                    1,
                    1),
                new DeltaPolynomial(
                    (xi) => 35.0 / 32.0 + (-105.0 / 32.0 + (105.0 / 32.0 - 35.0 / 32.0 * xi * xi) * xi * xi) * xi * xi,
                    1,
                    2),
                new DeltaPolynomial(
                    (xi) => 315.0 / 256.0 + (-315.0 / 64.0 + (945.0 / 128.0 + (-315.0 / 64.0 + 315.0 / 256.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    1,
                    3)
            },
            { // Three vanishing moments
                new DeltaPolynomial(
                    (xi) => 45.0 / 32.0 + (-75.0 / 16.0 + 105.0 / 32.0 * xi * xi) * xi * xi,
                    3,
                    0),
                new DeltaPolynomial(
                    (xi) => 105.0 / 64.0 + (-525.0 / 64.0 + (735.0 / 64.0 - 315.0 / 64.0 * xi * xi) * xi * xi) * xi * xi,
                    3,
                    1),
                new DeltaPolynomial(
                    (xi) => 945.0 / 512.0 + (-1575.0 / 128.0 + (6615.0 / 256.0 + (-2835.0 / 128.0 + 3465.0 / 512.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    3,
                    2),
                new DeltaPolynomial(
                    (xi) => 2079.0 / 1024.0 + (-17325.0 / 1024.0 + (24255.0 / 512.0 + (-31185.0 / 512.0 + (38115.0 / 1024.0 - 9009.0 / 1024.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    3,
                    3)
            },
            { // Five vanishing moments
                new DeltaPolynomial(
                    (xi) => 525.0 / 256.0 + (-3675.0 / 256.0 + (6615.0 / 256.0 - 3465.0 / 256.0 * xi * xi) * xi * xi) * xi * xi,
                    5,
                    0),
                new DeltaPolynomial(
                    (xi) => 4725.0 / 2048.0 + (-11025.0 / 512.0 + (59535.0 / 1024.0 + (-31185.0 / 512.0 + 45045.0 / 2048.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    5,
                    1),
                new DeltaPolynomial(
                    (xi) => 10395.0 / 4096.0 + (-121275.0 / 4096.0 + (218295.0 / 2048.0 + (-343035.0 / 2048.0 + (495495.0 / 4096.0 - 135135.0 / 4096.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    5,
                    2),
                new DeltaPolynomial(
                    (xi) => 45045.0 / 16384.0 + (-315315.0 / 8192.0 + (2837835.0 / 16384.0 + (-1486485.0 / 4096.0 + (6441435.0 / 16384.0 + (-1756755.0 / 8192.0 + 765765.0 / 16384.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi,
                    5,
                    3)
            }
        };

        /// <summary>
        /// Constructs a new polynomial
        /// </summary>
        /// <param name="polynomial"></param>
        /// <param name="numberOfVanishingMoments"></param>
        /// <param name="numberOfContinuousDerivatives"></param>
        private DeltaPolynomial(Func<double, double> polynomial, int numberOfVanishingMoments, int numberOfContinuousDerivatives)
            : base(polynomial, numberOfVanishingMoments, numberOfContinuousDerivatives) {
        }

        /// <summary>
        /// The order of the polynomial.
        /// </summary>
        public override int Order {
            get {
                return 2 * (NumberOfVanishingMoments / 2 + NumberOfContinuousDerivatives + 1);
            }
        }

        /// <summary>
        /// Evaluates the regularized polynomial
        /// </summary>
        /// <param name="distance">
        /// The signed distance from the zero iso-contour
        /// </param>
        /// <param name="width">
        /// The width of the regularization region
        /// </param>
        /// <returns>
        /// <list type="bullet">
        ///     <item>
        ///     Zero, if |<paramref name="distance"/>| &gt; <paramref name="width"/>
        ///     </item>
        ///     <item>
        ///     Otherwise, the actual value of the polynomial.
        ///     </item>
        /// </list>
        /// </returns>
        public override double Evaluate(double distance, double width) {
            if (Math.Abs(distance) > width) {
                return 0.0;
            } else {
                return polynomial(distance / width) / width;
            }
        }

        /// <summary>
        /// Retrieves the polynomial with the selected properties
        /// </summary>
        /// <param name="numberOfVanishingMoments">
        /// The number of vanishing moments the polynomial should have.
        /// </param>
        /// <param name="numberOfContinuousDerivatives">
        /// The number of continuous derivatives of the polynomial
        /// </param>
        /// <returns>
        /// A regularization polynomial
        /// </returns>
        public static DeltaPolynomial GetPolynomial(int numberOfVanishingMoments, int numberOfContinuousDerivatives) {
            return polynomials[numberOfVanishingMoments / 2, numberOfContinuousDerivatives];
        }
    }
}
