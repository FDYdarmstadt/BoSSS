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
    /// A polynomial regularization of the Heaviside function.
    /// </summary>
    /// <remarks>
    /// For the definition and the meaning of the parameters, see Tornberg2000.
    /// </remarks>
    public class HeavisidePolynomial : RegularizationPolynomoial {

        /// <summary>
        /// All the polynomials.
        /// </summary>
        private static HeavisidePolynomial[,] polynomials = new HeavisidePolynomial[,] {
            { // Two vanishing moments
                new HeavisidePolynomial(
                    (xi) => 0.5 + (9.0 / 8.0 - 5.0 / 8.0 * xi * xi)*xi,
                    2,
                    0),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (45.0 / 32.0 + (-25.0 / 16.0 + 21.0 / 32.0 * xi * xi) * xi * xi) * xi,
                    2,
                    1),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (105.0 / 64.0 + (-175.0 / 64.0 + (147.0 / 64.0 - 45.0 / 64.0 * xi * xi) * xi * xi) * xi * xi) * xi,
                    2,
                    2),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (945.0 / 512.0 + (-525.0 / 128.0 + (1323.0 / 256.0 + (-405.0 / 128.0 + 385.0 / 512.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi ) * xi,
                    2,
                    3),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (2079.0 / 1024.0 + (-5775.0 / 1024.0 + (4851.0 / 512.0 + (-4455.0 / 512.0 + (4235.0 / 1024.0 - 819.0 / 1024.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    2,
                    4),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (9009.0 / 4096.0 + (-15015.0 / 2048.0 + (63063.0 / 4096.0 + (-19305.0 / 1024.0 + (55055.0 / 4096.0 + (-10647.0 / 2048.0 + 3465.0 / 4096.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    2,
                    5)
            },
            { // Four vanishing moments
                new HeavisidePolynomial(
                    (xi) => 0.5 + xi * (225.0 + xi * xi * (-350.0 + 189.0 * xi * xi)) / 128.0,
                    4,
                    0),
                new HeavisidePolynomial(
                    (xi) => 0.5 + xi * (525.0 + xi * xi * (-1225.0 + xi * xi * (1323.0 - 495.0 * xi * xi))) / 256.0,
                    4,
                    1),
                new HeavisidePolynomial(
                    (xi) => 0.5 + xi * (4725.0 + xi * xi * (-14700.0 + xi * xi * (23814.0 + xi * xi * (-17820.0 + 5005.0 * xi * xi)))) / 2048.0,
                    4,
                    2),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (10395.0 / 4096.0 + (-40425.0 / 4096.0 + (43659.0 / 2048.0 + (-49005.0 / 2048.0 + (55055.0 / 4096.0 - 12285.0 / 4096.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    4,
                    3),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (45045.0 / 16384.0 + (-105105.0 / 8192.0 + (567567.0 / 16384.0 + (-212355.0 / 4096.0 + (715715.0 / 16384.0 + (-159705.0 / 8192.0 + 58905.0 / 16384.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    4,
                    4),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (96525.0 / 32768.0 + (-525525.0 / 32768.0 + (1702701.0 / 32768.0 + (-3185325.0 / 32768.0 + (3578575.0 / 32768.0 + (-2395575.0 / 32768.0 + (883575.0 / 32768.0 - 138567.0 / 32768.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    4,
                    5)
            },
            { // Six vanishing moments
                new HeavisidePolynomial(
                    (xi) => 0.5 + (1225.0/512.0 + (-3675.0 / 512.0 + (4851.0 / 512.0 - 2145.0 / 512.0 * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    0),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (11025.0 / 4096.0 + (-11025.0 / 1024.0 + (43659.0 / 2048.0 + (-19305.0 / 1024.0 + 25025.0 / 4096.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    1),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (24255.0 / 8192.0 + (-121275.0 / 8192.0 + (160083.0 / 4096.0 + (-212355.0 / 4096.0 + (275275.0 / 8192.0 - 69615.0 / 8192.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    2),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (105105.0 / 32768.0 + (-315315.0 / 16384.0 + (2081079.0 / 32768.0 + (-920205.0 / 8192.0 + (3578575.0 / 32768.0 + (-904995.0 / 16384.0 + 373065.0 / 32768.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    3),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (225225.0 / 65536.0 + (-1576575.9 / 65536.0 + (6243237.0 / 65536.0 + (-13803075.0 / 65536.0 + (17892875.0 / 65536.0 + (-13574925.0 / 65536.0 + (5595975.0 / 65536.0 - 969969.0 / 65536.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    4),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (3828825.0 / 1048576.0 + (-3828825.0 / 131072.0 + (35378343.0 / 262144.0 + (-46930455.0 / 131072.0 + (304178875.0 / 524288.0 + (-76924575.0 / 131072.0 + (95131575.0 / 262144.0 + (-16489473.0 / 131072.0 + 19684665.0 / 1048576.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    6,
                    5)
            },
            { // Eight vanishing moments
                new HeavisidePolynomial(
                    (xi) => 0.5 + (99225.0 / 32768.0 + (-121275.0 / 8192.0 + (567567.0 / 16384.0 + (-289575.0 / 8192.0 + 425425.0 / 32768.0 * xi * xi)* xi * xi) * xi * xi) * xi * xi)*xi,
                    8,
                    0),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (218295.0 / 65536.0 + (-1334025.0 / 65536.0 + (2081079.0 / 32768.0 + (-3185325.0 / 32768.0 + (4679675.0 / 65536.0 - 1322685.0 / 65536.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    8,
                    1),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (945945.0 / 262144.0 + (-3468465.0 / 131072.0 + (27054027.0 / 262144.0 + (-13803075.0 / 65536.0 + (60835775.0 / 262144.0 + (-17194905.0 / 131072.0 + 7834365.0 / 262144.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    8,
                    2),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (2027025.0 / 524288.0 + (-17342325.0 / 524288.0 + (81162081.0 / 524288.0 + (-207046125.0 / 524288.0 + (304178875.0 / 524288.0 + (-257923575.0 / 524288.0 + (117515475.0 / 524288.0 - 22309287.0 / 524288.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    8,
                    3),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (34459425.0 / 8388608.0 + (-42117075.0 / 1048576.0 + (459918459.0 / 2097152.0 + (-703956825.0 / 1048576.0 + (5171040875.0 / 4194304.0 + (-1461566925.0 / 1048576.0 + (1997763075.0 / 2097152.0 + (-379257879.0 / 1048576.0 + 492116625.0 / 8388608.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    8,
                    4),
                new HeavisidePolynomial(
                    (xi) => 0.5 + (72747675.0 / 16777216.0 + (-800224425.0 / 16777216.0 + (1248350103.0 / 4194304.0 + (-4458393225.0 / 4194304.0 + (19649955325.0 / 8388608.0 + (-27769771575.0 / 8388608.0 + (12652499475.0 / 4194304.0 + (-7205899701.0 / 4194304.0 + (9350215875.0 / 16777216.0 - 1320944625.0 / 16777216.0 * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi * xi) * xi,
                    8,
                    5)
            }
        };

        /// <summary>
        /// Constructs a new polynomial
        /// </summary>
        /// <param name="polynomial"></param>
        /// <param name="numberOfVanishingMoments"></param>
        /// <param name="numberOfContinuousDerivatives"></param>
        private HeavisidePolynomial(Func<double, double> polynomial, int numberOfVanishingMoments, int numberOfContinuousDerivatives)
            : base(polynomial, numberOfVanishingMoments, numberOfContinuousDerivatives) {
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
        public static HeavisidePolynomial GetPolynomial(int numberOfVanishingMoments, int numberOfContinuousDerivatives) {
            return polynomials[(numberOfVanishingMoments - 1) / 2, numberOfContinuousDerivatives];
        }

        /// <summary>
        /// The order of the polynomial.
        /// </summary>
        public override int Order {
            get {
                return 1 + 2 * ((NumberOfVanishingMoments + 1) / 2 + NumberOfContinuousDerivatives);
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
        ///     Zero, if <paramref name="distance"/> &lt; -<paramref name="width"/>
        ///     </item>
        ///     <item>
        ///     One, if <paramref name="distance"/> &gt; <paramref name="width"/>
        ///     </item>
        ///     <item>
        ///     Otherwise, the actual value of the polynomial.
        ///     </item>
        /// </list>
        /// </returns>
        public override double Evaluate(double distance, double width) {
            if (distance <= -width) {
                return 0.0;
            } else if (distance >= width) {
                return 1.0;
            } else {
                return polynomial(distance / width);
            }
        }
    }
}
