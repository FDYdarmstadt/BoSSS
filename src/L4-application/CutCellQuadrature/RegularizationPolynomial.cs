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
using BoSSS.Platform;
using ilPSP;

namespace CutCellQuadrature {

    /// <summary>
    /// Base class for regularizations of the Heaviside function and the delta
    /// distribution.
    /// </summary>
    /// <remarks>
    /// For the definition and the meaning of the parameters, see Tornberg2000.
    /// </remarks>
    public abstract class RegularizationPolynomoial {

        /// <summary>
        /// The actual polynomial being represented
        /// </summary>
        protected Func<double, double> polynomial;

        /// <summary>
        /// Constructs a new polynomial with the given properties.
        /// </summary>
        /// <param name="polynomial"></param>
        /// <param name="numberOfVanishingMoments"></param>
        /// <param name="numberOfContinuousDerivatives"></param>
        protected RegularizationPolynomoial(Func<double, double> polynomial, int numberOfVanishingMoments, int numberOfContinuousDerivatives) {
            this.polynomial = polynomial;
            this.NumberOfVanishingMoments = numberOfVanishingMoments;
            this.NumberOfContinuousDerivatives = numberOfContinuousDerivatives;
        }

        /// <summary>
        /// The number of vanishing moments the polynomial has
        /// </summary>
        public int NumberOfVanishingMoments {
            get;
            private set;
        }

        /// <summary>
        /// The number of continuous derivatives the polynomial has.
        /// </summary>
        public int NumberOfContinuousDerivatives {
            get;
            private set;
        }

        /// <summary>
        /// The order of the polynomial
        /// </summary>
        public abstract int Order {
            get;
        }

        /// <summary>
        /// Modulates the given function <paramref name="values"/> by multiplying
        /// the specific values with the value of regularization polynomial.
        /// </summary>
        /// <param name="regularizationWidth">
        /// The width of the regularization region.
        /// </param>
        /// <param name="values">
        /// The values to be modulated.
        /// </param>
        /// <param name="distances">
        /// The distances of the evaluation points (cf.
        /// <paramref name="values"/>) from the zero iso-contour.
        /// </param>
        public void Regularize(double regularizationWidth, MultidimensionalArray values, MultidimensionalArray distances) {
            int noOfCells = values.GetLength(0);
            int noOfNodesPerCell = values.GetLength(1);

            for (int i = 0; i < noOfCells; i++) {
                for (int j = 0; j < noOfNodesPerCell; j++) {
                    values[i, j] *= this.Evaluate(distances[i, j], regularizationWidth);
                }
            }
        }

        /// <summary>
        /// Variant of
        /// <see cref="Regularized(double, MultidimensionalArray, MultidimensionalArray)"/>
        /// that accounts for the gradient of the level set function
        /// </summary>
        /// <param name="regularizationWidth">
        /// <see cref="Regularized(double, MultidimensionalArray, MultidimensionalArray)"/>
        /// </param>
        /// <param name="values">
        /// <see cref="Regularized(double, MultidimensionalArray, MultidimensionalArray)"/>
        /// </param>
        /// <param name="distances">
        /// <see cref="Regularized(double, MultidimensionalArray, MultidimensionalArray)"/>
        /// </param>
        /// <param name="gradients">
        /// The values of the gradient of the level set function at the
        /// evaluation nodes
        /// </param>
        /// <remarks>
        /// See Tornberg2000
        /// </remarks>
        public void Regularize(double regularizationWidth, MultidimensionalArray values, MultidimensionalArray distances, MultidimensionalArray gradients) {
            int noOfCells = values.GetLength(0);
            int noOfNodesPerCell = values.GetLength(1);
            int spatialDimension = gradients.GetLength(2);

            for (int i = 0; i < noOfCells; i++) {
                for (int j = 0; j < noOfNodesPerCell; j++) {
                    double modifiedWidth = regularizationWidth;

                    double gradientL1Norm = 0.0;
                    double gradientL2Norm = 0.0;
                    for (int k = 0; k < spatialDimension; k++) {
                        double c = gradients[i, j, k];
                        gradientL1Norm += Math.Abs(c);
                        gradientL2Norm += c * c;
                    }
                    modifiedWidth *= gradientL1Norm / Math.Sqrt(gradientL2Norm);

                    values[i, j] *= this.Evaluate(distances[i, j], modifiedWidth);
                }
            }
        }

        /// <summary>
        /// Implement this method in order to account for the domain of
        /// definition of the regularization polynomial
        /// </summary>
        /// <param name="distance">
        /// The distance from the zero iso-contour
        /// </param>
        /// <param name="width">
        /// The width of the regularization region
        /// </param>
        /// <returns>
        /// The value of the polynomial.
        /// </returns>
        public abstract double Evaluate(double distance, double width);
    }
}
