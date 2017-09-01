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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using NUnit.Framework;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;

namespace QuadratureAndProjectionTest {

    /// <summary>
    /// The base class for quadrature and projection tests. Derive from this
    /// class and implement both abstract methods in order to test another
    /// element type
    /// </summary>
    [TestFixture]
    public abstract class QuadratueAndProjectionTest {

        /// <summary>
        /// The BoSSS context which we need in order to compute a DG projection
        /// </summary>
        private GridData context;

        /// <summary>
        /// The test application to be run
        /// </summary>
        private TestApplication testApplication;

        /// <summary>
        /// Sets up the current test application and the corresponding context
        /// </summary>
        [TestFixtureSetUp]
        public void SetUp() {
            testApplication = new TestApplication(this);
            context = testApplication.GridData;
        }

        /// <summary>
        /// Tests whether the nominal order of quadrature rules returned by
        /// <see cref="RefElement.GetQuadratureRule"/> is less than or equal to
        /// the requested order.
        /// </summary>
        [Test]
        public void TestNominalOrder() {
            var simplex = GetSimplex();

            int highestKnownOrder = simplex.HighestKnownOrder;

            for (int i = 0; i <= highestKnownOrder; i++) {
                QuadRule rule = simplex.GetQuadratureRule(i);
                int realOrder = rule.OrderOfPrecision;
                Assert.IsTrue(
                    realOrder >= i,
                    "Requested order < real order (" + realOrder + ")");
            }
        }

        /// <summary>
        /// Tests whether the sum of weights for each rule (approximately)
        /// equals the volume of the reference element.
        /// </summary>
        [Test]
        public void TestSumOfWeights() {
            var simplex = GetSimplex();

            int highestKnownOrder = simplex.HighestKnownOrder;

            QuadRule lastQuadRule = null;
            for (int i = 0; i <= highestKnownOrder; i++) {
                QuadRule rule = simplex.GetQuadratureRule(i);
                if (lastQuadRule == rule) {
                    // Same rule, nothing new to test
                    continue;
                }

                TestSumOfWeights(i, rule);
                lastQuadRule = rule;
            }
        }

        /// <summary>
        /// Tests whether the mass matrix following from
        /// <see cref="RefElement.OrthonormalPolynomials"/> is equal to the
        /// identity matrix
        /// </summary>
        [Test]
        public void TestMassMatrix() {
            var simplex = GetSimplex();

            // Limit order for run-time reasons
            int highestKnownOrder = Math.Min(simplex.HighestKnownOrder, 21);
             
            QuadRule lastQuadRule = null;
            for (int i = 0; i <= highestKnownOrder; i++) {
                QuadRule rule = simplex.GetQuadratureRule(i);
                if (lastQuadRule == rule) {
                    // Same rule, nothing new to test
                    continue;
                }

                TestMassMatrix(i, rule, simplex.GetOrthonormalPolynomials(Math.Min(highestKnownOrder, simplex.HighestSupportedPolynomialDegree)));
                lastQuadRule = rule;
            }
        }

        /// <summary>
        /// Tests the error of the projection of the orthonormal polynomials
        /// onto the reference element (i.e., the error should be zero)
        /// </summary>
        [Test]
        public void TestProjectionError() {
            var simplex = GetSimplex();

            int maxOrder = Math.Min(
                simplex.HighestKnownOrder,
                simplex.HighestSupportedPolynomialDegree);

            //In a field of order i, the a field automatically uses a
            //quadrature rule of order 2i+1 to compute the DG projection.
            //Thus, we can only (and need only to) orders i for which it is
            //guaranteed that a rule of order 2i-1 exists, too
            for (int i = 0; 2 * i + 1 <= maxOrder; i++) {
                QuadRule rule = simplex.GetQuadratureRule(i);
                TestProjectionError(i, simplex.GetOrthonormalPolynomials(maxOrder), rule);
            }
        }

        /// <summary>
        /// Implement this method by simply returning an instance of the
        /// simplex whose quadrature you want to test.
        /// </summary>
        /// <returns>
        /// A specific instance of a simplex, e.g. a <see cref="Line"/>
        /// </returns>
        protected abstract RefElement GetSimplex();

        /// <summary>
        /// Implement this method by generating a <b>one-cell</b> grid with
        /// correct dimensions (e.g. a <see cref="Grid1D"/> in case of a
        /// <see cref="Line"/>).
        /// </summary>
        /// <returns>A one-cell grid</returns>
        public abstract GridCommons GetSingleCellGrid();

        /// <summary>
        /// Sums up all weights associated with a quadrature rule and reports
        /// an error if the sum is not equal to the volume
        /// (<see cref="RefElement.Volume"/>) of the simplex.
        /// </summary>
        /// <param name="order">The <b>requested</b> order</param>
        /// <param name="rule">The quadrature rule</param>
        private void TestSumOfWeights(int order, QuadRule rule) {
            double error = Math.Abs(rule.Weights.Sum() - GetSimplex().Volume);

            // This error should be much smaller (and was much smaller on the
            // old master branch); Should be investigated.
            Assert.That(
                error <= 1e-9,
                "Sum of weights != Simplex volume (Deviance: " + String.Format("{0:e}", error) + ")");
        }

        /// <summary>
        /// Computes the matrix matrix
        /// V_ij = sum_k{p_i(nodes(k) * p_j(nodes(k)) * weights(k)}
        /// for a given quadrature rule (where p_i and p_j are the orthonormal
        /// polynomials supplied by the quadrature rule) and reports an error
        /// if this matrix differs from the identity matrix
        /// </summary>
        /// <param name="order">The <b>requested</b> order</param>
        /// <param name="rule">The quadrature rule</param>
        /// <param name="polynomials">The orthonormal polynomials</param>
        private void TestMassMatrix(int order, QuadRule rule, PolynomialList polynomials) {
            MultidimensionalArray[] values = new MultidimensionalArray[polynomials.Count];
            for (int i = 0; i < polynomials.Count; i++) {
                values[i] = MultidimensionalArray.Create(rule.NoOfNodes);
                polynomials[i].Evaluate(values[i], rule.Nodes);
            }

            //For every exactly integrable polynomial i
            for(int i = 0; i < polynomials.Count; i++) {
                // Check if quadrature rule is suitable
                if (polynomials[i].AbsoluteDegree > order) {
                    continue;
                }

                //For every exactly integrable co-polynomial j
                for (int j = 0; j <= i; j++) {
                    // Check if quadrature rule is suitable
                    if (polynomials[i].AbsoluteDegree + polynomials[j].AbsoluteDegree > order) {
                        continue;
                    }

                    double result = 0;
                    for (int k = 0; k < rule.NoOfNodes; k++) {
                        result += values[i][k] * values[j][k] * rule.Weights[k];
                    }

                    //Elements should be zero except of diagonal elements
                    double expectedResult = 0.0;
                    if (i == j) {
                        expectedResult = 1.0;
                    }

                    //Store maximum absolute error for the report we want to create
                    double error = Math.Abs(result - expectedResult);
                    //Console.WriteLine("MassMatrixError: {0}", error);

                    Assert.That(
                        error <= 1e-9,
                        String.Format(
                            "MassMatrix[" + i + "," + j + "] should be {0} but is off by {1:e} for order {2}",
                            expectedResult,
                            error,
                            order));
                }
            }
        }

        /// <summary>
        /// Projects every basis polynomial for the quadrature rule of the
        /// given order onto a single cell grid supplied by the sub-class and
        /// approximates the LInfinity error of the projection.
        /// </summary>
        /// <remarks>
        /// The DG projection uses a quadrature rule of 2*n+1 for a basis of
        /// order n. Thus, make sure <paramref name="order"/> that a quadrature
        /// rule of order 2*order+1 exists before calling this method</remarks>
        /// <param name="order">The <b>requested</b> order</param>
        /// <param name="polynomials">The orthonormal polynomials</param>
        /// <param name="rule">The tested quadrature rule</param>
        private void TestProjectionError(int order, PolynomialList polynomials, QuadRule rule) {
            Basis basis = new Basis(context, order);

            for (int i = 0; i < polynomials.Count; i++) {
                if (2 * polynomials[i].AbsoluteDegree + 1 > order) {
                    // Not integrable _exactly_ with this rule
                    continue;
                }

                DGField field = new SinglePhaseField(basis);
                ScalarFunction function = delegate(MultidimensionalArray input, MultidimensionalArray output) {
                    polynomials[i].Evaluate(output, new NodeSet(rule.RefElement, input));
                };
                field.ProjectField(function);

                double LInfError = 0;
                for (int j = 0; j < rule.NoOfNodes; j++) {
                    MultidimensionalArray result = MultidimensionalArray.Create(1, 1);
                    NodeSet point = new NodeSet(context.Cells.GetRefElement(0), 1, rule.SpatialDim);
                    for (int k = 0; k < rule.SpatialDim; k++) {
                        point[0, k] = rule.Nodes[j, k];
                    }
                    point.LockForever();

                    //Evaluate pointwise because we don't want the accumulated
                    //result
                    field.Evaluate(0, 1, point, result, 0.0);
                    
                    MultidimensionalArray exactResult = MultidimensionalArray.Create(1);
                    polynomials[i].Evaluate(exactResult, point);

                    LInfError = Math.Max(Math.Abs(result[0, 0] - exactResult[0]), LInfError);
                }
                //Console.WriteLine("ProjectionError: {0}", LInfError);

                Assert.IsTrue(
                    LInfError <= 1e-13,
                    String.Format(
                        "Projection error too high for polynomial {0}. Error: {1:e}",
                        i,
                        LInfError));
            }
        }
    }
}