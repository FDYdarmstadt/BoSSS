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

using BoSSS.Solution;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace CNS.Tests {

    static class TestUtils {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="table"></param>
        /// <param name="spatialDimension"></param>
        public static void CheckConvergenceRates(BoSSS.Foundation.IO.ISessionInfo table, int spatialDimension) {
            double eps = 0.12;

            foreach (var pair in ExtractOrders(table, "densityError", "dgDegree")) {
                Console.WriteLine("Density convergence rate for p={0}: {1:E16}", pair.Key, pair.Value);

                int pPlusOne = pair.Key + 1;
                Assert.That(
                    pair.Value >= pPlusOne - eps,
                    String.Format(
                        "Density convergence rate for p={0} too low (Should be {1} but is {2})",
                        pair.Key,
                        pPlusOne,
                        pair.Value));
            }

            foreach (var pair in ExtractOrders(table, "momentumError", "dgDegree")) {
                Console.WriteLine("Momentum convergence rate for p={0}: {1:E16}", pair.Key, pair.Value);

                int pPlusOne = pair.Key + 1;
                Assert.That(
                    pair.Value >= pPlusOne - eps,
                    String.Format(
                        "Momentum convergence rate for p={0} too low (Should be {1} but is {2})",
                        pair.Key,
                        pPlusOne,
                        pair.Value));
            }

            foreach (var pair in ExtractOrders(table, "energyError", "dgDegree")) {
                Console.WriteLine("Energy convergence rate for p={0}: {1:E16}", pair.Key, pair.Value);

                int pPlusOne = pair.Key + 1;
                Assert.That(
                    pair.Value >= pPlusOne - eps,
                    String.Format(
                        "Energy convergence rate for p={0} too low (Should be {1} but is {2})",
                        pair.Key,
                        pPlusOne,
                        pair.Value));
            }
        }

        private static IEnumerable<KeyValuePair<int, double>> ExtractOrders(BoSSS.Foundation.IO.ISessionInfo table, string errorKey, string degreeKey) {
            return table.KeysAndQueries[errorKey].
                GroupBy(keyResultPair => (int)keyResultPair.Key[degreeKey]).
                Select(g => new KeyValuePair<int, double>(g.Key,
                    Regression(g.OrderBy(v => (int)v.Key["divisions"]).Select(v => (double)v.Value))));
        }

        private static double Regression(IEnumerable<double> values) {
            double[] yValues = values.Select(y => Math.Log10(y)).ToArray();
            double yAvg = yValues.Average();

            IEnumerable<double> abscissae = Enumerable.Range(0, yValues.Length).Select(i => Math.Pow(2.0, -i));
            double[] xValues = abscissae.Select(x => Math.Log10(x)).ToArray();
            double xAvg = xValues.Average();

            double v1 = 0.0;
            double v2 = 0.0;

            for (int i = 0; i < yValues.Length; i++) {
                v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                v2 += Math.Pow(xValues[i] - xAvg, 2);
            }

            return v1 / v2;
        }
    }
}
