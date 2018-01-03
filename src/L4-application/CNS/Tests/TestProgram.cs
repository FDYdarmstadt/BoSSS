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
using System.Globalization;
using System.Linq;
using System.Threading;
using BoSSS.Solution;
using ilPSP.Utils;
using NUnit.Framework;

namespace CNS.Tests {

    /// <summary>
    /// Abstract base class for tests.
    /// </summary>
    [TestFixture]
    public abstract class TestProgram<T> : Program<T>
        where T : CNSControl, new() {

        /// <summary>
        /// Performs bootstrapping.
        /// </summary>
        [TestFixtureSetUp]
        public static void SetUp() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                Application.GetBoSSSInstallDir(),
                out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        /// <summary>
        /// Checks the errors for density, pressure and entropy against predefined
        /// thresholds.
        /// </summary>
        /// <param name="resultsTable"></param>
        /// <param name="thresholdsTable"></param>
        protected static void CheckErrorThresholds(Dictionary<string, object> resultsTable, params Tuple<string, double>[] thresholdsTable) {
            List<Action> assertions = new List<Action>();
            foreach (var valueThresholdPair in thresholdsTable) {
                string queryName = valueThresholdPair.Item1;
                double threshold = valueThresholdPair.Item2;

                double error = (double)resultsTable[queryName];
                string message = String.Format(
                    "{0}: {1} (Threshold is {2})",
                    //"{0:F16}: {1:F16} (Threshold is {2:F16})",
                    queryName,
                    error,
                    threshold);
                Console.WriteLine(message);

                assertions.Add(() => Assert.IsTrue(error < threshold, message));
            }

            assertions.ForEach(a => a());
        }

        protected static int GetTimeStepNumber(Program solver) {
            solver.QueryResultTable.FormatTable(out string[] KeyColumnNames, out string[] ValueColumnNames, out object[,] KeyTable, out object[,] ValueTable);
            int iCol = Array.IndexOf(KeyColumnNames, "Timestep");
            int iTimeStep = (int)(KeyTable.GetColumn(iCol).Last());

            return iTimeStep;
        }
    }
}
