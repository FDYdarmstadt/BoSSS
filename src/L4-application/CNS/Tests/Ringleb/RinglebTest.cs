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
using BoSSS.Foundation.Grid;
using BoSSS.Solution;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using NUnit.Framework;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace CNS.Tests.Ringleb {

    /// <summary>
    /// Contains tests related to the Ringleb test case with the exact solution
    /// given by <see cref="RinglebExactSolution"/>.
    /// </summary>
    public class RinglebTest : TestProgram<RinglebControl> {

        ///// <summary>
        ///// Alternative entry point allowing us to conduct Ringleb tests
        ///// easily.
        ///// </summary>
        ///// <param name="args">
        ///// Command line arguments
        ///// </param>
        //public static void Main(string[] args) {
        //    RinglebStiffenedTest();
        //}
        
        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        protected override BoundaryConditionMap GetBoundaryConditionMap() {
            return new RinglebBoundaryConditionMap(GridData, Control);
        }

        /// <summary>
        /// Tests the error for an ideal gas.
        /// </summary>
        [Test]
        public static void RinglebIdealGasTest() {
            Program<RinglebControl> p = null;
            Application<RinglebControl>._Main(
                new string[] { @"-c cs:CNS.Tests.Ringleb.ControlFiles.RinglebIdealGasTest()" },
                false,
                delegate() {
                    p = new RinglebTest();
                    return p;
                });

            double maxErrorEntropy = 2e-5;
            double maxErrorDensity = 5e-5;
            double maxErrorPressure = 4e-5;

            double entropyError = (double)p.QueryHandler.QueryResults["L2ErrorEntropy"];
            double densityError = (double)p.QueryHandler.QueryResults["L2ErrorDensity"];
            double pressureError = (double)p.QueryHandler.QueryResults["L2ErrorPressure"];

            Console.WriteLine(
                "Entropy L2 Error: {0} (Threshold is {1})", entropyError, maxErrorEntropy);
            Console.WriteLine(
                "Density L2 Error: {0} (Threshold is {1})", densityError, maxErrorDensity);
            Console.WriteLine(
                "Pressure L2 Error: {0} (Threshold is {1})", pressureError, maxErrorPressure);

            Assert.IsTrue(entropyError < maxErrorEntropy);
            Assert.IsTrue(densityError < maxErrorDensity);
            Assert.IsTrue(pressureError < maxErrorPressure);
        }

        /// <summary>
        /// Tests the error a stiffened gas
        /// </summary>
        [Test]
        public static void RinglebStiffenedTest() {
            Program<RinglebControl> p = null;
            Application<RinglebControl>._Main(
                new string[] { @"-c cs:CNS.Tests.Ringleb.ControlFiles.RinglebStiffenedGasTest()" },
                false,
                delegate() {
                    p = new RinglebTest();
                    return p;
                });

            double maxErrorEntropy = 3.016e-10;
            double maxErrorDensity = 6e-4;
            double maxErrorPressure = 5e-3;

            double entropyError = (double)p.QueryHandler.QueryResults["L2ErrorEntropy"];
            double densityError = (double)p.QueryHandler.QueryResults["L2ErrorDensity"];
            double pressureError = (double)p.QueryHandler.QueryResults["L2ErrorPressure"];

            Console.WriteLine(
                "Entropy L2 Error: {0} (Threshold is {1})", entropyError, maxErrorEntropy);
            Console.WriteLine(
                "Density L2 Error: {0} (Threshold is {1})", densityError, maxErrorDensity);
            Console.WriteLine(
                "Pressure L2 Error: {0} (Threshold is {1})", pressureError, maxErrorPressure);

            Assert.IsTrue(entropyError < maxErrorEntropy);
            Assert.IsTrue(densityError < maxErrorDensity);
            Assert.IsTrue(pressureError < maxErrorPressure);
        }

        /// <summary>
        /// Dummy variant of <see cref="BoundaryConditionMap"/> that always
        /// returns <see cref="ExactRinglebBoundaryState"/>.
        /// </summary>
        private class RinglebBoundaryConditionMap : BoundaryConditionMap {

            private RinglebControl ringlebControl;

            /// <summary>
            /// Constructs a new map
            /// </summary>
            /// <param name="gridData"></param>
            /// <param name="control"></param>
            public RinglebBoundaryConditionMap(IGridData gridData, RinglebControl control)
                : base(gridData, control, control.GetMaterial()) {
                this.ringlebControl = control;
            }

            /// <summary>
            /// Always returns an instance of <see cref="ExactRinglebBoundaryState"/>
            /// </summary>
            /// <param name="edgeTagName">Irrelevant</param>
            /// <returns>
            /// An appropriate instance of
            /// <see cref="ExactRinglebBoundaryState"/>.
            /// </returns>
            public override BoundaryCondition GetBoundaryCondition(string edgeTagName) {
                return new ExactRinglebBoundaryState(ringlebControl);
            }
        }
    }
}
