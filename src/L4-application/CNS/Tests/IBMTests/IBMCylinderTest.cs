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
using CNS.IBM;
using NUnit.Framework;
using System;

namespace CNS.Tests.IBMTests {

    /// <summary>
    /// Tests for the inviscid flow around a cylinder that is given by an
    /// immersed boundary. All tests restart a converged solution on a grid a
    /// 64x32 grid on the upper half of the problem domain [-40, 40] x [0, 40].
    /// </summary>
    [TestFixture]
    public class IBMCylinderTest : TestProgram<IBMControl> {

        /// <summary>
        /// Test using zeroth order DG
        /// </summary>
        [Test]
        public static void IBMCylinder0th() {
            IBMCylinderTest p = null;
            Application<IBMControl>._Main(
                new string[] {
                    @"-c ../../Tests/IBMTests/IBMCylinderControl.cs",
                    "-p 1" },
                false,
                "",
                delegate () {
                    p = new IBMCylinderTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2.7E-05),
                Tuple.Create("L2ErrorXMomentum", 4.5E-05),
                Tuple.Create("L2ErrorYMomentum", 1.1E-05),
                Tuple.Create("L2ErrorEnergy", 9.8E-06),
                Tuple.Create("L2ErrorEntropy", 0.098));
        }

        /// <summary>
        /// Test using first order DG
        /// </summary>
        [Test]
        public static void IBMCylinder1st() {
            IBMCylinderTest p = null;
            Application<IBMControl>._Main(
                new string[] {
                    @"-c ../../Tests/IBMTests/IBMCylinderControl.cs",
                    "-p 2" },
                false,
                "",
                delegate () {
                    p = new IBMCylinderTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 5.6E-08),
                Tuple.Create("L2ErrorXMomentum", 2.1E-07),
                Tuple.Create("L2ErrorYMomentum", 6.1E-07),
                Tuple.Create("L2ErrorEnergy", 1.8E-07),
                Tuple.Create("L2ErrorEntropy", 0.019));
        }

        /// <summary>
        /// Test using second order DG
        /// </summary>
        //[Test]
        public static void IBMCylinder2nd() {
            IBMCylinderTest p = null;
            Application<IBMControl>._Main(
                new string[] {
                    @"-c ../../Tests/IBMTests/IBMCylinderControl.cs",
                    "-p 3" },
                false,
                "",
                delegate () {
                    p = new IBMCylinderTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 1.2E-07),
                Tuple.Create("L2ErrorXMomentum", 3.7E-07),
                Tuple.Create("L2ErrorYMomentum", 1.7E-06),
                Tuple.Create("L2ErrorEnergy", 3.7E-07),
                Tuple.Create("L2ErrorEntropy", 0.0012));
        }

        /// <summary>
        /// Test using third order DG
        /// </summary>
        //[Test]
        public static void IBMCylinder3rd() {
            IBMCylinderTest p = null;
            Application<IBMControl>._Main(
                new string[] {
                    @"-c ../../Tests/IBMTests/IBMCylinderControl.cs",
                    "-p 4" },
                false,
                "",
                delegate () {
                    p = new IBMCylinderTest();
                    return p;
                });
            
            CheckErrorThresholds(
                p.QueryHandler.QueryResults,
                Tuple.Create("L2ErrorDensity", 2.1E-06),
                Tuple.Create("L2ErrorXMomentum", 2.1E-06),
                Tuple.Create("L2ErrorYMomentum", 5.8E-06),
                Tuple.Create("L2ErrorEnergy", 7.1E-06),
                Tuple.Create("L2ErrorEntropy", 0.00014));
        }
    }
}
