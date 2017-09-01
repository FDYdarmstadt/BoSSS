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

namespace CNS.Tests.BoundaryConditions {

    /// <summary>
    /// 
    /// </summary>
    [TestFixture]
    public class EulerBoundaryConditionTest : TestProgram<CNSControl> {

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSupersonicInletCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSupersonicInlet1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSubsonicOutletBoundaryCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSubsonicOutlet1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSubsonicInletBoundaryCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSubsonicInlet1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSubsonicInletAndOutletBoundaryCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSubsonicInletAndOutlet1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSubsonicPressureInletBoundaryCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSubsonicPressureInletTest1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestSubsonicPressureInletAndOutletBoundaryCondition1D() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.BoundaryConditions.ControlFiles.EulerSubsonicPressureInletAndOutletTest1D()" },
                false,
                "",
                delegate () {
                    p = new EulerBoundaryConditionTest();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }
    }
}