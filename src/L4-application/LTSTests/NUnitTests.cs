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
using System.Globalization;
using System.Threading;
using BoSSS.Solution;
using NUnit.Framework;

namespace LTSTests {
    /// <summary>
    /// NUNIT test class for LTS
    /// </summary>
    [TestFixture]
    class NUnitTests : Program {

        public static void Main(string[] args) {
            Application._Main(
                args,
                true,
                () => new NUnitTests());
        }

        [TestFixtureSetUp]
        static public void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out dummy);
            Thread.CurrentThread.CurrentCulture = CultureInfo.CurrentCulture;
        }

        [TestFixtureTearDown]
        static public void Cleanup() {
            //Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        public static void LTS2order_dtCoarse(bool LTS = true, bool ALTS = false, int numOfSubgrids = 3) {
            Program test = null;

            double dt = 2E-3;
            int order = 2;
            Application._Main(
                new string[0],
                true,
                delegate () {
                    test = new Program() {
                        dt_input = dt,
                        ABorder = order,
                        LTS = LTS,
                        ALTS = ALTS,
                        numOfSubgrids = numOfSubgrids
                    };
                    return test;
                });

            double L2error = (double)test.L2error;
            double maxL2error = double.MaxValue;
            if (LTS)
                //maxL2error = 1.94E-05;   // old (Stephan)
                //maxL2error = 1.72669556339186E-05 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input 
                //maxL2error = 1.7285117532218E-05 + 1e-14; // dt=dt_input
                maxL2error = 1.7285117532218E-05 + 1e-16; // dt=dt_input
            else if (ALTS)
                //maxL2error = 1.8958107158594E-05 + 1e-14; //if (TimestepNo < 3) dt = dt_input / 3, else dt = dt_input
                //maxL2error = 1.89762864341145E-05 + 1e-14; // dt=dt_input
                maxL2error = 1.89762864341145E-05 + 1e-16; // dt=dt_input

            Console.WriteLine(
                "LTS@dt=" + dt + "and order=" + order + " L2 Error: " + L2error + " (Threshold is " + maxL2error + ")");

            Assert.IsTrue(L2error < maxL2error);
        }

        public static void LTS2order_dtFine(bool LTS = true, bool ALTS = false, int numOfSubgrids = 3) {
            Program test = null;

            double dt = 1E-3;
            int order = 2;
            Application._Main(
                new string[0],
                true,
                delegate () {
                    test = new Program() {
                        dt_input = dt,
                        ABorder = order,
                        LTS = LTS,
                        ALTS = ALTS,
                        numOfSubgrids = numOfSubgrids
                    };
                    return test;
                });

            double L2error = (double)test.L2error;
            double maxL2error = double.MaxValue;
            if (LTS)
                //maxL2error = 4.94E-05;   // old (Stephan)
                //maxL2error = 4.3207686040932E-06 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input 
                //maxL2error = 4.32303875071948E-06 + 1e-14; // dt=dt_input 
                maxL2error = 4.32303875071948E-06 + 1e-16; // dt=dt_input 
            else if (ALTS)
                //maxL2error = 4.72677420970687E-06 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input
                //maxL2error = 4.72904629811328E-06 + 1e-14; // dt=dt_input
                maxL2error = 4.72904629811328E-06 + 1e-16; // dt=dt_input

            Console.WriteLine(
                "LTS@dt=" + dt + "and order=" + order + " L2 Error: " + L2error + " (Threshold is " + maxL2error + ")");

            Assert.IsTrue(L2error < maxL2error);
        }

        public static void LTS3order_dtCoarse(bool LTS = true, bool ALTS = false, int numOfSubgrids = 3) {
            Program test = null;

            double dt = 1E-3;
            int order = 3;
            Application._Main(
                new string[0],
                true,
                delegate () {
                    test = new Program() {
                        dt_input = dt,
                        ABorder = order,
                        LTS = LTS,
                        ALTS = ALTS,
                        numOfSubgrids = numOfSubgrids
                    };
                    return test;
                });

            double L2error = (double)test.L2error;
            double maxL2error = double.MaxValue;
            if (LTS)
                //maxL2error = 1.26E-08;   // old (Stephan)
                //maxL2error = 1.11948408784892E-08 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input 
                //maxL2error = 1.11970099038565E-08 + 1e-14; // dt=dt_input 
                maxL2error = 1.11970099038565E-08 + 1e-16; // dt=dt_input 
            else if (ALTS)
                //maxL2error = 1.23576348037153E-08 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input 
                //maxL2error = 1.23595770364496E-08 + 1e-14; // dt=dt_input 
                maxL2error = 1.23595770364496E-08 + 1e-16; // dt=dt_input 

            Console.WriteLine(
                "LTS@dt=" + dt + "and order=" + order + " L2 Error: " + L2error + " (Threshold is " + maxL2error + ")");

            Assert.IsTrue(L2error < maxL2error);
        }

        public static void LTS3order_dtFine(bool LTS = true, bool ALTS = false, int numOfSubgrids = 3) {
            Program test = null;

            double dt = 5E-4;
            int order = 3;
            Application._Main(
                new string[0],
                true,
                delegate () {
                    test = new Program() {
                        dt_input = dt,
                        ABorder = order,
                        LTS = LTS,
                        ALTS = ALTS,
                        numOfSubgrids = numOfSubgrids
                    };
                    return test;
                });

            double L2error = (double)test.L2error;
            double maxL2error = double.MaxValue;
            if (LTS)
                //maxL2error = 1.57E-09;   // old (Stephan)
                //maxL2error = 1.40012629740322E-09 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input
                //maxL2error = 1.40026137470915E-09 + 1e-14; // dt=dt_input 
                maxL2error = 1.40026137470915E-09 + 1e-16; // dt=dt_input 
            else if (ALTS)
                //maxL2error = 1.5438019656137E-09 + 1e-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input
                //maxL2error = 1.54550662596205E-09 + 1e-14; //dt=dt_input
                maxL2error = 1.54550662596205E-09 + 1e-16; //dt=dt_input

            Console.WriteLine(
                "LTS@dt=" + dt + "and order=" + order + " L2 Error: " + L2error + " (Threshold is " + maxL2error + ")");

            Assert.IsTrue(L2error < maxL2error);
        }

        public static void CheckTimestepper_order2_subgrid1(bool LTS = false, bool ALTS = false, int numOfSubgrids = 1, double dt = 2E-3 / 4) {
            Program test = null;

            int order = 2;
            Application._Main(
                new string[0],
                true,
                delegate () {
                    test = new Program() {
                        dt_input = dt,
                        ABorder = order,
                        LTS = LTS,
                        ALTS = ALTS,
                        numOfSubgrids = numOfSubgrids
                    };
                    return test;
                });

            double L2error = (double)test.L2error;
            //double maxL2error = 2.04207522617079E-06 + 1E-14; // if (TimestepNo < 3) dt=dt_input/3, else dt=dt_input
            //double maxL2error = 2.04241551595267E-06 + 1E-14; // dt=dt_input
            double maxL2error = 2.04241551595267E-06 + 1e-16; // dt=dt_input

            Console.WriteLine(
                "Check AB, LTS, A-LTS@dt=" + dt + "and order=" + order + " L2 Error: " + L2error + " (Threshold is " + maxL2error + ")");

            Assert.IsTrue(L2error < maxL2error);
        }

        // Call LTS tests
        [Test]
        public static void LTS2order_dtCoarse() {
            LTS2order_dtCoarse(LTS: true, ALTS: false, numOfSubgrids: 3);
        }

        [Test]
        public static void LTS2order_dtFine() {
            LTS2order_dtFine(LTS: true, ALTS: false, numOfSubgrids: 3);
        }

        [Test]
        public static void LTS3order_dtCoarse() {
            LTS3order_dtCoarse(LTS: true, ALTS: false, numOfSubgrids: 3);
        }

        [Test]
        public static void LTS3order_dtFine() {
            LTS2order_dtFine(LTS: true, ALTS: false, numOfSubgrids: 3);
        }

        // Call adaptive LTS tests
        // A-LTS gives same results as LTS for temporaly fixed clustering
        // when no flux correction for achieving conservativity is done
        // Here, the A-LTS test run without flux correction
        [Test]
        public static void ALTS2order_dtCoarse() {
            LTS2order_dtCoarse(LTS: false, ALTS: true, numOfSubgrids: 3);
        }

        [Test]
        public static void ALTS2order_dtFine() {
            LTS2order_dtFine(LTS: false, ALTS: true, numOfSubgrids: 3);
        }

        [Test]
        public static void ALTS3order_dtCoarse() {
            LTS3order_dtCoarse(LTS: false, ALTS: true, numOfSubgrids: 3);
        }

        [Test]
        public static void ALTS3order_dtFine() {
            LTS2order_dtFine(LTS: false, ALTS: true, numOfSubgrids: 3);
        }

        // Call tests for checking Adams-Bashforth, LTS and A-LTS for one subgrid
        // All time steppers should give the same results as pure Adams-Bashforth
        [Test]
        public static void Comparison_AB_to_AB_order2_subgrid1() {
            CheckTimestepper_order2_subgrid1(LTS: false, ALTS: false, numOfSubgrids: 1, dt: 2E-3 / 4);
        }

        [Test]
        public static void Comparison_LTS_to_AB_order2_subgrid1() {
            CheckTimestepper_order2_subgrid1(LTS: true, ALTS: false, numOfSubgrids: 1, dt: 2E-3 / 4);
        }

        [Test]
        public static void Comparison_ALTS_to_AB_order2_subgrid1() {
            CheckTimestepper_order2_subgrid1(LTS: false, ALTS: true, numOfSubgrids: 1, dt: 2E-3 / 4);
        }
    }
}
