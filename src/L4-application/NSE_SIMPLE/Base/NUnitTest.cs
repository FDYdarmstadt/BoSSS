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

namespace NSE_SIMPLE {

    /// <summary>
    /// An all-up NUnit test for the NSE_SIMPLE application.
    /// </summary>
    [TestFixture]
    static public class NUnitTest {

        [TestFixtureSetUp]
        static public void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out dummy);
        }

        [TestFixtureTearDown]
        static public void Cleanup() {
            //Console.Out.Dispose();
            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// Tests the steady 2D-Channel flow using the 'Steady_SIMPLE' algorithm.
        /// </summary>
        [Test]
        public static void IncompressibleSteadyPoiseuilleFlowTest() {
            NSE_SIMPLEMain p = null;
            Application<SIMPLEControl>._Main(new string[] { "-c cs:NSE_SIMPLE.Incompressible.ControlExamples.PoiseuilleFlow()" }, false, delegate () {
                p = new NSE_SIMPLEMain();
                return p;
            });

            double err_u = (double)p.QueryHandler.QueryResults["SolL2err_u"];
            double err_v = (double)p.QueryHandler.QueryResults["SolL2err_v"];
            double err_p = (double)p.QueryHandler.QueryResults["SolL2err_p"];
            double thres_u = 5.1e-6;
            double thres_v = 2.8e-6;
            double thres_p = 1.6e-5;


            Console.WriteLine("Number of SIMPLE iterations: " + p.SIMPLEStatus.SIMPLEStepNo + ". Expected number of iterations is less than 133.");
            Assert.IsTrue(p.SIMPLEStatus.SIMPLEStepNo < 133);

            Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
            Assert.IsTrue(err_u < thres_u);
            Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
            Assert.IsTrue(err_v < thres_v);
            Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            Assert.IsTrue(err_p < thres_p);
        }

        /// <summary>
        /// Tests the Taylor vortex flow using the 'Unsteady_SIMPLE' algorithm.
        /// </summary>
        [Test]
        public static void IncompressibleUnsteadyTaylorVortexTest() {
            NSE_SIMPLEMain p = null;
            Application<SIMPLEControl>._Main(new string[] { "-c cs:NSE_SIMPLE.Incompressible.ControlExamples.UnsteadyTaylorVortex()" }, false, delegate () {
                p = new NSE_SIMPLEMain();
                return p;
            });

            double err_u = (double)p.QueryHandler.QueryResults["SolL2err_u"];
            double err_v = (double)p.QueryHandler.QueryResults["SolL2err_v"];
            double err_p = (double)p.QueryHandler.QueryResults["SolL2err_p"];
            double thres_vel = 8.2e-3;
            double thres_p = 1.5e-2;

            if (p.SIMPLEStatus.CntMaxNoSIMPLEsteps > 0) {
                Console.WriteLine("WARNING: For this NUnitTest all time steps are expected to converge!");
            }
            Assert.IsTrue(p.SIMPLEStatus.CntMaxNoSIMPLEsteps == 0);

            Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_vel + ")");
            Assert.IsTrue(err_u < thres_vel);
            Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_vel + ")");
            Assert.IsTrue(err_v < thres_vel);
            Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            Assert.IsTrue(err_p < thres_p);
        }
        
        /// <summary>
        /// Tests the unsteady smooth interface solver.
        /// </summary>
        [Test]
        public static void MultiphaseUnsteadyWaveTest() {
            NSE_SIMPLEMain p = null;
            Application<SIMPLEControl>._Main(new string[] { "-c cs:NSE_SIMPLE.Multiphase.ControlExamples.UnsteadyMultiphaseWave()" }, false, delegate () {
                p = new NSE_SIMPLEMain();
                return p;
            });

            double err_u = (double)p.QueryHandler.QueryResults["SolL2err_u"];
            double err_p = (double)p.QueryHandler.QueryResults["SolL2err_p"];
            double err_phi = (double)p.QueryHandler.QueryResults["SolL2err_phi"];
            double err_rho = (double)p.QueryHandler.QueryResults["SolL2err_Rho"];
            double thres_u = 6.1e-6;
            double thres_p = 4.6e-3;
            double thres_phi = 2.0e-3;
            double thres_rho = 2.0;

            if (p.SIMPLEStatus.CntMaxNoSIMPLEsteps > 0) {
                Console.WriteLine("WARNING: For this NUnitTest all time steps are expected to converge!");
            }
            Assert.IsTrue(p.SIMPLEStatus.CntMaxNoSIMPLEsteps == 0);

            Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
            Assert.IsTrue(err_u < thres_u);
            Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            Assert.IsTrue(err_p < thres_p);
            Console.WriteLine("L2 Error of solution phi: " + err_phi + " (threshold is " + thres_phi + ")");
            Assert.IsTrue(err_phi < thres_phi);
            Console.WriteLine("L2 Error of solution rho: " + err_rho + " (threshold is " + thres_rho + ")");
            Assert.IsTrue(err_rho < thres_rho);
        }

        /// <summary>
        /// Tests the steady low-Mach solver for Couette flow with temperature gradient.
        /// </summary>
        [Test]
        public static void LowMachSteadyCouetteWithTemperatureGradientTest() {
            NSE_SIMPLEMain p = null;
            Application<SIMPLEControl>._Main(new string[] { "-c cs:NSE_SIMPLE.LowMach.ControlExamples.SteadyCouetteFlowWithTemperatureGradient()" }, false, delegate () {
                p = new NSE_SIMPLEMain();
                return p;
            });

            double err_u = (double)p.QueryHandler.QueryResults["SolL2err_u"];
            double err_v = (double)p.QueryHandler.QueryResults["SolL2err_v"];
            double err_p = (double)p.QueryHandler.QueryResults["SolL2err_p"];
            double err_T = (double)p.QueryHandler.QueryResults["SolL2err_T"];
            double thres_u = 3.7e-5;
            double thres_v = 2.3e-5;
            double thres_p = 5.1e-4;
            double thres_T = 2.0e-5;

            Console.WriteLine("Number of SIMPLE iterations: " + p.SIMPLEStatus.SIMPLEStepNo + ". Expected number of iterations is less than 430.");
            Assert.IsTrue(p.SIMPLEStatus.SIMPLEStepNo < 430);

            Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
            Assert.IsTrue(err_u < thres_u);
            Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
            Assert.IsTrue(err_v < thres_v);
            Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            Assert.IsTrue(err_p < thres_p);
            Console.WriteLine("L2 Error of solution T: " + err_T + " (threshold is " + thres_T + ")");
            Assert.IsTrue(err_T < thres_T);
        }
    }
}