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
using System.IO;
using BoSSS.Foundation;
using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.SipPoisson.Tests {

    [TestFixture]
    class TestProgram {

        [TestFixtureSetUp]
        public static void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application<SipControl>.GetBoSSSInstallDir(),
                out dummy);
        }

        [TestFixtureTearDown]
        public static void Cleanup() {
            //Console.Out.Dispose();
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        public static void TestCartesian() {

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian1();
                p.Init(ctrl);
                p.RunSolverMode();


                //Application<SipControl>._Main(new string[] {
                //        "--control", "cs:ipPoisson.ippHardcodedControl.TestCartesian1()"
                //    },
                //    false,
                //    "",
                //    delegate() {
                //        p = new SipPoissonMain();
                //        return p;
                //    });


                double err = (double)p.QueryHandler.QueryResults["SolL2err"];
                double thres = 5.0e-9;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
         
        }

               

        [Test]
        public static void TestCurved() {

            using(SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCurved();
                p.Init(ctrl);
                p.RunSolverMode();
                //Application<SipControl>._Main(new string[] {
                //        "--control", "cs:ipPoisson.ippHardcodedControl.TestCurved()"
                //    },
                //    false,
                //    "",
                //    delegate() {
                //        p = new SipPoissonMain();
                //        return p;
                //    });

                double err = (double)p.QueryHandler.QueryResults["SolL2err"];
                double thres = 2.0e-3;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
        }


        [Test]
        public static void TestSolver() {

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian2(8, 3, SolverCodes.exp_softpcg_mg, 3);
                p.Init(ctrl);
                p.RunSolverMode();


                //Application<SipControl>._Main(new string[] {
                //        "--control", "cs:ipPoisson.ippHardcodedControl.TestCartesian1()"
                //    },
                //    false,
                //    "",
                //    delegate() {
                //        p = new SipPoissonMain();
                //        return p;
                //    });


                double err = (double)p.QueryHandler.QueryResults["SolL2err"];
                double thres = 5.0e-9;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
         
        }
    }
}
