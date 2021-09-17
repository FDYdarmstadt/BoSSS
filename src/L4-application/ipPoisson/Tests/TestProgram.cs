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
using System.IO;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Gnuplot;
using MPI.Wrappers;
using NUnit.Framework;
using SolverCodes = BoSSS.Solution.Control.LinearSolverCode;

namespace BoSSS.Application.SipPoisson.Tests {
    



    /// <summary>
    /// NUnit tests
    /// </summary>
    [TestFixture]
    static class TestProgram {



        [Test]
        public static void TestCartesian() {
            //BoSSS.Application.SipPoisson.Tests.TestProgram.TestCartesian
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
                double thres = 2.5e-3;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
        }


        /// <summary>
        /// Cartesian problems with iterative solvers
        /// </summary>
        [Test]
        public static void TestIterativeSolver(
#if DEBUG            
            [Values(2)]int dgDeg,
            [Values(40)]int res,
            [Values(2)]int dim,
            [Values(SolverCodes.exp_gmres_levelpmg)] SolverCodes solver
#else
            [Values(3,4)]int dgDeg,
            [Values(8)]int res,
            [Values(3)]int dim,
            [Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz)] SolverCodes solver
#endif
            ) {

            //SolverCodes solver = SolverCodes.exp_gmres_levelpmg;

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian2(res, dim, solver, dgDeg);
                //ctrl.TracingNamespaces = "*";
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
                double h = ((Foundation.Grid.Classic.GridData)(p.GridData)).Cells.h_maxGlobal;
                double thres = 0.01 * Math.Pow(h, dgDeg);

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }

        }

        /*
        /// <summary>
        /// Cartesian problems with iterative solvers
        /// </summary>
        [Test]
        public static void TestIterativeSolver_Kcycle_schwarz(
#if DEBUG            
            [Values(2)]int dgDeg,
            [Values(40)]int res,
            [Values(2)]int dim
            //[Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz)] SolverCodes solver
#else
            [Values(3,4)]int dgDeg,
            [Values(8)]int res,
            [Values(3)]int dim
            //[Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz)] SolverCodes solver
#endif
            ) {

            SolverCodes solver = SolverCodes.exp_Kcycle_schwarz;

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian2(res, dim, solver, dgDeg);
                //ctrl.TracingNamespaces = "*";
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
                double h = ((Foundation.Grid.Classic.GridData)(p.GridData)).Cells.h_maxGlobal;
                double thres = 0.01 * Math.Pow(h, dgDeg);

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }

        }
        */

#if !DEBUG
        /// <summary>
        /// operator condition number scaling, 2D, for p=1 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling2D_p1() {
            TestOperatorScaling2D(1);
        }
        /// <summary>
        /// operator condition number scaling, 3D, for p=2 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling2D_p2() {
            TestOperatorScaling2D(2);
        }
        /// <summary>
        /// operator condition number scaling, 3D, for p=3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling2D_p3() {
            TestOperatorScaling2D(3);
        }
#endif

        /// <summary>
        /// operator condition number scaling, 2D
        /// </summary>
        public static void TestOperatorScaling2D(int dgDeg) {

            

            var Controls = new List<SipControl>();
            {
                int[] ResS = null;

                switch(dgDeg) {
                    case 1: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 2: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 3: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 4: ResS = new int[] { 8, 16, 32 }; break;
                    default: throw new NotImplementedException();
                }

                foreach(int res in ResS) {
                    var C = SipHardcodedControl.TestCartesian2(res, 2, solver_name: SolverCodes.classic_pardiso, deg: dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }

            }

            ConditionNumberScalingTest.Perform(Controls, false);
        }

#if !DEBUG
        /// <summary>
        /// operator condition number scaling, 3D, for p=1 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling3D_p1() {
            TestOperatorScaling3D(1);
        }
        /// <summary>
        /// operator condition number scaling, 3D, for p=2 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling3D_p2() {
            TestOperatorScaling3D(1);
        }
        /// <summary>
        /// operator condition number scaling, 3D, for p=3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void TestOperatorScaling3D_p3() {
            TestOperatorScaling3D(1);
        }
#endif

        /// <summary>
        /// operator condition number scaling
        /// </summary>
        public static void TestOperatorScaling3D(int dgDeg) {

  
            var Controls = new List<SipControl>();
            {
                int[] ResS = null;

                switch(dgDeg) {
                    case 1: ResS = new int[] { 4, 8, 16 }; break;
                    case 2: ResS = new int[] { 4, 8, 16 }; break;
                    case 3: ResS = new int[] { 4, 8 }; break;
                    case 4: ResS = new int[] { 4, 8 }; break;
                    default: throw new NotImplementedException();
                }
                foreach(int res in ResS) {
                    var C = SipHardcodedControl.TestCartesian2(res, 3, solver_name: SolverCodes.classic_pardiso, deg: dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }
            }

            ConditionNumberScalingTest.Perform(Controls);

        }
    }
}
