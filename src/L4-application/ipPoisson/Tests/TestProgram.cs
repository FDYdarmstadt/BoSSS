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
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Gnuplot;
using ilPSP;
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
            [Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz_CoarseMesh, SolverCodes.exp_Kcycle_schwarz_PerProcess)] SolverCodes solver
#endif
            ) {

            //SolverCodes solver = SolverCodes.exp_gmres_levelpmg;

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian2(res, dim, solver, dgDeg);
                //ctrl.TracingNamespaces = "*";

                if(ctrl.LinearSolver is OrthoMGSchwarzConfig omgsc) {
                    omgsc.CoarseKickIn = 10000;
                }

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
                    //var C = SipHardcodedControl.TestCartesian2(res, 2, solver_name: SolverCodes.classic_pardiso, deg: dgDeg);
                    var C = SipHardcodedControl.Square(res, res, deg: dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }

            }

            ConditionNumberScalingTest.Perform(Controls, new ConditionNumberScalingTest.Config() { plot = true, title = "CondnumberSlopes" });
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
                    var C = SipHardcodedControl.TestCartesian2(res, 3, solver_name: SolverCodes.direct_pardiso, deg: dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }
            }

            ConditionNumberScalingTest.Perform(Controls, new ConditionNumberScalingTest.Config() { plot = true, title = "SIP-TestOperatorScaling3D-p" + dgDeg });

        }

#if !DEBUG
        /// <summary>
        /// operator condition number scaling
        /// </summary>
        [Test()]
        public static void TestOperatorConvergence3D([Values(1,2,3)] int dgDeg) {
            //BoSSS.Application.SipPoisson.Tests.TestProgram.TestOperatorConvergence3D(2)

  
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
                    var C = SipHardcodedControl.TestCartesian2(res, 3, solver_name: SolverCodes.direct_pardiso, deg: dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }
            }

            SipSolverConvergenceTest(Controls, true, new double[] { dgDeg + 0.5 });

        }
#endif



        /// <summary>
        /// operator condition number scaling
        /// </summary>
        [Test()]
        public static void TestOperatorConvergence2D([Values(2,3)] int dgDeg) {

  
            var Controls = new List<SipControl>();
            {
                int[] ResS = null;

                switch(dgDeg) {
                    case 1: ResS = new int[] {  8, 16, 32, 64 }; break;
                    case 2: ResS = new int[] {  8, 16, 32, 64 }; break;
                    case 3: ResS = new int[] { 4, 8, 16 }; break;
                    case 4: ResS = new int[] { 4, 8, 16 }; break;
                    default: throw new NotImplementedException();
                }
                foreach(int res in ResS) {
                    var C = SipHardcodedControl.Square(res, res, dgDeg);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }
            }

            SipSolverConvergenceTest(Controls, true, new double[] { dgDeg + 0.5 });
        }


        private static void SipSolverConvergenceTest(IEnumerable<SipControl> _CS, bool useExactSolution, double[] ExpectedSlopes) {
            var CS = _CS.ToArray();
            int NoOfMeshes = CS.Length;

            double[] hS = new double[NoOfMeshes];
            MultidimensionalArray errorS = MultidimensionalArray.Create(NoOfMeshes, 1);
            string[] Names = new[] { "u" };

            SipPoissonMain[] solvers = new SipPoissonMain[NoOfMeshes];
            if (useExactSolution) {

                if (NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

                for (int k = 0; k < CS.Length; k++) {

                    var C = CS[k];
                    //using(var solver = new XNSE()) {
                    var solver = new SipPoissonMain();
                    solvers[k] = solver;
                    {
                        //Console.WriteLine("Warning! - enabled immediate plotting");
                        //C.ImmediatePlotPeriod = 1;
                        //C.SuperSampling = 3;

                        solver.Init(C);
                        solver.RunSolverMode();

                        //-------------------Evaluate Error ---------------------------------------- 

                        errorS[k,0] = solver.last_L2_ERR;
                        hS[k] = solver.GridData.iGeomCells.h_min.Min();
                    }

                }
            } else {
                if (NoOfMeshes < 3)
                    throw new ArgumentException("At least three meshes required for convergence if finest solution is assumed to be exact.");
                throw new NotImplementedException("todo");
            }


            //hS = hS.Take(hS.Length - 1).ToArray();

            


            for (int i = 0; i < errorS.GetLength(1); i++) {
                var slope = hS.LogLogRegressionSlope( errorS.GetColumn(i));

                Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{slope}\t(Expecting: {ExpectedSlopes[i]})");
            }

            for (int i = 0; i < errorS.GetLength(1); i++) {
                var slope = hS.LogLogRegressionSlope( errorS.GetColumn(i));
                Assert.IsTrue(slope >= ExpectedSlopes[i], $"Convergence Slope of {Names[i]} is degenerate.");
            }

            foreach (var s in solvers) {
                s.Dispose();
            }
        }
    }
}
