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

        /// <summary>
        /// MPI init
        /// </summary>
        [TestFixtureSetUp]
        public static void Init() {
            bool dummy;
            ilPSP.Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application<SipControl>.GetBoSSSInstallDir(),
                out dummy);
        }

        /// <summary>
        /// MPI teardown
        /// </summary>
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


        /// <summary>
        /// Cartesian problems with iterative solvers
        /// </summary>
        [Test]
        public static void TestIterativeSolver(
#if DEBUG            
            [Values(2)]int dgDeg,
            [Values(40)]int res,
            [Values(2)]int dim,
            [Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz)] SolverCodes solver
#else
            [Values(3,4)]int dgDeg,
            [Values(8)]int res,
            [Values(3)]int dim,
            [Values(SolverCodes.exp_gmres_levelpmg, SolverCodes.exp_Kcycle_schwarz)] SolverCodes solver

#endif
            ) {

            using(SipPoisson.SipPoissonMain p = new SipPoissonMain()) {
                var ctrl = SipHardcodedControl.TestCartesian2(res, dim, solver, dgDeg);
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


        /// <summary>
        /// operator condition number scaling
        /// </summary>
        [Test]
        public static void TestOperatorScaling(
#if DEBUG            
            [Values(1)]int dgDeg,
            [Values(2,3)]int dim
#else
            [Values(1,2,3,4)]int dgDeg,
            [Values(2,3)]int dim

#endif
            ) {

            var Controls = new List<SipControl>();

            int[] ResS = null;
            if(dim == 2) {
                switch(dgDeg) {
                    case 1: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 2: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 3: ResS = new int[] { 8, 16, 32, 64 }; break;
                    case 4: ResS = new int[] { 8, 16, 32 }; break;
                    default: throw new NotImplementedException();
                }
            } else if(dim == 3) {
                switch(dgDeg) {
                    case 1: ResS = new int[] { 4, 8, 16 }; break;
                    case 2: ResS = new int[] { 4, 8, 16 }; break;
                    case 3: ResS = new int[] { 4, 8 }; break;
                    case 4: ResS = new int[] { 4, 8 }; break;
                    default: throw new NotImplementedException();
                }
            } else {
                throw new ArgumentOutOfRangeException("un-supported spatial dimension");
            }
            
            foreach (int res in ResS) {
                var C = SipHardcodedControl.TestCartesian2(res, dim, solver_name:SolverCodes.classic_pardiso, deg: dgDeg);
                //C.TracingNamespaces = "*";
                C.savetodb = false;
                Controls.Add(C);
            }
            
            
            var Data = Solution.OpAnalysisBase.RunAndLog(Controls);

            /*
            using(var gp = new Gnuplot()) {

                var xVals = Data[OpAnalysisBase.XAxisDesignation.Grid_1Dres.ToString()];
                var yVal1 = Data["StencilCondNo-innerUncut-Var0"];
                var yVal2 = Data["TotCondNo-Var0"];

                gp.PlotLogXLogY(xVals, yVal1, "StencilCondNo", new PlotFormat("-ro"));
                gp.PlotLogXLogY(xVals, yVal2, "TotCondNo", new PlotFormat("-ro"));


                gp.Execute();

                Console.WriteLine("Press any key to continue...");
                Console.ReadKey();
            }

            /*
            Für p = 1:
            Regression of condition number slopes:
              slope of TotCondNo-Var0: 1.99104373245431
              slope of InnerCondNo-Var0: 2.10942711463848
              slope of StencilCondNo-innerUncut-Var0: 0.00980723109676584
            
            */

            
            var ExpectedSlopes = new List<ValueTuple<Solution.OpAnalysisBase.XAxisDesignation, string, double>>();

            ExpectedSlopes.Add((Solution.OpAnalysisBase.XAxisDesignation.Grid_1Dres, "TotCondNo-Var0", 2.2));
            ExpectedSlopes.Add((Solution.OpAnalysisBase.XAxisDesignation.Grid_1Dres, "StencilCondNo-innerUncut-Var0", 0.5));

            Solution.OpAnalysisBase.TestSlopes(Controls, ExpectedSlopes);
            //*/
        }


    }
}
