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

using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// An all-up NUnit test for the LowMachCombustion application.
    /// </summary>
    [TestFixture]
    public static partial class NUnitTest {
        //[OneTimeSetUp]
        //public static void SetUp() {
        //    BoSSS.Solution.Application.InitMPI();
        //}

        ///// <summary>
        ///// MPI shutdown.
        ///// </summary>
        //[OneTimeTearDown]
        //public static void TestFixtureTearDown() {
        //    csMPI.Raw.mpiFinalize();
        //}

        /// <summary>
        /// Tests the steady 2D-Channel flow using the 'Steady_SIMPLE' algorithm.***
        /// </summary>
        [Test]
        public static void IncompressibleSteadyPoiseuilleFlowTest() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.ChannelFlowTest_NUnit(false);
                p.Init(c);
                p.RunSolverMode();
                //p.OperatorAnalysis();
                // p.CheckJacobian();
                double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX];
                double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY];
                double err_p = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Pressure];
                double thres_u = 5.1e-6;
                double thres_v = 2.8e-6;
                double thres_p = 1.6e-5;

                Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");

                Assert.Less(err_u, thres_u, "L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Assert.Less(err_v, thres_v, "L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            }
        }

        /// <summary>
        /// Tests the steady 2D-Channel flow using the 'Steady_SIMPLE' algorithm.***
        /// </summary>
        [Test]
        public static void IncompressibleSteadyPoiseuilleFlowTest_WithImmersedBoundary() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.ChannelFlowTest_NUnit(true);     
                p.Init(c);
                p.RunSolverMode();
                //p.OperatorAnalysis();
                // p.CheckJacobian();
                double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX];
                double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY];
                double err_p = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Pressure];
                double thres_u = 5.1e-6;
                double thres_v = 2.8e-6;
                double thres_p = 1.6e-5;

                Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");

                Assert.Less(err_u, thres_u, "L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Assert.Less(err_v, thres_v, "L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            }
        }



   

        /// <summary>
        /// 
        /// </summary>
        //[Test]
        public static void ImmersedBoundaryTest_MixtureFraction_DropletCombustion() {

            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            //using (var p = new XNSEC_MixtureFraction()) {
            //    var c =  BoSSS.Application.XNSEC.FullNSEControlExamples.NuNit_Droplet_ImmersedBoundary_MixtureFraction( db.Path);                
            //    c.ImmediatePlotPeriod = 1;
            //    c.SuperSampling = 2;
            //    p.Init(c);
            //    p.RunSolverMode();
            //}

            Console.WriteLine("Flame sheet calculation done.");
            //using (var p = new XNSEC()) {
            //    var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_DropletFlame(2, 10,  db.Path);
            //    c.ImmediatePlotPeriod = 1;
            //    c.SuperSampling = 2;
            //    p.Init(c);
            //    p.RunSolverMode();

            //    var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
            //    var temp = temperatureXdg.ProjectToSinglePhaseField(4);
            //    double minT; double maxT;
            //    temp.GetExtremalValues(out minT, out maxT);
            //    Console.WriteLine("Maximum reached temperature is {0}K", maxT);
            //}
            //Console.WriteLine("Full calculation done.");

        }


        /// <summary>
        /// Tests the steady 2D-Channel flow using the 'Steady_SIMPLE' algorithm.***
        /// </summary>
        //[Test]
        public static void TwoPhaseIncompressibleSteadyPoiseuilleFlowTest() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.TwoPhaseChannelFlowTest_NUnit();
                p.Init(c);
                p.RunSolverMode();
                //// p.CheckJacobian();
                //double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX ];
                //double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY ];
                //double err_p = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Pressure ];
                //double thres_u = 5.1e-6;
                //double thres_v = 2.8e-6;
                //double thres_p = 1.6e-5;

                //Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                //Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                //Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");

                //Assert.Less(err_u, thres_u, "L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                //Assert.Less(err_v, thres_v, "L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                //Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            }
        }

        /// <summary>
        /// Tests the steady combustion solver with a manufactured solution.
        /// Current manufactured solutions used is T = cos(x*y), Y0 = 0.3 cos(x*y), Y1 = 0.6 cos(x*y), Y2 = 0.1 cos(x*y), u = -cos(x), v = -cos(y), p = sin(x*y).
        /// </summary>
        [Test]
        public static void ManufacturedSolutionLowMachCombustionTest() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.NUnitTestManuSol_3();
                p.Init(c);
                p.RunSolverMode();
                //p.OperatorAnalysis();
                //int[] varGroup_convDiff = new int[] { 0, 1 }; // u,v
                //int[] varGroup_Stokes = new int[] { 0, 1, 2 }; // u,v,p
                //int[] varGroup_onlyp = new int[] { 1 }; // u
                //int[] varGroup_Constitutive = new int[] { 3, 4 }; // T, Y
                //int[] varGroup_all = new int[] { 0, 1, 2, 3, 4 }; //all
                //p.Timestepping.TimesteppingBase.OperatorAnalysis(new[] { varGroup_convDiff, varGroup_Stokes, varGroup_onlyp, varGroup_Constitutive, varGroup_all }, true);

                double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX];
                double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY];
                double err_p = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Pressure];
                double err_T = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Temperature];
                double err_Y0 = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.MassFraction0];
                double err_Y1 = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.MassFraction1];
                double err_Y2 = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.MassFraction2];

                double factor = 1.2;
                double thres_u = 0.002 * factor;
                double thres_v = 0.002 * factor;
                double thres_p = 0.08 * factor;
                double thres_T = 0.002 * factor;
                double thres_Y0 = 0.0004 * factor;
                double thres_Y1 = 0.0004 * factor;
                double thres_Y2 = 0.0005 * factor;

                Console.WriteLine("L2 Error of solution v: " + err_u + " (threshold is " + thres_u + ")");
                Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
                Console.WriteLine("L2 Error of solution T: " + err_T + " (threshold is " + thres_T + ")");
                Console.WriteLine("L2 Error of solution Y0: " + err_Y0 + " (threshold is " + thres_Y0 + ")");
                Console.WriteLine("L2 Error of solution Y1: " + err_Y1 + " (threshold is " + thres_Y1 + ")");
                Console.WriteLine("L2 Error of solution Y2: " + err_Y2 + " (threshold is " + thres_Y2 + ")");

                Assert.Less(err_u, thres_u, "L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Assert.Less(err_v, thres_v, "L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
                Assert.Less(err_T, thres_T, "L2 Error of solution T: " + err_T + " (threshold is " + thres_T + ")");
                Assert.Less(err_Y0, thres_Y0, "L2 Error of solution Y0: " + err_Y0 + " (threshold is " + thres_Y0 + ")");
                Assert.Less(err_Y1, thres_Y1, "L2 Error of solution Y1: " + err_Y1 + " (threshold is " + thres_Y1 + ")");
                Assert.Less(err_Y2, thres_Y2, "L2 Error of solution Y2: " + err_Y2 + " (threshold is " + thres_Y2 + ")");

                //Console.WriteLine("Number fix point iterations: " + p.NumIterations + ". Expected number of iterations is less than 17.");
                //Assert.Less(p.NumIterations, 20 * factor);
            }
        }

        /// <summary>
        /// Tests the Taylor vortex flow
        /// </summary>
        [Test]
        public static void IncompressibleUnsteadyTaylorVortexTest() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.NUnitUnsteadyTaylorVortex();
                //c.ImmediatePlotPeriod = 1;
                p.Init(c);
                p.RunSolverMode();
                double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX];
                double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY];
                double err_p = (double)p.QueryHandler.QueryResults["SolL2err_p"];
                double thres_vel = 0.02;
                double thres_p = 0.55;

                Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_vel + ")");
                Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_vel + ")");
                Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");

                Assert.Less(err_u, thres_vel, "L2 Error of solution u: " + err_u + " (threshold is " + thres_vel + ")");
                Assert.Less(err_v, thres_vel, "L2 Error of solution v: " + err_v + " (threshold is " + thres_vel + ")");
                Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
            }
        }

        /// <summary>
        /// Tests the steady low-Mach solver for Couette flow with temperature gradient.
        /// Dg grad = 3 and nCells = 16x16 ***
        /// </summary>
        [Test]
        public static void LowMachSteadyCouetteWithTemperatureGradientTest() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.NUnitSteadyCouetteFlowWithTemperatureGradient();
                p.Init(c);
                p.RunSolverMode();
                //p.OperatorAna<lysis();
                double err_u = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityX];
                double err_v = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.VelocityY];
                double err_p = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Pressure];
                double err_T = (double)p.QueryHandler.QueryResults["Err_" + VariableNames.Temperature];
                double thres_u = 2e-5;
                double thres_v = 3e-5;
                double thres_p = 0.09;
                double thres_T = 6e-5;

                Console.WriteLine("L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Console.WriteLine("L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Console.WriteLine("L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
                Console.WriteLine("L2 Error of solution T: " + err_T + " (threshold is " + thres_T + ")");

                Assert.Less(err_u, thres_u, "L2 Error of solution u: " + err_u + " (threshold is " + thres_u + ")");
                Assert.Less(err_v, thres_v, "L2 Error of solution v: " + err_v + " (threshold is " + thres_v + ")");
                Assert.Less(err_p, thres_p, "L2 Error of solution p: " + err_p + " (threshold is " + thres_p + ")");
                Assert.Less(err_T, thres_T, "L2 Error of solution T: " + err_T + " (threshold is " + thres_T + ")");
            }
        }

        /// <summary>
        /// ***
        /// </summary>
        [Test]
        public static void CavityNaturalConvection() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.NaturalConvectionSquareCavityTest();
                p.Init(c);
                p.RunSolverMode();
                //p.OperatorAnalysis();
                var p0_solutions = new[] {
                    (1e2, 0.9573),
                    (1e3, 0.9381),
                    (1e4, 0.9146),
                    (1e5, 0.9220),
                    (1e6, 0.9245),
                    (1e7, 0.9226)
                };

                double Ra = c.Rayleigh;
                double p0Reference = -1;
                foreach (var sol in p0_solutions) {
                    if (Math.Abs(sol.Item1 - Ra) < 1e-3) {
                        p0Reference = sol.Item2;
                    }
                }
                if (p0Reference == -1)
                    throw new NotImplementedException();
                var thermoPressure = p.Parameters.Where(f => f.Identification == VariableNames.ThermodynamicPressure).FirstOrDefault();
                double ThermPressureCalculated = thermoPressure.GetMeanValueTotal(null);

                Console.WriteLine("The calculated thermodynamic pressure is  " + ThermPressureCalculated + " (and the reference value is " + p0Reference + ")");
                Console.WriteLine("aaaaaaaaaaaa"+Math.Abs(ThermPressureCalculated - p0Reference));
                if (Math.Abs(ThermPressureCalculated - p0Reference) > 1e-2) { 
                    throw new Exception("Error on calculation of the thermodynamic pressure. End value is not the correct one");

                }
                Console.WriteLine("The test passed! ");
            }
        }

        /// <summary>
        /// ***
        /// </summary>
        [Test]
        public static void CavityNaturalConvectionTest_Homotopy() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.NaturalConvectionSquareCavityTest_Homotopy();
                p.Init(c);
                p.RunSolverMode();
                var p0_solutions = new[] {
                    (1e2, 0.9573),
                    (1e3, 0.9381),
                    (1e4, 0.9146),
                    (1e5, 0.9220),
                    (1e6, 0.9245),
                    (1e7, 0.9226)
                };

                double Ra = c.Rayleigh;
                double p0Reference = -1;
                foreach (var sol in p0_solutions) {
                    if (Math.Abs(sol.Item1 - Ra) < 1e-3) {
                        p0Reference = sol.Item2;
                    }
                }
                if (p0Reference == -1)
                    throw new NotImplementedException();
                var thermoPressure = p.Parameters.Where(f => f.Identification == VariableNames.ThermodynamicPressure).FirstOrDefault();
                double ThermPressureCalculated = thermoPressure.GetMeanValueTotal(null);

                Console.WriteLine("The calculated thermodynamic pressure is  " + ThermPressureCalculated + " (and the reference value is " + p0Reference + ")");

                if (Math.Abs(ThermPressureCalculated - p0Reference) > 1e-2)
                    throw new Exception("Error on calculation of the thermodynamic pressure. End value is not the correct one");

                Console.WriteLine("The test passed! ");
            }
        }

        /// <summary>
        /// Testcase where Adaptive Mesh Refinement is used after each newton iteration
        /// </summary>
        //[Test]
        public static void CavityNaturalConvection_AMR_eachNewtonIteration() {
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.CavityNaturalConvection_AMR_eachNewtonIteration();
                p.Init(c);
                p.RunSolverMode();
                var p0_solutions = new[] {
                    (1e2, 0.9573),
                    (1e3, 0.9381),
                    (1e4, 0.9146),
                    (1e5, 0.9220),
                    (1e6, 0.9245),
                    (1e7, 0.9226)
                };

                double Ra = c.Rayleigh;
                double p0Reference = -1;
                foreach (var sol in p0_solutions) { // search for the corresponding reference p0
                    if (Math.Abs(sol.Item1 - Ra) < 1e-3) {
                        p0Reference = sol.Item2;
                    }
                }
                if (p0Reference == -1)
                    throw new NotImplementedException();
                var thermoPressure = p.Parameters.Where(f => f.Identification == VariableNames.ThermodynamicPressure).FirstOrDefault();
                double ThermPressureCalculated = thermoPressure.GetMeanValueTotal(null);

                Console.WriteLine("The calculated thermodynamic pressure is  " + ThermPressureCalculated + " (and the reference value is " + p0Reference + ")");

                if (Math.Abs(ThermPressureCalculated - p0Reference) > 1e-2)
                    throw new Exception("Error on calculation of the thermodynamic pressure. End value is not the correct one");

                Console.WriteLine("The test passed! ");
            }
        }

        //#if !DEBUG
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
        //[Test]
        public static void TestOperatorScaling2D_p3() {
            TestOperatorScaling2D(3);
        }

        /// <summary>
        /// operator condition number scaling, 2D
        /// </summary>
        public static void TestOperatorScaling2D(int dgDeg) {
            var Controls = new List<XNSEC_Control>();
            {
                int[] ResS = null;

                switch (dgDeg) {
                    //case 1: ResS = new int[] { 4, 5, 6, 7 }; break;
                    case 1: ResS = new int[] { 3, 4, 5, 6 }; break;
                    case 2: ResS = new int[] { 3, 4, 5 }; break; // more than 6 leads to out of memory in matlab
                    case 3: ResS = new int[] { 4, 5, 6, 7 }; break;
                    case 4: ResS = new int[] { 4, 5, 6 }; break;
                    default: throw new NotImplementedException();
                }

                foreach (int res in ResS) {
                    var C = BoSSS.Application.XNSEC.FullNSEControlExamples.ControlManuSolLowMachCombustion(dgDeg, res);
                    //C.TracingNamespaces = "*";
                    C.savetodb = false;
                    Controls.Add(C);
                }
            }

            ConditionNumberScalingTest.Perform(Controls, new ConditionNumberScalingTest.Config() { plot = true, title = "XNSEC-TestOperatorScaling2D-p" + dgDeg });
        }

        private static void XNSECSolverTest(IXNSECTest Tst, XNSEC_Control C) {
            using (var solver = new XNSEC()) {
                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();
                //int[] varGroup_convDiff = new int[] { 0, 1 }; // u,v
                //int[] varGroup_Stokes = new int[] { 0, 1, 2 }; // u,v,p
                //int[] varGroup_Constitutive = new int[] { 3, 4 }; // T, Y
                //int[] varGroup_all = new int[] { 0, 1, 2, 3, 4 }; //all
                //solver.Timestepping.TimesteppingBase.OperatorAnalysis(new[] { varGroup_convDiff, varGroup_Stokes,  varGroup_Constitutive, varGroup_all }, true);

                //-------------------Evaluate Error ----------------------------------------

                var evaluator = new XNSEErrorEvaluator<XNSEC_Control>(solver);
                double[] XNSE_Errors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);
                var combustionEvaluator = new CombustionErrorEvaluator<XNSEC_Control>(solver);
                double[] CombustionErrors = combustionEvaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                List<double> errors = new List<double>();
                errors.AddRange(XNSE_Errors);
                errors.AddRange(CombustionErrors);
                double[] AllErrors = errors.ToArray();

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (AllErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++) {
                    bool ok = AllErrors[i] <= ErrThresh[i];
                    Console.Write("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], AllErrors[i]);

                    if (ok)
                        Console.WriteLine("   (ok)");
                    else
                        Console.WriteLine("   Above Threshold (" + ErrThresh[i] + ")");
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    bool ok = ResNorms[i] <= ResThresh[i];
                    Console.Write("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);

                    if (ok)
                        Console.WriteLine("   (ok)");
                    else
                        Console.WriteLine("   Above Threshold (" + ResThresh[i] + ")");
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(AllErrors[i], ErrThresh[i], $"Error {solver.CurrentState.Fields[i].Identification} above threshold.");

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i], $"Residual {solver.CurrentResidual.Fields[i].Identification} above threshold.");
            }
        }

        private static XNSEC_Control TstObj2CtrlObj(IXNSECTest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            bool constantDensity,
            int GridResolution = 1, LinearSolverCode solvercode = LinearSolverCode.direct_pardiso) {
            XNSEC_Control C = new XNSEC_Control();
            int D = tst.SpatialDimension;
            int NoChemSpc = tst.NumberOfChemicalComponents;
            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSEC/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            C.NumberOfChemicalSpecies = tst.NumberOfChemicalComponents;
            C.SetDGdegree(FlowSolverDegree);

            // grid
            // ====

            C.GridFunc = () => tst.CreateGrid(GridResolution);

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================
            C.NumberOfChemicalSpecies = tst.NumberOfChemicalComponents;
            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;

            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;

            C.prescribedMassflux_Evaluator = (tst is IPrescribedMass ? (tst as IPrescribedMass).GetPrescribedMassflux_Evaluator() : null);
            C.ThermalParameters.hVap = (tst is IPrescribedMass ? 1.0 : 0.0);
            // initial values and exact solution
            // =================================

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionMassFractions = new Dictionary<string, Func<double[], double, double>[]>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));
                C.ExactSolutionMassFractions.Add(spc, NoChemSpc.ForLoop(q => tst.GetMassFractions(spc, q)));

                C.ExactSolutionTemperature.Add(spc, tst.GetTemperature(spc));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetTemperature(spc).Convert_Xt2X(0.0));

                C.InitialValues_Evaluators.Add(VariableNames.MassFraction0 + "#" + spc, X => 1.0);
            }
            if (tst.TestImmersedBoundary) {
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators_TimeDep.Add(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.Velocity_d(d)), tst.GetPhi2U(d));
                }
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCG, tst.GetPhi());

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            if (D == 3 && SurfTensionMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                Console.WriteLine($"Reminder: {SurfTensionMode} changed to LaplaceBeltrami_ContactLine for 3D test.");
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            } else {
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            }
            C.CutCellQuadratureType = CutCellQuadratureType;

            // immersed boundary
            // =================

            C.UseImmersedBoundary = tst.TestImmersedBoundary;
            if (C.UseImmersedBoundary) {
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), tst.GetPhi2());
            }

            // timestepping and solver
            // =======================

            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;
                C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            //C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.NonLinearSolver.MaxSolverIterations = 3;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver = solvercode.GetConfig();
            C.GravityDirection = tst.GravityDirection;
            C.ChemicalReactionActive = tst.ChemicalReactionTermsActive;
            C.EnableMassFractions = tst.EnableMassFractions;
            C.EnableTemperature = tst.EnableTemperature;
            C.rhoOne = constantDensity;
            // return
            // ======
            Assert.AreEqual(C.UseImmersedBoundary, tst.TestImmersedBoundary);
            return C;
        }


        private static (string Name, double slope, double Intercept)[] XNSECSolverConvergenceTest(IXNSETest Tst, XNSEC_Control[] CS, bool useExactSolution, (string Name, double Slope, double intercept, double interceptTol)[] RegResults) {
            int D = Tst.SpatialDimension;
            int NoOfMeshes = CS.Length;
            if (RegResults.Length != D + 1)
                throw new ArgumentException("Expecting slopes for velocity and pressure.");

            var Ret = new List<(string Name, double slope, double Intercept)>();

            if (useExactSolution) {
                if (NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

                MultidimensionalArray errorS = null;
                string[] Names = null;

                double[] hS = new double[NoOfMeshes];
                XNSEC[] solvers = new XNSEC[NoOfMeshes];

                for (int k = 0; k < CS.Length; k++) {

                    var C = CS[k];
                    //using(var solver = new XNSEC()) {
                    var solver = new XNSEC();
                    solvers[k] = solver;
                    {
                        //Console.WriteLine("Warning! - enabled immediate plotting");
                        //C.ImmediatePlotPeriod = 1;
                        //C.SuperSampling = 3;

                        solver.Init(C);
                        solver.RunSolverMode();

                        //-------------------Evaluate Error ---------------------------------------- 
                        var evaluator = new XNSEErrorEvaluator<XNSEC_Control>(solver);
                        double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);
                        double[] ErrThresh = Tst.AcceptableL2Error;


                        if (k == 0) {
                            errorS = MultidimensionalArray.Create(NoOfMeshes, LastErrors.Length);
                            Names = new string[LastErrors.Length];
                            if (RegResults.Length != Names.Length)
                                throw new ArgumentOutOfRangeException();
                        } else {
                            if (LastErrors.Length != Names.Length)
                                throw new ApplicationException();
                        }

                        if (LastErrors.Length != ErrThresh.Length)
                            throw new ApplicationException();
                        for (int i = 0; i < ErrThresh.Length; i++) {
                            Console.WriteLine($"L2 error, '{solver.Operator.DomainVar[i]}': \t{LastErrors[i]}");
                            Names[i] = solver.Operator.DomainVar[i];
                        }

                        errorS.SetRow(k, LastErrors);
                        hS[k] = evaluator.GetGrid_h();
                    }

                }

                for (int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);
                    Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{RegModel.Slope}\tIntercept: \t{RegModel.Intercept}\t(Expecting: {RegRef.Slope}, {RegRef.intercept}+/-{RegRef.interceptTol})");
                    Ret.Add((Names[i], RegModel.Slope, RegModel.Intercept));
                }

                for (int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);


                    Assert.GreaterOrEqual(RegModel.Slope, RegRef.Slope, $"Convergence Slope of {Names[i]} is degenerate.");
                    Assert.GreaterOrEqual(RegModel.Intercept, RegRef.intercept - RegRef.interceptTol, $"Convergence Intercept of {Names[i]} is degenerate.");
                    Assert.LessOrEqual(RegModel.Intercept, RegRef.intercept + RegRef.interceptTol, $"Convergence Intercept of {Names[i]} overshoot.");
                }

                foreach (var s in solvers) {
                    s.Dispose();
                }
            } else {
                if (NoOfMeshes < 3)
                    throw new ArgumentException("At least three meshes required for convergence if finest solution is assumed to be exact.");



                BoSSS.Solution.Statistic.ConvergenceTest.SolverConvergenceTest_Experimental(
                    CS,
                    "Experimental Convergence",
                    (D + 1).ForLoop(iVar => (iVar < D ? VariableNames.Velocity_d(iVar) : VariableNames.Pressure,
                                             iVar < D ? BoSSS.Solution.Statistic.NormType.L2_approximate : BoSSS.Solution.Statistic.NormType.L2noMean_approximate,
                                             RegResults[iVar].Slope,
                                             RegResults[iVar].intercept, RegResults[iVar].interceptTol)));

            }



            return Ret.ToArray();


        }

        private static XNSEC_Control TstObj2CtrlObj(IXNSECTest_Heat tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
        XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
        SurfaceStressTensor_IsotropicMode SurfTensionMode,
        bool constantDensity,
        int GridResolution = 1, LinearSolverCode solvercode = LinearSolverCode.direct_pardiso) {
            XNSEC_Control C = new XNSEC_Control();
            int D = tst.SpatialDimension;
            int NoChemSpc = tst.NumberOfChemicalComponents;
            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSEC/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            C.NumberOfChemicalSpecies = tst.NumberOfChemicalComponents;
            C.SetDGdegree(FlowSolverDegree);

            // grid
            // ====

            C.GridFunc = () => tst.CreateGrid(GridResolution);

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================
            C.NumberOfChemicalSpecies = tst.NumberOfChemicalComponents;
            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;

            //C.PhysicalParametersCombustion.rho_A = tst.rho_A;
            //C.PhysicalParametersCombustion.rho_B = tst.rho_B;
            //C.PhysicalParametersCombustion.mu_A = tst.mu_A;
            //C.PhysicalParametersCombustion.mu_B = tst.mu_B;
            //C.PhysicalParametersCombustion.Sigma = tst.Sigma;
            //C.PhysicalParametersCombustion.IncludeConvection = tst.IncludeConvection;
            //C.PhysicalParametersCombustion.rhoD_A = tst.rhoD_A;
            //C.PhysicalParametersCombustion.rhoD_B = tst.rhoD_B;

            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;

            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;
            C.ThermalParameters.c_A = tst.c_A;
            C.ThermalParameters.c_B = tst.c_B;
            C.ThermalParameters.k_A = tst.k_A;
            C.ThermalParameters.k_B = tst.k_B;
            C.ThermalParameters.T_sat = tst.T_sat;
            C.ThermalParameters.IncludeConvection = tst.IncludeConvection;

            C.prescribedMassflux_Evaluator = (tst is IPrescribedMass ? (tst as IPrescribedMass).GetPrescribedMassflux_Evaluator() : null);
            C.ThermalParameters.hVap = (tst is IPrescribedMass ? 1.0 : tst.h_vap);

            // initial values and exact solution
            // =================================
            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionMassFractions = new Dictionary<string, Func<double[], double, double>[]>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));
                C.ExactSolutionMassFractions.Add(spc, NoChemSpc.ForLoop(q => tst.GetMassFractions(spc, q)));
                C.ExactSolutionTemperature.Add(spc, tst.GetTemperature(spc));

                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }
                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetTemperature(spc).Convert_Xt2X(0.0));
                for (int i = 0; i < tst.NumberOfChemicalComponents; i++) {
                    C.InitialValues_Evaluators.Add(VariableNames.MassFraction_n(i) + "#" + spc, tst.GetMassFractions(spc, i).Convert_Xt2X(0.0));
                }
            }
            if (tst.TestImmersedBoundary) {
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators_TimeDep.Add(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.Velocity_d(d)), tst.GetPhi2U(d));
                }
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCG, tst.GetPhi());

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            if (D == 3 && SurfTensionMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                Console.WriteLine($"Reminder: {SurfTensionMode} changed to LaplaceBeltrami_ContactLine for 3D test.");
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            } else {
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            }
            C.CutCellQuadratureType = CutCellQuadratureType;

            // immersed boundary
            // =================

            C.UseImmersedBoundary = tst.TestImmersedBoundary;
            if (C.UseImmersedBoundary) {
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), tst.GetPhi2());
            }

            // timestepping and solver
            // =======================

            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;
                C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.verbose = true;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            //C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.NonLinearSolver.MaxSolverIterations = 3;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver = solvercode.GetConfig();
            C.GravityDirection = tst.GravityDirection;
            C.ChemicalReactionActive = tst.ChemicalReactionTermsActive;
            C.EnableMassFractions = tst.EnableMassFractions;
            C.EnableTemperature = tst.EnableTemperature;
            C.rhoOne = constantDensity;
            // return
            // ======
            Assert.AreEqual(C.UseImmersedBoundary, tst.TestImmersedBoundary);
            return C;
        }

        public static void COMBUSTION_TEST() {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            using (var p = new XNSEC_MixtureFraction()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_CounterDiffusionFlame(2, 6, 1.0, db.Path);

                p.Init(c);
                p.RunSolverMode();
            }

            Console.WriteLine("Flame sheet calculation done.");
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_CounterDiffusionFlame(2, 6, 1.0, db.Path);
                p.Init(c);
                p.RunSolverMode();

                var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                double minT; double maxT;
                temp.GetExtremalValues(out minT, out maxT);
                Console.WriteLine("Maximum reached temperature is {0}K", maxT);
            }
            Console.WriteLine("Full calculation done.");
        }
        /// <summary>
        /// pseudo 1D combustion test
        /// </summary>
        /// <param name="dg"></param>
        /// <param name="ncells"></param>
        [Test]
        public static void XDG_PSEUDO1D_COMBUSTION_TEST([Values(3)] int dg = 3, [Values(4)] int ncells = 4) {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            // activate for combustion calculation, if deactivated, only evaporation
            using (var p = new XNSEC_MixtureFraction()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_XDG_pseudo2dCombustion(dg, ncells, db.Path);

                p.Init(c);
                p.RunSolverMode();
            }

            Console.WriteLine("Flame sheet calculation done.");
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_XDG_pseudo2dCombustion(dg, ncells, db.Path);
                p.Init(c);
                p.RunSolverMode();

                var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var massFraction0Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction0).SingleOrDefault());
                var massFraction1Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction1).SingleOrDefault());

                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                var mF0 = massFraction0Xdg.ProjectToSinglePhaseField(4);
                var mF1 = massFraction1Xdg.ProjectToSinglePhaseField(4);

                double minT; double maxT;
                double minMF0; double maxMF0;
                double minMF1; double maxMF1;
                temp.GetExtremalValues(out minT, out maxT);
                mF0.GetExtremalValues(out minMF0, out maxMF0);
                mF1.GetExtremalValues(out minMF1, out maxMF1);

                Console.WriteLine("Maximum reached temperature is {0}K", maxT);
                Console.WriteLine("Maximum reached massfraction0 is {0}", maxMF0);
                Console.WriteLine("Maximum reached massfraction1 is {0}", maxMF1);

                Assert.Greater(maxT, 2.6, "Maximum temperature is smaller than expected");
                Assert.Less(maxT, 2.75, "Maximum temperature is higher than expected");

                Assert.Greater(maxMF0, 0.95, "Maximum massfraction0 is smaller than expected");
                Assert.Less(maxMF0, 1.05, "Maximum massfraction0 is higher than expected");

                Assert.Greater(maxMF1, 0.95, "Maximum massfraction1 is smaller than expected");
                Assert.Less(maxMF1, 1.05, "Maximum massfraction1 is higher than expected");

            }
            Console.WriteLine("Full calculation done.");
        }
        /// <summary>
        /// 2D droplet combustion test with fixed mass flux
        /// </summary>
        /// <param name="dg"></param>
        /// <param name="ncells"></param>
        //[Test]
        public static void XDG_DROPLET_COMBUSTION_TEST(int dg = 2, int ncells = 6) {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            // activate for combustion calculation, if deactivated, only evaporation
            using (var p = new XNSEC_MixtureFraction()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_XDG_Droplet_2dCombustion(3, 7, db.Path);
                //c.SkipSolveAndEvaluateResidual = true;
                p.Init(c);
                p.RunSolverMode();
            }

            Console.WriteLine("Flame sheet calculation done.");
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_XDG_Droplet_2dCombustion(3, 7, db.Path);
                p.Init(c);
                p.RunSolverMode();


                var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var massFraction0Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction0).SingleOrDefault());
                var massFraction1Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction1).SingleOrDefault());

                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                var mF0 = massFraction0Xdg.ProjectToSinglePhaseField(4);
                var mF1 = massFraction1Xdg.ProjectToSinglePhaseField(4);

                double minT; double maxT;
                double minMF0; double maxMF0;
                double minMF1; double maxMF1;
                temp.GetExtremalValues(out minT, out maxT);
                mF0.GetExtremalValues(out minMF0, out maxMF0);
                mF1.GetExtremalValues(out minMF1, out maxMF1);

                Console.WriteLine("Maximum reached temperature is {0}K", maxT);
                Console.WriteLine("Maximum reached massfraction0 is {0}", maxMF0);
                Console.WriteLine("Maximum reached massfraction1 is {0}", maxMF1);

                Assert.Greater(maxT, 1.7, "Maximum temperature is smaller than expected");
                Assert.Less(maxT, 1.8, "Maximum temperature is higher than expected");

                Assert.Greater(maxMF0, 0.9, "Maximum massfraction0 is smaller than expected");
                Assert.Less(maxMF0, 1.1, "Maximum massfraction0 is higher than expected");

                Assert.Greater(maxMF1, 0.15, "Maximum massfraction1 is smaller than expected");
                Assert.Less(maxMF1, 0.3, "Maximum massfraction1 is higher than expected");

            }
            Console.WriteLine("Full calculation done.");
        }


        /// <summary>
        /// pseudo 1D evaporation test with fixed mass flux
        /// </summary>
        /// <param name="dg"></param>
        /// <param name="ncells"></param>
        [Test]
        public static void XDG_PSEUDO1D_EVAPORATION_TEST([Values(3)] int dg = 3, [Values(4)] int ncells = 4) {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            // activate for combustion calculation, if deactivated, only evaporation
            //using (var p = new XNSEC_MixtureFraction()) {
            //    var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_XDG_pseudo2dEvaporation(dg, ncells, db.Path);

            //    p.Init(c);
            //    p.RunSolverMode();
            //}

            Console.WriteLine("Flame sheet calculation done.");
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_XDG_pseudo2dEvaporation(dg, ncells, db.Path);
                p.Init(c);
                p.RunSolverMode();

                var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var massFraction0Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction0).SingleOrDefault());
                var massFraction1Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction1).SingleOrDefault());

                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                var mF0 = massFraction0Xdg.ProjectToSinglePhaseField(4);
                var mF1 = massFraction1Xdg.ProjectToSinglePhaseField(4);

                double minT; double maxT;
                double minMF0; double maxMF0;
                double minMF1; double maxMF1;
                temp.GetExtremalValues(out minT, out maxT);
                mF0.GetExtremalValues(out minMF0, out maxMF0);
                mF1.GetExtremalValues(out minMF1, out maxMF1);

                Console.WriteLine("Maximum reached temperature is {0}K", maxT);


                Assert.Greater(maxT, 0.95, "Maximum temperature is smaller than expected");
                Assert.Less(maxT, 1.05, "Maximum temperature is higher than expected");

                Assert.Greater(maxMF0, 0.95, "Maximum massfraction0 is smaller than expected");
                Assert.Less(maxMF0, 1.05, "Maximum massfraction0 is higher than expected");

                Assert.Greater(maxMF1, 0.95, "Maximum massfraction1 is smaller than expected");
                Assert.Less(maxMF1, 1.05, "Maximum massfraction1 is higher than expected");
            }
            Console.WriteLine("Full calculation done.");
        }

        /// <summary>
        /// 2D droplet evaporation test wit fixed mass flux
        /// </summary>
        /// <param name="dg"></param>
        /// <param name="ncells"></param>
        [Test]
        public static void XDG_DROPLET_EVAPORATION_TEST([Values(2)] int dg = 2, [Values(6)] int ncells = 6) {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            // activate for combustion calculation, if deactivated, only evaporation
            //using (var p = new XNSEC_MixtureFraction()) {
            //    var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_XDG_Droplet_2dEvaporation(3, 7, db.Path);
            //    //c.SkipSolveAndEvaluateResidual = true;
            //    p.Init(c);
            //    p.RunSolverMode();
            //}

            Console.WriteLine("Flame sheet calculation done.");
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_XDG_Droplet_2dEvaporation(3, 7, db.Path);
                p.Init(c);
                p.RunSolverMode();

                var temperatureXdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.Temperature).SingleOrDefault());
                var massFraction0Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction0).SingleOrDefault());
                var massFraction1Xdg = (XDGField)(p.CurrentStateVector.Fields.Where(f => f.Identification == VariableNames.MassFraction1).SingleOrDefault());

                var temp = temperatureXdg.ProjectToSinglePhaseField(4);
                var mF0 = massFraction0Xdg.ProjectToSinglePhaseField(4);
                var mF1 = massFraction1Xdg.ProjectToSinglePhaseField(4);

                double minT; double maxT;
                double minMF0; double maxMF0;
                double minMF1; double maxMF1;
                temp.GetExtremalValues(out minT, out maxT);
                mF0.GetExtremalValues(out minMF0, out maxMF0);
                mF1.GetExtremalValues(out minMF1, out maxMF1);

                Console.WriteLine("Maximum reached temperature is {0}K", maxT);

                Assert.Greater(maxT, 0.95, "Maximum temperature is smaller than expected");
                Assert.Less(maxT, 1.05, "Maximum temperature is higher than expected");

                Assert.Greater(maxMF0, 0.95, "Maximum massfraction0 is smaller than expected");
                Assert.Less(maxMF0, 1.05, "Maximum massfraction0 is higher than expected");

                Assert.Greater(maxMF1, 0.15, "Maximum massfraction1 is smaller than expected");
                Assert.Less(maxMF1, 0.3, "Maximum massfraction1 is higher than expected");

            }
            Console.WriteLine("Full calculation done.");
        }

        /// <summary>
        /// This test checks if the solution obtained with an homotopy strategy for increasing
        /// the velocity inlets (by introducing a multiplier as the homotopy-variable)
        /// is the same as the one if one wouldnt use homotopy
        /// Obviously a moderate velocity multiplier has to be used to obtain a solution without homotopy.
        /// </summary>
        //[Test]
        public static void CounterDiffFlameHomotopy_TEST() {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, @"default_bosss_db_comb23");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            // First, do a calculation without homotopy
            double desiredVelMultiplier = 4.0;
            XDGField VelocityX_NoHomotopy;
            XDGField VelocityX_WithHomotopy;

            List<double> L2NormsWithoutHomotopy;
            List<double> L2NormsWithHomotopy;
            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_CounterDiffusionFlame(2, 6, desiredVelMultiplier, db.Path);
                c.activeAMRlevelIndicators.Clear(); // no mesh refinement
                c.AdaptiveMeshRefinement = false;
                c.NoOfTimesteps = 1;
                p.Init(c);
                p.RunSolverMode();
            }

            // Now do same simulation, using homotopy

            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_CounterDiffusionFlame(2, 6, desiredVelMultiplier, db.Path); //
                c.activeAMRlevelIndicators.Clear(); // no mesh refinement
                c.AdaptiveMeshRefinement = false;
                c.NoOfTimesteps = 1;

                c.HomotopyApproach = XNSEC_Control.HomotopyType.Automatic;
                c.HomotopyVariable = XNSEC_Control.HomotopyVariableEnum.VelocityInletMultiplier;
                c.homotopieAimedValue = desiredVelMultiplier;

                p.Init(c);
                p.RunSolverMode();
            }

            var fieldsNoHomotopy = db.Sessions[1].Timesteps.Last().Fields;
            var fieldsWithHomotopy = db.Sessions[0].Timesteps.Last().Fields;

            foreach (var field in fieldsNoHomotopy) {
                var fieldWithHomotopy = fieldsWithHomotopy.Where(f => f.Identification == field.Identification).FirstOrDefault();
                //Console.WriteLine(field.Identification);
                //Console.WriteLine(field.L2Norm() - fieldWithHomotopy.L2Norm());
                Assert.LessOrEqual(field.L2Norm() - fieldWithHomotopy.L2Norm(), 1e-4, "Comparison of Norms for variable " + field.Identification + " show that they are different when using homotopy");
            }
        }

        public static void COMBUSTION_CoFlowFlame_TEST() {
            //string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            //if (basepath.IsEmptyOrWhite())
            //    basepath = System.Environment.GetEnvironmentVariable("HOME");
            string basepath = Directory.GetCurrentDirectory();
            string path = Path.Combine(basepath, "default_bosss_db_CoFlowcombustion");

            bool alreadyExists = Directory.Exists(path);
            var db = DatabaseInfo.CreateOrOpen(path);

            int rank;
            var comm = csMPI.Raw._COMM.WORLD;
            csMPI.Raw.Comm_Rank(comm, out rank);

            if (rank == 0 && alreadyExists) {
                db.Controller.ClearDatabase();
            }

            //using(var p = new XNSEC_MixtureFraction()) {
            //    var c = BoSSS.Application.XNSEC.FullNSEControlExamples.FS_CoFlowDiffusionFlame(2, 7, 0.5, db.Path);
            //    p.Init(c);
            //    p.RunSolverMode();
            //}

            //Console.WriteLine("Flame sheet calculation done.");

            using (var p = new XNSEC()) {
                var c = BoSSS.Application.XNSEC.FullNSEControlExamples.Full_CoFlowDiffusionFlame(2, 7, 1, 0.5, db.Path);
                p.Init(c);
                p.RunSolverMode();
            }

            Console.WriteLine("Full calculation done.");
        }
    }
}