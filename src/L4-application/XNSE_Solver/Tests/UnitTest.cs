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

using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using ilPSP;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// A collection of all-up NUnit tests for the XNSE solver.
    /// </summary>
    [TestFixture]
    static public partial class UnitTest {

        /// <summary>
        /// MPI finalize.
        /// </summary>
        [TestFixtureTearDown]
        static public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init.
        /// </summary>
        [TestFixtureSetUp]
        static public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);
            XQuadFactoryHelper.CheckQuadRules = true;
        }



        [Test]
        public static void ViscosityJumpTest(
#if DEBUG
            [Values(1)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode
#else
            [Values(1, 2, 3, 4)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode
#endif
            ) {

            var Tst = new ViscosityJumpTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            C.SkipSolveAndEvaluateResidual = C.AdvancedDiscretizationOptions.CellAgglomerationThreshold <= 1e-6;
                
            GenericTest(Tst, C);
        }

        [Test]
        public static void BcTest_PressureOutletTest(
            [Values(1)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true, false)] bool performsolve
            ) {
            //XNSE_ConsistencyTestMain p = null;
            //BoSSS.Solution.Application._Main(new string[0], true, "", delegate() {
            //    p = new XNSE_ConsistencyTestMain();
            //    p.Testcase = new BcTest_PressureOutlet();
            //    p.FlowSolverDegree = deg;
            //    p.m_dntParams.CellAgglomerationThreshold = AgglomerationTreshold;
            //    p.SolverMode = performsolve;
            //    p.m_dntParams.ViscosityMode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter
            //    return p;
            //});

            var Tst = new BcTest_PressureOutlet();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, ViscosityMode.Standard);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            GenericTest(Tst, C);
        }

        
        //[Test]
        public static void MovingDropletTest(
#if DEBUG
            [Values(1)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.8)] double Radius,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
#else
            [Values(2, 3)] int deg,
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
#endif
            ) {

            if (deg == 3 && AgglomerationTreshold <= 0.01)
                return;

            var Tst = new MovingDropletTest(Radius, includeConvection, bSteady);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, SurfTensionMode: stm);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            GenericTest(Tst, C);
        }

        [Test]
        public static void ChannelTest(
#if DEBUG
            [Values(2)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0)] double angle
#else
            [Values(3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0, 60.0 * Math.PI / 180.0)] double angle
#endif
            ) {

            var Tst = new ChannelTest(angle);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);

            GenericTest(Tst, C);
        }


        
        [Test]
        public static void TranspiratingChannelTest(
            [Values(2)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(0.0, 0.1)] double U2,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(true,false)] bool periodicity
            ) {

            var Tst = new TranspiratingChannelTest(U2, periodicity);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            //C.SkipSolveAndEvaluateResidual = true;
            C.Solver_MaxIterations = 100;
            GenericTest(Tst, C);
        }
        

        [Test]
        public static void PolynomialTestForConvectionTest(
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false)] bool SolverMode_performsolve
            ) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new PolynomialTestForConvection();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            GenericTest(Tst, C);
        }



        private static void GenericTest(ITest Tst, XNSE_Control C) {
            using (var solver = new XNSE_SolverMain()) {
                
                solver.Init(C);
                solver.RunSolverMode();

                // matrix analysis
                // ===============

                solver.SpatialOperatorMatrixAnalysis(true, !C.SkipSolveAndEvaluateResidual ? 1 : 0);

                // check residuals and errors
                // ==========================

                double[] LastErrors = solver.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt);
                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", solver.CurrentSolution.Mapping.Fields[i].Identification, LastErrors[i]);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Mapping.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Mapping.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Mapping.Fields[i].Identification, ResNorms[i]);
                }
                                
                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);
                
                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);
            }
        }

        static XNSE_Control TstObj2CtrlObj(ITest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode, 
            SurfaceStressTensor_IsotropicMode SurfTensionMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local) {
            XNSE_Control C = new XNSE_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            C.SetFieldOptions(FlowSolverDegree, tst.LevelsetPolynomialDegree);

            // grid
            // ====

            C.GridFunc = tst.CreateGrid;

            // boundary conditions
            // ===================

            foreach (var kv in tst.GetBoundaryConfig()) {
                C.BoundaryValues.Add(kv);
            }

            // Physical parameters
            // ====================

            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;

            // initial values and exact solution
            // =================================

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));

                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                    C.InitialValues_Evaluators.Add(VariableNames.Gravity_d(d) + "#" + spc, tst.GetF(spc, d));
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = AgglomerationTreshold;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;


            // timestepping and solver
            // =======================


            if (tst.steady) {
                C.CompMode = AppControl._CompMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.CompMode = AppControl._CompMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.Prescribed;
                C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }

            C.Solver_ConvergenceCriterion = 1e-9;

            // return
            // ======

            return C;

        }
    }
}
