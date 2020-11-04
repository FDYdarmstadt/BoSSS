﻿/* =======================================================================
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
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP.Connectors.Matlab;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// A collection of all-up NUnit tests for the XNSE solver.
    /// </summary>
    [TestFixture]
    static public partial class ASUnitTest {

        /// <summary>
        /// MPI init.
        /// </summary>
        [OneTimeSetUp]
        static public void OneTimeSetUp() {
            XQuadFactoryHelper.CheckQuadRules = true;
        }


        [Test]
        public static void ASViscosityJumpTest(
            [Values(1)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode)
        {
            var Tst = new ViscosityJumpTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            C.SkipSolveAndEvaluateResidual = C.AdvancedDiscretizationOptions.CellAgglomerationThreshold <= 1e-6;
            C.ImmediatePlotPeriod = 1;
            ApplicationWithSolverTest(Tst, C);
        }


#if !DEBUG
        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ASScalingViscosityJumpTest_p2(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode
            ) {
            ASScalingViscosityJumpTest(2, vmode);
        }

        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ASScalingViscosityJumpTest_p3(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode
            ) {
            ASScalingViscosityJumpTest(3, vmode);
        }

        /// <summary>
        /// <see cref="ViscosityJumpTest"/>
        /// </summary>
        public static void ASScalingViscosityJumpTest(

            [Values(2, 3)] int deg,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode
            ) {

            double AgglomerationTreshold = 0.1;

            var Tst = new ViscosityJumpTest();
            var LaLa = new List<XNSE_Control>();
            foreach(var Res in new[] { 4, 8, 16 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingViscosityJumpTest-p" + deg);
        }
#endif


#if !DEBUG        
        /// <summary>
        /// <see cref="SinglePhaseChannel"/>
        /// </summary>
        [Test]
        public static void ASScalingSinglePhaseChannelTest(
            [Values(1, 2, 3)] int deg,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode) //
        {

            double AgglomerationTreshold = 0.1;

            var Tst = new SinglePhaseChannel(0.0);
            var LaLa = new List<XNSE_Control>();
            foreach(var Res in new[] { 1, 2, 3, 4 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res);
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingSinglePhaseChannelTest-p" + deg);
        }
#endif      

        [Test]
        public static void ASBcTest_PressureOutletTest(
            [Values(1)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true, false)] bool performsolve
            )
        {
            var Tst = new BcTest_PressureOutlet();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, ViscosityMode.Standard);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            ApplicationWithSolverTest(Tst, C);
        }

        [Test]
        public static void ASMovingDropletTest(
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
            )
        {

            if (deg == 3 && AgglomerationTreshold <= 0.01)
                return;

            var Tst = new MovingDropletTest(Radius, includeConvection, bSteady);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, SurfTensionMode: stm);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            ApplicationWithSolverTest(Tst, C);
        }

        [Test]
        public static void ASChannelTest(
#if DEBUG
            [Values(2)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0)] double angle
#else
            [Values(2, 3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0, 60.0 * Math.PI / 180.0)] double angle
#endif
            ) {

            var Tst = new ChannelTest(angle);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);

            ApplicationWithSolverTest(Tst, C);
            if(deg < 3)
                ASScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg);

        }
      
        /// <summary>
        /// <see cref="Tests.TranspiratingChannelTest"/>
        /// </summary>
        [Test]
        public static void ASTranspiratingChannelTest(
            [Values(2)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(0.0, 0.1)] double U2,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(true, false)] bool periodicity
            ) {

            var Tst = new TranspiratingChannelTest(U2, periodicity);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MaxSolverIterations = 100;
            ApplicationWithSolverTest(Tst, C);
        }

        /// <summary>
        /// <see cref="Tests.PolynomialTestForConvection"/>
        /// </summary>
        [Test]
        public static void ASPolynomialTestForConvectionTest(
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false)] bool SolverMode_performsolve
            )
        {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new PolynomialTestForConvection();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            ApplicationWithSolverTest(Tst, C);
        }

        private static void ApplicationWithSolverTest(ITest Tst, XNSE_Control C) {
            using(var solver = new XNSE()) {

                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();

                //solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                ErrorEvaluator evaluator = new ErrorEvaluator(solver);
                double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++)
                {
                    Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++)
                {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);
            }
        }

        private static void ASScalingTest(ITest Tst, int[] ResolutionS, ViscosityMode vmode, int deg) {
#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.1;
                        
            var LaLa = new List<AS_XNSE_Control>();
            foreach(var Res in ResolutionS) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res);
                C.SkipSolveAndEvaluateResidual = false;
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: Name);
#endif
        }

        class AS_XNSE_Control : XNSE_Control
        {
            public override Type GetSolverType()
            {
                return typeof(XNSE);
            }
        }

        static AS_XNSE_Control TstObj2CtrlObj(ITest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode, 
            SurfaceStressTensor_IsotropicMode SurfTensionMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local,
            int GridResolution = 1) {
            AS_XNSE_Control C = new AS_XNSE_Control();
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

            C.GridFunc = () => tst.CreateGrid(GridResolution);

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

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            // return
            // ======

            return C;

        }
    }
}
