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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;


namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// An all-up NUnit test for the XNSEC application.
    /// </summary>
    [TestFixture]
    public static partial class NUnitTest {
        /// <summary>
        /// Simple Test for Evaporation of a straight interface with two chemical species
        /// </summary>
        //[Test]
        public static void SteadyStateEvaporationTest_XNSEC_MixtureFraction(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm
            ) {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTestXNSEC_MixtureFraction(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, true, 7);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            C.ImmediatePlotPeriod = 1;
            C.PlotNewtonIterations = true;
            C.PlotAdditionalParameters = false;
            //C.AgglomerationThreshold = 0.1;

            XNSECSolverTest(Tst, C);
        }

        private static XNSEC_Control TstObj2CtrlObj(IXNSECTest_MixtureFraction tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
                 XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
                 SurfaceStressTensor_IsotropicMode SurfTensionMode,
                 bool constantDensity,
                 int GridResolution = 1, LinearSolverCode solvercode = LinearSolverCode.direct_pardiso) {
            XNSEC_Control C = new XNSEC_Control();
            int D = tst.SpatialDimension;
            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSEC/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

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
            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;

            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;

            C.prescribedMassflux_Evaluator = (tst is IPrescribedMass ? (tst as IPrescribedMass).GetPrescribedMassflux_Evaluator() : null);
            C.ThermalParameters.hVap = (tst is IPrescribedMass ? 1.0 : 0.0);
            // initial values and exact solution
            // =================================

            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionMixtureFraction = new Dictionary<string, Func<double[], double, double>>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));

                C.ExactSolutionMixtureFraction.Add(spc, tst.GetMixtureFraction(spc));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));
                C.InitialValues_Evaluators.Add(VariableNames.MixtureFraction + "#" + spc, tst.GetMixtureFraction(spc).Convert_Xt2X(0.0));
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

            C.rhoOne = constantDensity;
            // return
            // ======
            Assert.AreEqual(C.UseImmersedBoundary, tst.TestImmersedBoundary);
            return C;
        }
        private static void XNSECSolverTest(IXNSECTest_MixtureFraction Tst, XNSEC_Control C) {
            using (var solver = new XNSEC_MixtureFraction()) {
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
                var MixtureFractionEvaluator = new MixtureFractionErrorEvaluator<XNSEC_Control>(solver);
                double[] MF_Errors= MixtureFractionEvaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                List<double> errors = new List<double>();
                errors.AddRange(XNSE_Errors);
                errors.AddRange(MF_Errors);
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


    }
}