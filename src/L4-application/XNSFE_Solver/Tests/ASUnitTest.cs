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
using BoSSS.Solution.AdvancedSolvers.Testing;
using ilPSP.Connectors.Matlab;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Application.XNSE_Solver.Tests;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Solution.Gnuplot;
using System.Diagnostics;
<<<<<<< HEAD
//using BoSSS.Foundation.Quadrature;
//using ilPSP.LinSolvers.MUMPS;
using static BoSSS.Solution.AdvancedSolvers.Testing.ConditionNumberScalingTest;
=======
using System.IO;
using BoSSS.Foundation.IO;

<<<<<<< HEAD
>>>>>>> 7c3a165f9c (third stage of parameterized LS implementation)
=======
using BoSSS.Solution.LevelSetTools.ParameterizedLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;

>>>>>>> 716f72fae4 (forth stage of parameterized LS implementation)

namespace BoSSS.Application.XNSFE_Solver.Tests {

    /// <summary>
    /// A collection of all-up NUnit tests for the new XNSE solver (<see cref="XNSE"/>).
    /// </summary>
    [TestFixture]
    static public partial class ASUnitTest {

        /// <summary>
        /// Scaling Test for XNSFE Operator Analysis in increasing complexity.
        /// Basically a:
        ///     0. Horizontal interface
        ///     1. Two-Phase Heat conduction
        ///     2. Horizontal interface + Heat conduction
        ///     3. Horizontal interface + Evaporation, 10/2021, scaling not archieved, leave this test out for now
        /// </summary>
        [Test]
        public static void XNSFEScalingTest(
            [Values(3)] int deg,
            [Values(0, 1, 2)] int Setup,
            [Values(true, false)] bool EqualFluids) {

            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new XNSFEScalingTest(Setup, EqualFluids);

            var LaLa = new List<XNSFE_Control>();
            foreach (var Res in new[] { 4, 8, 16 }) {
                var C = TstObj2CtrlObj(Tst, deg, 0.1,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.Saye);
                C.SkipSolveAndEvaluateResidual = true;

                C.InitSignedDistance = false;

                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = "XSNFEScalingTest-p" + deg+"-Setup" + Setup });
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSFE_Solver.Tests.ParameterizedLevelSetTest"/>
        /// </summary>
        [Test]
        public static void ParameterizedLevelSetTest(

            [Values(2)] int deg

            ) {

            string basepath = System.Environment.GetEnvironmentVariable("USERPROFILE");
            if (basepath.IsEmptyOrWhite())
                basepath = System.Environment.GetEnvironmentVariable("HOME");

            string path = Path.Combine(basepath, "CapillaryHeatTest");
            var db = DatabaseInfo.CreateOrOpen(path);

            double AgglomerationTreshold = 0.1;

            var Tst = new ParameterizedLevelSetTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, ViscosityMode.FullySymmetric, XQuadFactoryHelper.MomentFittingVariants.Saye, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, 16);

            C.ProjectName = "XNSFE_ParamLevelSetTestInCapillary";
            string jobName = C.ProjectName;
            C.SessionName = " " + jobName;

            C.savetodb = true;
            C.DbPath = db.Path;

            C.ImmediatePlotPeriod = 1;

            //C.solveCoupledHeatEquation = true;
            //C.IncludeRecoilPressure = true;

            //C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.ConstrainedDG;
            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;//ParameterizedLevelSet;
            //if (C.Option_LevelSetEvolution == LevelSetEvolution.ParameterizedLevelSet){
            //    C.ParameterizedLevelSetControl = new ParameterizedLevelSetControl(2 * 1e-04 / (3.0).Sqrt(), 2 * 1e-04 / (3.0).Sqrt(), 3e-04);
                
            //}
            //C.PhysicalParameters.theta_e = Math.PI / 6.0;
            //C.PhysicalParameters.betaS_A = 0;
            //C.PhysicalParameters.betaS_B = 0;
            //C.PhysicalParameters.betaL = 0.0;
            //C.PhysicalParameters.sliplength = 0.0;
            //C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            //C.SkipSolveAndEvaluateResidual = false;

            XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// Convergence Test for XNSFE using a manufactured solution.
        /// </summary>
        [Test]
        public static void XNSFEConvergenceTest([Values(2)] int deg) {

            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new XNSFEConvergenceTest();
            int[] GridResolutionS = new int[] { 3, 5, 9, 17 };
            //int[] GridResolutionS = new int[] { 8, 16, 32, 64 };
            var CS = new XNSFE_Control[GridResolutionS.Length];
            for(int i = 0; i< GridResolutionS.Length; i++) {
                CS[i] = TstObj2CtrlObj(Tst, deg, 0.1,
                    vmode: vmode,
                    GridResolution: GridResolutionS[i],
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.Saye);
            }

            XNSFESolverConvergenceTest(Tst, CS, true, new double[] { deg, deg, deg - 1, deg }); // be **very** generous with the expected slopes
        }        

        /// <summary>
        /// Convergence Test for XNSFE using a manufactured solution.
        /// 10/2021, currently scaling is off, test inactive
        /// </summary>
        public static void XNSFEConvergenceTestScaling([Values(2)] int deg) {

            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new XNSFEConvergenceTest();
            int[] GridResolutionS = new int[] { 3, 5, 9, 17 };
            var CS = new XNSFE_Control[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                CS[i] = TstObj2CtrlObj(Tst, deg, 0.1,
                    vmode: vmode,
                    GridResolution: GridResolutionS[i],
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.Saye);
                CS[i].SkipSolveAndEvaluateResidual = true;
            }

            ConditionNumberScalingTest.Perform(CS, new ConditionNumberScalingTest.Config());
        }

        /// <summary>
        /// Simple pure Heatconductivity Test with Temperature kink at the Interface
        /// </summary>
        [Test]
        public static void HeatConductivityTest(
            [Values(3)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new HeatConductivityTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XHeatSolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Heatconductivity Test to ensure Energy conservation of Heat equation
        /// </summary>
        [Test]
        public static void HeatDecayTest(
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.8598)] double r,
            [Values(-130, -50, 0.0, 10, 167)] double q,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm) {
            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new HeatDecayTest(r, q);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XHeatSolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Test for Evaporation of a straight interface
        /// </summary>
        [Test]
        public static void SteadyStateEvaporationTest(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver) // evaporation currently only implemented with use of newton solver
            {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTest(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 2, nonlinsolver: nonlinsolver);
            XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Test for Evaporation of a straight interface, Test Splitting / Moving Mesh
        /// Currently the Test would be run only one timestep, is this even meaningful?
        /// </summary>
        public static void TransientEvaporationTest(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(LevelSetHandling.LieSplitting, LevelSetHandling.Coupled_Once)] LevelSetHandling levelSetHandling) // evaporation currently only implemented with use of newton solver
            {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new TransientEvaporationTest(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 2, nonlinsolver: nonlinsolver, lsHandling: levelSetHandling);
            XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// Simple Test special "double cut" Quadrules, Test that the pressure jump and the Temperature profile are calculated correct
        /// </summary>
        public static void ZwoLsBasicTest(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver) // evaporation currently only implemented with use of newton solver
            {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter
             
            // TODO

            //var Tst = new SteadyStateEvaporationTest(rawangle * Math.PI / 180.0);
            //var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 2, nonlinsolver: nonlinsolver);
            //XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// A simple Shear flow with interfacial slip and evaporation
        /// <see cref="BoSSS.Application.XNSFE_Solver.Tests.InterfaceSlipTest"/>.
        /// 
        /// NOTE: something about this test is fishy; it seems to fail on certain machines, especially when multi-threading is activated.
        /// However, it also fails in non-multithreaded execution on certain computers,
        /// e.g., it passes on an i7-6700 and  i7-9700K, but fails on i7-11800H.
        /// </summary>
        [Test]
        public static void InterfaceSlipTestLin(
#if DEBUG
            [Values(3)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver,
            [Values(1.0, double.PositiveInfinity)] double slipI,
            [Values(0.143)] double viscosityratio,
            [Values(1.2)] double massflux
#else
            [Values(3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(0.0, 1.0, double.PositiveInfinity)] double slipI,
            [Values(1.0, 0.143)] double viscosityratio,
            [Values(0.0, 0.27, 1.2)] double massflux
#endif
            ) {
            var Tst = new InterfaceSlipTestLin(angle, slipI, viscosityratio, massflux);

            if (deg == 3 && slipI == 0) {
                //Note (fk,24jan24)
                //Ich habe mal die Fehler aus mehreren Runs gesammelt.
                //Die Gemeinsamkeit scheint zu sein, 
                //Param 1: es tritt nur bei DG-Grad 3 auftritt(und nicht bei 4);
                //                Param 8: es tritt nur bei slipI = 0 auf.
                //                (Param 2-- 7 werden im Test gar nicht nicht variiere, sondern sind immer[0.0d, FullySymmetric, 0.0d, Saye, Newton]).
                //                //                   1 2    3              4    5    7      7    8      9
                //                InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 1.2d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 0.27d)

                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 1.2d)

                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 0.27d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 1.2d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 0.27d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 1.2d)

                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 0.27d)

                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 0.27d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 1.0d, 1.2d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 0.27d)
                //InterfaceSlipTestLin(3, 0.0d, FullySymmetric, 0.0d, Saye, Newton, 0.0d, 0.143d, 1.2d)



                ilPSP.Environment.NumThreads = 1;


            }
            //Quadrature_Settings.ENABLE_MULTITHREAD_CHECKING = true;

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, nonlinsolver: nonlinsolver);
            C.PhysicalParameters.slipI = Tst.slipI;
            C.NonLinearSolver.MaxSolverIterations = 10;

            // clear initial values, such that not only consistency is checked
            C.InitialValues.Clear();
            C.InitialValues_Evaluators.Clear();

            C.Phi = Tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", Tst.GetPhi().Convert_Xt2X(0.0));

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// A simple Shear flow with interfacial slip and evaporation
        /// <see cref="BoSSS.Application.XNSFE_Solver.Tests.InterfaceSlipTest"/>
        /// </summary>
        [Test]
        public static void InterfaceSlipTestNonLin(
#if DEBUG
            [Values(3)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver,
            [Values(1.0, double.PositiveInfinity)] double slipI,
            [Values(0.143)] double viscosityratio,
            [Values(1.2)] double massflux
#else
            [Values(3)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(0.0, 1.0, double.PositiveInfinity)] double slipI,
            [Values(1.0, 0.143)] double viscosityratio,
            [Values(0.27, 1.2)] double massflux
#endif
            ) {
            var Tst = new InterfaceSlipTestNonLin(angle, slipI, viscosityratio, massflux);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, nonlinsolver: nonlinsolver, GridResolution: 3);
            C.PhysicalParameters.slipI = Tst.slipI;
            C.NonLinearSolver.MaxSolverIterations = 20;

            // clear initial values, such that not only consistency is checked
            C.InitialValues.Clear();
            C.InitialValues_Evaluators.Clear();

            C.Phi = Tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", Tst.GetPhi().Convert_Xt2X(0.0));

            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();

            XNSFESolverTest(Tst, C);
        }

        /// <summary>
        /// A simple Shear flow with interfacial slip and evaporation, test scaling
        /// <see cref="BoSSS.Application.XNSFE_Solver.Tests.InterfaceSlipTest"/>
        /// Different setups <paramref name="setup"/> possible:
        /// 0: linear (no convection no recoil pressure)
        /// 1: semilinear (no ceonvection)
        /// 2: nonlinear (full complexity) - currently disabled, the scaling is not achieved as in <see cref="XNSFEScalingTest(int, int, bool)"/>
        /// </summary>
        #if !DEBUG
        [Test]
        #endif
        public static void InterfaceSlipTestScaling(
            [Values(2, 3)] int deg,
            [Values(0, 1)] byte setup,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(1.0)] double slipI,
            [Values(1.21)] double viscosityratio,
            [Values(0.27)] double massflux
            ) {

            IXNSFETest Tst;
            switch (setup) {
                case 0:
                    Tst = new InterfaceSlipTestLin(angle, slipI, viscosityratio, massflux);
                    break;
                case 1:
                    Tst = new InterfaceSlipTestLin(angle, slipI, viscosityratio, massflux);
                    break;
                case 2:
                    Tst = new InterfaceSlipTestNonLin(angle, slipI, viscosityratio, massflux);
                    break;
                default:
                    throw new NotImplementedException();
            }

            var LaLa = new List<XNSFE_Control>();
            foreach (var Res in new[] { 1, 3, 5, 7, 9 }) {
                var C = TstObj2CtrlObj(Tst, deg, 0.1,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.Saye);

                if(setup == 0) {
                    C.IncludeRecoilPressure = false; // disable recoil, the exact coded solution is not analytical for that case, so the residuals wont be close to zero.
                }

                C.SkipSolveAndEvaluateResidual = true;
                C.PhysicalParameters.slipI = ((ISlipTest)Tst).slipI;

                LaLa.Add(C);
            }

            // somehow some stencil scaling are negative, but we allow this.
            var config = new ConditionNumberScalingTest.Config() { plot = true, title = "InterfaceSlipScalingTest-p" + deg + "-Setup" + setup };
            config.ExpectedSlopes[ConditionNumberScalingTest.Config.TotCondNo] = (XAxisDesignation.Grid_1Dres, 2.4, 1.5);
            config.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_innerUncut] = (XAxisDesignation.Grid_1Dres, 0.5, -1.0);
            config.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_innerCut] = (XAxisDesignation.Grid_1Dres, 0.5, -1.0);
            config.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_bndyUncut] = (XAxisDesignation.Grid_1Dres, 0.5, -1.0);
            config.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_bndyCut] = (XAxisDesignation.Grid_1Dres, 0.5, -1.0);

            ConditionNumberScalingTest.Perform(LaLa, config);
        }

        private static void XHeatSolverTest(IXHeatTest Tst, XNSFE_Control C) {
            using (var solver = new XHeat()) {
                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis(); // deavtivated; has only value for a series of meshes, but not for a single calc.

                //-------------------Evaluate Temperature Error ---------------------------------------- 
                var evaluator = new XHeatErrorEvaluator<XNSFE_Control>(solver);
                if (Tst.CheckT) {                    
                    double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                    double[] ErrThresh = Tst.AcceptableL2Error;
                    if (LastErrors.Length != ErrThresh.Length)
                        throw new ApplicationException();
                    for (int i = 0; i < ErrThresh.Length; i++) {
                        Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                    }
                    for (int i = 0; i < ErrThresh.Length; i++)
                        Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);
                }

                //-------------------Evaluate Energy Error ---------------------------------------- 
                if (Tst.CheckE) {
                    double LastErrors = evaluator.ComputeEnergyError(Tst.GetE(), Tst.dt);

                    double ErrThresh = Tst.AcceptableL2Error[0];   
                    
                    Console.WriteLine("relative error, '{0}': \t{1}", "total thermal Energy", LastErrors);

                    // less than one promille
                    Assert.LessOrEqual(LastErrors, 1e-3);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);
                }

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);
                
            }
        }
        private static void XNSFESolverTest(IXNSFETest Tst, XNSFE_Control C) {

            using (var solver = new XNSFE<XNSFE_Control>()) {

                solver.Init(C);
                solver.RunSolverMode();
                //solver.OperatorAnalysis(); // deavtivated; has only value for a series of meshes, but not for a single calc.

                //-------------------Evaluate Flow Error ---------------------------------------- 
                var flowevaluator = new XNSEErrorEvaluator<XNSFE_Control>(solver);
                double[] LastErrors = flowevaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length - 1) // Last Error thershold is for temperature
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length - 1; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                }
                for (int i = 0; i < LastErrors.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                //-------------------Evaluate Temperature Error ---------------------------------------- 
                var heatevaluator = new XHeatErrorEvaluator<XNSFE_Control>(solver);
                if (Tst.CheckT) {
                    LastErrors = LastErrors.Cat(heatevaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C) );

                    // Last Error belongs to temperature
                    for (int i = ErrThresh.Length - 1; i < ErrThresh.Length; i++) {
                        Console.WriteLine("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);
                    }
                    for (int i = ErrThresh.Length - 1; i < ErrThresh.Length; i++)
                        Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);
                }

                //-------------------Evaluate Energy Error ---------------------------------------- 
                if (Tst.CheckE) {
                    double LastError = heatevaluator.ComputeEnergyError(Tst.GetE(), Tst.dt);

                    Console.WriteLine("relative error, '{0}': \t{1}", "total thermal Energy", LastErrors);

                    // less than one promille
                    Assert.LessOrEqual(LastError, 1e-3);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);
                }

                for (int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);

            }
        }
        private static void XNSFESolverConvergenceTest(IXNSFETest Tst, XNSFE_Control[] CS, bool useExactSolution, double[] ExpectedSlopes) {
            int D = Tst.SpatialDimension;
            int NoOfMeshes = CS.Length;

            double[] hS = new double[NoOfMeshes];
            MultidimensionalArray errorS = null;
            string[] Names = null;

            XNSFE[] solvers = new XNSFE[NoOfMeshes];
            if (useExactSolution) {

                if (NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

                for (int k = 0; k < CS.Length; k++) {

                    Console.WriteLine("================================================================");
                    Console.WriteLine($"Convergence Test:  Run {k+1} of {CS.Length}");
                    Console.WriteLine("================================================================");

                    var C = CS[k];
                    //using(var solver = new XNSE()) {
                    var solver = new XNSFE();
                    solvers[k] = solver;
                    {
                        //Console.WriteLine("Warning! - enabled immediate plotting");
                        //C.ImmediatePlotPeriod = 1;
                        //C.SuperSampling = 3;

                        solver.Init(C);
                        solver.RunSolverMode();

                        //-------------------Evaluate Error ---------------------------------------- 
                        var evaluator = new XNSEErrorEvaluator<XNSFE_Control>(solver);
                        var heatevaluator = new XHeatErrorEvaluator<XNSFE_Control>(solver);
                        double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);
                        LastErrors = LastErrors.Cat(heatevaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C));
                        double[] ErrThresh = Tst.AcceptableL2Error;


                        if (k == 0) {
                            errorS = MultidimensionalArray.Create(NoOfMeshes, LastErrors.Length);
                            Names = new string[LastErrors.Length];
                            if (ExpectedSlopes.Length != Names.Length)
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
            } else {
                if (NoOfMeshes < 3)
                    throw new ArgumentException("At least three meshes required for convergence if finest solution is assumed to be exact.");
                throw new NotImplementedException("todo");
            }


            //hS = hS.Take(hS.Length - 1).ToArray();

            double LogLogRegression(IEnumerable<double> _xValues, IEnumerable<double> _yValues) {
                double[] xValues = _xValues.Select(x => Math.Log10(x)).ToArray();
                double[] yValues = _yValues.Select(y => Math.Log10(y)).ToArray();

                double xAvg = xValues.Average();
                double yAvg = yValues.Average();

                double v1 = 0.0;
                double v2 = 0.0;

                for (int i = 0; i < yValues.Length; i++) {
                    v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                    v2 += Math.Pow(xValues[i] - xAvg, 2);
                }

                double a = v1 / v2;
                double b = yAvg - a * xAvg;

                return a;
            }


            for (int i = 0; i < errorS.GetLength(1); i++) {
                var slope = LogLogRegression(hS, errorS.GetColumn(i));

                Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{slope}\t(Expecting: {ExpectedSlopes[i]})");
            }

            // set terminal, for automatic testing we do not want this output
            if (false) {
                int xRes = 1024;
                int yRes = 768;
                var gp = new Gnuplot();
                var fmt = new PlotFormat("rx-");
                int Kount = 1;
                for (int i = 0; i < errorS.GetLength(1); i++) {
                    gp.PlotXY(hS, errorS.GetColumn(i), logX: true, logY: true, title: Names[i], format: (fmt.WithLineColor(Kount).WithPointType(Kount)));
                    gp.SetXLabel("h-Grid");
                    Kount++;
                }
                gp.Execute();
                Console.WriteLine("plotting in interactive gnuplot session - press any key to continue...");
                Console.ReadKey();
            }

            for (int i = 0; i < errorS.GetLength(1); i++) {
                var slope = LogLogRegression(hS, errorS.GetColumn(i));
                Assert.IsTrue(slope >= ExpectedSlopes[i], $"Convergence Slope of {Names[i]} is degenerate.");
            }

            foreach (var s in solvers) {
                s.Dispose();
            }
        }

        class AS_XHeat_Control : XNSFE_Control {
            public override Type GetSolverType() {
                return typeof(XHeat);
            }
        }

        class AS_XNSFE_Control : XNSFE_Control {
            public override Type GetSolverType() {
                return typeof(XNSFE);
            }
        }

        static AS_XHeat_Control TstObj2CtrlObj(IXHeatTest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1) {
            AS_XHeat_Control C = new AS_XHeat_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XHEAT/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

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

            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;
            C.ThermalParameters.c_A = tst.c_A;
            C.ThermalParameters.c_B = tst.c_B;
            C.ThermalParameters.k_A = tst.k_A;
            C.ThermalParameters.k_B = tst.k_B;
            C.ThermalParameters.hVap = tst.h_vap;
            C.ThermalParameters.T_sat = tst.T_sat;
            C.ThermalParameters.IncludeConvection = tst.IncludeConvection;

            // initial values and exact solution
            // =================================

            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();

            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionTemperature.Add(spc, tst.GetT(spc));

                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetT(spc).Convert_Xt2X(0.0));
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            C.CutCellQuadratureType = CutCellQuadratureType;

            // timestepping and solver
            // =======================


            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }

            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();

            // return
            // ======

            return C;
        }

        static AS_XNSFE_Control TstObj2CtrlObj(IXNSFETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton,
            LevelSetHandling lsHandling = LevelSetHandling.LieSplitting) {
            AS_XNSFE_Control C = new AS_XNSFE_Control();
            int D = tst.SpatialDimension;


            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSFE/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========
            C.FieldOptions.Add("VelocityX", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("VelocityY", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Pressure", new FieldOpts() {
                Degree = FlowSolverDegree - 1,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Temperature", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatFlux*", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("PhiDG", new FieldOpts() {
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("Phi", new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add(VariableNames.Curvature, new FieldOpts() {
                Degree = tst.LevelsetPolynomialDegree ,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("GravityX", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            }); 
            C.FieldOptions.Add("GravityY", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });
            C.FieldOptions.Add("HeatSource", new FieldOpts() {
                Degree = FlowSolverDegree,
                SaveToDB = FieldOpts.SaveToDBOpt.TRUE
            });

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

            C.ThermalParameters.rho_A = tst.rho_A;
            C.ThermalParameters.rho_B = tst.rho_B;
            C.ThermalParameters.c_A = tst.c_A;
            C.ThermalParameters.c_B = tst.c_B;
            C.ThermalParameters.k_A = tst.k_A;
            C.ThermalParameters.k_B = tst.k_B;
            C.ThermalParameters.hVap = tst.h_vap;
            C.ThermalParameters.T_sat = tst.T_sat;
            C.ThermalParameters.IncludeConvection = tst.IncludeConvection;


            C.PhysicalParameters.rho_A = tst.rho_A;
            C.PhysicalParameters.rho_B = tst.rho_B;
            C.PhysicalParameters.mu_A = tst.mu_A;
            C.PhysicalParameters.mu_B = tst.mu_B;
            C.PhysicalParameters.Sigma = tst.Sigma;
            C.PhysicalParameters.IncludeConvection = tst.IncludeConvection;
            C.PhysicalParameters.theta_e = tst.theta_e;
            // initial values and exact solution
            // =================================

            C.ExactSolutionTemperature = new Dictionary<string, Func<double[], double, double>>();
            C.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            C.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();
                        
            foreach (var spc in new[] { "A", "B" }) {
                C.ExactSolutionTemperature.Add(spc, tst.GetT(spc));
                C.InitialValues_Evaluators.Add(VariableNames.Temperature + "#" + spc, tst.GetT(spc).Convert_Xt2X(0.0));

                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));

                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));
                for (int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators.Add(VariableNames.Velocity_d(d) + "#" + spc, tst.GetU(spc, d).Convert_Xt2X(0.0));
                }

                Func<double[], double, double>[] Gravity = new Func<double[], double, double>[D];
                for (int d = 0; d < D; d++) {
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }

                var Source = tst.GetQ(spc).Convert_X2Xt();
                C.SetHeatSource(spc, Source);
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            if (tst.TestImmersedBoundary) {
                C.FieldOptions.Add("Phi2DG", new FieldOpts() {
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
                C.FieldOptions.Add("Phi2", new FieldOpts() {
                    Degree = tst.LevelsetPolynomialDegree,
                    SaveToDB = FieldOpts.SaveToDBOpt.TRUE
                });
                C.UseImmersedBoundary = true;
                C.InitialValues_Evaluators.Add("Phi2", tst.GetPhi2().Convert_Xt2X(0.0));
                C.AdvancedDiscretizationOptions.DoubleCutSpecialQuadrature = true;
            }

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AgglomerationThreshold = AgglomerationTreshold;
            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            C.CutCellQuadratureType = CutCellQuadratureType;

            // timestepping and solver
            // =======================


            if (tst.steady) {
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                C.Option_LevelSetEvolution = LevelSetEvolution.None;
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            } else {
                C.TimesteppingMode = AppControl._TimesteppingMode.Transient;

                C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;
                C.Timestepper_LevelSetHandling = lsHandling;

                C.NoOfTimesteps = 1;
                C.dtFixed = tst.dt;
            }

            C.NonLinearSolver.ConvergenceCriterion = 1e-8;
            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.NonLinearSolver.SolverCode = nonlinsolver;

            // return
            // ======

            return C;
        }
    }
}
