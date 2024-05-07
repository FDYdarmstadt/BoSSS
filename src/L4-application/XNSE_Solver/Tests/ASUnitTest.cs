/* =======================================================================
Copyright 2017 
Technische Universitaet Darmstadt, 
Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

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
using BoSSS.Solution.Statistic;
using BoSSS.Solution.AdvancedSolvers;
using static BoSSS.Solution.AdvancedSolvers.Testing.ConditionNumberScalingTest;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// A collection of all-up NUnit tests for the new XNSE solver (<see cref="XNSE"/>).
    /// </summary>
    [TestFixture]
    static public partial class ASUnitTest {


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ViscosityJumpTest"/>
        /// </summary>
        [Test]
        public static void ViscosityJumpTest(
#if DEBUG
            [Values(2,3)] int spatialDimension,
            [Values(1)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode
#else
            [Values(2, 3)] int spatialDimension,
            [Values(1, 2, 3, 4)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode
#endif
            ) {

            var Tst = new ViscosityJumpTest(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfTensionMode);
            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.SkipSolveAndEvaluateResidual = C.AgglomerationThreshold <= 1e-6;
            XNSESolverTest(Tst, C);
        }

#if !DEBUG
        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ScalingViscosityJumpTest_p2(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType) {
            ScalingViscosityJumpTest(2, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ScalingViscosityJumpTest_p3(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) {
            ScalingViscosityJumpTest(3, vmode, CutCellQuadratureType);
        }
#endif

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ViscosityJumpTest"/>
        /// </summary>
        public static void ScalingViscosityJumpTest(int deg, ViscosityMode vmode, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType) {

            double AgglomerationTreshold = 0.1;
            int spatialDimension = 2;

            var Tst = new ViscosityJumpTest(spatialDimension);
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 4, 8, 16 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local,
                    CutCellQuadratureType: CutCellQuadratureType);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = "ScalingViscosityJumpTest-p" + deg });
        }

#if !DEBUG
        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p2_Standard_OneStepGaussAndStokes() //1
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p3_Standard_OneStepGaussAndStokes() //2
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p2_FullySymmetric_OneStepGaussAndStokes() //3
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p3_FullySymmetric_OneStepGaussAndStokes() //4
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p2_Standard_Saye() //5
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p3_Standard_Saye() //6
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p2_FullySymmetric_Saye() //7
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest_p3_FullySymmetric_Saye() //8
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            ScalingStaticDropletTest(deg, vmode, CutCellQuadratureType);
        }

#endif


        /// <summary>
        /// <see cref="XNSE_Solver.Tests.StaticDropletTest"/>
        /// </summary>
        [Test]
        public static void ScalingStaticDropletTest(
            [Values(1, 2, 3)] int deg,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) {

            double AgglomerationTreshold = 0.1;

            var Tst = new StaticDropletTest();
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 2, 4, 8 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: CutCellQuadratureType);
                //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
                //C.LinearSolver = LinearSolverCode.classic_pardiso.GetConfig();

                C.InitSignedDistance = false;

                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = "ScalingStaticDropletTest-p" + deg });
        }



#if !DEBUG        
        /// <summary>
        /// <see cref="SinglePhaseChannel"/>
        /// </summary>
        [Test]
        public static void ScalingSinglePhaseChannelTest(
            [Values(1, 2, 3)] int deg,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
            ) //

        {

            double AgglomerationTreshold = 0.1;

            var Tst = new SinglePhaseChannel(0.0);
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 1, 2, 3, 4 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res, SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, CutCellQuadratureType: CutCellQuadratureType, nonlinsolver: nonlinsolver);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = false, title = "ScalingSinglePhaseChannelTest-p" + deg });
        }
#endif      

        /// <summary>
        /// <see cref="BcTest_PressureOutlet"/>
        /// </summary>
        [Test]
        public static void BcTest_PressureOutletTest(
            [Values(2, 3)] int spatialDimension,
            [Values(1)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.Curvature_Projected, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode,
            [Values(true, false)] bool performsolve
            ) {

            var Tst = new BcTest_PressureOutlet(spatialDimension);
#if DEBUG
            const int GridRes = 4; // resolutions 1, 3, etc. place level-set at cell center
#else
            const int GridRes = 4; // resolutions 2, 4, etc. place level-set at cell boundary
#endif

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, ViscosityMode.Standard, GridResolution: GridRes, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            XNSESolverTest(Tst, C);
            if (spatialDimension == 2)
                ASScalingTest(Tst, new[] { 4, 8, 16 }, ViscosityMode.Standard, deg, CutCellQuadratureType, SurfTensionMode);

        }


#if DEBUG
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_debug(
            [Values(1)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.8)] double Radius,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) //
        {
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }
#else 
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p2_OneStepGaussAndStokes_Standard(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//1
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p3_OneStepGaussAndStokes_Standard(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//2
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p2_Saye_Standard(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//3
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p3_Saye_Standard(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//4
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p2_OneStepGaussAndStokes_FullySymmetric(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//5
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p3_OneStepGaussAndStokes_FullySymmetric(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//6
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p2_Saye_FullySymmetric(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//7
        {
            int deg = 2;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        [Test]
        public static void MovingDropletTest_rel_p3_Saye_FullySymmetric(
            [Values(0.01, 0.1, 0.3)] double AgglomerationTreshold,
            [Values(true)] bool performsolve,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(0.69711, 0.70611, 0.70711, 0.70811, 0.71711, 0.75468, 0.80226, 0.83984, 0.84884, 0.84984, 0.85084, 0.85984)] double Radius,
            [Values(true)] bool bSteady,
            [Values(false)] bool includeConvection
            )//8
        {
            int deg = 3;
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
            MovingDropletTest(deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }
#endif


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        public static void MovingDropletTest(
            int deg, double AgglomerationTreshold,
            bool performsolve,
            SurfaceStressTensor_IsotropicMode stm,
            double Radius,
            ViscosityMode vmode,
            bool bSteady,
            bool includeConvection,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType) //
       {

            if (deg == 3 && AgglomerationTreshold <= 0.01)
                return;

            var Tst = new MovingDropletTest(Radius, includeConvection, bSteady);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, SurfTensionMode: stm, CutCellQuadratureType: CutCellQuadratureType, GridResolution: 2);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            XNSESolverTest(Tst, C);
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.BasicThreePhase"/>
        /// </summary>
        [Test]
        public static void BasicThreePhaseTest(
            [Values(true, false)] bool bSteady = false,
            [Values(true, false)] bool performsolve = false,
            [Values(true, false)] bool bConvection = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux,
            [Values(2)] int spatDim = 2
            ) {
            double R = 0.8;
            int FlowSolverDegree = 2;
            double AgglomerationTreshold = 0.3;
            ViscosityMode vmode = ViscosityMode.Standard;
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

            int GridResolution = 1;

            var Tst = new BasicThreePhase(R, bConvection, bSteady, spatDim);

            var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, vmode, SurfTensionMode: stm, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, nonlinsolver: nonlinsolver);
            C.SkipSolveAndEvaluateResidual = !performsolve;


            //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!1   remove me !!!!!!!!!!!!!!!!!!!!!!1");
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 3;


            XNSESolverTest(Tst, C);

        }


#if !DEBUG
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_IBM_SchurComplOn(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            bool SchurCompl = true;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.TestIBM;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_IBM_SchurComplOff(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            bool SchurCompl = false;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.TestIBM;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }



        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_LaplaceBeltrami_Flux(
            [Values(2, 3)] int FlowSolverDegree,
            [Values(false, true)] bool SchurCompl,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver
            ) {
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>; Schur complement off
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_Curvature_Proj_Soff_p2(
            //[Values(2, 3)] int FlowSolverDegree = 3,
            //[Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {

            int FlowSolverDegree = 2;
            bool SchurCompl = false;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_Curvature_Proj_Son_p3(
            //[Values(2, 3)] int FlowSolverDegree = 3,
            //[Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver 
            ) {
            int FlowSolverDegree = 3;
            bool SchurCompl = true;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>; Schur complement off
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_Curvature_Proj_Soff_p3(
            //[Values(2, 3)] int FlowSolverDegree = 3,
            //[Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver 
            ) {

            int FlowSolverDegree = 3;
            bool SchurCompl = false;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_Curvature_Proj_Son_p2(
            //[Values(2, 3)] int FlowSolverDegree = 3,
            //[Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton/*, NonLinearSolverCode.Picard*/)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            int FlowSolverDegree = 2;
            bool SchurCompl = true;
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }
#endif

        /*
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_temp(
            [Values(2, 3)] int FlowSolverDegree,
            [Values(Tests.TaylorCouette.Mode.Test2Phase, Tests.TaylorCouette.Mode.TestIBM)] Tests.TaylorCouette.Mode modus
            ) //
        {
            // --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.TaylorCouetteConvergenceTest_temp
            SurfaceStressTensor_IsotropicMode stm = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux;
            bool SchurCompl = true;
            NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton;

            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, stm, SchurCompl, nonlinsolver);
        }
        */

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        public static void TaylorCouetteConvergenceTest(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(Tests.TaylorCouette.Mode.Test2Phase, Tests.TaylorCouette.Mode.TestIBM)] Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.TestIBM,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux,
            [Values(false, true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {

            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            ViscosityMode vmode;
            switch (modus) {
                case TaylorCouette.Mode.Test2Phase: vmode = ViscosityMode.FullySymmetric; break;
                case TaylorCouette.Mode.TestIBM: vmode = ViscosityMode.Standard; break; // for IBM, the FullySymmetric implementation is still missing
                default: throw new ArgumentOutOfRangeException();
            }

            int[] GridResolutionS = new int[] { 2, 4, 8, 16 };

            var Tst = new TaylorCouette(modus);

            var CS = new XNSE_Control[GridResolutionS.Length];
            for (int i = 0; i < GridResolutionS.Length; i++) {
                var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, vmode, SurfTensionMode: stm, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolutionS[i], nonlinsolver: nonlinsolver);
                C.SkipSolveAndEvaluateResidual = false;
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                C.NonLinearSolver.ConvergenceCriterion = 1e-10;
                C.UseSchurBlockPrec = SchurCompl;
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 3;
                
                CS[i] = C;
            }

            // wrongly scaled mu (see commit 4451965a7b620e7173949cabfd22c156aa9de23c, 16feb22, fk):
            //FlowSolverDegree: 2 modus: Test2Phase (VelocityX, 3.083483164355715, 0.7989661925786997) (VelocityY, 3.0834840895521585, 0.7989627364682361) (Pressure, 3.4906457282676495, 1.8723140951578137)
            //FlowSolverDegree: 2 modus: TestIBM (VelocityX, 2.515012543968083, 1.3838718070630835) (VelocityY, 2.515012543968087, 1.3838718070630929) (Pressure, 2.214610364766888, 0.3311383592014634)     
            //FlowSolverDegree: 3 modus: Test2Phase (VelocityX, 4.654926583305634, 1.0928175327396197) (VelocityY, 4.654927897587302, 1.0928285940487372) (Pressure, 3.130217351146254, 0.7776322436159941)  
            //FlowSolverDegree: 3 modus: TestIBM (VelocityX, 4.00386309467254, 0.8724039374232135) (VelocityY, 4.002991471367835, 0.8712665660677734) (Pressure, 2.6954476819767255, -0.2606463145099438)    

            // correctly scaled mu:
            //FlowSolverDegree: 2 modus: Test2Phase (VelocityX, 3.1254595326793955, 0.7404020602004735) (VelocityY, 3.1254601695618773, 0.7403994267116443) (Pressure, 3.437159541536483, 1.6811550756258642)
            //FlowSolverDegree: 2 modus: TestIBM (VelocityX, 3.6865441522172957, 0.855131845172898) (VelocityY, 3.686544152216214, 0.8551318451720138) (Pressure, 2.4618884145401307, -0.9380731392988562)
            //FlowSolverDegree: 3 modus: Test2Phase (VelocityX, 4.6639625753059795, 0.9897897821002215) (VelocityY, 4.6639506311723515, 0.9897856559892002) (Pressure, 3.195132673497341, 0.7451829927632292)
            //FlowSolverDegree: 3 modus: TestIBM (VelocityX, 4.453883287724354, -0.10906474316034576) (VelocityY, 4.453610133084815, -0.10508285176620813) (Pressure, 3.2604571163979905, -1.6391369361012353)
            var RegressionBounds = new Dictionary<(int Degree, TaylorCouette.Mode mode), (string Name, double Slope, double intercept, double interceptTol)[]>();
            RegressionBounds.Add((2, TaylorCouette.Mode.Test2Phase), new[] { ("VelocityX", 3.0, 0.740, 0.1), ("VelocityY", 3.0, 0.740, 0.1), ("Pressure", 2.0, 1.68, 0.1) });
            RegressionBounds.Add((2, TaylorCouette.Mode.TestIBM), new[] { ("VelocityX", 3.0, 0.855, 0.1), ("VelocityY", 3.0, 0.855, 0.1), ("Pressure", 2.0, -0.938, 0.1) });
            RegressionBounds.Add((3, TaylorCouette.Mode.Test2Phase), new[] { ("VelocityX", 4.0, 0.990, 0.1), ("VelocityY", 4.0, 0.990, 0.1), ("Pressure", 3.0, 0.745, 0.1) });
            RegressionBounds.Add((3, TaylorCouette.Mode.TestIBM), new[] { ("VelocityX", 4.0, -0.109, 0.1), ("VelocityY", 4.0, -0.105, 0.1), ("Pressure", 3.0, -1.639, 0.1) });
            
            XNSESolverConvergenceTest(Tst, CS, true, RegressionBounds[(FlowSolverDegree, modus)]);
            
            
            //writing the Reference Data:
            //var RegData = XNSESolverConvergenceTest(Tst, CS, true, RegressionBounds[(FlowSolverDegree, modus)]);
            //using(var dataFile = new System.IO.StreamWriter("RegData.txt", true)) {
            //    dataFile.Write("FlowSolverDegree: " + FlowSolverDegree);
            //    dataFile.Write(" modus: " + modus.ToString());
            //    dataFile.Write(" ");
            //    foreach(var ttt in RegData)
            //        dataFile.Write(ttt + " ");
            //    dataFile.WriteLine();
            //}


            // be **very** generous with the expected slopes
        }

        /// <summary>
        /// <see cref="TaylorCouette_CurvElm"/>
        /// </summary>
        public static void CurvedElementsTest(
            [Values(2, 3)] int FlowSolverDegree = 1,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {

            var Tst = new TaylorCouette_CurvElm();
            var C = TstObj2CtrlObj(Tst, FlowSolverDegree,
                AgglomerationTreshold: 0.1,
                vmode: ViscosityMode.Standard,
                CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes,
                SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                GridResolution: 1,
                nonlinsolver: nonlinsolver);
            C.NoOfMultigridLevels = 1;
            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 4;
            C.SkipSolveAndEvaluateResidual = true;
            XNSESolverTest(Tst, C);

        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.IBMChannel"/>
        /// </summary>
        [Test]
        public static void IBMChannelTest(
            [Values(1, 2, 3)] int FlowSolverDegree = 2,
            [Values(0, 30 * Math.PI / 180, 40 * Math.PI / 180, 45 * Math.PI / 180, 60 * Math.PI / 180, Math.PI / 2)] double angle = 0.0,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {

            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;


            int GridResolution = 1;

            var Tst = new IBMChannel(30 * Math.PI / 180, true);

            var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, ViscosityMode.Standard, SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, nonlinsolver: nonlinsolver);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;


            XNSESolverTest(Tst, C);

        }

        [Test]
        public static void IBMChannelSolverTest(
            [Values(1, 2, 3)] int FlowSolverDegree = 2,
            [Values(0)] double angle = 0.0,
            [Values(LinearSolverCode.direct_pardiso, LinearSolverCode.exp_Kcycle_schwarz, LinearSolverCode.exp_Kcycle_schwarz_CoarseMesh, LinearSolverCode.exp_Kcycle_schwarz_PerProcess, LinearSolverCode.exp_gmres_levelpmg)] LinearSolverCode solvercode = LinearSolverCode.direct_pardiso
            ) {
            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;


            int GridResolution = 3;

            var Tst = new IBMChannel(30 * Math.PI / 180, true);

            var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, ViscosityMode.Standard, SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, solvercode: solvercode);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            // set to transient...
            C.dtFixed = 0.1;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            // and reset to steady-state
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            Assert.IsTrue(C.Option_LevelSetEvolution == LevelSetEvolution.None, "resetting option for 'Option_LevelSetEvolution' for steady-state seems messed up");
            Assert.IsTrue(C.Option_LevelSetEvolution2 == LevelSetEvolution.None, "resetting option for 'Option_LevelSetEvolution2' for steady-state seems messed up");
            Assert.IsTrue(C.Timestepper_LevelSetHandling == LevelSetHandling.None, "resetting option for 'Timestepper_LevelSetHandling' for steady-state seems messed up");

            //var Ccomp = AppControl.Deserialize(System.IO.File.ReadAllText("control.obj"));

            XNSESolverTest(Tst, C);
        }


        [Test]
        public static void IBMChannelSolverTest_Transient(
            [Values(1, 2, 3)] int FlowSolverDegree = 2,
            [Values(0)] double angle = 0.0,
            [Values(LinearSolverCode.direct_pardiso, LinearSolverCode.exp_Kcycle_schwarz, LinearSolverCode.exp_Kcycle_schwarz_CoarseMesh, LinearSolverCode.exp_Kcycle_schwarz_PerProcess, LinearSolverCode.exp_gmres_levelpmg)] LinearSolverCode solvercode = LinearSolverCode.direct_pardiso
            ) {
            double AgglomerationTreshold = 0.3;

            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;


            int GridResolution = 3;

            var Tst = new IBMChannel(30 * Math.PI / 180, true);

            var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, ViscosityMode.Standard, SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolution, solvercode: solvercode);
            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.NonLinearSolver.ConvergenceCriterion = 1e-11;

            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 2;

            C.PhysicalParameters.IncludeConvection = false;

            // set to transient...
            C.dtFixed = 0.1;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            Assert.IsTrue(C.TimesteppingMode == AppControl._TimesteppingMode.Transient, "'TimesteppingMode' should be 'Transient' automatically when someone sets a small dt.");

            //var Ccomp = AppControl.Deserialize(System.IO.File.ReadAllText("control.obj"));

            XNSESolverTest(Tst, C);
        }



        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ChannelTest"/>
        /// </summary>
        [Test]
        public static void ChannelTest(
#if DEBUG
            [Values(2)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(true, false)] bool IncludeConvection, 
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
#else
            [Values(2, 3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0, 60.0 * Math.PI / 180.0)] double angle,
            [Values(true, false)] bool IncludeConvection, 
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
#endif      
            ) {
            var Tst = new ChannelTest(angle, IncludeConvection);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);

            XNSESolverTest(Tst, C);
            if (deg == 2)
                ASScalingTest(Tst, new[] { 2, 3, 4 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);
            if (deg == 3)
                ASScalingTest(Tst, new[] { 2, 3, 4, 5 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);
        }


        /// <summary>
        /// <see cref="TranspiratingChannelTest"/>
        /// </summary>
        [Test]
        public static void TranspiratingChannelTest(
            [Values(2)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(0.0, 0.1)] double U2,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(true, false)] bool periodicity,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
            ) {

            var Tst = new TranspiratingChannelTest(U2, periodicity);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver); // surface tension plays no role in this test, so ignore it
            //C.SkipSolveAndEvaluateResidual = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            //C.Solver_MaxIterations = 100;
            XNSESolverTest(Tst, C);
            //if(AgglomerationTreshold > 0) {
            //    ScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg);
            //}
        }

        /// <summary>
        /// <see cref="Tests.PolynomialTestForConvection"/>
        /// </summary>
        [Test]
        public static void PolynomialTestForConvectionTest(
            [Values(2, 3)] int spatialDimension,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
            ) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new PolynomialTestForConvection(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, nonlinsolver: nonlinsolver);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XNSESolverTest(Tst, C);
        }


        /// <summary>
        /// <see cref="SphericalHarmonicsTest"/>
        /// </summary>
        [Test]
        public static void SphericalHarmonicsPostprocessingTest(
            [Values(true, false)] bool OnQuarterDomain,
            [Values(3, 5)] int MaxLength,
            [Values(true, false)] bool RotationalSymmetric) {
            // --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.SphericalHarmonicsPostprocessingTest()

            if(MaxLength == 5 && (OnQuarterDomain || RotationalSymmetric))
                return;

            var Tst = new SphericalHarmonicsTest();
            Tst.ComputeOnQuarterDomain = OnQuarterDomain;
            Tst.IsRotationalSymmetric = RotationalSymmetric;
            var C = TstObj2CtrlObj(Tst, 2, 0.1,
                ViscosityMode.FullySymmetric, 
                XQuadFactoryHelper.MomentFittingVariants.Saye, 
                SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);

            var pp = new Logging.SphericalHarmonicsLogging() { MaxL = MaxLength, RotSymmetric = RotationalSymmetric};
            C.PostprocessingModules.Add(pp);
            C.SkipSolveAndEvaluateResidual = true;

            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            //C.AMR_startUpSweeps = 1;

            //C.ImmediatePlotPeriod = 1;
            //C.SuperSampling = 0;

            XNSESolverTest(Tst, C);

            double[] SH_modes = pp.LoggedValues.Last();
            double totErr = Tst.ComputeModeError(SH_modes);

            Console.WriteLine("Total mode error: " + totErr);

            Assert.LessOrEqual(totErr, 0.003, "Large Error on spherical mode decomposition.");

        }


#if !DEBUG
        /// <summary>
        /// Tests the combination of AMR and 
        /// The test is considered a success if no exception is thrown
        /// </summary>
        [Test]
        public static void AMRAndBDFTest([Values(LevelSetHandling.None, LevelSetHandling.LieSplitting)] LevelSetHandling lsh) {
            int _res = 10;
            var C = HardcodedControl.KarmanVortexStreet(k: 1, Res: _res, writeToDB: false, NoOfTimeSteps: 50);
            C.Timestepper_LevelSetHandling = lsh;

            Assert.IsTrue(C.AdaptiveMeshRefinement == true);
            Assert.IsTrue(C.activeAMRlevelIndicators.Count > 0);
            Assert.IsTrue(C.TimeSteppingScheme == TimeSteppingScheme.BDF2 || C.TimeSteppingScheme == TimeSteppingScheme.BDF3 || C.TimeSteppingScheme == TimeSteppingScheme.BDF4);

            long JEnd = 0;
            using (var solver = new XNSE()) {

                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();

                JEnd = solver.Grid.NumberOfCells;
            }

            Assert.IsTrue(JEnd != _res * _res, "AMR was probably not active.");
        }
#endif


        private static void XNSESolverTest(IXNSETest Tst, XNSE_Control C) {

            if (Tst.SpatialDimension == 3 && C.CutCellQuadratureType != XQuadFactoryHelper.MomentFittingVariants.Saye) {
                Console.WriteLine($"Reminder: only starting 3D tests with Saye's quadrature. All other 3D Tests are skipped automatically!");
                return;
            }

            using (var solver = new XNSE()) {

                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;

                solver.Init(C);
                solver.RunSolverMode();
                //if(C.TimesteppingMode == AppControl._TimesteppingMode.Steady) // deavtivated by FK; has only value for a series of meshes, but not for a single calc.
                //    solver.OperatorAnalysis();

                //-------------------Evaluate Error ---------------------------------------- 
                var evaluator = new XNSEErrorEvaluator<XNSE_Control>(solver);
                double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);

                double[] ErrThresh = Tst.AcceptableL2Error;
                if (LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ErrThresh.Length; i++) {
                    bool ok = LastErrors[i] <= ErrThresh[i];
                    Console.Write("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);

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
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i], $"Error {solver.CurrentState.Fields[i].Identification} above threshold.");

                for (int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i], $"Residual {solver.CurrentResidual.Fields[i].Identification} above threshold.");
            }
        }


        private static (string Name, double slope, double Intercept)[] XNSESolverConvergenceTest(IXNSETest Tst, XNSE_Control[] CS, bool useExactSolution,  (string Name, double Slope, double intercept, double interceptTol)[] RegResults) {
            int D = Tst.SpatialDimension;
            int NoOfMeshes = CS.Length;
            if(RegResults.Length != D + 1)
                throw new ArgumentException("Expecting slopes for velocity and pressure.");

            var Ret = new List<(string Name, double slope, double Intercept)>();

            if(useExactSolution) {
                if(NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

                MultidimensionalArray errorS = null;
                string[] Names = null;

                double[] hS = new double[NoOfMeshes];
                XNSE[] solvers = new XNSE[NoOfMeshes];

                for(int k = 0; k < CS.Length; k++) {

                    var C = CS[k];
                    //using(var solver = new XNSE()) {
                    var solver = new XNSE();
                    solvers[k] = solver;
                    {
                        //Console.WriteLine("Warning! - enabled immediate plotting");
                        //C.ImmediatePlotPeriod = 1;
                        //C.SuperSampling = 3;

                        solver.Init(C);
                        solver.RunSolverMode();

                        //-------------------Evaluate Error ---------------------------------------- 
                        var evaluator = new XNSEErrorEvaluator<XNSE_Control>(solver);
                        double[] LastErrors = evaluator.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt, C);
                        double[] ErrThresh = Tst.AcceptableL2Error;


                        if(k == 0) {
                            errorS = MultidimensionalArray.Create(NoOfMeshes, LastErrors.Length);
                            Names = new string[LastErrors.Length];
                            if(RegResults.Length != Names.Length)
                                throw new ArgumentOutOfRangeException();
                        } else {
                            if(LastErrors.Length != Names.Length)
                                throw new ApplicationException();
                        }

                        if(LastErrors.Length != ErrThresh.Length)
                            throw new ApplicationException();
                        for(int i = 0; i < ErrThresh.Length; i++) {
                            Console.WriteLine($"L2 error, '{solver.Operator.DomainVar[i]}': \t{LastErrors[i]}");
                            Names[i] = solver.Operator.DomainVar[i];
                        }

                        errorS.SetRow(k, LastErrors);
                        hS[k] = evaluator.GetGrid_h();
                    }

                }

                for(int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);
                    Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{RegModel.Slope}\tIntercept: \t{RegModel.Intercept}\t(Expecting: {RegRef.Slope}, {RegRef.intercept}+/-{RegRef.interceptTol})");
                    Ret.Add((Names[i], RegModel.Slope, RegModel.Intercept));
                }

                for(int i = 0; i < errorS.GetLength(1); i++) {
                    var RegModel = hS.LogLogRegression(errorS.GetColumn(i));
                    var RegRef = RegResults.Single(ttt => ttt.Name == Names[i]);


                    Assert.GreaterOrEqual(RegModel.Slope, RegRef.Slope, $"Convergence Slope of {Names[i]} is degenerate.");
                    Assert.GreaterOrEqual(RegModel.Intercept, RegRef.intercept - RegRef.interceptTol, $"Convergence Intercept of {Names[i]} is degenerate.");
                    Assert.LessOrEqual(RegModel.Intercept, RegRef.intercept + RegRef.interceptTol, $"Convergence Intercept of {Names[i]} overshoot.");
                }

                foreach(var s in solvers) {
                    s.Dispose();
                }
            } else {
                if(NoOfMeshes < 3)
                    throw new ArgumentException("At least three meshes required for convergence if finest solution is assumed to be exact.");



                BoSSS.Solution.Statistic.ConvergenceTest.SolverConvergenceTest_Experimental(
                    CS,
                    "Experimental Convergence",
                    (D + 1).ForLoop(iVar => (iVar < D ? VariableNames.Velocity_d(iVar) : VariableNames.Pressure,
                                             iVar < D ? NormType.L2_approximate : NormType.L2noMean_approximate,
                                             RegResults[iVar].Slope,
                                             RegResults[iVar].intercept, RegResults[iVar].interceptTol)));

            }



            return Ret.ToArray();

           
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.HardcodedControl.Rotating_Something_Unsteady"/>
        /// </summary>
        [Test]
        public static void ScalingRotCubeTests2D_p1([Values(0.1, 0.2, 0.3)] double AggTreshold  ) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotCubeTests2D_p1 --result=blabla.xml

            int deg = 1;
            ScalingRotCubeTest(deg, 2, AggTreshold);
        }

        [Test]
        public static void ScalingRotCubeTests2D_p2([Values(0.1, 0.2, 0.3)] double AggTreshold) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotCubeTests2D_p2 --result=blabla.xml

            int deg = 2;
            ScalingRotCubeTest(deg, 2, AggTreshold);
        }

        [Test]
        public static void ScalingRotCubeTests2D_p3([Values(0.1, 0.2, 0.3)] double AggTreshold) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotCubeTests2D_p3 --result=blabla.xml

            int deg = 3;
            ScalingRotCubeTest(deg, 2, AggTreshold);
        }

        public static void ScalingRotTorusTests2D_p1([Values(0.1, 0.2, 0.3)] double AggTreshold) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotTorusTests2D_p1 --result=blabla.xml

            int deg = 1;
            ScalingRotTorusTest(deg, 2, AggTreshold);
        }

        public static void ScalingRotTorusTests2D_p2([Values(0.1, 0.2, 0.3)] double AggTreshold) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotTorusTests2D_p2 --result=blabla.xml

            int deg = 2;
            ScalingRotTorusTest(deg, 2, AggTreshold);
        }

        public static void ScalingRotTorusTests2D_p3([Values(0.1, 0.2, 0.3)] double AggTreshold) {
            // PublicTestRunner.exe nunit3 'XNSE_Solver' --test=BoSSS.Application.XNSE_Solver.Tests.ASUnitTest.ScalingRotTorusTests2D_p3 --result=blabla.xml

            int deg = 3;
            ScalingRotTorusTest(deg, 2, AggTreshold);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.HardcodedControl.RotatingTiltedXRigid"/>
        /// </summary>
        public static void ScalingRotTorusTest(int deg = 1, int spatialDimension = 2, double AgglomerationTreshold = 0.1) {
            Console.WriteLine("Agg Threshold is set to " + AgglomerationTreshold);

            var Controls = new List<XNSE_Control>();

            foreach (var Res in new[] { 16, 32, 64 }) {
                var C = HardcodedControl.RotatingTiltedXRigid(k: deg, Res: Res , SpaceDim: spatialDimension, AMR: true, AMRLevel:1, shape: Shape.Torus, TiltAngle: 0.0, RotAxis: "z");
                C.Timestepper_LevelSetHandling = LevelSetHandling.None;
                C.PhysicalParameters.IncludeConvection = true;
                C.NonLinearSolver.ConvergenceCriterion = 10e-6;
                C.NoOfTimesteps = 1;
                C.AgglomerationThreshold = AgglomerationTreshold;
                Controls.Add(C);
            }

            ConditionNumberScalingTest.Perform(Controls, new ConditionNumberScalingTest.Config() { plot = true, title = $"Scaling{spatialDimension}DRotTorus-p{deg}-Agg{AgglomerationTreshold}" });
        }

        /// <summary>
        /// A simple Shear flow with interfacial slip
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.InterfaceSlipTest"/>
        /// </summary>
        [Test]
        public static void InterfaceSlipTest(
#if DEBUG
            [Values(2)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(1.0, double.PositiveInfinity)] double slipI,
            [Values(0.143)] double viscosityratio
#else
            [Values(2, 3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver,
            [Values(0.0, 1.0, double.PositiveInfinity)] double slipI,
            [Values(1.0, 0.143, 13.0)] double viscosityratio
#endif      
            ) {
            var Tst = new InterfaceSlipTest(angle, slipI, viscosityratio);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, nonlinsolver: nonlinsolver);
            C.PhysicalParameters.slipI = Tst.slipI;
            C.NonLinearSolver.MaxSolverIterations = 10;

            // clear initial values, such that not only consistency is checked
            C.InitialValues.Clear();
            C.InitialValues_Evaluators.Clear();

            C.Phi = Tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", Tst.GetPhi().Convert_Xt2X(0.0));

            XNSESolverTest(Tst, C);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.HardcodedControl.Rotating_Something_Unsteady"/>
        /// </summary>
        public static void ScalingRotCubeTest(int deg = 1, int spatialDimension = 2, double AgglomerationTreshold = 0.1) {
            Console.WriteLine("Agg Threshold is set to " + AgglomerationTreshold);

            var Controls = new List<XNSE_Control>();

            foreach (var Res in new[] { 2, 4, 8, 16}) {
                var C = HardcodedControl.Rotating_Something_Unsteady(k: deg, Res: Res * 5, SpaceDim: spatialDimension, useAMR: true, useLoadBal: true, Gshape: Shape.Cube);
                C.PhysicalParameters.IncludeConvection = true;
                C.NonLinearSolver.ConvergenceCriterion = 10e-6;
                C.NoOfTimesteps = 3;
                C.AgglomerationThreshold = AgglomerationTreshold;
                Controls.Add(C);
            }
            var conf = new ConditionNumberScalingTest.Config() { plot = true, title = $"Scaling{spatialDimension}DRotCube-p{deg}-Agg{AgglomerationTreshold}" };
            conf.ExpectedSlopes[ConditionNumberScalingTest.Config.TotCondNo] = (XAxisDesignation.Grid_1Dres, 2.4, 1.3);
            ConditionNumberScalingTest.Perform(Controls, conf);
        }

#if !DEBUG
        [Test]
        public static void RotatingCubeTest(
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants MFV,
            [Values(Newton.GlobalizationOption.None, Newton.GlobalizationOption.LineSearch, Newton.GlobalizationOption.Dogleg)] Newton.GlobalizationOption globalizationOption,
            [Values(IncompressibleBcType.Pressure_Outlet)] IncompressibleBcType __OuterBcType 
            ) {

            //this error does happen when number of processors =2 but not with np=4
            var C = HardcodedControl.Rotating_Something_Unsteady(k: 4, Res: 30, SpaceDim: 2, useAMR: false, Gshape: Shape.Cube, OuterBcType: __OuterBcType);

            //C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            //C.NonLinearSolver.ConvergenceCriterion = 1.0e-8;
            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Newton;
            C.NonLinearSolver.Globalization = globalizationOption;
            C.NonLinearSolver.ConvergenceCriterion = 0.0; // As accurate as possible
            C.NonLinearSolver.MinSolverIterations = 10;
            C.NonLinearSolver.MaxSolverIterations = 20;

            C.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            C.UseSchurBlockPrec = false;

            C.PhysicalParameters.IncludeConvection = true;
            C.NoOfTimesteps = 1;

            C.CutCellQuadratureType = MFV;

            C.TracingNamespaces = "BoSSS";


            XNSE.DeleteOldPlotFiles();
            using (var Slv = new XNSE()) {
                Slv.Init(C);
                Slv.RunSolverMode();

                var resNorms = new List<(string id, double val)>();
                foreach(XDGField f in Slv.CurrentResidual.Fields) {
                    double ResNorm = f.L2NormAllSpecies();
                    resNorms.Add((f.Identification, ResNorm));
                    Console.WriteLine($"Residual norm of {f.Identification}: {ResNorm}");
                }
                foreach(var tt in resNorms)
                    Assert.Less(tt.val, 1.0e-5, $"Residual norm of {tt.id} out of bounds.");
            }
            
        }
#endif

        private static void ASScalingTest(IXNSETest Tst, int[] ResolutionS, ViscosityMode vmode, int deg, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType, SurfaceStressTensor_IsotropicMode SurfTensionMode, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard) {
#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.1;

            var LaLa = new List<XNSE_Control>();
            foreach (var Res in ResolutionS) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode, GridResolution: Res, nonlinsolver: nonlinsolver);
                C.SkipSolveAndEvaluateResidual = false;
                LaLa.Add(C);
            }

            var conf = new ConditionNumberScalingTest.Config() { plot = true, title = Name };
            conf.ExpectedSlopes[ConditionNumberScalingTest.Config.StencilCondNo_innerCut] = (XAxisDesignation.Grid_1Dres, 0.5, -0.7);

            ConditionNumberScalingTest.Perform(LaLa, conf);
#endif
        }

        internal static XNSE_Control TstObj2CtrlObj(IXNSETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1,
            NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard,
            LinearSolverCode solvercode = LinearSolverCode.direct_pardiso) {
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
                    var Gravity_d = tst.GetF(spc, d).Convert_X2Xt();
                    C.SetGravity(spc, d, Gravity_d);
                }

                C.InitialValues_Evaluators.Add(VariableNames.Pressure + "#" + spc, tst.GetPress(spc).Convert_Xt2X(0.0));

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

            //C.Solver_ConvergenceCriterion = 1e-9;

            C.NonLinearSolver.SolverCode = nonlinsolver;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            
            C.LinearSolver = solvercode.GetConfig();
            if(C.LinearSolver is Solution.AdvancedSolvers.IterativeSolverConfig isc) {
                isc.ConvergenceCriterion = 1e-9;
            }

            // return
            // ======
            Assert.AreEqual(C.UseImmersedBoundary, tst.TestImmersedBoundary);
            return C;
        }
    }
}
