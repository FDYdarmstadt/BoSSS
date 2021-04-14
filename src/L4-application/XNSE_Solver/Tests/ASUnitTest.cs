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
            [Values(2,3)] int spatialDimension,
            [Values(1, 2, 3, 4)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode
#endif
            )
        {

            var Tst = new ViscosityJumpTest(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfTensionMode);
            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;

            C.SkipSolveAndEvaluateResidual = C.AdvancedDiscretizationOptions.CellAgglomerationThreshold <= 1e-6;

            XNSESolverTest(Tst, C);
        }

#if !DEBUG
        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ScalingViscosityJumpTest_p2(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType)
        {
            ScalingViscosityJumpTest(2, vmode, CutCellQuadratureType);
        }

        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ScalingViscosityJumpTest_p3(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            )
        {
            ScalingViscosityJumpTest(3, vmode, CutCellQuadratureType);
        }
#endif

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ViscosityJumpTest"/>
        /// </summary>
        public static void ScalingViscosityJumpTest(int deg, ViscosityMode vmode, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType)
        {

            double AgglomerationTreshold = 0.1;
            int spatialDimension = 2;

            var Tst = new ViscosityJumpTest(spatialDimension);
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 4, 8, 16 })
            {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local,
                    CutCellQuadratureType: CutCellQuadratureType);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingViscosityJumpTest-p" + deg);
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
        public static void ScalingStaticDropletTest(
            int deg,
            ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            )
        {

            double AgglomerationTreshold = 0.1;

            var Tst = new StaticDropletTest();
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 2, 4, 8 })
            {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold,
                    vmode: vmode,
                    GridResolution: Res,
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType: CutCellQuadratureType);
                //C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;
                //C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

                C.InitSignedDistance = false;

                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingStaticDropletTest-p" + deg);
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
            foreach (var Res in new[] { 1, 2, 3, 4 })
            {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res, SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, CutCellQuadratureType: CutCellQuadratureType, nonlinsolver: nonlinsolver);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: false, title: "ScalingSinglePhaseChannelTest-p" + deg);
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
            )
        {

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

            if(deg == 3 && AgglomerationTreshold <= 0.01)
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
        public static void TaylorCouetteConvergenceTest_IBM(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.TestIBM;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_LaplaceBeltrami_Flux(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SchurCompl, nonlinsolver: nonlinsolver);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        [Test]
        public static void TaylorCouetteConvergenceTest_2Phase_Curvature_Projected(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {
            Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.Test2Phase;
            TaylorCouetteConvergenceTest(FlowSolverDegree, modus, SurfaceStressTensor_IsotropicMode.Curvature_Projected, SchurCompl, nonlinsolver: nonlinsolver);
        }

#endif

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
        /// </summary>
        public static void TaylorCouetteConvergenceTest(
            [Values(2, 3)] int FlowSolverDegree = 3,
            [Values(Tests.TaylorCouette.Mode.Test2Phase, Tests.TaylorCouette.Mode.TestIBM)] Tests.TaylorCouette.Mode modus = Tests.TaylorCouette.Mode.TestIBM,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux,
            [Values(false,true)] bool SchurCompl = true,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard
            ) {

            double AgglomerationTreshold = 0.3;
            
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;

            ViscosityMode vmode;
            switch(modus) {
                case TaylorCouette.Mode.Test2Phase: vmode = ViscosityMode.FullySymmetric; break;
                case TaylorCouette.Mode.TestIBM: vmode = ViscosityMode.Standard; break; // for IBM, the FullySymmetric implementation is still missing
                default: throw new ArgumentOutOfRangeException();
            }

            int[] GridResolutionS = new int[] { 2, 4, 8, 16 };

            var Tst = new TaylorCouette(modus);

            var CS = new XNSE_Control[GridResolutionS.Length];
            for(int i = 0; i < GridResolutionS.Length; i++) {
                var C = TstObj2CtrlObj(Tst, FlowSolverDegree, AgglomerationTreshold, vmode, SurfTensionMode: stm, CutCellQuadratureType: CutCellQuadratureType, GridResolution: GridResolutionS[i], nonlinsolver: nonlinsolver);
                C.SkipSolveAndEvaluateResidual = false;
                C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
                C.NonLinearSolver.ConvergenceCriterion = 1e-10;
                C.UseSchurBlockPrec = SchurCompl;
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;
                CS[i] = C;

                //Console.WriteLine("!!!!!!!!!!!!!!!!!!!!!!!!1   remove me !!!!!!!!!!!!!!!!!!!!!!1");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 3;
                //C.SkipSolveAndEvaluateResidual = false;
                //C.UseSchurBlockPrec = true;
                //XNSESolverTest(Tst, C);
                //break;
            }

            XNSESolverConvergenceTest(Tst, CS, true, new double[] { FlowSolverDegree, FlowSolverDegree, FlowSolverDegree - 1 } ); // be **very** generous with the expected slopes
        }

        /// <summary>
        /// <see cref="TaylorCouette_CurvElm"/>
        /// </summary>
        public static void CurvedElementsTest(
            [Values(2, 3)] int FlowSolverDegree = 1,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver= NonLinearSolverCode.Picard
            ) {

            var Tst = new TaylorCouette_CurvElm();
            var C = TstObj2CtrlObj(Tst, FlowSolverDegree, 
                AgglomerationTreshold:0.1, 
                vmode: ViscosityMode.Standard, 
                CutCellQuadratureType: XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, 
                SurfTensionMode: SurfaceStressTensor_IsotropicMode.Curvature_Projected, 
                GridResolution: 1, 
                nonlinsolver: nonlinsolver);
            C.NoOfMultigridLevels = 1;
            C.ImmediatePlotPeriod = 1;
            C.SuperSampling = 4;
            C.SkipSolveAndEvaluateResidual = true;
            XNSESolverTest(Tst, C);
            
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.TaylorCouette"/>
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
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
#else
            [Values(2, 3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(0.0, 60.0 * Math.PI / 180.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver
#endif      
            )
        {
            var Tst = new ChannelTest(angle);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);

            XNSESolverTest(Tst, C);
            if (deg == 2)
                ASScalingTest(Tst, new[] { 2, 3, 4 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);
            if (deg == 3)
                ASScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver);
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
            )
        {

            var Tst = new TranspiratingChannelTest(U2, periodicity);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, nonlinsolver: nonlinsolver); // surface tension plays no role in this test, so ignore it
            //C.SkipSolveAndEvaluateResidual = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MaxSolverIterations = 100;
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
            )
        {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new PolynomialTestForConvection(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, nonlinsolver: nonlinsolver);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            XNSESolverTest(Tst, C);
        }

        /// <summary>
        /// Simple pure Heatconductivity Test with Temperature kink at the Interface
        /// </summary>
        [Test]
        public static void HeatConductivityTest(
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
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
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm,
            [Values(NonLinearSolverCode.Newton, NonLinearSolverCode.Picard)] NonLinearSolverCode nonlinsolver) {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTest(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, 2, nonlinsolver: nonlinsolver);
            XNSFESolverTest(Tst, C);
        }


        private static void XNSESolverTest(IXNSETest Tst, XNSE_Control C) {
            
            if(Tst.SpatialDimension == 3) {
                Console.WriteLine($"Reminder: skipping 3D test for now...");
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
                for(int i = 0; i < ErrThresh.Length; i++) {
                    bool ok = LastErrors[i] <= ErrThresh[i];
                    Console.Write("L2 error, '{0}': \t{1}", solver.Operator.DomainVar[i], LastErrors[i]);

                    if(ok)
                        Console.WriteLine("   (ok)");
                    else
                        Console.WriteLine("   Above Threshold (" + ErrThresh[i] + ")");
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if (solver.CurrentResidual.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for (int i = 0; i < ResNorms.Length; i++)
                {
                    ResNorms[i] = solver.CurrentResidual.Fields[i].L2Norm();
                    bool ok = ResNorms[i] <= ResThresh[i];
                    Console.Write("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Fields[i].Identification, ResNorms[i]);

                    if(ok)
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


        private static void XNSESolverConvergenceTest(IXNSETest Tst, XNSE_Control[] CS, bool useExactSolution, double[] ExpectedSlopes) {
            int D = Tst.SpatialDimension;
            int NoOfMeshes = CS.Length;

            double[] hS = new double[NoOfMeshes];
            MultidimensionalArray errorS = null;
            string[] Names = null;

            XNSE[] solvers = new XNSE[NoOfMeshes];
            if(useExactSolution) {

                if(NoOfMeshes < 2)
                    throw new ArgumentException("At least two meshes required for convergence against exact solution.");

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
                            if(ExpectedSlopes.Length != Names.Length)
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
            } else {
                if(NoOfMeshes < 3)
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

                for(int i = 0; i < yValues.Length; i++) {
                    v1 += (xValues[i] - xAvg) * (yValues[i] - yAvg);
                    v2 += Math.Pow(xValues[i] - xAvg, 2);
                }

                double a = v1 / v2;
                double b = yAvg - a * xAvg;

                return a;
            }


            for(int i = 0; i < errorS.GetLength(1); i++) {
                var slope = LogLogRegression(hS, errorS.GetColumn(i));

                Console.WriteLine($"Convergence slope for Error of '{Names[i]}': \t{slope}\t(Expecting: {ExpectedSlopes[i]})");
            }

            for(int i = 0; i < errorS.GetLength(1); i++) {
                var slope = LogLogRegression(hS, errorS.GetColumn(i));
                Assert.IsTrue(slope >= ExpectedSlopes[i], $"Convergence Slope of {Names[i]} is degenerate.");
            }

            foreach(var s in solvers) {
                s.Dispose();
            }
        }


        private static void XHeatSolverTest(IXHeatTest Tst, XNSE_Control C) {
            using (var solver = new XHeat()) {

                solver.Init(C);
                solver.RunSolverMode();
                solver.OperatorAnalysis();

                //-------------------Evaluate Temperature Error ---------------------------------------- 
                var evaluator = new XHeatErrorEvaluator<XNSE_Control>(solver);
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
        private static void XNSFESolverTest(IXNSFETest Tst, XNSE_Control C) {

            using (var solver = new XNSFE()) {

                solver.Init(C);
                solver.RunSolverMode();
                solver.OperatorAnalysis();

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

        private static void ASScalingTest(IXNSETest Tst, int[] ResolutionS, ViscosityMode vmode, int deg, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType, SurfaceStressTensor_IsotropicMode SurfTensionMode, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard)
        {
#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.1;

            var LaLa = new List<XNSE_Control>();
            foreach (var Res in ResolutionS)
            {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode, GridResolution: Res, nonlinsolver: nonlinsolver);
                C.SkipSolveAndEvaluateResidual = false;
                LaLa.Add(C);
            }
            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: Name);
#endif
        }

        static XNSE_Control TstObj2CtrlObj(IXNSETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1,
            NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Picard) {
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

            foreach(var kv in tst.GetBoundaryConfig()) {
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
            if(tst.TestImmersedBoundary) {
                for(int d = 0; d < D; d++) {
                    C.InitialValues_Evaluators_TimeDep.Add(VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(1), VariableNames.Velocity_d(d)), tst.GetPhi2U(d));
                }
            }


            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCG, tst.GetPhi());

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = AgglomerationTreshold;
            if(D == 3 && SurfTensionMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                Console.WriteLine($"Reminder: {SurfTensionMode} changed to LaplaceBeltrami_ContactLine for 3D test.");
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            } else {
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            }
            C.CutCellQuadratureType = CutCellQuadratureType;

            // immersed boundary
            // =================

            C.UseImmersedBoundary = tst.TestImmersedBoundary;
            if(C.UseImmersedBoundary) {
                C.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), tst.GetPhi2());
            }

            // timestepping and solver
            // =======================


            if(tst.steady) {
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
            C.NonLinearSolver.SolverCode = nonlinsolver;

            // return
            // ======
            Assert.AreEqual(C.UseImmersedBoundary, tst.TestImmersedBoundary);
            return C;
        }


        class AS_XHeat_Control : XNSFE_Control {
            public override Type GetSolverType() {
                return typeof(XHeat);
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
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = AgglomerationTreshold;
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
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;

            // return
            // ======

            return C;
        }

        static AS_XHeat_Control TstObj2CtrlObj(IXNSFETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton) {
            AS_XHeat_Control C = new AS_XHeat_Control();
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
            }

            C.Phi = tst.GetPhi();
            C.InitialValues_Evaluators.Add("Phi", tst.GetPhi().Convert_Xt2X(0.0));

            // advanced spatial discretization settings
            // ========================================

            C.AdvancedDiscretizationOptions.ViscosityMode = vmode;
            C.AdvancedDiscretizationOptions.CellAgglomerationThreshold = AgglomerationTreshold;
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
            C.LinearSolver.ConvergenceCriterion = 1e-9;
            //C.Solver_ConvergenceCriterion = 1e-9;

            C.LinearSolver.SolverCode = LinearSolverCode.classic_pardiso;
            C.NonLinearSolver.SolverCode = nonlinsolver;

            // return
            // ======

            return C;
        }
    }
}
