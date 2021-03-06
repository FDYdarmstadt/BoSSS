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
using BoSSS.Solution.LevelSetTools;
using BoSSS.Application.XNSE_Solver.Legacy;
using BoSSS.Application.XNSE_Solver.Tests;

namespace BoSSS.Application.XNSE_Solver.Legacy.LegacyTests {

    /// <summary>
    /// A collection of all-up NUnit tests for the old XNSE solver <see cref="XNSE_SolverMain"/>.
    /// </summary>
    [TestFixture]
    static public partial class UnitTest {

        /*

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ViscosityJumpTest"/>
        /// </summary>
        [Test]
        public static void ViscosityJumpTest(
#if DEBUG
            [Values(2, 3)] int spatialDimension,
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

            C.SkipSolveAndEvaluateResidual = C.AdvancedDiscretizationOptions.CellAgglomerationThreshold <= 1e-6;

            GenericTest(Tst, C);
        }
/*
#if !DEBUG
        /// <summary>
        /// scaling of condition number for polynomial order 3 (polynomial order parameter is unwrapped for better parallelism of test execution)
        /// </summary>
        [Test]
        public static void ScalingViscosityJumpTest_p2(
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType           
            ) {

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
            int[] resolutions = new[] { 4, 8, 16 };
            //if (spatialDimension == 3)
            //    resolutions = new[] { 1, 2, 4 }; 

            var LaLa = new List<XNSE_Control>();

            foreach (var Res in resolutions) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, 
                    vmode: vmode, 
                    GridResolution: Res, 
                    SurfTensionMode:SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, 
                    CutCellQuadratureType:CutCellQuadratureType);
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingViscosityJumpTest-p" + deg);
        }


#if !DEBUG
        /// <summary>
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
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
        /// <see cref="ViscosityJumpTest"/>
        /// </summary>
        public static void ScalingStaticDropletTest(
            [Values(2, 3)] int deg,
            [Values(ViscosityMode.Standard, ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) {

            double AgglomerationTreshold = 0.1;

            var Tst = new StaticDropletTest();
            var LaLa = new List<XNSE_Control>();
            foreach (var Res in new[] { 2, 4, 8 }) { 
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, 
                    vmode: vmode, 
                    GridResolution: Res, 
                    SurfTensionMode: SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine,
                    CutCellQuadratureType:CutCellQuadratureType);
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
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) //
        {

            double AgglomerationTreshold = 0.1;

            var Tst = new SinglePhaseChannel(0.0);
            var LaLa = new List<XNSE_Control>();
            foreach(var Res in new[] { 1, 2, 3, 4 }) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, GridResolution: Res, SurfTensionMode:SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, CutCellQuadratureType:CutCellQuadratureType);
                C.ImmediatePlotPeriod = 1;
                C.SuperSampling = 2;
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: "ScalingSinglePhaseChannelTest-p" + deg);
        }
#endif      
*/

        /// <summary>
        /// <see cref="BcTest_PressureOutlet"/>
        /// </summary>
        //[Test]
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

            GenericTest(Tst, C);
            //if(spatialDimension == 2)
            //    ScalingTest(Tst, new[] { 4, 8, 16 }, ViscosityMode.Standard, deg, CutCellQuadratureType, SurfTensionMode);

        }


#if DEBUG
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_debug(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }
#else 
        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p2_OneStepGaussAndStokes_Standard(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p3_OneStepGaussAndStokes_Standard(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p2_Saye_Standard(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p3_Saye_Standard(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p2_OneStepGaussAndStokes_FullySymmetric(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p3_OneStepGaussAndStokes_FullySymmetric(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p2_Saye_FullySymmetric(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);

        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        //[Test]
        public static void MovingDropletTest_rel_p3_Saye_FullySymmetric(
            [Values(2, 3)] int spatialDimension,
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
            MovingDropletTest(spatialDimension, deg, AgglomerationTreshold, performsolve, stm, Radius, vmode, bSteady, includeConvection, CutCellQuadratureType);
        }
#endif


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.MovingDropletTest"/>
        /// </summary>
        public static void MovingDropletTest(
            int spatialDimension,
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

            var Tst = new MovingDropletTest(Radius, includeConvection, bSteady, spatialDimension);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, SurfTensionMode: stm,
                CutCellQuadratureType: CutCellQuadratureType, GridResolution: (spatialDimension == 2) ? 2 : 1);
            C.SkipSolveAndEvaluateResidual = !performsolve;

            GenericTest(Tst, C);
        }


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ChannelTest"/>
        /// </summary>
        //[Test]
        public static void ChannelTest(
#if DEBUG
            [Values(2)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
#else
            [Values(2, 3, 4)] int deg,
            [Values(0.0)] double AgglomerationTreshold,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(0.0, 60.0 * Math.PI / 180.0)] double angle,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
#endif
            ) {

            var Tst = new ChannelTest(angle);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local);

            GenericTest(Tst, C);
            //if(deg == 2)
            //    ScalingTest(Tst, new[] { 2, 3, 4 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local);
            //if(deg == 3)
            //    ScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local);
        }


        /// <summary>
        /// <see cref="TranspiratingChannelTest"/>
        /// </summary>
        //[Test]
        public static void TranspiratingChannelTest(
            [Values(2, 3)] int spatialDimension,
            [Values(2)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(0.0, 0.1)] double U2,
            [Values(ViscosityMode.Standard)] ViscosityMode vmode,
            [Values(true, false)] bool periodicity,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType
            ) {

            if(spatialDimension == 3 && periodicity) {
                Console.WriteLine($"Reminder: check constained DG projection in 3D for periodic BC");
                return;
            }

            var Tst = new TranspiratingChannelTest(U2, periodicity, spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local); // surface tension plays no role in this test, so ignore it

            //C.SkipSolveAndEvaluateResidual = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.LinearSolver.MaxSolverIterations = 100;
            //C.Solver_MaxIterations = 100;
            GenericTest(Tst, C);
            //if(AgglomerationTreshold > 0) {
            //    ScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg);
            //}
        }

        /// <summary>
        /// <see cref="Tests.PolynomialTestForConvection"/>
        /// </summary>
        //[Test]
        public static void PolynomialTestForConvectionTest(
            [Values(2, 3)] int spatialDimension,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm
            ) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new PolynomialTestForConvection(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm);
            //if (spatialDimension == 3)
            //    C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            GenericTest(Tst, C);
        }

        private static void GenericTest(IXNSETest Tst, XNSE_Control C) {
            if(Tst.SpatialDimension == 3) {
                Console.WriteLine($"Reminder: skipping 3D test for now...");
                return;
                //if (C.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.Saye) {
                //    Console.WriteLine($"Reminder: skipping 3D test of {C.CutCellQuadratureType} for now...");
                //    return;
                //}
                //if (C.CutCellQuadratureType == XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes) {
                //    Console.WriteLine($"Reminder: {C.CutCellQuadratureType} changed to classic for 3D test.");
                //    C.CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Classic;
                //}
            }


            using(var solver = new XNSE_SolverMain()) {

                var ErrorPoPo = new L2ErrorLogger();
                C.PostprocessingModules.Add(ErrorPoPo);


                //Console.WriteLine("Warning! - enabled immediate plotting");
                //C.ImmediatePlotPeriod = 1;
                //C.SuperSampling = 4;

                solver.Init(C);
                solver.RunSolverMode();

                // matrix analysis
                // ===============

                solver.SpatialOperatorMatrixAnalysis(true, !C.SkipSolveAndEvaluateResidual ? 1 : 0);

                // check residuals and errors
                // ==========================

                double[] LastErrors = ErrorPoPo.ComputeL2Error(Tst.steady ? 0.0 : Tst.dt);
                double[] ErrThresh = Tst.AcceptableL2Error;
                if(LastErrors.Length != ErrThresh.Length)
                    throw new ApplicationException();
                for(int i = 0; i < ErrThresh.Length; i++) {
                    Console.WriteLine("L2 error, '{0}': \t{1}", solver.CurrentSolution.Mapping.Fields[i].Identification, LastErrors[i]);
                }

                double[] ResThresh = Tst.AcceptableResidual;
                double[] ResNorms = new double[ResThresh.Length];
                if(solver.CurrentResidual.Mapping.Fields.Count != ResThresh.Length)
                    throw new ApplicationException();
                for(int i = 0; i < ResNorms.Length; i++) {
                    ResNorms[i] = solver.CurrentResidual.Mapping.Fields[i].L2Norm();
                    Console.WriteLine("L2 norm, '{0}': \t{1}", solver.CurrentResidual.Mapping.Fields[i].Identification, ResNorms[i]);
                }

                for(int i = 0; i < ErrThresh.Length; i++)
                    Assert.LessOrEqual(LastErrors[i], ErrThresh[i]);

                for(int i = 0; i < ResNorms.Length; i++)
                    Assert.LessOrEqual(ResNorms[i], ResThresh[i]);
            }
        }

        /*
         * For the legacy solver, we don't want to run the expensive Scaling tests
         * 

        private static void ScalingTest(IXNSETest Tst, int[] ResolutionS, ViscosityMode vmode, int deg, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType, SurfaceStressTensor_IsotropicMode SurfTensionMode) {
            Console.WriteLine("Remenber do deactivate");
#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.1;

            var LaLa = new List<XNSE_Control>();
            foreach(var Res in ResolutionS) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode, GridResolution: Res);
                C.SkipSolveAndEvaluateResidual = false;
                LaLa.Add(C);
            }

            ConditionNumberScalingTest.Perform(LaLa, plotAndWait: true, title: Name);
#endif
        }
*/

        static XNSE_Control TstObj2CtrlObj(IXNSETest tst, int FlowSolverDegree, double AgglomerationTreshold, ViscosityMode vmode,
            XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            SurfaceStressTensor_IsotropicMode SurfTensionMode,
            int GridResolution = 1) {

            XNSE_Control C = new XNSE_Control();
            int D = tst.SpatialDimension;

            if(tst.TestImmersedBoundary)
                throw new NotSupportedException("Immersed boundary is not supported for the old XNSE solver.");

            // database setup
            // ==============

            C.DbPath = null;
            C.savetodb = false;
            C.ProjectName = "XNSE/" + tst.GetType().Name;
            C.ProjectDescription = "Test";

            // DG degree
            // =========

            //C.SetFieldOptions(FlowSolverDegree, tst.LevelsetPolynomialDegree);
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

            foreach(var spc in new[] { "A", "B" }) {
                C.ExactSolutionPressure.Add(spc, tst.GetPress(spc));
                C.ExactSolutionVelocity.Add(spc, D.ForLoop(d => tst.GetU(spc, d)));

                for(int d = 0; d < D; d++) {
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
            C.AgglomerationThreshold = AgglomerationTreshold;
            if(D == 3 && SurfTensionMode != SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine) {
                Console.WriteLine($"Reminder: {SurfTensionMode} changed to LaplaceBeltrami_ContactLine for 3D test.");
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            } else {
                C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfTensionMode;
            }
            C.CutCellQuadratureType = CutCellQuadratureType;


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

            // return
            // ======

            return C;

        }
    }

}
