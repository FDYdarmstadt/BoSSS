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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.AdvancedSolvers.Testing;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace BoSSS.Application.XNSEC {

    /// <summary>
    /// An all-up NUnit test for the XNSEC application.
    /// </summary>
    [TestFixture]
    static public partial class NUnitTest {

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ViscosityJumpTest"/>
        /// </summary>
        [Test]
        public static void ViscosityJumpTest(
#if DEBUG
            [Values(2)] int spatialDimension,
            [Values(1)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode
#else
            [Values(2)] int spatialDimension,
            [Values(1, 2, 3, 4)] int deg,
            [Values(0.0, 0.1)] double AgglomerationTreshold,
            [Values(ViscosityMode.FullySymmetric)] ViscosityMode vmode,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode
#endif
            ) {
            var Tst = new ViscosityJumpTest(spatialDimension);

            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfTensionMode, constantDensity: true, GridResolution: 1);
            C.LSContiProjectionMethod = Solution.LevelSetTools.ContinuityProjectionOption.None;
            C.SkipSolveAndEvaluateResidual =/* false;//*/ C.AgglomerationThreshold <= 1e-6;

            XNSECSolverTest(Tst, C);
        }

        /// <summary>
        /// <see cref="BcTest_PressureOutlet"/>
        /// </summary>
        [Test]
        public static void BcTest_PressureOutletTest(
            [Values(2)] int spatialDimension,
            [Values(1)] int deg,
            [Values(0.1)] double AgglomerationTreshold,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.Curvature_Projected, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local)] SurfaceStressTensor_IsotropicMode SurfTensionMode,
            [Values(true, false)] bool performsolve
            ) {
            var Tst = new BcTest_PressureOutlet(spatialDimension);
#if DEBUG
            const int GridRes = 2; // resolutions 1, 3, etc. place level-set at cell center
#else
            const int GridRes = 3; // resolutions 2, 4, etc. place level-set at cell boundary
#endif
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, ViscosityMode.FullySymmetric, constantDensity: true, GridResolution: GridRes, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode);
            C.SkipSolveAndEvaluateResidual = performsolve;
            
            XNSECSolverTest(Tst, C);
            if (spatialDimension == 2) // not working?...
                ASScalingTest(Tst, new[] { 4, 8, 16 }, ViscosityMode.FullySymmetric, deg, CutCellQuadratureType, SurfTensionMode);
        }

        private static void ASScalingTest(IXNSECTest Tst, int[] ResolutionS, ViscosityMode vmode, int deg, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType, SurfaceStressTensor_IsotropicMode SurfTensionMode, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton) {
#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.0;

            var LaLa = new List<XNSEC_Control>();
            foreach (var Res in ResolutionS) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, constantDensity: true, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode, GridResolution: Res);
                C.SkipSolveAndEvaluateResidual = true;

                LaLa.Add(C);
            }
            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = Name });
#endif
        }



        private static void ASScalingTest(IXNSECTest_Heat Tst, int[] ResolutionS, ViscosityMode vmode, int deg, XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType, SurfaceStressTensor_IsotropicMode SurfTensionMode, NonLinearSolverCode nonlinsolver = NonLinearSolverCode.Newton) {
//#if !DEBUG
            string Name = "Scaling" + Tst.GetType().Name + "-" + vmode + "-p" + deg;

            double AgglomerationTreshold = 0.0;

            var LaLa = new List<XNSEC_Control>();
            foreach (var Res in ResolutionS) {
                var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode: vmode, constantDensity: true, CutCellQuadratureType: CutCellQuadratureType, SurfTensionMode: SurfTensionMode, GridResolution: Res);
                C.SkipSolveAndEvaluateResidual = true;

                LaLa.Add(C);
            }
            ConditionNumberScalingTest.Perform(LaLa, new ConditionNumberScalingTest.Config() { plot = true, title = Name });
//#endif
        }

        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.ChannelTest"/>
        /// </summary>
        [Test]
        public static void NuNit_ChannelTest(
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
            BoSSS.Solution.Application.InitMPI();

            var Tst = new ChannelTest(angle);
            int gridResolution = 1;
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, constantDensity: true, gridResolution);
            C.EnableTemperature = false;
            C.EnableMassFractions = false;
            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;
            C.NonLinearSolver.verbose = true;
            XNSECSolverTest(Tst, C);
        }

        /// <summary>
        /// <see cref="Tests.PolynomialTestForConvection"/>
        /// </summary>
        [Test]
        public static void PolynomialTestForConvectionTest(
            [Values(2)] int spatialDimension,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm
            ) {
            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            int resolution = 2;

            var Tst = new BoSSS.Application.XNSEC.FullNSEControlExamples.PolynomialTestForConvection(spatialDimension);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, constantDensity: true, resolution);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            C.EnableTemperature = false;
            C.EnableMassFractions = false;
            //C.NonLinearSolver.verbose = true;
            XNSECSolverTest(Tst, C);
        }

        /// <summary>
        /// Tests a fixed level set in a constant velocity field for each phase
        /// Mass flux is predefined (not calculated from temperature gradients)
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="AgglomerationTreshold"></param>
        /// <param name="SolverMode_performsolve"></param>
        /// <param name="CutCellQuadratureType"></param>
        /// <param name="stm"></param>
        //[Test]
        public static void PseudoTwoDimensionalTwoPhaseFlow(
            [Values(2)] int deg,
            [Values(0, 0.1)] double AgglomerationTreshold,
            [Values(false, true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(false, true)] bool differentFluids,
            [Values(false, true)] bool TopBC_PressureOutlet,
            [Values(false, true)] bool BotBC_PressureOutlet
            ) {
            BoSSS.Solution.Application.InitMPI();
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter
            int resolution =5;
            bool RecoilPressure = true;
            var Tst = new BoSSS.Application.XNSEC.FullNSEControlExamples.PseudoTwoDimensional_TwoPhaseFlow(differentFluids, false, TopBC_PressureOutlet, BotBC_PressureOutlet, _Prescribed_MassFlux:1.0, recoilPressure: RecoilPressure);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, constantDensity: true, resolution);
            C.IncludeRecoilPressure = RecoilPressure;
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            C.ImmediatePlotPeriod = 1;
            C.rhoOne = true;
            C.ThermalParameters.T_sat = 1;
            C.PlotAdditionalParameters = false;
            XNSECSolverTest(Tst, C);
             }

        /// <summary>
        /// Tests a fixed level set in a constant velocity field for each phase
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="AgglomerationTreshold"></param>
        /// <param name="SolverMode_performsolve"></param>
        /// <param name="CutCellQuadratureType"></param>
        /// <param name="stm"></param>
        //[Test]
        public static void PseudoTwoDimensionalTwoPhaseFlow_ScalingTest(
            [Values(2)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(false, true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(true)] bool differentFluids,
            [Values(false, true)] bool RightBC_PressureOutlet
            ) {
            BoSSS.Solution.Application.InitMPI();
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter
            var Tst = new BoSSS.Application.XNSEC.FullNSEControlExamples.PseudoTwoDimensional_TwoPhaseFlow(differentFluids, false, RightBC_PressureOutlet, true);
            ASScalingTest(Tst, new[] { 8, 16, 32, 64 }, vmode, deg, CutCellQuadratureType, stm);
        }

        /// <summary>
        /// Tests a fixed level set in a constant velocity field
        /// </summary>
        /// <param name="deg"></param>
        /// <param name="AgglomerationTreshold"></param>
        /// <param name="SolverMode_performsolve"></param>
        /// <param name="CutCellQuadratureType"></param>
        /// <param name="stm"></param>
        //[Test]
        public static void PseudoTwoDimensionalTwoPhaseFlow_withviscosity(
            [Values(2)] int deg,
            [Values(0, 0.1)] double AgglomerationTreshold,
            [Values(false, true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm,
            [Values(false, true)] bool differentFluids,
                        [Values(false, true)] bool RightBC_PressureOutlet

            ) {
            BoSSS.Solution.Application.InitMPI();
            ViscosityMode vmode = ViscosityMode.FullySymmetric;
            int resolution = 5;
            var Tst = new BoSSS.Application.XNSEC.FullNSEControlExamples.PseudoTwoDimensional_TwoPhaseFlow(differentFluids, ViscosityActive: true, RightBC_PressureOutlet, LeftPressureOutlet :true);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, constantDensity: true, resolution);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;

            XNSECSolverTest(Tst, C);

            //ASScalingTest(Tst, new[] { 8, 16, 32,64 }, ViscosityMode.FullySymmetric, deg, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux);
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
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm) {

            ViscosityMode vmode = ViscosityMode.Standard; // viscosity is 0.0 => this selection does not matter

            var Tst = new HeatConductivityTest();
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, true, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;

            XNSECSolverTest(Tst, C);
        }




        /// <summary>
        /// Simple Test for Evaporation of a straight interface
        /// </summary>
        [Test]
        public static void SteadyStateEvaporationTestXNSEC(
            [Values(0.0, 15.0, 45.0, 73.1264)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux)] SurfaceStressTensor_IsotropicMode stm
            ) {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTestXNSEC(rawangle * Math.PI / 180.0, false);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, true, 3);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            //C.AgglomerationThreshold = 0.1;
           
            XNSECSolverTest(Tst, C);
        }


        /// <summary>
        /// Simple Test for Evaporation of a straight interface with two chemical species
        /// </summary>
        //[Test]
        public static void SteadyStateEvaporationTest_TwoSpecies_XNSEC(
            [Values(0.0, 15.0, 45.0, 73.1264, 90.0)] double rawangle,
            [Values(3)] int deg,
            [Values(0)] double AgglomerationTreshold,
            [Values(true)] bool SolverMode_performsolve,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType,
            [Values(SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Flux, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine, SurfaceStressTensor_IsotropicMode.Curvature_Projected)] SurfaceStressTensor_IsotropicMode stm
            ) {
            ViscosityMode vmode = ViscosityMode.FullySymmetric; // viscosity is 0.0 => this selection does not matter

            var Tst = new SteadyStateEvaporationTestXNSEC_TwoChemicalComponents(rawangle * Math.PI / 180.0);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, stm, true, 17);
            C.SkipSolveAndEvaluateResidual = !SolverMode_performsolve;
            //C.ImmediatePlotPeriod = 1;
            C.PlotAdditionalParameters = false;
            //C.AgglomerationThreshold = 0.1;

            XNSECSolverTest(Tst, C);
        }



 


        /// <summary>
        /// <see cref="BoSSS.Application.XNSE_Solver.Tests.LevelSetAdvectionTest"/>
        /// </summary>

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
            [Values(NonLinearSolverCode.Newton)] NonLinearSolverCode nonlinsolver
            ) {
            var Tst = new TranspiratingChannelTestXNSEC(U2, periodicity);
            var C = TstObj2CtrlObj(Tst, deg, AgglomerationTreshold, vmode, CutCellQuadratureType, SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_Local, constantDensity: true, GridResolution: 1);
            //C.SkipSolveAndEvaluateResidual = true;
            C.NonLinearSolver.MaxSolverIterations = 100;
            C.EnableTemperature = false;
            C.EnableMassFractions = false;

            C.SkipSolveAndEvaluateResidual = false;
            C.rhoOne = true;
            //C.savetodb = true;
            //C.DbPath = @"C:\Databases\BoSSS_DB";

            //C.Solver_MaxIterations = 100;
            XNSECSolverTest(Tst, C);
            //if(AgglomerationTreshold > 0) {
            //    ScalingTest(Tst, new[] { 1, 2, 3 }, vmode, deg);
            //}
        }

        public class TranspiratingChannelTestXNSEC : TranspiratingChannelTest, IXNSECTest {

            public TranspiratingChannelTestXNSEC(double _U2, bool periodic = false, int spatDim = 2) : base(_U2, periodic, spatDim) {
            }

            public bool EnableMassFractions => false;
            public bool EnableTemperature => false;
            public int NumberOfChemicalComponents => 1;

            public bool ChemicalReactionTermsActive => false;

            public double[] GravityDirection => new double[] { 0, 0, 0 };

            /// <summary>
            ///
            /// </summary>
            public new double[] AcceptableL2Error {
                get {
                    return (base.SpatialDimension == 2) ?
                        new double[] { 5.0e-2, 5.0e-2, 5.0e-1, 5.0e-2, 5.0e-2 } : new double[] { 5.0e-2, 5.0e-2, 5.0e-2, 5.0e-1, 5.0e-2, 5.0e-2 };
                }
            }

            /// <summary>
            ///
            /// </summary>
            public new double[] AcceptableResidual {
                get {
                    double[] Resi = (SpatialDimension == 2) ?
                        new double[] { 5e-8, 5e-8, 5e-8, 5e-8, 5e-8 } : new double[] { 5e-8, 5e-8, 5e-8, 5e-8, 5e-8, 5e-8 };

                    if (base.m_periodic)
                        Resi.ScaleV(5);
                    return Resi;
                }
            }

            public Func<double[], double, double> GetMassFractions(string species, int q) {
                switch (base.SpatialDimension) {
                    case 2:
                        if (q == 0) {
                            return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();
                        } else {
                            throw new ArgumentOutOfRangeException();
                        }
                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }

            public Func<double[], double, double> GetTemperature(string species) {
                switch (base.SpatialDimension) {
                    case 2:
                        return ((_3D)((t, x, y) => 1.0)).Convert_txy2Xt();

                    default:
                        throw new ArgumentOutOfRangeException();
                }
            }
        }

        /// <summary>
        /// a periodic channel flow, with water on bottom and air on top
        /// </summary>
        public class ChannelTest : IXNSECTest {
            public bool TestImmersedBoundary => false;

            /// <summary>
            /// nix
            /// </summary>
            public Func<double[], double, double> GetPhi2() {
                throw new NotImplementedException(); // will never be called, as long as 'TestImmersedBoundary' == false;
            }

            public Func<double[], double, double> GetPhi2U(int d) {
                throw new NotImplementedException();
            }

            public bool Material {
                get {
                    return true;
                }
            }

            public bool steady {
                get {
                    return true;
                }
            }

            public bool EnableMassFractions => false;
            public bool EnableTemperature => false;

            public bool IncludeConvection {
                get {
                    return true;
                }
            }

            /// <summary>
            /// the zero-level-set is identical to the x-axis
            /// </summary>
            /// <param name="time"></param>
            /// <returns></returns>
            public Func<double[], double, double> GetPhi() {
                return delegate (double[] X, double t) {
                    var Coord = ROTinv.Transform(X[0], X[1]);
                    double x = Coord[0];
                    double y = Coord[1];

                    return -y;
                };
            }

            public int LevelsetPolynomialDegree {
                get {
                    return 1;
                }
            }

            public ChannelTest(double angle) {
                //double angle = 0.0;
                //double angle = 60.0 * Math.PI / 180.0;
                ROT = AffineTrafo.Some2DRotation(angle);
                ROTinv = ROT.Invert();
            }

            private AffineTrafo ROT;

            private AffineTrafo ROTinv;

            private bool periodic = true;

            public Func<double[], double> GetF(string species, int d) {
                double rho = double.NaN;
                switch (species) {
                    case "A": rho = rho_A; break;
                    case "B":
                        rho = rho_B; break;
                        throw new ArgumentException();
                }

                if (mu_A == 0.0 && mu_B == 0) {
                    return (X => 0.0);
                } else {
                    double sc = Math.Min(this.mu_A, this.mu_B);
                    double[] Fvec = new double[] { (1.0) * sc, 0 };
                    var FvecT = ROT.Transform(Fvec);

                    return (X => FvecT[d] / rho);
                }
            }

            public Func<double[], double, double> GetU(string species, int d) {
                double a2, a1, a0;

                if (species == "A") {
                    ParabolaCoeffs_A(out a2, out a1, out a0);
                } else if (species == "B") {
                    ParabolaCoeffs_B(out a2, out a1, out a0);
                } else
                    throw new ArgumentException();

                return ((_2D)(delegate (double _x, double _y) {
                    var Coord = ROTinv.Transform(_x, _y);
                    double y = Coord[1];

                    Debug.Assert(Coord[0] >= -2);
                    Debug.Assert(Coord[0] <= +2 * 4);
                    Debug.Assert(Coord[1] >= -1);
                    Debug.Assert(Coord[1] <= +1);

                    double u = (a0 + a1 * y + a2 * y * y);
                    //double u = 1 - y*y;
                    var UT = ROT.Transform(u, 0.0);

                    return UT[d];
                })).Convert_xy2X().Convert_X2Xt();
            }

            private void ParabolaCoeffs_B(out double a2, out double a1, out double a0) {
                double muA = this.mu_A;
                double muB = this.mu_B;

                if (muA <= 0.0 && muB <= 0.0) {
                    muA = 0.001;
                    muB = 0.001;
                }
                if (muA <= 0.0 != muB <= 0.0)
                    throw new NotSupportedException();

                double sc = Math.Min(muA, muB);

                a0 = 1.0 / (muA + muB);
                a1 = (muB - muA) / (2.0 * (muB + muA) * muB);
                a2 = (-1.0) / (2 * muB);

                a0 *= sc;
                a1 *= sc;
                a2 *= sc;
            }

            private void ParabolaCoeffs_A(out double a2, out double a1, out double a0) {
                double muA = this.mu_A;
                double muB = this.mu_B;

                if (muA <= 0.0 && muB <= 0.0) {
                    muA = 0.001;
                    muB = 0.001;
                }
                if (muA <= 0.0 != muB <= 0.0)
                    throw new NotSupportedException();

                double sc = Math.Min(muA, muB);

                a0 = 1.0 / (muA + muB);
                a1 = (muB - muA) / (2.0 * (muB + muA) * muA);
                a2 = (-1.0) / (2 * muA);

                a0 *= sc;
                a1 *= sc;
                a2 *= sc;
            }

            public double dt {
                get {
                    return 1.0;
                }
            }

            public GridCommons CreateGrid(int Resolution) {
                if (Resolution < 1)
                    throw new ArgumentException();

                //var yNodes = GenericBlas.Linspace(-1, 1, 9);

                //var yNodes1 = yNodes.GetSubVector(0, 4);
                //var yNodes2 = yNodes.GetSubVector(5, 4);
                //var _yNodes = ArrayTools.Cat(yNodes1, yNodes2);
                //var _yNodes = yNodes;

                var __yNodes = new double[] { -1, -0.8, -0.6, 0.6, 1.0 };
                var _yNodes = new double[(__yNodes.Length - 1) * Resolution + 1];
                for (int i = 0; i < __yNodes.Length - 1; i++) {
                    double[] part = GenericBlas.Linspace(__yNodes[i], __yNodes[i + 1], Resolution + 1);
                    _yNodes.SetSubVector(part, i * Resolution, part.Length);
                }

                //var _yNodes = GenericBlas.Linspace(-1, 1, 6);

                var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-2, 2 * 4, 3 * Resolution + 1), _yNodes, periodicX: periodic);
                if (periodic) {
                    grd.EdgeTagNames.Add(1, "wall_top");
                    grd.EdgeTagNames.Add(2, "wall_bottom");

                    grd.DefineEdgeTags(delegate (double[] _X) {
                        var X = _X;
                        double x = X[0];
                        double y = X[1];

                        if (Math.Abs(y - (-1)) < 1.0e-6)
                            // bottom wall
                            return 2;

                        if (Math.Abs(y - (+1)) < 1.0e-6)
                            // top wall
                            return 1;

                        throw new ArgumentOutOfRangeException();
                        //return 1;
                    });

                    Console.WriteLine("ChannelTest, periodic.");
                } else {
                    grd.EdgeTagNames.Add(1, "wall_top");
                    grd.EdgeTagNames.Add(2, "wall_bottom");
                    grd.EdgeTagNames.Add(3, "Velocity_Inlet");
                    grd.EdgeTagNames.Add(4, "Pressure_Outlet");

                    grd.DefineEdgeTags(delegate (double[] _X) {
                        var X = _X;
                        double x = X[0];
                        double y = X[1];

                        if (Math.Abs(y - (-1)) < 1.0e-6)
                            // bottom wall
                            return 2;

                        if (Math.Abs(y - (+1)) < 1.0e-6)
                            // top wall
                            return 1;

                        if (Math.Abs(x - (-2)) < 1.0e-6)
                            // inlet
                            return 3;

                        if (Math.Abs(x - (2 * 4)) < 1.0e-6)
                            // outlet
                            return 4;

                        throw new ArgumentOutOfRangeException();
                        //return 1;
                    });
                }
                var grdT = grd.Transform(ROT);

                return grdT;
            }

            public IDictionary<string, AppControl.BoundaryValueCollection> GetBoundaryConfig() {
                var config = new Dictionary<string, AppControl.BoundaryValueCollection>();
                config.Add("wall_top", new AppControl.BoundaryValueCollection());
                config.Add("wall_bottom", new AppControl.BoundaryValueCollection());

                config["wall_top"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                config["wall_top"].Evaluators.Add(VariableNames.Temperature + "#B", (X, t) => 1.0);
                config["wall_top"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                config["wall_top"].Evaluators.Add(VariableNames.MassFraction0 + "#B", (X, t) => 1.0);

                config["wall_bottom"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                config["wall_bottom"].Evaluators.Add(VariableNames.Temperature + "#B", (X, t) => 1.0);
                config["wall_bottom"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                config["wall_bottom"].Evaluators.Add(VariableNames.MassFraction0 + "#B", (X, t) => 1.0);
                if (!periodic) {
                    if (!this.ROT.ApproximateEquals(AffineTrafo.Some2DRotation(0.0)))
                        throw new NotSupportedException();

                    double A_a0, A_a1, A_a2, B_a0, B_a1, B_a2;
                    this.ParabolaCoeffs_A(out A_a2, out A_a1, out A_a0);
                    this.ParabolaCoeffs_B(out B_a2, out B_a1, out B_a0);

                    config.Add("velocity_inlet", new AppControl.BoundaryValueCollection());
                    config["velocity_inlet"].Evaluators.Add(
                        VariableNames.Velocity_d(0) + "#A",
                        (X, t) => A_a0 + A_a1 * X[1] + A_a2 * X[1] * X[1]);
                    config["velocity_inlet"].Evaluators.Add(
                        VariableNames.Velocity_d(0) + "#B",
                        (X, t) => B_a0 + B_a1 * X[1] + B_a2 * X[1] * X[1]);
                    config["velocity_inlet"].Evaluators.Add(VariableNames.Temperature + "#A", (X, t) => 1.0);
                    config["velocity_inlet"].Evaluators.Add(VariableNames.Temperature + "#B", (X, t) => 1.0);
                    config["velocity_inlet"].Evaluators.Add(VariableNames.MassFraction0 + "#A", (X, t) => 1.0);
                    config["velocity_inlet"].Evaluators.Add(VariableNames.MassFraction0 + "#B", (X, t) => 1.0);

                    config.Add("Pressure_Outlet", new AppControl.BoundaryValueCollection());
                }

                return config;
            }

            public Func<double[], double, double> GetPress(string species) {
                return (X, t) => 0.0;
            }

            public Func<double[], double, double> GetTemperature(string species) {
                return (X, t) => 1.0;
            }

            public Func<double[], double, double> GetMassFractions(string species, int comp) {
                return (X, t) => 1.0;
            }

            /// <summary>
            /// specific weight, air
            /// </summary>
            public double rho_B {
                get {
                    return 1.2;
                }
            }

            /// <summary>
            /// specific weight, water
            /// </summary>
            public double rho_A {
                get {
                    return 1000;
                }
            }

            /// <summary>
            /// dynamic viscosity, air
            /// </summary>
            public double mu_B {
                get {
                    return 17.1e-3;
                }
            }

            /// <summary>
            /// dynamic viscosity, water
            /// </summary>
            public double mu_A {
                get {
                    return 1.0;
                }
            }

            /// <summary>
            /// surface tension of water (surface tension has no effect due to the planar interface)
            /// </summary>
            public double Sigma {
                get {
                    //return 0.0;
                    return 72.75e-3;
                }
            }

            public double[] AcceptableL2Error {
                get {
                    return new double[] { 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7 };
                }
            }

            public double[] AcceptableResidual {
                get {
                    return new double[] { 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7, 1.0e-7 };
                }
            }

            public int SpatialDimension {
                get {
                    return 2;
                }
            }

            public int NumberOfChemicalComponents => 1;

            public bool ChemicalReactionTermsActive => false;

            public double[] GravityDirection => new double[] { 0, 0, 0 };
        }
    }
}