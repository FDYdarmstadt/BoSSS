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
using System.IO;
using BoSSS.Foundation;
using BoSSS.Solution;
using MPI.Wrappers;
using NUnit.Framework;
using System.Diagnostics;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.Control;

namespace BoSSS.Application.XdgTimesteppingTest {

    /// <summary>
    /// Nunit testing.
    /// </summary>
    [TestFixture]
    class TestProgram {

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_MovingInterface_SingleInitLowOrder_BDF_dt02(
#if DEBUG
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF2)]
#else
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs
            ) {
            TestConvection_MovingInterface_SingleInitLowOrder(tsc, 0.2, NoOfTs);
        }
        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_MovingInterface_SingleInitLowOrder_RK_dt02(
#if DEBUG
            [Values(TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#else
            [Values(TimeSteppingScheme.RK1, TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK4, TimeSteppingScheme.RK_ImplicitEuler, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs
            ) {
            TestConvection_MovingInterface_SingleInitLowOrder(tsc, 0.2, NoOfTs);
        }
        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_MovingInterface_SingleInitLowOrder_BDF_dt023(
#if DEBUG
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF2)]
#else
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs
            ) {
            TestConvection_MovingInterface_SingleInitLowOrder(tsc, 0.23, NoOfTs);
        }
        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_MovingInterface_SingleInitLowOrder_RK_dt023(
#if DEBUG
            [Values(TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#else
            [Values(TimeSteppingScheme.RK1, TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK4, TimeSteppingScheme.RK_ImplicitEuler, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs
            ) {
            TestConvection_MovingInterface_SingleInitLowOrder(tsc, 0.23, NoOfTs);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/>
        /// as well as the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        public static void TestConvection_MovingInterface_SingleInitLowOrder(
            TimeSteppingScheme tsc,
            double TimestepSize,
            int NoOfTs
            ) {

            // set up
            // ------------------------------------------

            XdgTimesteppingTestControl ctrl = HardCodedControl.Gerade(angle: 0, degree: 0, GridResolutionFactor: 1);
            ctrl.NoOfTimesteps = NoOfTs;
            ctrl.dtFixed = TimestepSize; 
            ctrl.Endtime = ctrl.dtFixed * ctrl.NoOfTimesteps;
            ctrl.MultiStepInit = false;
            ctrl.TimeSteppingScheme = tsc;
            ctrl.InterfaceMode = InterfaceMode.MovingInterface;

            ctrl.NonLinearSolver.SolverCode = Solution.Control.NonLinearSolverCode.Newton;
            
            // run
            // ------------------------------------------

            XdgTimesteppingMain p = new XdgTimesteppingMain();
            p.Init(ctrl);
            p.RunSolverMode();

            // evaluate/check
            // ------------------------------------------
            double thres = 5.0e-13;
            double uA_Err = (double)p.QueryHandler.QueryResults["uA_Err"];
            double uB_Err = (double)p.QueryHandler.QueryResults["uB_Err"];
            double JmpL2Err = (double)p.QueryHandler.QueryResults["uJmp_Err"];
            Console.WriteLine("L2 Error of solution (A/B/jmp): {0}/{1}/{2} (threshold is {3}).", uA_Err, uB_Err, JmpL2Err, thres);
            Assert.LessOrEqual(uA_Err, thres);
            Assert.LessOrEqual(uB_Err, thres);
            Assert.LessOrEqual(JmpL2Err, thres);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.MultiInit(double, double, Action{int, double, DGField[]})"/>.
        /// </summary>
        [Test]
        public static void TestConvection_MovingInterface_MultiinitHighOrder(
#if DEBUG
            [Values(1)] int PolyOrder,
            [Values(0.2, 0.23)] double TimestepSize
#else
            [Values(1, 2, 3)] int PolyOrder,
            [Values(0.2, 0.23)] double TimestepSize
#endif
            ) {

            // set up
            // ------------------------------------------

            TimeSteppingScheme tsc;
            switch (PolyOrder) {
                case 0: tsc = TimeSteppingScheme.ImplicitEuler; break;
                case 1: tsc = TimeSteppingScheme.BDF2; break;
                case 2: tsc = TimeSteppingScheme.BDF3; break;
                case 3: tsc = TimeSteppingScheme.BDF4; break;
                default: throw new ArgumentOutOfRangeException();

            }


            XdgTimesteppingTestControl ctrl = HardCodedControl.Gerade(angle: 0, degree: PolyOrder, GridResolutionFactor: 1);
            ctrl.NoOfTimesteps = 8;
            ctrl.dtFixed = TimestepSize;
            ctrl.Endtime = ctrl.dtFixed * ctrl.NoOfTimesteps;
            ctrl.MultiStepInit = true;
            ctrl.TimeSteppingScheme = tsc;
            ctrl.InterfaceMode = InterfaceMode.MovingInterface;

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 3;

            // run
            // ------------------------------------------

            XdgTimesteppingMain p = new XdgTimesteppingMain();
            p.Init(ctrl);
            p.RunSolverMode();

            // evaluate/check
            // ------------------------------------------

            double thres = 5.0e-11;
            double uA_Err = (double)p.QueryHandler.QueryResults["uA_Err"];
            double uB_Err = (double)p.QueryHandler.QueryResults["uB_Err"];
            double JmpL2Err = (double)p.QueryHandler.QueryResults["uJmp_Err"];
            Console.WriteLine("L2 Error of solution (A/B/jmp): {0}/{1}/{2} (threshold is {3}).", uA_Err, uB_Err, JmpL2Err, thres);
            Assert.LessOrEqual(uA_Err, thres);
            Assert.LessOrEqual(uB_Err, thres);
            Assert.LessOrEqual(JmpL2Err, thres);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_Splitting_LowOrder_BDF_t02(
#if DEBUG
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF2)]
#else
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs,
            [Values(0.0)] double TimeOffest
            ) {
            TestConvection_Splitting_LowOrder(tsc, 0.2, NoOfTs, TimeOffest);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_Splitting_LowOrder_RK_t02(
#if DEBUG
            [Values(TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#else
            [Values(TimeSteppingScheme.RK1, TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK4, TimeSteppingScheme.RK_ImplicitEuler, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs,
            [Values(0.0)] double TimeOffest
            ) {
            TestConvection_Splitting_LowOrder(tsc, 0.2, NoOfTs, TimeOffest);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_Splitting_LowOrder_BDF_t023(
#if DEBUG
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.BDF2)]
#else
            [Values(TimeSteppingScheme.ExplicitEuler, TimeSteppingScheme.CrankNicolson, TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3, TimeSteppingScheme.BDF4)]
#endif
            TimeSteppingScheme tsc,

            [Values(8)] int NoOfTs,
            [Values(0.0)] double TimeOffest
            ) {
            TestConvection_Splitting_LowOrder(tsc, 0.23, NoOfTs, TimeOffest);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        [Test]
        public static void TestConvection_Splitting_LowOrder_RK_t023(
#if DEBUG
            [Values(TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#else
            [Values(TimeSteppingScheme.RK1, TimeSteppingScheme.RK1u1, TimeSteppingScheme.RK3, TimeSteppingScheme.RK4, TimeSteppingScheme.RK_ImplicitEuler, TimeSteppingScheme.RK_CrankNic, TimeSteppingScheme.RK_IMEX3)]
#endif
            TimeSteppingScheme tsc,
            [Values(8)] int NoOfTs,
            [Values(0.0)] double TimeOffest
            ) {
            TestConvection_Splitting_LowOrder(tsc, 0.23, NoOfTs, TimeOffest);
        }

        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/>
        /// as well as the <see cref="BoSSS.Solution.XdgTimestepping.XdgRKTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.SingleInit"/>.
        /// </summary>
        public static void TestConvection_Splitting_LowOrder(
            TimeSteppingScheme tsc,
            double TimestepSize,
            int NoOfTs,
            double TimeOffest
            ) {

            // set up
            // ------------------------------------------

            XdgTimesteppingTestControl ctrl = HardCodedControl.Gerade(angle: 0, degree: 0, GridResolutionFactor: 1);
            ctrl.NoOfTimesteps = NoOfTs;
            ctrl.dtFixed = TimestepSize;
            ctrl.Endtime = ctrl.dtFixed * ctrl.NoOfTimesteps;
            ctrl.MultiStepInit = false;
            ctrl.TimeSteppingScheme = tsc;
            ctrl.InterfaceMode = InterfaceMode.Splitting;

            

            // run
            // ------------------------------------------

            XdgTimesteppingMain p = new XdgTimesteppingMain();
            p.Init(ctrl);
            p.RunSolverMode();

            // evaluate/check
            // ------------------------------------------

            double thres = 5.0e-13;
            double uA_Err = (double)p.QueryHandler.QueryResults["uA_Err"];
            double uB_Err = (double)p.QueryHandler.QueryResults["uB_Err"];
            double JmpL2Err = (double)p.QueryHandler.QueryResults["uJmp_Err"];
            Console.WriteLine("L2 Error of solution (A/B/jmp): {0}/{1}/{2} (threshold is {3}).", uA_Err, uB_Err, JmpL2Err, thres);
            Assert.LessOrEqual(uA_Err, thres);

            double uB_Min = (double)p.QueryHandler.QueryResults["uB_Min"];
            double uB_Max = (double)p.QueryHandler.QueryResults["uB_Max"];

            Console.WriteLine("Min/Max of uB: {0} / {1}", uB_Min, uB_Max);

            Assert.GreaterOrEqual(uB_Min, ctrl.uA_Ex(new double[2], 0)*0.99999);
            Assert.LessOrEqual(uB_Max, ctrl.uB_Ex(new double[2], 0) * 1.00001);
            //Assert.LessOrEqual(uB_Err, thres);
            //Assert.LessOrEqual(JmpL2Err, thres);
        }


        /// <summary>
        /// Tests the <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping"/> time-stepper at 
        /// polynomial order 0 with single-value init, see <see cref="BoSSS.Solution.XdgTimestepping.XdgBDFTimestepping.MultiInit(double, double, Action{int, double, DGField[]})"/>.
        /// </summary>
        [Test, Sequential]
        public static void TestBurgers_HighOrder(
            [Values(0, 1, 2, 3, 0, 1)] int PolyOrder,
            [Values(0.08, 0.08, 0.08, 0.08, 0.08, 0.08)] double TimestepSize,
            [Values("bdf", "bdf", "bdf", "bdf", "rk", "rk")] string Timestepper,
            [Values(8, 8, 8, 8, 8, 8)] int NoOfTs
            ) {

            // set up
            // ------------------------------------------

            XdgTimesteppingTestControl ctrl = HardCodedControl.Burgers(angle: 0, degree: PolyOrder, GridResolutionFactor: 1, tsm: InterfaceMode.MovingInterface);
            if (Timestepper == "bdf") {
                switch (PolyOrder) {
                    case 0: ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler; break;
                    case 1: ctrl.TimeSteppingScheme = TimeSteppingScheme.CrankNicolson; break;
                    case 2: ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF3; break;
                    case 3: ctrl.TimeSteppingScheme = TimeSteppingScheme.BDF4; break;
                    default: throw new ArgumentOutOfRangeException();
                }
                ctrl.MultiStepInit = true;
            } else if (Timestepper == "rk") {
                switch (PolyOrder) {
                    case 0: ctrl.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler; break;
                    case 1: ctrl.TimeSteppingScheme = TimeSteppingScheme.CrankNicolson; break;
                    //case 2: ctrl.TimeSteppingScheme = TimeSteppingScheme.RK4; break;
                    //case 3: ctrl.TimeSteppingScheme = TimeSteppingScheme.RK4; break;
                    default: throw new ArgumentOutOfRangeException();
                }
                ctrl.MultiStepInit = false;
            } else {
                throw new ArgumentOutOfRangeException();
            }

            Console.WriteLine("Polyorder = {0}, timestepper = {1}", PolyOrder, ctrl.TimeSteppingScheme);

            ctrl.NoOfTimesteps = NoOfTs;
            ctrl.dtFixed = TimestepSize;
            ctrl.Endtime = ctrl.dtFixed * ctrl.NoOfTimesteps;

            ctrl.NonLinearSolver.SolverCode = Solution.Control.NonLinearSolverCode.Picard;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1.0e-8; // note: strangely, the test fails for **lower** thresholds; something is weired here.

            // run
            // ------------------------------------------

            XdgTimesteppingMain p = new XdgTimesteppingMain();
            p.Init(ctrl);
            p.RunSolverMode();

            // evaluate/check
            // ------------------------------------------

            double thres = 5.0e-7;
            double uA_Err = (double)p.QueryHandler.QueryResults["uA_Err"];
            double uB_Err = (double)p.QueryHandler.QueryResults["uB_Err"];
            double JmpL2Err = (double)p.QueryHandler.QueryResults["uJmp_Err"];
            Console.WriteLine("L2 Error of solution (A/B/jmp): {0}/{1}/{2} (threshold is {3}).", uA_Err, uB_Err, JmpL2Err, thres);
            Assert.LessOrEqual(uA_Err, thres);
            Assert.LessOrEqual(uB_Err, thres);
            Assert.LessOrEqual(JmpL2Err, thres);
        }

        /// <summary>
        /// A Resetting of the time-stepper is required in order for doing temporal refinement/adaptive timestepping
        /// </summary>
        [Test]
        public static void TestTimestepperReset() {
            
            // the following is one setup which has shown to be critical w.r.t. Tracker Update, etc., 
            // when the time-stepper is resetted:
            var C = BoSSS.Application.XdgTimesteppingTest.HardCodedControl2.Rarefaction(degree: 1, GridResolutionFactor: 4);
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.NoOfTimesteps = 1;
            C.dtFixed = 0.0005;
            C.Endtime = C.dtFixed * C.NoOfTimesteps;
            C.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();
            C.CutCellQuadratureType = Foundation.XDG.XQuadFactoryHelper.MomentFittingVariants.Saye;


            using(XdgTimesteppingMain p = new XdgTimesteppingMain()) {
                p.Init(C);
                p.RunSolverMode();

                p.RunFromReset(C.Endtime, C.NoOfTimesteps);
            }
        }

    }
}
