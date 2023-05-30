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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;
using System;

namespace IntersectingLevelSetTest {

    [TestFixture]
    static public class AllUpTest {

        [OneTimeSetUp]
        static public void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }

        [Test]
        // Test two LS-line in one single cell. For same shape, but two different combinations of line. 
        // For first time step, the error magnitude is E-15, it's straight horizontal line + curve line. 
        // For second time step, the error magnitude is E-3, it's straight vertical line + curve line. Might be something wrong. 
        public static void ParabolaTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => -(y + (0 + t) * x * x);
            Func<double, double, double, double> levelSet1 = (x, y, t) => (x - (1 - t) * y * y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 2;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }

        [Test]
        // Test two LS-line in one single cell.
        // For one straight horizontal line + the other moving straight line with slope line.
        // Set NoOfTimesteps = 5;
        // The error magnitude is E-15 for good cases, E-4 for bad cases. 
        public static void TwoStraightTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => (x - (Math.Tan((t + 382) / 90 * Math.PI) * y) + (t + 382) / 9000 * Math.PI * Math.Sin((t + 382) / 9 * Math.PI));
            Func<double, double, double, double> levelSet1 = (x, y, t) => -(y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 5;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }

        [Test]
        // Test two LS-line in 5 x 5 cells.
        // For one straight horizontal line + curve line.
        // The largest error might not caused by the central cell, but by the cells around. 
        // The error magnitude is E-15 for good cases, E-7 for bad cases. 
        public static void TransformTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => (x + 1 * y * y * 10 * Math.Sin((t + 92) / 90.1 * Math.PI));
            Func<double, double, double, double> levelSet1 = (x, y, t) => -(y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 6; //number of nodes per line
            C.NoOfTimesteps = 5;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }

        //[Test]
        //// Test two straight LS-line in 5 x 5 cells.
        //// translate + rotate => "Random" movement of lines. 
        //public static void RandomTest(
        //    [Values(1, 2, 3)] int DGdegree) {
        //    BoSSS.Solution.Application.InitMPI();
        //    //BoSSS.Solution.Application.DeleteOldPlotFiles();
        //    Func<double, double, double, double> levelSet0 = (x, y, t) => (x - 1 * (Math.Tan(t / 90.1 * Math.PI) * y) + t / 90.1 * Math.PI * 0.01 * Math.Sin(t / 90.1 * Math.PI * 10));
        //    Func<double, double, double, double> levelSet1 = (x, y, t) => (1 * (Math.Tan(t / 90.1 * Math.PI) * x) - (y) + t / 90.1 * Math.PI * 0.01 * Math.Cos(t / 90.1 * Math.PI * 10));
        //    var C = new TestControl(levelSet0, levelSet1);
        //    C.Resolution = 6; //number of nodes per line
        //    C.NoOfTimesteps = 1300;
        //    C.ErrorThreshold = 1e-6;
        //    var p = new ZwoLsSolver<TestControl>();
        //    p.DEGREE = DGdegree;
        //    p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //    p.Init(C);
        //    p.RunSolverMode();
        //}

        //[Test]
        //// Test two curve LS-lines in 5 x 5 cells.
        //// The shape of each line changes with time. 
        //public static void RealTransformTest(
        //    [Values(1, 2, 3)] int DGdegree) {
        //    BoSSS.Solution.Application.InitMPI();
        //    //BoSSS.Solution.Application.DeleteOldPlotFiles();
        //    Func<double, double, double, double> levelSet0 = (x, y, t) => (x + y * y * 10 * Math.Sin(t / 90.1 * Math.PI));
        //    Func<double, double, double, double> levelSet1 = (x, y, t) => -(y + x * x * 10 * Math.Cos(t / 90.1 * Math.PI));
        //    var C = new TestControl(levelSet0, levelSet1);
        //    C.Resolution = 6; //number of nodes per line
        //    C.NoOfTimesteps = 181;
        //    C.ErrorThreshold = 1e-6;
        //    var p = new ZwoLsSolver<TestControl>();
        //    p.DEGREE = DGdegree;
        //    p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //    p.Init(C);
        //    p.RunSolverMode();
        //}

        [Test]
        // This is the 3D test of two rotating planes in one single cell. 
        // The errors are quite low, E-15
        // However, might throw ex: 'Root not found' at timesetp 24.
        // However, if the C.Resolution = 3; calculation will be faster, and throw ex will not happened. 
        public static void Rotation3DTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => (x - 1 * (Math.Tan(t / 90.1 * Math.PI) * y) + t / 90.1 * Math.PI * 0.01 * Math.Sin(t / 90.1 * Math.PI * 10));
            Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
            var C = new TestControl();
            C.Dimension = 3;
            C.LevelSet3D_0 = levelSet3D_0;
            C.LevelSet3D_1 = levelSet3D_1;
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 50;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }


        [Test]
        // This is the 3D test of a flat horizontal plane + a curved plane in one single cell. 
        // The errors are quite low, E-15
        // However, might throw ex: 'Root not found' at timesetp 7.
        // However, if the C.Resolution = 3; calculation will be faster, and throw ex will not happened. But the error is not small after timestep 7. 
        public static void Transform3DTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => (x + 1 * y * y * 10 * Math.Sin(t / 90.1 * Math.PI));
            Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
            var C = new TestControl();
            C.Dimension = 3;
            C.LevelSet3D_0 = levelSet3D_0;
            C.LevelSet3D_1 = levelSet3D_1;
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 10;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }


        [Test]
        // This is the 3D test of a horizontal planes and a moving sphere in one single cell. 
        // The error is E-1 to E-4. 
        public static void MovingSphere3DTest(
            [Values(1, 2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => -(0.09 - (x - 0.5 + t / 90.1 * Math.PI) * (x - 0.5 + t / 90.1 * Math.PI) - y * y - z * z);
            Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
            var C = new TestControl();
            C.Dimension = 3;
            C.LevelSet3D_0 = levelSet3D_0;
            C.LevelSet3D_1 = levelSet3D_1;
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 38;
            C.ErrorThreshold = 1e-6;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }


        //[Test]
        //static public void AllUp(
        //    [Values(1, 2, 3)] int DGdegree,
        //    [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants quadVariant) {
        //    var C = new PlotControl();
        //    ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl> p = null;
        //    BoSSS.Solution.Application._Main(
        //        new string[0],
        //        true,
        //        delegate () {
        //            p = new ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl>();
        //            p.DEGREE = DGdegree;
        //            p.MomentFittingVariant = quadVariant;
        //            return p;
        //        });
        //}

        //static public void LocalTestWithPlotting(
        //    [Values(1, 2, 3)] int DGdegree,
        //    [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants quadVariant) {
        //    BoSSS.Solution.Application.InitMPI();
        //    BoSSS.Solution.Application.DeleteOldPlotFiles();
        //    var C = new PlotControl();
        //    var p = new ZwoLsSolver<PlotControl>();
        //    p.DEGREE = DGdegree;
        //    p.MomentFittingVariant = quadVariant;
        //    C.SuperSampling = 5;
        //    p.Init(C);
        //    p.RunSolverMode();
        //}
        //static public void OneCellTwoIntersectingLevelSets() {
        //    var C = new PlotControl();
        //    C.SetDGdegree(0);

        //    C.SetGrid(Grid2D.Cartesian2DGrid(new double[] { 0, 1 }, new double[] { 0, 2 }));

        //    ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl> p = null;
        //    BoSSS.Solution.Application._Main(
        //        new string[0],
        //        true,
        //        delegate () {
        //            p = new ZwoLsSolver<BoSSS.Solution.Application.EmptyAppControl>();
        //            p.DEGREE = 0;
        //            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //            return p;
        //        });
        //}
    }
}
