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
using System.Linq;

namespace IntersectingLevelSetTest {

    [TestFixture]
    static public class AllUpTest {

        [OneTimeSetUp]
        static public void SetUp() {
            BoSSS.Solution.Application.InitMPI();
        }

        [Test]
        // Test two LS-line in one single cell. For same shape, but two different combinations of line. 

        // When u.ProjectField((x, y) => x * x), which means du_dx_Exact.ProjectField((x, y) => 2 * x), and DGdegree = 2: 
        // For first time step, the error magnitude is E-15, it's straight horizontal line + curve line. 
        // For second time step, the error magnitude is E-3, it's straight vertical line + curve line. Might be something wrong. 

        // When u.ProjectField((x, y) => Math.Sin(x) * Math.Cos(y)) 
        // This issue no longer appears. 

        public static void ParabolaTest(
            [Values(3, 4, 5)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => -(y + (0 + t) * x * x);
            Func<double, double, double, double> levelSet1 = (x, y, t) => (x - (1 - t) * y * y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 2;
            //C.ErrorThreshold = 1e-6;
            C.ErrorThreshold = 5e-3;

            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            //C.SuperSampling = 5;
            p.Init(C);
            p.RunSolverMode();
        }
        
        [Test]
        // Test two LS-line in one single cell. 
        // For one straight horizontal line + the other moving straight line with slope line.
        // Set NoOfTimesteps = 5;

        // When u.ProjectField((x, y) => x * x), which means du_dx_Exact.ProjectField((x, y) => 2 * x), and DGdegree = 2: 
        // The error magnitude is E-15 for good cases, E-4 for bad cases. 

        // When u.ProjectField((x, y) => Math.Sin(x) * Math.Cos(y)) 
        // This issue no longer appears. 
        public static void TwoStraightTest(
            [Values(3, 4, 5)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => (x - (Math.Tan((t + 382) / 90 * Math.PI) * y) + (t + 382) / 9000 * Math.PI * Math.Sin((t + 382) / 9 * Math.PI));
            Func<double, double, double, double> levelSet1 = (x, y, t) => -(y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 5;
            //C.ErrorThreshold = 1e-6;
            C.ErrorThreshold = 5e-3;
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

        // When u.ProjectField((x, y) => x * x), which means du_dx_Exact.ProjectField((x, y) => 2 * x), and DGdegree = 2: 
        // The error magnitude is E-15 for good cases, E-7 for bad cases. 

        // When u.ProjectField((x, y) => Math.Sin(x) * Math.Cos(y)) 
        // This issue no longer appears. 
        public static void TransformTest(
            [Values(2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double> levelSet0 = (x, y, t) => (x + 1 * y * y * 10 * Math.Sin((t + 92) / 90.1 * Math.PI));
            Func<double, double, double, double> levelSet1 = (x, y, t) => -(y);
            var C = new TestControl(levelSet0, levelSet1);
            C.Resolution = 6; //number of nodes per line
            C.NoOfTimesteps = 5;
            //C.ErrorThreshold = 1e-6;
            C.ErrorThreshold = 5e-3;
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
        //    Func<double, double, double, double> levelSet0 = (x, y, t) =>  (x + 0 * y * y * 10 * Math.Sin(t / 90.1 * Math.PI));
        //    Func<double, double, double, double> levelSet1 = (x, y, t) => -(y + 1 * x * x * 10 * Math.Cos(t / 90.1 * Math.PI));
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
        // This is the 2D Convergence test of two straight line. 
        public static void Convergence2DTest(
            [Values(1, 2)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            int[] numbersOfCells = new int[] {2, 4, 8};
            double[] errorList = new double[numbersOfCells.Length];
            int i = 0;
            foreach (int cell in numbersOfCells) {
                Func<double, double, double, double> levelSet0 = (x, y, t) => (y + 0.3 * Math.Cos(x * 0.5 * Math.PI) - 0.4);
                Func<double, double, double, double> levelSet1 = (x, y, t) => (x - 0.3 * Math.Cos(y * 0.5 * Math.PI) + 0.4);
                var C = new TestControl(levelSet0, levelSet1);
                C.Dimension = 2;
                C.Resolution = cell + 1; //number of nodes per line
                C.NoOfTimesteps = 1;
                C.ErrorThreshold = 1e-0;
                var p = new ZwoLsSolver<TestControl>();
                p.DEGREE = DGdegree;
                p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
                p.Init(C);
                p.RunSolverMode();
                errorList[i] = p.errorBack; //get the error from solver
                i++;
            }

            //starting to calculate the slope of loglog regression: 
            numbersOfCells.Select(x => x).ToList();
            errorList.Select(x => x).ToList();

            double[] xValues = numbersOfCells.Select(x => Math.Log10(x)).ToArray();
            double[] yValues = errorList.Select(y => Math.Log10(y)).ToArray();


            double xAvg = xValues.Average();
            double yAvg = yValues.Average();

            double v1 = 0.0;
            double v2 = 0.0;

            for (int j = 0; j < yValues.Length; j++) {
                v1 += (xValues[j] - xAvg) * (yValues[j] - yAvg);
                v2 += Math.Pow(xValues[j] - xAvg, 2);
            }

            double a = - v1 / v2; //The slope

            // Output the summary
            for (int j = 0;j < yValues.Length; j++) {
                Console.WriteLine("error for " + numbersOfCells[j] + " cells: " + errorList[j]);
            }
            Console.WriteLine("Slope: " + a);

            Assert.IsTrue(Math.Abs(a - DGdegree) < 0.1, "The regression slope is not correct!");
            Console.WriteLine("Convergence2DTest for DG degree " + DGdegree + " PASSED!");
        }


        [Test]
        // This is the 3D test of two rotating planes in one single cell. 
        // The errors are quite low, E-15
        // However, might throw ex: 'Root not found' at timesetp 24.
        // However, if the C.Resolution = 3; calculation will be faster, and throw ex will not happened. 
        public static void Rotation3DTest(
            [Values(2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => (x - 1 * (Math.Tan(t / 90.1 * Math.PI) * y) + t / 90.1 * Math.PI * 0.01 * Math.Sin(t / 90.1 * Math.PI * 10));
            Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
            var C = new TestControl();
            C.Dimension = 3;
            C.LevelSet3D_0 = levelSet3D_0;
            C.LevelSet3D_1 = levelSet3D_1;
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 6;
            C.ErrorThreshold = 1e-10;
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
            [Values(2, 3)] int DGdegree) {
            BoSSS.Solution.Application.InitMPI();
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => (x + 1 * y * y * 10 * Math.Sin(t / 90.1 * Math.PI));
            Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
            var C = new TestControl();
            C.Dimension = 3;
            C.LevelSet3D_0 = levelSet3D_0;
            C.LevelSet3D_1 = levelSet3D_1;
            C.Resolution = 2; //number of nodes per line
            C.NoOfTimesteps = 6;
            C.ErrorThreshold = 1e-10;
            var p = new ZwoLsSolver<TestControl>();
            p.DEGREE = DGdegree;
            p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
            p.Init(C);
            p.RunSolverMode();
        }


        //[Test]
        //// This is the 3D test of a horizontal planes and a moving sphere in one single cell. 
        //// The error is E-1 to E-4. 
        //public static void MovingSphere3DTest(
        //    [Values(2, 3)] int DGdegree) {
        //    BoSSS.Solution.Application.InitMPI();
        //    //BoSSS.Solution.Application.DeleteOldPlotFiles();
        //    Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => -(0.09 - (x - 0.5 + t / 90.1 * Math.PI) * (x - 0.5 + t / 90.1 * Math.PI) - y * y - z * z);
        //    Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => -(y);
        //    var C = new TestControl();
        //    C.Dimension = 3;
        //    C.LevelSet3D_0 = levelSet3D_0;
        //    C.LevelSet3D_1 = levelSet3D_1;
        //    C.Resolution = 2; //number of nodes per line
        //    C.NoOfTimesteps = 38;
        //    C.ErrorThreshold = 1e-6;
        //    var p = new ZwoLsSolver<TestControl>();
        //    p.DEGREE = DGdegree;
        //    p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //    p.Init(C);
        //    p.RunSolverMode();
        //}


        //[Test]
        //// This is the 3D Convergence test of a horizontal planes and a sphere. 
        //public static void Convergence3DTest(
        //    [Values(1, 2, 3)] int DGdegree) {
        //    BoSSS.Solution.Application.InitMPI();
        //    BoSSS.Solution.Application.DeleteOldPlotFiles();
        //    for (int i = 0; i < 5; i++) {
        //        //Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => -(0.09 - x * x - y * y - z * z);
        //        Func<double, double, double, double, double> levelSet3D_0 = (x, y, z, t) => (z - 0.2 * Math.Cos(20 / 11 * Math.PI * x));
        //        Func<double, double, double, double, double> levelSet3D_1 = (x, y, z, t) => (-z);
        //        var C = new TestControl();
        //        C.Dimension = 3;
        //        C.LevelSet3D_0 = levelSet3D_0;
        //        C.LevelSet3D_1 = levelSet3D_1;
        //        C.Resolution = 3 + 2 * i; //number of nodes per line
        //        C.NoOfTimesteps = 1;
        //        C.ErrorThreshold = 1e-6;
        //        var p = new ZwoLsSolver<TestControl>();
        //        p.DEGREE = DGdegree;
        //        //p.MPISize = 2;
        //        p.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Saye;
        //        p.Init(C);
        //        p.RunSolverMode();
        //    }
        //}

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
