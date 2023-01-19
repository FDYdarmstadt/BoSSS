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
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution;
using CutCellQuadrature.TestCases;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;

namespace CutCellQuadrature {

    [TestFixture]
    public partial class Program : Application {

        [Test]
        public static void Test3DSphereVolume() {
            ITestCase testCase = new SphereVolume3DTestCase(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.SayeGaussRules,
                    7,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-2,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test3DSphereSurface() {
            ITestCase testCase = new ConstantIntgreandSphereSurfaceIntegral3DTestCase(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.SayeGaussRules,
                    7,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-2,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }



        //Due to the sharp kinks of cube shape, the analytical results are not accurate.
        //
        //[Test]
        //public static void Test3DCubeSurface() {
        //    ITestCase testCase = new SingleCubeCubeSurfaceTestCase(GridSizes.Tiny, GridTypes.Structured);
        //    testCase.ScaleShifts(0.5 * testCase.GridSpacing);

        //    Program app = new Program(testCase);
        //    app.Init(null);
        //    app.SetUpEnvironment();
        //    app.SetInitial(0);

        //    int i = 1;
        //    while (testCase.ProceedToNextShift()) {
        //        double referenceValue = app.SetUpConfiguration();
        //        var result = app.PerformConfiguration(
        //            Modes.SayeGaussRules,
        //            7,
        //            rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));

        //        double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;
        //        //Console.WriteLine($"result.Item1: {result.Item1} , referenceValue: {referenceValue} , relError: {relError}");
        //        Assert.That(
        //            relError < 1e-2,
        //            "Relative error too large for shift number " + i);
        //        i++;
        //    }
        //}

        //[Test]
        //public static void Test3DCubeVolume() {
        //    Console.WriteLine("Hello from Test3DCubeVolume...");
        //    ITestCase testCase = new SingleCubeVolumeTestCase(GridSizes.Tiny, GridTypes.Structured);
        //    testCase.ScaleShifts(0.5 * testCase.GridSpacing);

        //    Program app = new Program(testCase);
        //    app.Init(null);
        //    app.SetUpEnvironment();
        //    app.SetInitial(0);
        //    int i = 1;
        //    while (testCase.ProceedToNextShift()) {
        //        double referenceValue = app.SetUpConfiguration();
        //        var result = app.PerformConfiguration(
        //            Modes.SayeGaussRules,
        //            7,
        //            rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
        //        double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;
        //        Console.WriteLine($"result.Item1: {result.Item1} , referenceValue: {referenceValue} , relError: {relError} , shitfI: {i}");
        //        Assert.That(relError < 1e-2, "Relative error too large for shift number " + i);
        //        i++;
        //    }
        //}

    }
}