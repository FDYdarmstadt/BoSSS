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

using System;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution;
using CutCellQuadrature.TestCases;
using MPI.Wrappers;
using NUnit.Framework;
using ilPSP.Utils;

namespace CutCellQuadrature {

    [TestFixture]
    public partial class Program : Application {

        [Test]
        public static void Test2DSurfaceHighOrderRobustnessStructured() {
            ITestCase testCase = new Smereka2EllipseArcLength(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.HMFClassic,
                    7,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-4,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DSurfaceHighOrderRobustnessStructuredSaye() {
            ITestCase testCase = new Smereka2EllipseArcLength(GridSizes.Tiny, GridTypes.Structured);
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
                    8,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-4,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DSurfaceHighOrderRobustnessStructuredAlgoim() {
            ITestCase testCase = new Smereka2EllipseArcLength(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.Algoim,
                    8);
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;
                
                Assert.That(
                    relError < 1e-4,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DSurfaceConvergenceStructured() {
            int[] orders = Enumerable.Range(0, 10).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Tiny, GridSizes.Small, GridSizes.Normal };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new Smereka2EllipseArcLength(sizes[i], GridTypes.Structured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.HMFClassic,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Assert.That(
                    eoc > orders[j] + 1,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        [Test]
        public static void Test2DSurfaceConvergenceStructuredSaye() {
            int[] orders = Enumerable.Range(0, 10).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Tiny, GridSizes.Small, GridSizes.Normal };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new Smereka2EllipseArcLength(sizes[i], GridTypes.Structured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.SayeGaussRules,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Assert.That(
                    eoc > orders[j] / 2,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        /// <summary>
        /// [Test] not activated yet
        /// </summary>
        public static void Test2DSurfaceConvergenceStructuredAlgoim() {
            int[] orders = Enumerable.Range(3, 5).ToArray();

            GridSizes[] sizes = new GridSizes[] { GridSizes.Ultra, GridSizes.Mega, GridSizes.Giga };
            double[,] results = new double[sizes.Length, orders.Length];

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new Saye2022EllipsePerimeter(sizes[i]);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.Algoim,
                        orders[j]);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Console.WriteLine(orders[j] + ": " + eoc);
                Assert.That(
                    eoc > orders[j] + 1,
                    "Convergence order (" + eoc + ") too low for order " + orders[j] );
            }
        }

        [Test]
        public static void Test2DVolumeHighOrderRobustnessStructured() {
            ITestCase testCase = new MinGibou1EllipseArea(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.HMFClassic,
                    7,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-5,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DVolumeHighOrderRobustnessStructuredSaye() {
            ITestCase testCase = new MinGibou1EllipseArea(GridSizes.Tiny, GridTypes.Structured);
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
                    relError < 1e-4,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DVolumeHighOrderRobustnessStructuredAlgoim() {
            ITestCase testCase = new MinGibou1EllipseArea(GridSizes.Tiny, GridTypes.Structured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.Algoim,
                    7);
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-5,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }

        [Test]
        public static void Test2DVolumeConvergenceStructured() {
            int[] orders = Enumerable.Range(0, 9).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Tiny, GridSizes.Small, GridSizes.Normal };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new MinGibou1EllipseArea(sizes[i], GridTypes.Structured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.HMFClassic,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Console.WriteLine(orders[j] + ": " + eoc);
                Assert.That(
                    eoc > orders[j] + 1,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        [Test]
        public static void Test2DVolumeConvergenceStructuredSaye() {
            int[] orders = Enumerable.Range(0, 9).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Tiny, GridSizes.Small, GridSizes.Normal };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new MinGibou1EllipseArea(sizes[i], GridTypes.Structured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.SayeGaussRules,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Console.WriteLine(orders[j] + ": " + eoc);
                Assert.That(
                    eoc > orders[j] - 1,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        [Test]
        public static void Test2DVolumeHConvergenceStructuredAlgoim() {
            int[] orders = Enumerable.Range(3,10).ToArray();
            //int[] orders = new int[] { 1 };

            GridSizes[] sizes = new GridSizes[] { GridSizes.Ultra, GridSizes.Mega, GridSizes.Giga };
            double[,] results = new double[sizes.Length, orders.Length];

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new Saye2022EllipseArea(sizes[i]);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.Algoim,
                        orders[j]);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Console.WriteLine(orders[j] + ": " + eoc);
                Assert.That(
                    eoc > orders[j] + 1 ,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        [Test]
        public void Test2DSurfaceHighOrderRobustnessUnstructured() {
            ITestCase testCase = new Smereka2EllipseArcLength(GridSizes.Tiny, GridTypes.PseudoStructured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.HMFClassic,
                    7,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Assert.That(
                    relError < 1e-4,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }


        [Test]
        public void Test2DSurfaceConvergenceUnstructured() {
            int[] orders = Enumerable.Range(0, 7).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Small, GridSizes.Normal, GridSizes.Large };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new Smereka2EllipseArcLength(sizes[i], GridTypes.PseudoStructured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.HMFClassic,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Assert.That(
                    eoc > orders[j] + 1,
                    "Convergence order too low for order " + orders[j]);
            }
        }

        [Test]
        public void Test2DVolumeHighOrderRobustnessUnstructured() {
            ITestCase testCase = new MinGibou1EllipseArea(GridSizes.Tiny, GridTypes.PseudoStructured);
            testCase.ScaleShifts(0.5 * testCase.GridSpacing);

            Program app = new Program(testCase);
            app.Init(null);
            app.SetUpEnvironment();
            app.SetInitial(0);

            int i = 1;
            while (testCase.ProceedToNextShift()) {
                double referenceValue = app.SetUpConfiguration();
                var result = app.PerformConfiguration(
                    Modes.HMFClassic,
                    6,
                    rootFindingAlgorithm: new LineSegment.SafeGuardedNewtonMethod(1e-14));
                double relError = Math.Abs(result.Item1 - referenceValue) / testCase.Solution;

                Console.WriteLine(relError);
                Assert.That(
                    relError < 1e-6,
                    "Relative error too large for shift number " + i);
                i++;
            }
        }


        [Test]
        public void Test2DVolumeConvergenceUnstructured() {
            int[] orders = Enumerable.Range(0, 6).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Small, GridSizes.Normal, GridSizes.Large };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new MinGibou1EllipseArea(sizes[i], GridTypes.PseudoStructured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.HMFClassic,
                        orders[j],
                        rootFindingAlgorithm: rootFindingAlgorithm);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                }
        }

        [Test]
        public void Test2DVolumeConvergenceUnstructuredAlgoim() {
            int[] orders = Enumerable.Range(0, 6).ToArray();
            GridSizes[] sizes = new GridSizes[] { GridSizes.Small, GridSizes.Normal, GridSizes.Large };
            double[,] results = new double[sizes.Length, orders.Length];
            var rootFindingAlgorithm = new LineSegment.SafeGuardedNewtonMethod(1e-14);

            for (int i = 0; i < sizes.Length; i++) {
                ITestCase testCase = new MinGibou1EllipseArea(sizes[i], GridTypes.PseudoStructured);
                testCase.ScaleShifts(0.5 * testCase.GridSpacing);

                Program app = new Program(testCase);
                app.Init(null);
                app.SetUpEnvironment();
                app.SetInitial(0);
                testCase.ProceedToNextShift();
                double referenceValue = app.SetUpConfiguration();

                for (int j = 0; j < orders.Length; j++) {
                    var result = app.PerformConfiguration(
                        Modes.Algoim,
                        orders[j]);
                    results[i, j] = Math.Abs(result.Item1 - referenceValue);
                }
            }

            double[] xValues = sizes.Select(s => -Math.Log(2.0) * (int)s).ToArray();
            for (int j = 0; j < orders.Length; j++) {
                double[] yValues = new double[sizes.Length];

                for (int i = 0; i < sizes.Length; i++) {
                    yValues[i] = Math.Log(results[i, j]);
                }

                double eoc = Regression(xValues, yValues);

                Console.WriteLine(eoc);
                Assert.That(
                    eoc > orders[j] + 1,
                    "Convergence order too low for order " + orders[j]);
            }
        }
    }
}