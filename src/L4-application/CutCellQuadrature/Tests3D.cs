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
using BoSSS.Foundation.XDG.Quadrature.HMF;
using BoSSS.Solution;
using CutCellQuadrature.TestCases;
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
                    Modes.HMFClassic,
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
                    Modes.HMFClassic,
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
        public static void Test3DSurfaceSubdivision() {
            
        }

        [Test]
        public static void Test3DVolumeSubdivision() { 
        
        }
        }
}