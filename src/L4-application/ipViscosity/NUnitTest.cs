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
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using BoSSS.Solution;
using System.IO;
using BoSSS.Solution.NSECommon;
using BoSSS.Platform;
using MPI.Wrappers;
using ilPSP;

namespace BoSSS.Application.ipViscosity {

    /// <summary>
    /// An all-up NUnit test for the ipPoisson application
    /// </summary>
    [TestFixture]
    public class _Test {

        [TestFixtureTearDown]
        public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        [TestFixtureSetUp]
        public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);
        }

        static TestSolution[] solutions = new TestSolution[] { new Polynomial2D_ConstantVisc(), new Polynomial2D_VariableVisc(), new Transcendent2D() };

        [Test]
        public static void ConsistencyTest(
            [Values(Terms.T1, Terms.T2, Terms.T3)] Terms t,
            [Values(0, 1)] int iSol
            ) {
            ipViscosityMain p = null;
            BoSSS.Solution.Application._Main(new string[0], true, delegate() {
                p = new ipViscosityMain();
                p.mode = TestMode.CheckResidual;
                p.solution = solutions[iSol];
                p.PolynomialDegree = 4;
                p.grid  = new MixedBcGrid();
                p.whichTerms = t;
                return p;
            });

            p.L2ResidualNorm.ForEach(res => Assert.LessOrEqual(res, 1.0e-8, " Residual L2 Norm "));
        }


        [Test]
        public static void solverTest(
            [Values(Terms.T1, Terms.T1 | Terms.T2, Terms.T1 | Terms.T2 | Terms.T3)] Terms t,
            [Values(0, 1, 2)] int iSol,
            [Values(2,3,4)] int deg
            ) {
            ipViscosityMain p = null;
            BoSSS.Solution.Application._Main(new string[0], true, delegate() {
                p = new ipViscosityMain();
                p.mode = TestMode.Solve;
                p.solution = solutions[iSol];
                p.PolynomialDegree = deg;
                p.grid  = new MixedBcGrid();
                p.whichTerms = t;
                return p;
            });


            p.L2ResidualNorm.ForEach(res => Assert.LessOrEqual(res, 1.0e-8));
            p.L2ErrorNorm.ForEach(err => Assert.LessOrEqual(err, 0.5));
        }
    }
}
