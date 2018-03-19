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

using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.SpecFEM {

    /// <summary>
    /// Unit test for the continuous Galerkin method.
    /// </summary>
    [TestFixture]
    public class AllUpTest {

        /// <summary>
        /// MPI Finalization.
        /// </summary>
        [TestFixtureTearDown]
        public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init.
        /// </summary>
        [TestFixtureSetUp]
        public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);
        }

        /// <summary>
        /// not the smartest way to define such a test...
        /// </summary>
        [Test]
        public void AllUp([Values(false, true)] bool perX,
                          [Values(false, true)] bool perY
            ) {

            SpecFEMMain p = null;

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new SpecFEMMain();
                p.m_periodicX = perX;
                p.m_periodicY = perY;
                return p;
            });

            Assert.IsTrue(p.Passed);
        }

    }
}
