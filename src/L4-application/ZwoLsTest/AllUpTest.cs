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

namespace BoSSS.Application.ZwoLsTest {

    [TestFixture]
    static public class AllUpTest {

        [TestFixtureSetUp]
        public static void SetUp() {
            bool MpiInit;
            Environment.Bootstrap(
                new string[0],
                BoSSS.Solution.Application.GetBoSSSInstallDir(),
                out MpiInit);
        }

        [TestFixtureTearDown]
        public static void Teardown() {
            csMPI.Raw.mpiFinalize();
        }

        [Test]
        static public void AllUp([Values(0.0, 0.3)] double AggTresh, [Values(1, 2, 3)] int DGdegree, [Values(false, true)] bool DynamicBalance) {
            ZwoLsTestMain p = null;
            if(AggTresh <= 0.001 && DGdegree > 1)
                // this combination is not supposed to work
                return;

            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                null,
                delegate() {
                    p = new ZwoLsTestMain();
                    p.THRESHOLD = AggTresh;
                    p.DEGREE = DGdegree;
                    p.DYNAMIC_BALANCE = DynamicBalance;
                    return p;
                });
        }
    }
}
