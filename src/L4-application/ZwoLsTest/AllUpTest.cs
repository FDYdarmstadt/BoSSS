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

using BoSSS.Foundation.XDG;
using ilPSP;
using MPI.Wrappers;
using NUnit.Framework;

namespace BoSSS.Application.ZwoLsTest {

    [TestFixture]
    static public class AllUpTest {

        [OneTimeSetUp]
        static public void SetUp()
        {
            BoSSS.Solution.Application.InitMPI();
        }

        [Test]
        static public void AllUp([Values(0.0, 0.3)] double AggTresh,
#if DEBUG            
            [Values(1)] int DGdegree,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)] XQuadFactoryHelper.MomentFittingVariants quadVariant,
            [Values(false)] bool DynamicBalance)
#else
            [Values(1, 2, 3)] int DGdegree,
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)] XQuadFactoryHelper.MomentFittingVariants quadVariant,
            [Values(false, true)] bool DynamicBalance)            
#endif
         {
            ZwoLsTestMain p = null;
            if(AggTresh <= 0.001 && DGdegree > 1)
                // this combination is not supposed to work
                return;

            BoSSS.Solution.Application._Main(
                new string[0],
                true,
                delegate() {
                    p = new ZwoLsTestMain();
                    p.THRESHOLD = AggTresh;
                    p.DEGREE = DGdegree;
                    p.DYNAMIC_BALANCE = DynamicBalance;
                    p.MomentFittingVariant = quadVariant;
                    return p;
                });
        }
    }
}
