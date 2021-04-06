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
using ilPSP;
using NUnit.Framework;
using BoSSS.Foundation.XDG;
using MPI.Wrappers;

namespace BoSSS.Application.XdgNastyLevsetLocationTest {

    /// <summary>
    /// Nunit entry point
    /// </summary>
    [TestFixture]
    public class AllUpTest { 
        /// <summary>
        /// Level-Set is parallel resp. close-to-parallel ot a cell edge
        /// </summary>
        [Test]
        public static void ParalleTest_2D(
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)]
            XQuadFactoryHelper.MomentFittingVariants variant) {

            TestTemplate(variant, new Parallel(XdgNastyLevsetLocationTest.GetTestRange(), XdgNastyLevsetLocationTest.GetTestRange()));
        }

        /// <summary>
        /// Level-Set passes through a corner 
        /// </summary>
        [Test]
        public static void CornerTest_2D(
            [Values(XQuadFactoryHelper.MomentFittingVariants.OneStepGauss, XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes, XQuadFactoryHelper.MomentFittingVariants.Saye)]
            XQuadFactoryHelper.MomentFittingVariants variant) {

            TestTemplate(variant, new Schraeg(XdgNastyLevsetLocationTest.GetTestRange(), XdgNastyLevsetLocationTest.GetTestRange()));
        }


        static void TestTemplate(XQuadFactoryHelper.MomentFittingVariants variant, ITest tst) {
            XdgNastyLevsetLocationTest p = null;

            tst.ResetTest();

            BoSSS.Solution.Application._Main(new string[0], true, delegate () {
                p = new XdgNastyLevsetLocationTest();
                p.test = tst;
                p.momentFittingVariant = variant;
                return p;
            });

            Assert.IsTrue(p.IsPassed);
        }
    }
}

