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

namespace BoSSS.Application.XdgNastyLevsetLocationTest {


    [TestFixture]
    public class AllUpTest {


        /// <summary>
        /// not the smartest way to define such a test...
        /// </summary>
        [Test]
        public void AllUp() {
            //static void Main(string[] args) {

            bool MpiInit;
            ilPSP.Environment.Bootstrap(new string[0], BoSSS.Solution.Application.GetBoSSSInstallDir(), out MpiInit);

            var Tests = new ITest[] { new Schraeg(XdgNastyLevsetLocationTest.GetTestRange(), XdgNastyLevsetLocationTest.GetTestRange()),
                new Parallel(XdgNastyLevsetLocationTest.GetTestRange(), XdgNastyLevsetLocationTest.GetTestRange()) };

            XQuadFactoryHelper.MomentFittingVariants[] Variants = new[] {
                XQuadFactoryHelper.MomentFittingVariants.OneStepGauss,
                XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes };


            foreach (var tst in Tests) {
                foreach (var variant in Variants) {

                    XdgNastyLevsetLocationTest p = null;

                    tst.ResetTest();

                    BoSSS.Solution.Application._Main(new string[0], true, null, delegate() {
                        p = new XdgNastyLevsetLocationTest();
                        p.test = tst;
                        p.momentFittingVariant = variant;
                        return p;
                    });

                    Assert.IsTrue(p.IsPassed);
                }
            }

            MPI.Wrappers.csMPI.Raw.mpiFinalize();
        }

    }
}
