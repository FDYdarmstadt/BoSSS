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
using System.IO;
using MPI.Wrappers;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// An all-up NUnit for XDG
    /// </summary>
    [TestFixture]
    public class UnitTest {

        [Test]
        public static void AllUp() {
            XDGTestMain p = null;
            BoSSS.Solution.Application._Main(new string[0], true, "", delegate() {
                p = new XDGTestMain();
                return p;
            });


            double err = p.AutoExtrapolationErr;
            double thres = 1.0e-10;

            Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
            Assert.IsTrue(err < thres);
        }
    }


}
