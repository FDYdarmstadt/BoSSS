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

using BoSSS.Solution;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CNS.Tests.MMS {
    /// <summary>
    /// 
    /// </summary>
    //[TestFixture]
    class MMS_Test: TestProgram<CNSControl> {
        //[Test]
        public static void TestMMS2D_unsteadyCNS() {
            Program<CNSControl> p = null;
            Application<CNSControl>._Main(
                new string[] { @"-c cs:CNS.Tests.MMS.MMS_unsteady.Gassner2DStudy_conserved(3,3)" },
                false,
                delegate () {
                    p = new MMS_Test();
                    return p;
                });

            TestUtils.CheckConvergenceRates(p.QueryResultTable, p.Grid.SpatialDimension);
        }
    }
}
