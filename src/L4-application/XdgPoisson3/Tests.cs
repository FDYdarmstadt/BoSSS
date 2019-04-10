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
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using Code = BoSSS.Solution.Control.LinearSolverConfig.Code;

namespace BoSSS.Application.XdgPoisson3 {
    [TestFixture]
    static public class Tests {


        /// <summary>
        /// MPI finalize.
        /// </summary>
        [TestFixtureTearDown]
        static public void TestFixtureTearDown() {
            csMPI.Raw.mpiFinalize();
        }

        /// <summary>
        /// MPI init.
        /// </summary>
        [TestFixtureSetUp]
        static public void TestFixtureSetUp() {
            BoSSS.Solution.Application.InitMPI(new string[0]);
            XQuadFactoryHelper.CheckQuadRules = true;
        }


        [Test]
        public static void SolverTest([Values(Code.exp_softpcg_mg
            //, "pcg+schwarz", "pcg+mg+schwarz", "gmres+mg+schwarz", "ono+mg+schwarz"
            )] Code SolverName) {
        using (var solver = new XdgPoisson3Main()) {
            var C = HardCodedControl.Circle(solver: SolverName);


            solver.Init(C);
            solver.RunSolverMode();
        }
    }


    }
}
