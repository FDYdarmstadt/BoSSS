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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid.Classic;

namespace QuadratureAndProjectionTest {

    class TestApplication : Application {

        private static bool MPIInialized = false;

        private QuadratueAndProjectionTest test;

        public TestApplication(QuadratueAndProjectionTest test) {
            this.test = test;
            if (!MPIInialized) {
                Application.InitMPI(new string[0]);
                MPIInialized = true;
            }
            Init(null, m_opt, "");
            SetUpEnvironment();
        }

        public static void Main(string[] args) {
            //LineTest bla = new LineTest();
            //bla.SetUp();
            //bla.TestNominalOrder();
            //bla.TestSumOfWeights();
            //bla.TestMassMatrix();
            //bla.TestProjectionError();
            TetraTest blubb = new TetraTest();
            blubb.SetUp();
            blubb.TestNominalOrder();
            blubb.TestMassMatrix();
            blubb.TestSumOfWeights();
            blubb.TestProjectionError();
        }

        protected override GridCommons CreateOrLoadGrid() {
            return test.GetSingleCellGrid();
        }

        protected override void CreateEquationsAndSolvers(LoadBalancingData L) {
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
        }
    }
}
