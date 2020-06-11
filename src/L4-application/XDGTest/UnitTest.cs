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
using BoSSS.Foundation.IO;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XDGTest {

    /// <summary>
    /// An all-up NUnit for XDG
    /// </summary>
    [TestFixture]
    public class UnitTest {


        static public XDGTestControl AllUpTestControl() {
            var ctrl = new XDGTestControl();

            ctrl.SetDGdegree(2);

            ctrl.InitialValues_Evaluators.Add("Phi", X => ((X[0] - 0.83) / 0.8).Pow2() + (X[1] / 0.8).Pow2() - 1.0);
            ctrl.InitialValues_Evaluators.Add("Pressure#A", X => 2 + 0.3 * X[0] * X[1]);
            ctrl.InitialValues_Evaluators.Add("Pressure#B", X => 1 - X[1].Pow2());

            ctrl.GridFunc = delegate () {
                var xNodes = GenericBlas.Linspace(-0.33333, 0.666667, 7);
                var yNodes = GenericBlas.Linspace(-1, 1, 13);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                return grd;
            };

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 3;

            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 1.0;
            ctrl.NoOfTimesteps = 1;


            return ctrl;
        }



        [Test]
        public static void AllUp() {
            XDGTestMain p = null;
            XDGTestMain._Main(new string[] { "--control", "cs:BoSSS.Application.XDGTest.UnitTest.AllUpTestControl()" /*"--delplt", "--implt", "1", "-u4"*/}, false, delegate() {
                p = new XDGTestMain();
                return p;
            });


            double err = p.AutoExtrapolationErr;
            double thres = 1.0e-10;

            Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
            Assert.LessOrEqual(err, thres);
        }


         static public XDGTestControl RestartTest_FirstControl(Guid gridGuid, IDatabaseInfo db) {
            var ctrl = new XDGTestControl();

            ctrl.SetDGdegree(2);

            ctrl.InitialValues_Evaluators.Add("Phi", X => ((X[0] - 0.83) / 0.8).Pow2() + (X[1] / 0.8).Pow2() - 1.0);
            ctrl.InitialValues_Evaluators.Add("Pressure#A", X => 1 - X[1].Pow2());
            ctrl.InitialValues_Evaluators.Add("Pressure#B", X => 1 - X[1].Pow2());

            ctrl.GridFunc = delegate () {
                var xNodes = GenericBlas.Linspace(-0.33333, 0.666667, 7);
                var yNodes = GenericBlas.Linspace(-1, 1, 13);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                return grd;
            };


            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 1.0;
            ctrl.NoOfTimesteps = 5;


            return ctrl;
        }


        [Test]
        public static void RestartTest() {
            string TestDbDir = "testdb_" + DateTime.Now.ToString("MMMdd_HHmm");

            var TestDb = DatabaseInfo.CreateOrOpen(Path.Combine(Directory.GetCurrentDirectory(), TestDbDir));

            Guid gridGuid;
            {
                var xNodes = GenericBlas.Linspace(-1.0/3.0, 10.0/3.0, 11*3 + 1);
                var yNodes = GenericBlas.Linspace(-1, 1, 6*2 +1);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                Guid GridGuid = TestDb.Controller.DBDriver.SaveGrid(grd, TestDb);
            };



            using(var FirstRun = new XDGTestMain()) {
                var ctrl = new XDGTestControl();

                //FirstRun.Init()
                //FirstRun.RunSolverMode();
            }


        }

    }
}
