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

            ctrl.InitialValues_Evaluators.Add("Phi", XDGTestMain.Phi0);
            ctrl.InitialValues_Evaluators.Add("Pressure#A", XDGTestMain.PressureExactA);
            ctrl.InitialValues_Evaluators.Add("Pressure#B", XDGTestMain.PressureExactB);

            ctrl.GridFunc = delegate () {
                var xNodes = GenericBlas.Linspace(-0.33333, 0.666667, 7);
                var yNodes = GenericBlas.Linspace(-1, 1, 13);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                return grd;
            };

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 3;

            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 1.0;
            ctrl.NoOfTimesteps = 1;


            return ctrl;
        }


        static public XDGTestControl PeriodicBoundaryTest() {
            var ctrl = new XDGTestControl();

            ctrl.SetDGdegree(2);

            ctrl.InitialValues_Evaluators.Add("Phi", XDGTestMain.Phi0);
            ctrl.InitialValues_Evaluators.Add("Pressure#A", XDGTestMain.PressureExactA);
            ctrl.InitialValues_Evaluators.Add("Pressure#B", XDGTestMain.PressureExactB);

            ctrl.GridFunc = delegate () {
                var xNodes = GenericBlas.Linspace(-1, 2, 25);
                xNodes = xNodes.Take(24).ToArray();
                var yNodes = GenericBlas.Linspace(-1, 1, 13);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX:true, periodicY: true);
                return grd;
            };

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 3;

            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 1.0;
            ctrl.NoOfTimesteps = 1;


            return ctrl;
        }




        [Test]
        public static void AllUp() {
            //XDGTestMain p = null;
            //XDGTestMain._Main(new string[] { "--control", "cs:BoSSS.Application.XDGTest.UnitTest.AllUpTestControl()" /*"--delplt", "--implt", "1", "-u4"*/}, false, delegate() {
            //    p = new XDGTestMain();
            //    return p;
            //});


            var ctrl = BoSSS.Application.XDGTest.UnitTest.AllUpTestControl();
            using(var p = new XDGTestMain()) {
                p.Init(ctrl);
                p.RunSolverMode();

                double err = p.AutoExtrapolationErr;
                double thres = 1.0e-10;

                Console.WriteLine("L2 Error of solution: " + err + " (threshold is " + thres + ")");
                Assert.LessOrEqual(err, thres);
            }
        }


        static public XDGTestControl RestartTest_FirstControl(string DbPath, out int[] ExpectedTimeSteps) {
            
            var ctrl = new XDGTestControl();

            ctrl.SetDGdegree(2);

            ctrl.InitialValues_Evaluators.Add("Phi", XDGTestMain.Phi0);
            ctrl.InitialValues_Evaluators.Add("Pressure#A", XDGTestMain.PressureExactA);
            ctrl.InitialValues_Evaluators.Add("Pressure#B", XDGTestMain.PressureExactB);

            ctrl.GridFunc = delegate ()  {
                var xNodes = GenericBlas.Linspace(-1.0/3.0, 10.0/3.0, 11*3 + 1);
                var yNodes = GenericBlas.Linspace(-1, 1, 6 * 2 + 1);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                //Guid GridGuid = TestDb.Controller.DBDriver.SaveGrid(grd, TestDb);
                return grd;
            };

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 3;

            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 0.1;
            ctrl.NoOfTimesteps = 10;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = 5;
            ctrl.rollingSaves = false;
            ExpectedTimeSteps = new int[] { 0, 4, 5, 9, 10 };

            return ctrl;
        }
        

        static public XDGTestControl RestartTest_SecondControl(string DbPath, Guid RestartSession, out int[] ExpectedTimeSteps) {
            
            var ctrl = new XDGTestControl();

            ctrl.SetDGdegree(2);

            ctrl.RestartInfo = Tuple.Create(RestartSession, default(TimestepNumber));
           
            ctrl.TimesteppingMode = Solution.Control.AppControl._TimesteppingMode.Transient;
            ctrl.dtFixed = 0.1;
            ctrl.NoOfTimesteps = 20;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = 50;
            ctrl.rollingSaves = true;
            ExpectedTimeSteps = new int[] { 10, 29, 30 };


            return ctrl;
        }


        [Test]
        public static void RestartTest() {
            string TestDbDir = "testdb_" + DateTime.Now.ToString("MMMdd_HHmm");
            string TestDbFullPath = Path.Combine(Directory.GetCurrentDirectory(), TestDbDir);

            {
                var TestDb = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                DatabaseInfo.Close(TestDb);
            }
            /*
            Guid gridGuid;
            {
                var xNodes = GenericBlas.Linspace(-1.0/3.0, 10.0/3.0, 11*3 + 1);
                var yNodes = GenericBlas.Linspace(-1, 1, 6 * 2 + 1);
                var grd = Grid2D.Cartesian2DGrid(xNodes, yNodes);

                Guid GridGuid = TestDb.Controller.DBDriver.SaveGrid(grd, TestDb);
            };
            */

            var ctrl1 = RestartTest_FirstControl(TestDbFullPath, out var ExpectedTs1stRun);
            using(var FirstRun = new XDGTestMain()) {

                FirstRun.Init(ctrl1);
                FirstRun.RunSolverMode();
            }


            Guid RestartSession;
            {
                var TestDb2 = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                Assert.IsTrue(TestDb2.Grids.Count() == 1, "Number of grids seems to be wrong.");
                Assert.IsTrue(TestDb2.Sessions.Count() == 1, "Number of sessions seems to be wrong.");


                var si = TestDb2.Sessions.Single();
                int[] tsiNumbers = si.Timesteps.Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                Assert.IsTrue(ExpectedTs1stRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");

                //var tend = si.Timesteps.Last();
                RestartSession = si.ID;
                
                DatabaseInfo.Close(TestDb2);
            }


            var ctrl2 = RestartTest_SecondControl(TestDbFullPath, RestartSession, out var ExpectedTs2ndRun);
            using(var SecondRun = new XDGTestMain()) {

                SecondRun.Init(ctrl2);
                SecondRun.RunSolverMode();
            }


            {
                var TestDb3 = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                Assert.IsTrue(TestDb3.Grids.Count() == 1, "Number of grids seems to be wrong.");
                Assert.IsTrue(TestDb3.Sessions.Count() == 2, "Number of sessions seems to be wrong.");


                var si = TestDb3.Sessions.First(); 
                int[] tsiNumbers = si.Timesteps.Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                Assert.IsTrue(ExpectedTs2ndRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");

                   
                DatabaseInfo.Close(TestDb3);
            }

            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            if(rank == 0) {
                Console.WriteLine($"Deleting test database at {TestDbFullPath} ...");
                Directory.Delete(TestDbFullPath, true);
                Console.WriteLine("done.");
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
        }


    }
}
