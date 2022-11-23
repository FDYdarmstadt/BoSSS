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
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Unit tests for the XNSE restart procedures
    /// </summary>
    [TestFixture]
    public class RestartTests {


        /// <summary>
        /// reference simulation (rising bubble testcase 1) for restart test 
        /// </summary>
        /// <param name="DbPath"></param>
        /// <param name="ExpectedTimeSteps"></param>
        /// <returns></returns>
        static public XNSE_Control RestartTest_ReferenceControl(string DbPath, bool transient, TimeSteppingScheme timeStepScheme, bool AMRon, out int[] ExpectedTimeSteps) {

            var ctrl = new XNSE_Control();

            ctrl.SetDGdegree(2);

            ctrl.PhysicalParameters.rho_A = 100;
            ctrl.PhysicalParameters.rho_B = 1000;
            ctrl.PhysicalParameters.mu_A = 1;
            ctrl.PhysicalParameters.mu_B = 10;
            ctrl.PhysicalParameters.Sigma = 24.5;
            ctrl.PhysicalParameters.IncludeConvection = false;

            int kelem = 10;
            ctrl.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, 1.0, kelem + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 1.0, kelem + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: false);

                grd.DefineEdgeTags(delegate(Vector X) {
                   string ret = null;
                    if(X.x.Abs() <= 1e-8 || X.x.Abs() - 1.0 <= 1.0e-8)
                        ret = IncompressibleBcType.FreeSlip.ToString();
                    if (X.y.Abs() <= 1e-8 || X.y.Abs() - 1.0 <= 1.0e-8)
                        ret = IncompressibleBcType.Wall.ToString();
                    return ret;
                });  

                return grd;
            };
            ctrl.AdaptiveMeshRefinement = AMRon;
            ctrl.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });

            double[] center = new double[] { 0.5, 0.5 };
            double radius = 0.25;
            Func<double[], double> PhiFunc = (X => Math.Sqrt((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - radius); 
            ctrl.InitialValues_Evaluators.Add("Phi", PhiFunc);

            ctrl.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            ctrl.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-9;


            ctrl.TimesteppingMode = transient ? AppControl._TimesteppingMode.Transient : AppControl._TimesteppingMode.Steady;
            ctrl.TimeSteppingScheme = timeStepScheme;
            ctrl.Timestepper_LevelSetHandling = transient ? LevelSetHandling.Coupled_Once : LevelSetHandling.None;
            ctrl.Option_LevelSetEvolution = transient ? Solution.LevelSetTools.LevelSetEvolution.FastMarching : Solution.LevelSetTools.LevelSetEvolution.None;


            ctrl.dtFixed = 0.02;
            ctrl.NoOfTimesteps = 6;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = 3;
            ctrl.rollingSaves = false;
            ctrl.BurstSave = 1;
            ctrl.MultiStepInit = false;

            ExpectedTimeSteps = new int[] { 0, 1, 2, 3, 4, 5, 6 };
            //switch (timeStepScheme) {
            //    case TimeSteppingScheme.ImplicitEuler:
            //        ExpectedTimeSteps = new int[] { 0, 3, 6 };
            //        break;
            //    case TimeSteppingScheme.BDF2:
            //        ExpectedTimeSteps = new int[] { 0, 2, 3, 5, 6 };
            //        break;
            //    case TimeSteppingScheme.BDF3:
            //        ExpectedTimeSteps = new int[] { 0, 1, 2, 3, 4, 5, 6 };
            //        break;
            //    default:
            //        throw new ArgumentException("Chosen timestepping scheme not supported for current test setting");
            //}

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 2;


            return ctrl;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="DbPath"></param>
        /// <param name="RestartSession"></param>
        /// <param name="ExpectedTimeSteps"></param>
        /// <returns></returns>
        static public XNSE_Control RestartTest_RestartControl(string DbPath, Guid RestartSession, bool transient, TimeSteppingScheme timeStepScheme, bool AMRon, out int[] ExpectedTimeSteps) {

            var ctrl = new XNSE_Control();

            ctrl.SetDGdegree(2);

            ctrl.PhysicalParameters.rho_A = 100;
            ctrl.PhysicalParameters.rho_B = 1000;
            ctrl.PhysicalParameters.mu_A = 1;
            ctrl.PhysicalParameters.mu_B = 10;
            ctrl.PhysicalParameters.Sigma = 24.5;
            ctrl.PhysicalParameters.IncludeConvection = false;

            ctrl.RestartInfo = Tuple.Create(RestartSession, new TimestepNumber(6));

            ctrl.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            ctrl.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);
            //ctrl.SetVolumeForce("A", 1, (X, t) => -9.81e-1);
            //ctrl.SetVolumeForce("B", 1, (X, t) => -9.81e-1);


            ctrl.AdaptiveMeshRefinement = AMRon;
            ctrl.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });


            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-9;


            ctrl.TimesteppingMode = transient ? AppControl._TimesteppingMode.Transient : AppControl._TimesteppingMode.Steady;
            ctrl.TimeSteppingScheme = timeStepScheme;
            ctrl.Timestepper_LevelSetHandling = transient ? LevelSetHandling.Coupled_Once : LevelSetHandling.None;
            ctrl.Option_LevelSetEvolution = transient ? Solution.LevelSetTools.LevelSetEvolution.FastMarching : Solution.LevelSetTools.LevelSetEvolution.None;

            ctrl.dtFixed = 0.02;
            ctrl.NoOfTimesteps = 3;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = 9;
            ctrl.rollingSaves = true;

            switch (timeStepScheme) {
                case TimeSteppingScheme.ImplicitEuler:
                    ctrl.MultiStepInit = false;
                    ExpectedTimeSteps = new int[] { 3, 6 };
                    break;
                case TimeSteppingScheme.BDF2:
                    ctrl.MultiStepInit = true;
                    ExpectedTimeSteps = new int[] { 3, 5, 6 };
                    break;
                case TimeSteppingScheme.BDF3:
                    ctrl.MultiStepInit = true;
                    ExpectedTimeSteps = new int[] { 3, 4, 5, 6 };
                    break;
                default:
                    throw new ArgumentException("Chosen timestepping scheme not supported for current test setting");
            }

            ctrl.ImmediatePlotPeriod = 1;
            ctrl.SuperSampling = 2;


            return ctrl;
        }




        /// <summary>
        /// Start a XNSE simulation with full complexity, calculate a few timesteps, save and load from db.
        /// Checks that all specified fields are correctly stored and loaded.
        /// </summary>
        [Test]
        public static void RestartTest(
              //[Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timestepScheme,
              //[Values(false, true)] bool AMRon) {
                [Values(false)] bool transient,
                [Values(TimeSteppingScheme.BDF3)] TimeSteppingScheme timestepScheme,
                [Values(false)] bool AMRon) {

            string TestDbDir = "testdb_" + DateTime.Now.ToString("MMMdd_HHmm");
            string TestDbFullPath = Path.Combine(Directory.GetCurrentDirectory(), TestDbDir);

            {
                var TestDb = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                DatabaseInfo.Close(TestDb);
            }

            var ctrl1 = RestartTest_ReferenceControl(TestDbFullPath, transient, timestepScheme, AMRon, out var ExpectedTs1stRun);
            //List<DGField> FieldsToCompare = new List<DGField>();
            using (var FirstRun = new XNSE()) {
                FirstRun.Init(ctrl1);
                FirstRun.RunSolverMode();
                // get current state
                //foreach (var f in FirstRun.IOFields) {
                //    FieldsToCompare.Add(f);
                //}
            }
            //Assert.IsTrue(ctrl1.FieldOptions.Where(s => s.Value.SaveToDB == FieldOpts.SaveToDBOpt.TRUE).ToArray().Count() == FieldsToCompare.Count, "Number of IO_Fields mismatch");

            Guid RestartSession;
            {
                var TestDb2 = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                int nGrids = 1;
                if (AMRon) nGrids++;
                if (AMRon && transient) nGrids++;
                Assert.IsTrue(TestDb2.Grids.Count() == nGrids, "Number of grids seems to be wrong.");
                Assert.IsTrue(TestDb2.Sessions.Count() == 1, "Number of sessions seems to be wrong.");

                var si = TestDb2.Sessions.Last(); // Single();
                int[] tsiNumbers = si.Timesteps.Skip(AMRon ? 2 : 0).Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                Assert.IsTrue(ExpectedTs1stRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");

                //var tend = si.Timesteps.Last();
                RestartSession = si.ID;

                DatabaseInfo.Close(TestDb2);
            }


            var ctrl2 = RestartTest_RestartControl(TestDbFullPath, RestartSession, transient, timestepScheme, AMRon, out var ExpectedTs2ndRun);
            using (var SecondRun = new XNSE()) {

                SecondRun.Init(ctrl2);
                SecondRun.RunSolverMode();
                // check current state
                //foreach (var f in FieldsToCompare) {
                //    var s = SecondRun.IOFields.Single(t => t.Identification == f.Identification);
                //    Assert.IsTrue(s != null);
                //    s.Coordinates.Acc(-1.0, f.Coordinates);
                //    //Assert.IsTrue(s.L2Norm() <= 1e-12, $"mismatch between reference and restarted computation for field {f.Identification}: L2Norm = {s.L2Norm()}");
                //    Console.WriteLine($"field {f.Identification}: L2-norm = {s.L2Norm()}");
                //}
            }


            {
                var TestDb3 = DatabaseInfo.CreateOrOpen(TestDbFullPath);
                int nGrids = 1;
                if (AMRon) nGrids++;
                if (AMRon && transient) nGrids += 2;
                //Assert.IsTrue(TestDb3.Grids.Count() == nGrids, "Number of grids seems to be wrong.");
                //Assert.IsTrue(TestDb3.Sessions.Count() == 2, "Number of sessions seems to be wrong.");


                var siRestart = TestDb3.Sessions.First();
                int[] tsiNumbers = siRestart.Timesteps.Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                //Assert.IsTrue(ExpectedTs2ndRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");


                var siRef = TestDb3.Sessions.Where(s => s.ID.Equals(siRestart.RestartedFrom)).Single();
                Assert.AreNotEqual(siRestart.ID, siRef.ID);
                var tsiRestart = siRestart.Timesteps; 
                foreach (var tsi in tsiRestart) {
                    Console.WriteLine($"========== timestep {tsi.TimeStepNumber.MajorNumber} ==========");
                    var tsiComparison = siRef.Timesteps.Single(t => t.TimeStepNumber.MajorNumber == tsi.TimeStepNumber.MajorNumber);
                    foreach (var f in tsi.Fields) {
                        var s = tsiComparison.Fields.Single(fRe => fRe.Identification == f.Identification);
                        Assert.IsTrue(s != null);
                        s.Coordinates.Acc(-1.0, f.Coordinates);
                        Console.WriteLine($"field {f.Identification}: L2-norm = {s.L2Norm()}");
                    }
                }

                DatabaseInfo.Close(TestDb3);
            }


            //csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            //if (rank == 0) {
            //    Console.WriteLine($"Deleting test database at {TestDbFullPath} ...");
            //    Directory.Delete(TestDbFullPath, true);
            //    Console.WriteLine("done.");
            //}
            //csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);
        }

    }
}
