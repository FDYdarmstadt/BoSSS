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
using BoSSS.Solution.Tecplot;
using Newtonsoft.Json.Linq;

namespace BoSSS.Application.XNSE_Solver.Tests {

    /// <summary>
    /// Unit tests for the XNSE restart procedures
    /// </summary>
    [TestFixture]
    public class RestartTest {


        /// <summary>
        /// reference simulation (rising bubble testcase 1) for restart test 
        /// </summary>
        /// <param name="DbPath"></param>
        /// <param name="ExpectedTimeSteps"></param>
        /// <returns></returns>
        static public XNSE_Control RestartTest_ReferenceControl(string DbPath, bool transient, LevelSetHandling LevSetHandling, TimeSteppingScheme timeStepScheme, bool AMRon, int savePeriod, out int[] ExpectedTimeSteps) {

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
            ctrl.DynamicLoadBalancing_RedistributeAtStartup = false;

            double[] center = new double[] { 0.5, 0.5 };
            double radius = 0.25;
            Func<double[], double> PhiFunc = (X => Math.Sqrt((X[0] - center[0]).Pow2() + (X[1] - center[1]).Pow2()) - radius); 
            ctrl.InitialValues_Evaluators.Add("Phi", PhiFunc);

            ctrl.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            ctrl.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);


            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-15;
            ctrl.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            //ctrl.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();


            ctrl.TimesteppingMode = transient ? AppControl._TimesteppingMode.Transient : AppControl._TimesteppingMode.Steady;
            ctrl.TimeSteppingScheme = timeStepScheme;
            ctrl.Timestepper_LevelSetHandling = transient ? LevSetHandling : LevelSetHandling.None;
            ctrl.Option_LevelSetEvolution = transient ? Solution.LevelSetTools.LevelSetEvolution.FastMarching : Solution.LevelSetTools.LevelSetEvolution.None;


            ctrl.dtFixed = 0.02;
            ctrl.NoOfTimesteps = 6;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = savePeriod;
            ctrl.rollingSaves = true;
            ctrl.MultiStepInit = false;


            ExpectedTimeSteps = new int[] { }; // { 0, 3, 4, 3, 5, 8, 9, 10 };

            switch (timeStepScheme) {
                case TimeSteppingScheme.ImplicitEuler:
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 0, 3, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 0, 4, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 0, 5, 6 };
                    break;
                case TimeSteppingScheme.BDF2:
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 0, 2, 3, 5, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 0, 3, 4, 3, 5, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 0, 4, 5, 6 };
                    break;
                case TimeSteppingScheme.BDF3:
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 0, 1, 2, 3, 4, 5, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 0, 2, 3, 4, 3, 2, 5, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 0, 3, 4, 3, 5, 6 };
                    break;
                default:
                    throw new ArgumentException("Chosen timestepping scheme not supported for current test setting");
            }

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 3;

            return ctrl;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="DbPath"></param>
        /// <param name="RestartSession"></param>
        /// <param name="ExpectedTimeSteps"></param>
        /// <returns></returns>
        static public XNSE_Control RestartTest_RestartControl(string DbPath, Guid RestartSession, bool transient, LevelSetHandling LevSetHandling, TimeSteppingScheme timeStepScheme, bool AMRon, int savePeriod,  out int[] ExpectedTimeSteps) {

            var ctrl = new XNSE_Control();

            ctrl.SetDGdegree(2);

            ctrl.PhysicalParameters.rho_A = 100;
            ctrl.PhysicalParameters.rho_B = 1000;
            ctrl.PhysicalParameters.mu_A = 1;
            ctrl.PhysicalParameters.mu_B = 10;
            ctrl.PhysicalParameters.Sigma = 24.5;
            ctrl.PhysicalParameters.IncludeConvection = false;

            ctrl.RestartInfo = Tuple.Create(RestartSession, new TimestepNumber(savePeriod));

            ctrl.InitialValues_Evaluators.Add("GravityY#A", X => -9.81e-1);
            ctrl.InitialValues_Evaluators.Add("GravityY#B", X => -9.81e-1);
            //ctrl.SetVolumeForce("A", 1, (X, t) => -9.81e-1);
            //ctrl.SetVolumeForce("B", 1, (X, t) => -9.81e-1);


            ctrl.AdaptiveMeshRefinement = AMRon;
            ctrl.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            ctrl.DynamicLoadBalancing_RedistributeAtStartup = false;


            ctrl.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            ctrl.NonLinearSolver.ConvergenceCriterion = 1e-15;
            ctrl.LinearSolver = LinearSolverCode.direct_pardiso.GetConfig();
            //ctrl.LinearSolver = LinearSolverCode.direct_mumps.GetConfig();


            ctrl.TimesteppingMode = transient ? AppControl._TimesteppingMode.Transient : AppControl._TimesteppingMode.Steady;
            ctrl.TimeSteppingScheme = timeStepScheme;
            ctrl.Timestepper_LevelSetHandling = transient ? LevSetHandling : LevelSetHandling.None;
            ctrl.Option_LevelSetEvolution = transient ? Solution.LevelSetTools.LevelSetEvolution.FastMarching : Solution.LevelSetTools.LevelSetEvolution.None;

            ctrl.dtFixed = 0.02;
            ctrl.NoOfTimesteps = 6 - savePeriod;

            ctrl.DbPath = DbPath;
            ctrl.saveperiod = 6;
            ctrl.rollingSaves = true;


            ExpectedTimeSteps = new int[] { }; // { 3, 4, 5, 6 };

            switch (timeStepScheme) {
                case TimeSteppingScheme.ImplicitEuler:
                    ctrl.MultiStepInit = false;
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 3, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 4, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 5, 6 };

                    break;
                case TimeSteppingScheme.BDF2:
                    ctrl.MultiStepInit = true;
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 3, 5, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 4, 5, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 5, 6 };
                    break;
                case TimeSteppingScheme.BDF3:
                    ctrl.MultiStepInit = true;
                    if (savePeriod == 3)
                        ExpectedTimeSteps = new int[] { 3, 4, 5, 6 };
                    if (savePeriod == 4)
                        ExpectedTimeSteps = new int[] { 4, 5, 6 };
                    if (savePeriod == 5)
                        ExpectedTimeSteps = new int[] { 5, 6 };
                    break;
                default:
                    throw new ArgumentException("Chosen timestepping scheme not supported for current test setting");
            }

            //ctrl.ImmediatePlotPeriod = 1;
            //ctrl.SuperSampling = 2;


            return ctrl;
        }


        /// <summary>
        /// Collection of fields to be checked with corresponding max allowed comparison error
        /// </summary>
        static Dictionary<string, double> AllowedErrors = new Dictionary<string, double>()
        {
            { "Phi", 1e-15 },
            { "PhiDG", 1e-15 },
            { "VelocityX", 1e-13 },
            { "VelocityY", 1e-13 },
            { "Pressure", 1e-10 },
            { "Residual-MomentumX", 1e-9 },
            { "Residual-MomentumY", 1e-9 },
            { "Residual-ContiEq", 1e-14 },
            { "Velocity0X_Mean", 1e-15 },
            { "Velocity0Y_Mean", 1e-15 },
            { "VelocityX@Phi", 1e-13 },
            { "VelocityY@Phi", 1e-13 },
        };


        /// <summary>
        /// Start an XNSE simulation with full complexity, calculate a few timesteps, save and load from db.
        /// Checks that all specified fields are correctly stored and loaded.
        /// </summary>
        [Test]
        public static void Run_RestartTests(
                [Values(false, true)] bool transient,
                [Values(LevelSetHandling.LieSplitting, LevelSetHandling.Coupled_Once)] LevelSetHandling LevSetHandling,
                [Values(TimeSteppingScheme.ImplicitEuler, TimeSteppingScheme.BDF2, TimeSteppingScheme.BDF3)] TimeSteppingScheme timestepScheme,
                [Values(false, true)] bool AMRon,
                [Values(3, 4, 5)] int savePeriod){

            if (savePeriod > 3 && (!AMRon || !transient))
                return;

            string restartDB = "restartDB_" + transient + "_" + timestepScheme.ToString() + "_" + AMRon + "_" + savePeriod + DateTime.Now.ToString("MMMdd_HHmm");
            string restartDBfullPath = Path.Combine(Directory.GetCurrentDirectory(), restartDB);

            {
                var TestDb = DatabaseInfo.CreateOrOpen(restartDBfullPath);
                DatabaseInfo.Close(TestDb);
            }

            var ctrl1 = RestartTest_ReferenceControl(restartDBfullPath, transient, LevSetHandling, timestepScheme, AMRon, savePeriod, out var ExpectedTs1stRun);
            using (var FirstRun = new XNSE()) {

                FirstRun.Init(ctrl1);
                FirstRun.RunSolverMode();
            }


            Guid RestartSession;
            {
                var TestDb2 = DatabaseInfo.CreateOrOpen(restartDBfullPath);
                int nGrids = 1;
                if (AMRon) nGrids++;
                if (AMRon && transient) nGrids = 3;
                Assert.IsTrue(TestDb2.Grids.Count() == nGrids, "Number of grids seems to be wrong.");
                Assert.IsTrue(TestDb2.Sessions.Count() == 1, "Number of sessions seems to be wrong.");

                var si = TestDb2.Sessions.Single();
                int[] tsiNumbers = si.Timesteps.Skip(AMRon ? 2 : 0).Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                string output = "tsiNumbers = ";
                foreach (int n in tsiNumbers)
                    output += n + " "; 
                Console.WriteLine(output);
                Assert.IsTrue(ExpectedTs1stRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");

                RestartSession = si.ID;

                DatabaseInfo.Close(TestDb2);
            }


            var ctrl2 = RestartTest_RestartControl(restartDBfullPath, RestartSession, transient, LevSetHandling, timestepScheme, AMRon, savePeriod, out var ExpectedTs2ndRun);
            using (var SecondRun = new XNSE()) {

                SecondRun.Init(ctrl2);
                SecondRun.RunSolverMode();

            }

           
            {
                var TestDb3 = DatabaseInfo.CreateOrOpen(restartDBfullPath);
                int nGrids = 1;
                if (AMRon) nGrids++;
                if (AMRon && transient) nGrids = (savePeriod == 3) ? 4 : 3;
                Assert.IsTrue(TestDb3.Grids.Count() == nGrids, "Number of grids seems to be wrong.");
                Assert.IsTrue(TestDb3.Sessions.Count() == 2, "Number of sessions seems to be wrong.");


                var siRestart = TestDb3.Sessions.First();
                int[] tsiNumbers = siRestart.Timesteps.Select(tsi => tsi.TimeStepNumber.MajorNumber).ToArray();
                string output2 = "tsiNumbers = ";
                foreach (int n in tsiNumbers)
                    output2 += n + " ";
                Console.WriteLine(output2);
                Assert.IsTrue(ExpectedTs2ndRun.ListEquals(tsiNumbers), "mismatch between saved time-steps in test database and expected saves.");

                bool comparisonFailed = false;
                var siRef = TestDb3.Sessions.Where(s => s.ID.Equals(siRestart.RestartedFrom)).Single();
                Assert.AreNotEqual(siRestart.ID, siRef.ID);
                var tsiRestart = siRestart.Timesteps;
                foreach (var tsi in tsiRestart) {
                    //Console.WriteLine($"========== timestep {tsi.TimeStepNumber.MajorNumber} ==========");
                    var tsiComparison = siRef.Timesteps.Single(t => t.TimeStepNumber.Equals(tsi.TimeStepNumber));
                    foreach (var f in tsi.Fields) {
                        var s = tsiComparison.Fields.Single(fRe => fRe.Identification == f.Identification);
                        Assert.IsTrue(s != null);
                        s.Coordinates.Acc(-1.0, f.Coordinates);

                        //Console.WriteLine($"timestep {tsi.TimeStepNumber.MajorNumber}: field {f.Identification} L2-norm = {s.L2Norm()}");
                        if (tsi.TimeStepNumber.MajorNumber == savePeriod) {
                            if (s.L2Norm() > 0.0) { // loaded data needs to be exact
                                Console.WriteLine($"Loaded data at timestep {tsi.TimeStepNumber.MajorNumber} for field {f.Identification} are not exact: L2-norm = {s.L2Norm()}");
                                Tecplot.PlotFields(new List<DGField>() { s }, $"{f.Identification}_errorField", 0.0, 3);
                                comparisonFailed = true;
                            }
                        } else if (AllowedErrors.TryGetValue(f.Identification, out double value)) {
                            if (s.L2Norm() > value) {
                                Console.WriteLine($"Allowed error exceeded at timestep {tsi.TimeStepNumber.MajorNumber} for field {f.Identification}: L2-norm = {s.L2Norm()}");
                                Tecplot.PlotFields(new List<DGField>() { s }, $"{f.Identification}_errorField", 0.0, 3);
                                comparisonFailed = true;
                            }
                        }

                    }
                }
                DatabaseInfo.Close(TestDb3);

                Assert.IsFalse(comparisonFailed, "comparison between reference solution and restart solution not equal");
            }


            csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out int rank);
            if (rank == 0) {
                Console.WriteLine($"Deleting test database at {restartDBfullPath} ...");
                Directory.Delete(restartDBfullPath, true);
                Console.WriteLine("done.");
            }
            csMPI.Raw.Barrier(csMPI.Raw._COMM.WORLD);

        }

    }
}
