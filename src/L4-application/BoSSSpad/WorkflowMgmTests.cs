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


using BoSSS.Application.BoSSSpad;
using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.ComponentModel.DataAnnotations;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading;

namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// Test for the workflow management; should actually be part of BoSSSpad, but here for dependency reasons
    /// </summary>
    [TestFixture]
    [NUnitNumThreads(1)]
    static public class WorkflowMgmTests {


        /// <summary>
        /// test of job persistence, i.e., that a second try to submit a job actually finds the run which is already in place.
        /// </summary>
        [Test]
        static public void Run__WorkFlowManagerFunctionalityCheck_1(
            [Values(1, 2)] int SessionEqualityComparisonVariant
            ) {
            BoSSS.Solution.Application.InitMPI(num_threads: 1);

            string basename = "WorkFlowManagerFunctionalityCheck_1" + SessionEqualityComparisonVariant;

            GridCommons GridFunction() {
                var xNodes = GenericBlas.Linspace(0, 10, 41);
                var yNodes = GenericBlas.Linspace(-1, 1, 9);
                GridCommons grid = Grid2D.Cartesian2DGrid(xNodes, yNodes);
                IGrid_Extensions.DefineEdgeTags(grid, delegate (double[] X) {
                    double x = X[0];
                    double y = X[1];
                    if(Math.Abs(y - (-1)) <= 1.0e-8)
                        return "wall"; // lower wall
                    if(Math.Abs(y - (+1)) <= 1.0e-8)
                        return "wall"; // upper wall
                    if(Math.Abs(x - (0.0)) <= 1.0e-8)
                        return "Velocity_Inlet"; // inlet
                    if(Math.Abs(x - (+10.0)) <= 1.0e-8)
                        return "Pressure_Outlet"; // outlet
                    throw new ArgumentOutOfRangeException("unknown domain");
                });
                return grid;
            }

            XNSE_Control ControlFunc() {
                var c = new XNSE_Control();

                // general description:
                int k = 1;
                string desc = "Steady state, channel, k" + k;
                c.SessionName = "SteadyStateChannel";
                c.ProjectDescription = desc;
                c.savetodb = true;
                c.Tags.Add("k" + k);
                // setting the grid:
                c.SetGrid(GridFunction());
                // DG polynomial degree
                c.SetDGdegree(k);
                // Physical parameters:
                double reynolds = 20;
                c.PhysicalParameters.rho_A = 1;
                c.PhysicalParameters.mu_A = 1.0 / reynolds;
                // Timestepping properties:
                c.TimesteppingMode = AppControl._TimesteppingMode.Steady;

                c.SessionName = "Channel_40x8k" + k;

                c.AddBoundaryValue("Velocity_Inlet", "VelocityX", new Formula("X => 1 - X[0]*X[0]", false));

                return c;
            }


            Mutex TestMutex = new Mutex(false, basename);
            try {
                TestMutex.WaitOne();
                
                // 0.) Clean leftover of previous runs
                // -----------------------------------
                {
                    WorkflowMgmTestUtils.DeleteDatabase(basename);
                    WorkflowMgmTestUtils.DeleteDeployments(basename + "*");
                    
                    BoSSSshell.WorkflowMgm.Init(basename);

                    switch(SessionEqualityComparisonVariant) {
                        case 1:
                        BoSSSshell.WorkflowMgm.SetEqualityBasedSessionJobControlCorrelation();
                        break;

                        case 2:
                        BoSSSshell.WorkflowMgm.SetNameBasedSessionJobControlCorrelation();
                        break;

                        default:
                        throw new ArgumentOutOfRangeException(nameof(SessionEqualityComparisonVariant) + " is expected to be either 1 or 2");
                    }



                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.Zero(allDepl.Length, $"expecting 0 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));
                }

                // 1.) First run
                // -----------------------------------
                Guid sessionGuid;
                {
                    Console.WriteLine("========================================");
                    Console.WriteLine("  first job run...");
                    Console.WriteLine("========================================");

                    var c0 = ControlFunc();
                    var J = c0.CreateJob();
                    Assert.AreEqual(J.Status, JobStatus.PreActivation, "unexpected job status");
                    J.NumberOfMPIProcs = 2;
                    J.Activate();

                    BoSSSshell.WorkflowMgm.BlockUntilAllJobsTerminate(1000);

                    Assert.IsTrue(J.Control.Equals(BoSSSshell.WorkflowMgm.Sessions.First().GetControl()), "equality check for control objects does not seem to work");

                    Assert.AreEqual(J.Status, JobStatus.FinishedSuccessful, "unexpected job status");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, basename, "un-expected project name");
                    Assert.AreEqual(J.AllDeployments.Count(), 1, $"expecting 1 deployment");
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(allDepl.Length, 1, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));


                    sessionGuid = J.LatestSession.ID;
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Console.WriteLine("done run 1. =============================");
                    Console.WriteLine();
                    Console.WriteLine();
                }

                // 2.) Persistence test: Test that, if the Worksheet is already executed, the job will **not** be re-submitted
                // --------------------------------------------------------------------------------------------------------------------------
                {
                    BoSSSshell.WorkflowMgm.ResetProject();

                    Console.WriteLine("========================================");
                    Console.WriteLine("  second job run...");
                    Console.WriteLine("========================================");


                    // execute the worksheet a second time!
                    var c1 = ControlFunc();
                    var J = c1.CreateJob();
                    Assert.AreEqual(J.Status, JobStatus.PreActivation, "unexpected job status");
                    J.NumberOfMPIProcs = 2;
                    J.Activate();

                    BoSSSshell.WorkflowMgm.BlockUntilAllJobsTerminate(1000);

                    Assert.IsTrue(J.Control.Equals(BoSSSshell.WorkflowMgm.Sessions.First().GetControl()), "equality check for control objects does not seem to work");

                    Assert.AreEqual(J.Status, JobStatus.FinishedSuccessful, "unexpected job status");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, basename, "un-expected project name");
                    Assert.AreEqual(J.AllDeployments.Count(), 1, $"expecting 1 deployment");
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(allDepl.Length, 1, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                    Assert.AreEqual(sessionGuid, J.LatestSession.ID, "expecting session ID known from first run");
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Console.WriteLine("done run 2. =============================");
                    Console.WriteLine();
                    Console.WriteLine();                
                }


                // 3.) delete deployments and test a third time:
                // -----------------------------------------------
                {
                    BoSSSshell.WorkflowMgm.ResetProject(ResetJobs:true, deleteDeployments:false);
                    WorkflowMgmTestUtils.DeleteDeployments(basename + "*");

                    Console.WriteLine("========================================");
                    Console.WriteLine("  third job run...");
                    Console.WriteLine("========================================");


                    // execute the worksheet a third time!
                    var c3 = ControlFunc();
                    var J = c3.CreateJob();
                    Assert.AreEqual(J.Status, JobStatus.PreActivation, "unexpected job status");
                    J.NumberOfMPIProcs = 1;
                    J.Activate();

                    Assert.IsTrue(J.Control.Equals(BoSSSshell.WorkflowMgm.Sessions.First().GetControl()), "equality check for control objects does not seem to work");

                    Assert.AreEqual(J.Status, JobStatus.FinishedSuccessful, "unexpected job status");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, basename, "un-expected project name");
                    Assert.AreEqual(J.AllDeployments.Count(), 1, $"expecting 1 deployment");

                    // note: job should have one deployment, if we can find an existing session,
                    // but there will be no more deployment directory.

                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(allDepl.Length, 0, $"expecting 0 deployment directories, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                    Assert.AreEqual(sessionGuid, J.LatestSession.ID, "expecting session ID known from first run");
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Console.WriteLine("done run 3. =============================");
                    Console.WriteLine();
                    Console.WriteLine();
                }
            } catch(Exception) {
                throw;
            } finally {
                TestMutex.ReleaseMutex();
            }
            Console.WriteLine("success!");


        }
    }
}