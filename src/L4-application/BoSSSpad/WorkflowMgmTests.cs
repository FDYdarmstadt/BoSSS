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
        /// Tests batch-processor database path translation and generated remote roots.
        /// </summary>
        [Test]
        static public void AlternateDatabasePathTranslationTest() {
            string localRoot = Path.Combine(Path.GetTempPath(), "bosss-local-root");
            string localDbPath = Path.Combine(localRoot, @"case\database");

            var client = new AlternateDatabasePathTestBatchProcessor();
            client.AllowedDatabasesPaths.Clear();
            client.AllowedDatabasesPaths.Add(new AllowedDatabasesPair(localRoot, "/remote/root"));

            Assert.IsTrue(client.TryTranslateDatabasePath(localDbPath, out string remoteDbPath));
            Assert.AreEqual("/remote/root/case/database", remoteDbPath);

            var remoteRoots = client.GetAllowedDatabaseRemotePaths().ToArray();
            Assert.AreEqual(1, remoteRoots.Length);
            Assert.AreEqual("/remote/root", remoteRoots[0].Item1);
            Assert.IsNull(remoteRoots[0].Item2);
        }

        class AlternateDatabasePathTestBatchProcessor : BatchProcessorClient {
            public override void Cancel(string idToken, string message) {
                throw new NotImplementedException();
            }

            public override (JobStatus, int? ExitCode) EvaluateStatus(string idToken, object optInfo, string DeployDir) {
                throw new NotImplementedException();
            }

            public override string GetStderrFile(string idToken, string DeployDir) {
                throw new NotImplementedException();
            }

            public override string GetStdoutFile(string idToken, string DeployDir) {
                throw new NotImplementedException();
            }

            public override (string id, object optJobObj) Submit(Job myJob, string DeploymentDirectory) {
                throw new NotImplementedException();
            }
        }


        /// <summary>
        /// Session name set by <see cref="ControlFunc"/>; for name-based correlation of sessions and jobs, 
        /// the job must have the same name.
        /// </summary>
        public const string ChannelSessionName = WorkflowTestControls.ChannelSessionName;

        /// <summary>
        /// Control object for the channel test, see <see cref="WorkflowTestControls.ChannelControl"/>.
        /// 
        /// Note: to be used as a control statement (see <see cref="Job.SetControlStatement(string, IDatabaseInfo, bool?)"/>) on a batch system, 
        /// the statement must refer to the solver assembly, see <see cref="WorkflowTestControls"/>.
        /// </summary>
        static public XNSE_Control ControlFunc() {
            return WorkflowTestControls.ChannelControl();
        }

        /// <summary>
        /// test of job persistence, i.e., that a second try to submit a job actually finds the run which is already in place.
        /// </summary>
        [Test]
        static public void Run__WorkFlowManagerFunctionalityCheck_1(
            [Values(1, 2)] int SessionEqualityComparisonVariant
            ) {
            BoSSS.Solution.Application.InitMPI(num_threads: 1);

            string basename = "WorkFlowManagerFunctionalityCheck_1" + SessionEqualityComparisonVariant;

            Mutex TestMutex = new Mutex(false, basename);
            try {
                TestMutex.WaitOne();
                
                // 0.) Clean leftover of previous runs
                // -----------------------------------
                {
                    WorkflowMgmTestUtils.DeleteDatabase(basename);
                    WorkflowMgmTestUtils.DeleteDeployments(basename + "*");

                    BoSSSshell.Reset();
                    BoSSSshell.WorkflowMgm.Init(basename);
                    BoSSSshell.WorkflowMgm.ResetProject(true, true, true, true);

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

                    Console.WriteLine("Waiting for jobs to finish...");
                    BoSSSshell.WorkflowMgm.BlockUntilAllJobsTerminate(1000);
                    Console.WriteLine("All jobs finished.");
                    
                    /*
                    var si = BoSSSshell.WorkflowMgm.Sessions.First();
                    Console.WriteLine("First session: " +  si);
                    while(si.SuccessfulTermination == false) {
                        Console.WriteLine("   Triggering reload...");
                        if(si is SessionProxy sp)
                            sp.TriggerReload = true;
                        Thread.Sleep(1000*3);
                        Console.WriteLine("    First session: " + si);
                    }
                    */

                    Assert.IsTrue(J.Control.Equals(BoSSSshell.WorkflowMgm.Sessions.First().GetControl()), "equality check for control objects does not seem to work");

                    var s = J.GetStatus(true);
                    Assert.AreEqual(JobStatus.FinishedSuccessful, s, "unexpected job status: " + s);
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Length, 1, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, basename, "un-expected project name");
                    Assert.AreEqual(J.AllDeployments.Count(), 1, $"expecting 1 deployment");
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(allDepl.Length, 1, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));


                    sessionGuid = J.LatestSession.ID;
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Console.WriteLine("Session runtime: " + J.LatestSession.GetRunTime());
                    Console.WriteLine("Job runtime: " + J.LatestDeployment.RunTime);

                    //Assert.IsTrue(J.LatestDeployment.RunTime.TotalSeconds > 0, "expecting a non-zero runtime");

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

                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Grids.Length, "expecting exactly one grid");

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
                    J.NumberOfMPIProcs = 1; // different number of procs should habe no effect
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

                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Grids.Length, "expecting exactly one grid");

                    Console.WriteLine("done run 3. =============================");
                    Console.WriteLine();
                    Console.WriteLine();
                }
            } catch(Exception) {
                throw;
            } finally {
                //WorkflowMgmTestUtils.DeleteDeployments(basename + "*");
                //WorkflowMgmTestUtils.DeleteDatabase(basename);

                TestMutex.ReleaseMutex();

            }
            Console.WriteLine("success!");


        }

        /// <summary>
        /// Same as <see cref="Run__WorkFlowManagerFunctionalityCheck_1"/>, but the job is configured with a 
        /// control statement (<see cref="Job.SetControlStatement(string, IDatabaseInfo, bool?)"/>) instead of a control object;
        /// the correlation of job and session can therefore only be name-based
        /// (<see cref="WorkflowMgm.SetNameBasedSessionJobControlCorrelation"/>).
        /// </summary>
        [Test]
        static public void Run__WorkFlowManagerFunctionalityCheck_ControlStatement() {
            BoSSS.Solution.Application.InitMPI(num_threads: 1);

            string basename = "WorkFlowManagerFunctionalityCheck_CS";

            // Job name must be equal to the session name set by the control statement, 
            // otherwise, session and job are not correlated.
            Job CreateJobFromStatement(int NoOfProcs) {
                var J = new Job(ChannelSessionName, typeof(BoSSS.Application.XNSE_Solver.XNSE));
                // no database specified: the job management sets it (see `Job.FiddleControlStatement`)
                J.SetControlStatement("cs:BoSSS.Application.XNSE_Solver.WorkflowTestControls.ChannelControl(true)");
                J.NumberOfMPIProcs = NoOfProcs;
                return J;
            }

            Mutex TestMutex = new Mutex(false, basename);
            try {
                TestMutex.WaitOne();

                // 0.) Clean leftover of previous runs
                // -----------------------------------
                {
                    WorkflowMgmTestUtils.DeleteDatabase(basename);
                    WorkflowMgmTestUtils.DeleteDeployments(basename + "*");

                    BoSSSshell.Reset();
                    BoSSSshell.WorkflowMgm.Init(basename);
                    BoSSSshell.WorkflowMgm.ResetProject(true, true, true, true);
                    BoSSSshell.WorkflowMgm.SetNameBasedSessionJobControlCorrelation();

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

                    var J = CreateJobFromStatement(2);
                    Assert.AreEqual(JobStatus.PreActivation, J.Status, "unexpected job status");
                    J.Activate();

                    Console.WriteLine("Waiting for jobs to finish...");
                    BoSSSshell.WorkflowMgm.BlockUntilAllJobsTerminate(1000);
                    Console.WriteLine("All jobs finished.");

                    var s = J.GetStatus(true);
                    Assert.AreEqual(JobStatus.FinishedSuccessful, s, "unexpected job status: " + s);
                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Sessions.Length, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(basename, BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, "un-expected project name");
                    Assert.AreEqual(1, J.AllDeployments.Count(), "expecting 1 deployment");
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(1, allDepl.Length, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                    sessionGuid = J.LatestSession.ID;
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Console.WriteLine("done run 1. =============================");
                    Console.WriteLine();
                    Console.WriteLine();
                }

                // 2.) Persistence test: if the Worksheet is already executed, the job will **not** be re-submitted
                // ------------------------------------------------------------------------------------------------
                {
                    BoSSSshell.WorkflowMgm.ResetProject();

                    Console.WriteLine("========================================");
                    Console.WriteLine("  second job run...");
                    Console.WriteLine("========================================");

                    var J = CreateJobFromStatement(2);
                    Assert.AreEqual(JobStatus.PreActivation, J.Status, "unexpected job status");
                    J.Activate();

                    BoSSSshell.WorkflowMgm.BlockUntilAllJobsTerminate(1000);

                    Assert.AreEqual(JobStatus.FinishedSuccessful, J.Status, "unexpected job status");
                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Sessions.Length, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(basename, BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, "un-expected project name");
                    Assert.AreEqual(1, J.AllDeployments.Count(), "expecting 1 deployment");
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(1, allDepl.Length, $"expecting 1 deployments, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                    Assert.AreEqual(sessionGuid, J.LatestSession.ID, "expecting session ID known from first run");
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Grids.Length, "expecting exactly one grid");

                    Console.WriteLine("done run 2. =============================");
                    Console.WriteLine();
                    Console.WriteLine();
                }

                // 3.) delete deployments and test a third time:
                // -----------------------------------------------
                {
                    BoSSSshell.WorkflowMgm.ResetProject(ResetJobs: true, deleteDeployments: false);
                    WorkflowMgmTestUtils.DeleteDeployments(basename + "*");

                    Console.WriteLine("========================================");
                    Console.WriteLine("  third job run...");
                    Console.WriteLine("========================================");

                    var J = CreateJobFromStatement(1); // different number of procs should have no effect
                    Assert.AreEqual(JobStatus.PreActivation, J.Status, "unexpected job status");
                    J.Activate();

                    Assert.AreEqual(JobStatus.FinishedSuccessful, J.Status, "unexpected job status");
                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Sessions.Length, "expecting exactly one session after the workseet has been executed");
                    Assert.AreEqual(basename, BoSSSshell.WorkflowMgm.Sessions.Single().ProjectName, "un-expected project name");
                    Assert.AreEqual(1, J.AllDeployments.Count(), "expecting 1 deployment");

                    // note: job should have one deployment, if we can find an existing session,
                    // but there will be no more deployment directory.
                    var allDepl = WorkflowMgmTestUtils.GetAllDeployments(basename + "*");
                    Assert.AreEqual(0, allDepl.Length, $"expecting 0 deployment directories, but got: {allDepl.Length}: " + allDepl.ToConcatString("", ", ", ";"));

                    Assert.AreEqual(sessionGuid, J.LatestSession.ID, "expecting session ID known from first run");
                    Assert.IsTrue(J.LatestSession.SuccessfulTermination, "session is expected to be successful");

                    Assert.AreEqual(1, BoSSSshell.WorkflowMgm.Grids.Length, "expecting exactly one grid");

                    Console.WriteLine("done run 3. =============================");
                }
            } finally {
                TestMutex.ReleaseMutex();
            }
            Console.WriteLine("success!");
        }
    }
}
