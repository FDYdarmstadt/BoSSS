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

using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Reflection;
using BoSSS.Platform;

namespace BoSSS.Application.BoSSSpad {
    
    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation using the mini batch processor, i.e. the local computer,
    /// see <see cref="MiniBatchProcessor.Client"/>.
    /// </summary>
    public class MiniBatchProcessorClient : BatchProcessorClient {

        /// <summary>
        /// Path to standard output file, if present - otherwise null.
        /// </summary>
        public override string GetStdoutFile(Job myJob) {
            var AllProblems = FilterJobData(myJob);

            if (AllProblems.Length > 0) {
                return MiniBatchProcessor.Client.GetStdoutFile(AllProblems[0].ID);
            } else {
                return null;
            }
        }

        /// <summary>
        /// Path to standard error file, if present - otherwise null.
        /// </summary>
        public override string GetStderrFile(Job myJob) {
            var AllProblems = FilterJobData(myJob);

            if (AllProblems.Length > 0) {
                return MiniBatchProcessor.Client.GetStderrFile(AllProblems[0].ID);
            } else {
                return null;
            }
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="DeployDir">
        /// If null, a default choice is made.
        /// </param>
        public MiniBatchProcessorClient(string DeployDir = null) {
            var userDir = BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath();
            if (userDir == null || userDir.Length <= 0 || !Directory.Exists(userDir)) {
                throw new ApplicationException("Unable to create local machine batch, user settings path ('.BoSSS' - directory) does not exist or unable to find.");
            }

            //base.DeployDirectory = Path.Combine(userDir, "batch");

            if (string.IsNullOrWhiteSpace(DeployDir)) {
                string localAppData = System.Environment.GetEnvironmentVariable("LOCALAPPDATA")
                    ?? System.Environment.GetEnvironmentVariable("HOME");

                this.DeploymentBaseDirectory = Path.Combine(localAppData, "BoSSS-LocalJobs");

            } else {
                this.DeploymentBaseDirectory = DeployDir;
            }

            if (!Directory.Exists(this.DeploymentBaseDirectory))
                throw new IOException("Deploy directory '" + this.DeploymentBaseDirectory + "' does not exist.");
        }

        
        private string GetFullJobName(Job myJob) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(InteractiveShell.WorkflowMgm.CurrentProject)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }
            return PrjName + "__" + myJob.Name;
        }

        /// <summary>
        /// See <see cref="BatchProcessorClient.EvaluateStatus(Job, out int, out bool, out bool, out bool, out string)"/>.  
        /// </summary>
        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            if (!object.ReferenceEquals(this, myJob.AssignedBatchProc))
                throw new ArgumentException("Why you ask me?");
            string FullName = GetFullJobName(myJob);
            MiniBatchProcessor.JobData[] AllProblems = FilterJobData(myJob);
            MiniBatchProcessor.JobData JD = null;
            if (AllProblems.Length > 0) {
                if (myJob.BatchProcessorIdentifierToken == null) {
                    JD = AllProblems.ElementAtMax(jd => jd.SubmitTime);
                } else {
                    int idSearch = (int)(myJob.BatchProcessorIdentifierToken);
                    JD = AllProblems.SingleOrDefault(jobDat => jobDat.ID == idSearch);
                }
            }
            SubmitCount = AllProblems.Length;
            if (AllProblems.Length <= 0 || JD == null) {
                // we know nothing
                isRunning = false;
                isFailed = false;
                wasSuccessful = false;
                DeployDir = null;
                return;
            }
                
            int ExitCode;
            var mbpStatus = MiniBatchProcessor.ClientAndServer.GetStatusFromID(JD.ID, out ExitCode);
            DeployDir = JD.ExeDir;

            switch(mbpStatus) {
                case MiniBatchProcessor.JobStatus.Queued:
                // we know nothing
                isRunning = false;
                isFailed = false;
                wasSuccessful = false;
                return;

                case MiniBatchProcessor.JobStatus.Finished:
                // we know nothing
                isRunning = false;
                wasSuccessful = (ExitCode == 0);
                isFailed = !wasSuccessful;
                return;

                case MiniBatchProcessor.JobStatus.Working:
                // we know nothing
                isRunning = true;
                isFailed = false;
                wasSuccessful = false;
                return;

                case MiniBatchProcessor.JobStatus.Undefined:
                // we know nothing
                isRunning = false;
                isFailed = false;
                wasSuccessful = false;
                return;

                default:
                throw new NotImplementedException();
            }



        }

        private MiniBatchProcessor.JobData[] FilterJobData(Job myJob) {
            string FullName = GetFullJobName(myJob);
            var all = new List<MiniBatchProcessor.JobData>();
            int idSearch = -1;
            if(myJob.BatchProcessorIdentifierToken != null) {
                try {
                    idSearch = (int)myJob.BatchProcessorIdentifierToken;
                } catch(Exception) {

                }
            }

            //var All = MiniBatchProcessor.ClientAndServer.AllJobs.Where(mbpJob => mbpJob.Name.Equals(FullName)).ToArray();
            foreach(var mbpJob in MiniBatchProcessor.ClientAndServer.AllJobs) {
                if (!mbpJob.Name.Equals(FullName))
                    continue;

                if(idSearch >= 0) {
                    if (mbpJob.ID != idSearch)
                        continue;
                }

                all.Add(mbpJob);
            }

            return all.OrderByDescending(job => job.SubmitTime).ToArray();
        }


        /// <summary>
        /// See <see cref="MiniBatchProcessorClient.Submit(Job)"/>.
        /// </summary>
        public override object Submit(Job myJob) {
            string FullName = GetFullJobName(myJob);
            //var AllProblems = FilterJobData(myJob);
            //if (AllProblems.Length > 0) {
            //    throw new ApplicationException("There are already " + AllProblems.Length + " jobs with the name '" + FullName + "' in the MiniBatchProcessor. Since the job name must be unique, we cannot submit - try another project name.");
            //}

            var JD = new MiniBatchProcessor.JobData() {
                Name = FullName,
                NoOfProcs = myJob.NumberOfMPIProcs,
                ExeDir = myJob.DeploymentDirectory,
                exefile = Path.GetFileName(myJob.EntryAssembly.Location),
                Arguments = myJob.CommandLineArguments,
                EnvVars = myJob.EnvironmentVars.Select(kv => new Tuple<string, string>(kv.Key, kv.Value)).ToArray(),
                UseComputeNodesExclusive = myJob.UseComputeNodesExclusive
            };

            int id = MiniBatchProcessor.Client.SubmitJob(JD);
            return id;
        }
    }
}
