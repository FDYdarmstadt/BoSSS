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
using System.Runtime.Serialization;
using ilPSP.Tracing;
using ilPSP.Utils;

namespace BoSSS.Application.BoSSSpad {
    
    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation using the mini batch processor, i.e. the local computer,
    /// see <see cref="MiniBatchProcessor.Client"/>.
    /// </summary>
    [DataContract]
    public class MiniBatchProcessorClient : BatchProcessorClient {


        /// <summary>
        /// Optional override for <see cref="MiniBatchProcessor.Configuration.BatchInstructionDir"/>
        /// </summary>
        [DataMember]
        public string BatchInstructionDir;

        /// <summary>
        /// Empty constructor for de-serialization
        /// </summary>
        private MiniBatchProcessorClient() : base() {
        }

        [NonSerialized]
        MiniBatchProcessor.Client m_Clint;

        /// <summary>
        /// %
        /// </summary>
        MiniBatchProcessor.Client Clint {
            get {
                if(m_Clint == null) {
                    m_Clint = new MiniBatchProcessor.Client(BatchInstructionDir);
                }
                return m_Clint;
            }
        }


        /// <summary>
        /// Path to standard output file, if present - otherwise null.
        /// </summary>
        public override string GetStdoutFile(string idToken, string DeployDir) {
            return Clint.GetStdoutFile(int.Parse(idToken));
        }

        /// <summary>
        /// Path to standard error file, if present - otherwise null.
        /// </summary>
        public override string GetStderrFile(string idToken, string DeployDir) {
            return Clint.GetStderrFile(int.Parse(idToken));
        }

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="DeployDir">
        /// If null, a default choice is made.
        /// </param>
        public MiniBatchProcessorClient(string DeployDir = null) : base() {
            var userDir = BoSSS.Foundation.IO.Utils.GetBoSSSUserSettingsPath();
            if(userDir == null || userDir.Length <= 0 || !Directory.Exists(userDir)) {
                throw new ApplicationException("Unable to create local machine batch, user settings path ('.BoSSS' - directory) does not exist or unable to find.");
            }

            if(System.OperatingSystem.IsWindows())
                base.RuntimeLocation = "win\\amd64";
            else
                base.RuntimeLocation = "linux\\amd64-openmpi";

            //base.DeployDirectory = Path.Combine(userDir, "batch");

            if (string.IsNullOrWhiteSpace(DeployDir)) {
                string localAppData = System.Environment.GetEnvironmentVariable("LOCALAPPDATA")
                    ?? System.Environment.GetEnvironmentVariable("HOME");

                this.DeploymentBaseDirectory = Path.Combine(localAppData, "BoSSS-LocalJobs");
                if(!Directory.Exists(this.DeploymentBaseDirectory)) {
                    Directory.CreateDirectory(this.DeploymentBaseDirectory);
                }
            } else {
                this.DeploymentBaseDirectory = DeployDir;
            }

            if(!Directory.Exists(this.DeploymentBaseDirectory))
                throw new IOException("Deploy directory '" + this.DeploymentBaseDirectory + "' does not exist.");

            {
                string localUserDir = System.Environment.GetEnvironmentVariable("USERPROFILE") ?? System.Environment.GetEnvironmentVariable("HOME");
                if(localUserDir != null)
                    base.AllowedDatabasesPaths.Add(new AllowedDatabasesPair(localUserDir, null));
                if(Path.DirectorySeparatorChar == '\\')
                    base.AllowedDatabasesPaths.Add(new AllowedDatabasesPair("C:\\", null));
                else
                    base.AllowedDatabasesPaths.Add(new AllowedDatabasesPair("/", null));
            }
        }

        
        private string GetFullJobName(Job myJob) {
            string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
            if (string.IsNullOrWhiteSpace(BoSSSshell.WorkflowMgm.CurrentProject)) {
                throw new NotSupportedException("Project management not initialized - set project name (try e.g. 'WorkflowMgm.CurrentProject = \"BlaBla\"').");
            }
            return PrjName + "__" + myJob.Name;
        }

        /// <summary>
        /// See <see cref="BatchProcessorClient.EvaluateStatus"/>.  
        /// </summary>
        public override (BoSSSpad.JobStatus,int? ExitCode) EvaluateStatus(string idToken, object optInfo, string DeployDir) { 
        //public override void EvaluateStatus(string idToken, object optInfo, string DeployDir, out bool isRunning, out bool isTerminated, out int ExitCode) {
            using (new FuncTrace()) {
                

                int ID = int.Parse(idToken);
                var mbpStatus = Clint.GetStatusFromID(ID);
                int ExitCode = mbpStatus.ExitCode;


                switch (mbpStatus.stat) {
                    case MiniBatchProcessor.JobStatus.Queued:
                        // we know nothing
                        return (BoSSSpad.JobStatus.PendingInExecutionQueue, null);

                    case MiniBatchProcessor.JobStatus.Finished:
                        // we know nothing
                        return (ExitCode == 0 ? BoSSSpad.JobStatus.FinishedSuccessful : BoSSSpad.JobStatus.FailedOrCanceled, ExitCode);

                    case MiniBatchProcessor.JobStatus.Working:
                        // we know nothing
                        return (BoSSSpad.JobStatus.InProgress, null);

                    case MiniBatchProcessor.JobStatus.Undefined:
                        // we know nothing
                        return (BoSSSpad.JobStatus.Unknown, null);

                    default:
                        throw new NotImplementedException();
                }

            }

        }

        /*
        private MiniBatchProcessor.JobData FilterJobData(Job myJob) {
            int idSearch;
            try {
                idSearch = int.Parse(myJob.BatchProcessorIdentifierToken);
            } catch(Exception) {
                return null;
            }

            return Clint.AllJobs.FirstOrDefault(jd => jd.ID == idSearch);
        }
        */

        /// <summary>
        /// See <see cref="BatchProcessorClient.Submit"/>. 
        /// </summary>
        public override (string id, object optJobObj) Submit(Job myJob, string DeploymentDirectory) {
            var started = MiniBatchProcessor.Server.StartIfNotRunning(RunExternal: true);
            if(started) {
                Console.WriteLine("Warning: MiniBatchProcessor server was not running, started by job activation; it might be beneficial to start `MiniBatchProcessor.dll` externally, for the future.");
            }


            string FullName = GetFullJobName(myJob);
            //var AllProblems = FilterJobData(myJob);
            //if (AllProblems.Length > 0) {
            //    throw new ApplicationException("There are already " + AllProblems.Length + " jobs with the name '" + FullName + "' in the MiniBatchProcessor. Since the job name must be unique, we cannot submit - try another project name.");
            //}

            var JD = new MiniBatchProcessor.JobData() {
                Name = FullName,
                NoOfProcs = myJob.NumberOfMPIProcs,
                ExeDir = DeploymentDirectory,
                exefile = base.DotnetRuntime,
                Arguments = ArrayTools.Cat(new[] { myJob.EntryAssemblyName }, myJob.CommandLineArguments),
                EnvVars = myJob.EnvironmentVars.Select(kv => new Tuple<string, string>(kv.Key, kv.Value)).ToArray(),
                UseComputeNodesExclusive = myJob.UseComputeNodesExclusive
            };

            int id = Clint.SubmitJob(JD);
            return (id.ToString(), JD);
        }

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            string NameString = "";
            if(!base.Name.IsEmptyOrWhite())
                NameString = " " + base.Name + " ";

            return $"MiniBatchProcessor client {NameString}@{this.DeploymentBaseDirectory}";
        }

    }
}
