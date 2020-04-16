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
using Microsoft.Hpc.Scheduler;
using System.IO;
using Microsoft.Hpc.Scheduler.Properties;
using BoSSS.Platform;
using ilPSP;
using System.Runtime.Serialization;

namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation which uses a Microsoft HPC 2012 server.
    /// </summary>
    [DataContract]
    [Serializable]
    public class MsHPC2012Client : BatchProcessorClient {
        /*
        /// <summary>
        /// Configuration options specific to the <see cref="MiniBatchProcessorClient"/>
        /// </summary>
        [Serializable]
        public new class Config : BatchProcessorClient.Config {

            /// <summary>
            /// %
            /// </summary>
            public string ServerName;
            
            /// <summary>
            /// %
            /// </summary>
            public string Username;
            
            /// <summary>
            /// %
            /// </summary>
            public string Password;
            
            /// <summary>
            /// %
            /// </summary>
            public string[] ComputeNodes = null;

            /// <summary>
            /// %
            /// </summary>
            public override BatchProcessorClient Instance() {
                return new MsHPC2012Client(
                    base.DeploymentBaseDirectory,
                    ServerName,
                    Username,
                    Password,
                    ComputeNodes,
                    base.DeployRuntime);
            }
        }
        
        /// <summary>
        /// .
        /// </summary>
        public override BatchProcessorClient.Config GetConfig() {
            return new MsHPC2012Client.Config() {
                DeploymentBaseDirectory = this.DeploymentBaseDirectory,
                DeployRuntime = this.DeployRuntime,
                ComputeNodes = this.m_ComputeNodes.CloneAs(),
                Username = this.m_Username,
                Password = this.m_Password,
                ServerName = this.m_ServerName
            };
        }
        */

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="DeploymentBaseDirectory">
        /// A directory location which must be accessible from both, the HPC server as well as the local machine.
        /// </param>
        /// <param name="ServerName">
        /// Name of the HPC server.
        /// </param>
        /// <param name="Username">
        /// Can be null for the local user.
        /// </param>
        /// <param name="Password">
        /// Password for user <paramref name="Username"/>, can be null if the user is the local user.
        /// </param>
        /// <param name="ComputeNodes">
        /// </param>
        /// <param name="DeployRuntime">
        /// See <see cref="BatchProcessorClient.DeployRuntime"/>.
        /// </param>
        public MsHPC2012Client(string DeploymentBaseDirectory, string ServerName, string Username = null, string Password = null, string[] ComputeNodes = null, bool DeployRuntime = true) {
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            base.DeployRuntime = DeployRuntime;


            this.Username = Username;
            this.Password = Password;
            this.ComputeNodes = ComputeNodes;
            this.ServerName = ServerName;

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            if (this.Username == null)
                this.Username = System.Security.Principal.WindowsIdentity.GetCurrent().Name;

           
        }

        [NonSerialized]
        IScheduler m__scheduler;

        [DataMember]
        string Username;

        [DataMember]
        string Password;

        [DataMember]
        string ServerName;

        [DataMember]
        string[] ComputeNodes;

        /// <summary>
        /// Access to the Microsoft HPC job scheduler interface.
        /// </summary>
        IScheduler Scheduler {
            get {
                if (m__scheduler == null) {
                    m__scheduler = new Scheduler();
                    m__scheduler.Connect(ServerName);
                }
                return m__scheduler;
            }
        }


        /// <summary>
        /// Job status.
        /// </summary>
        public override void EvaluateStatus(string idToken, string DeployDir, out bool isRunning, out bool isTerminated, out int ExitCode) {

            int id = int.Parse(idToken);

            

            List<SchedulerJob> allFoundJobs = new List<SchedulerJob>();

            ISchedulerCollection allJobs = Scheduler.GetJobList(null, null);
            foreach (SchedulerJob sJob in allJobs) {
                if(sJob.Id != id)
                    continue;
                allFoundJobs.Add(sJob);
            }

            if (allFoundJobs.Count <= 0) {
                // some weird state
                isRunning = false;
                isTerminated = false;
                ExitCode = int.MinValue;
                return;
            }

            SchedulerJob JD = allFoundJobs.ElementAtMax(MsHpcJob => MsHpcJob.SubmitTime);
            
            ISchedulerCollection tasks = JD.GetTaskList(null, null, false);
            ExitCode = int.MinValue;
            foreach (ISchedulerTask t in tasks) {
                DeployDir = t.WorkDirectory;
                ExitCode = t.ExitCode;
            }
            
            switch (JD.State) {
                case JobState.Configuring:
                case JobState.Submitted:
                case JobState.Validating:
                case JobState.ExternalValidation:
                case JobState.Queued:
                isRunning = false;
                isTerminated = false;
                break;

                case JobState.Running:
                case JobState.Finishing:
                isRunning = true;
                isTerminated = false;
                break;

                case JobState.Finished:
                isRunning = false;
                isTerminated = true;
                break;

                case JobState.Failed:
                case JobState.Canceled:
                case JobState.Canceling:
                isRunning = false;
                isTerminated = true;
                break;

                default:
                throw new NotImplementedException("Unknown job state: " + JD.State);
            }

            
        }


        /// <summary>
        /// Path to standard error file.
        /// </summary>
        public override string GetStderrFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stderr.txt");
            return fp;
        }
        /// <summary>
        /// Path to standard output file.
        /// </summary>
        public override string GetStdoutFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stdout.txt");
            return fp;
            
        }

        /// <summary>
        /// Submit the job to the Microsoft HPC server.
        /// </summary>
        public override string Submit(Job myJob) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            
            //Console.WriteLine("MsHPC2012Client: submitting job...");
            //Console.Out.Flush();


            ISchedulerJob job = null;
            ISchedulerTask task = null;

            // Create a job and add a task to the job.
            job = Scheduler.CreateJob();
            
            job.Name = myJob.Name;
            job.Project = PrjName;
            job.MaximumNumberOfCores = myJob.NumberOfMPIProcs;
            job.MinimumNumberOfCores = myJob.NumberOfMPIProcs;

            job.UserName = Username;
            
            task = job.CreateTask();
            task.MaximumNumberOfCores = myJob.NumberOfMPIProcs;
            task.MinimumNumberOfCores = myJob.NumberOfMPIProcs;

            task.WorkDirectory = myJob.DeploymentDirectory;

            using(var str = new StringWriter()) {
                str.Write("mpiexec ");
                str.Write(Path.GetFileName(myJob.EntryAssembly.Location));
                foreach (string arg in myJob.CommandLineArguments) {
                    str.Write(" ");
                    str.Write(arg);
                }

                task.CommandLine = str.ToString();
            }
            foreach( var kv in myJob.EnvironmentVars) {
                string name = kv.Key;
                string valu = kv.Value;
                task.SetEnvironmentVariable(name, valu);
            }

            task.StdOutFilePath = Path.Combine(myJob.DeploymentDirectory, "stdout.txt");
            task.StdErrFilePath = Path.Combine(myJob.DeploymentDirectory, "stderr.txt");

            if(ComputeNodes != null) {
                foreach(string node in ComputeNodes)
                    job.RequestedNodes.Add(node);
            }


            job.AddTask(task);

            // Start the job.
            Scheduler.SubmitJob(job, Username != null ? Username : null, Password);

            return job.Id.ToString();
        }
    }
}
