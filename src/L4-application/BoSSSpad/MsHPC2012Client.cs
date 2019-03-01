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


namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation which uses a Microsoft HPC 2012 server.
    /// </summary>
    public class MsHPC2012Client : BatchProcessorClient {

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

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            m_Username = Username;
            m_Password = Password;
            m_ComputeNodes = ComputeNodes;

            if (m_Username == null)
                m_Username = System.Security.Principal.WindowsIdentity.GetCurrent().Name;

            m_scheduler = new Scheduler();
            m_scheduler.Connect(ServerName);
        }

        IScheduler m_scheduler;
        string m_Username;
        string m_Password;
        string[] m_ComputeNodes;

        /// <summary>
        /// Access to the Microsoft HPC job scheduler interface.
        /// </summary>
        IScheduler Scheduler {
            get {
                return m_scheduler;
            }
        }


        /// <summary>
        /// Job status.
        /// </summary>
        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            DeployDir = null;

            List<SchedulerJob> allFoundJobs = new List<SchedulerJob>();
            ISchedulerCollection allJobs = m_scheduler.GetJobList(null, null);
            foreach (SchedulerJob sJob in allJobs) {
                if (!sJob.Project.Equals(PrjName))
                    continue;
                if (!sJob.Name.Equals(myJob.Name))
                    continue;
                if (!sJob.UserName.Equals(m_Username, StringComparison.OrdinalIgnoreCase) && !sJob.Owner.Equals(m_Username, StringComparison.OrdinalIgnoreCase))
                    continue;
                if (!sJob.UserName.Equals(sJob.Owner, StringComparison.OrdinalIgnoreCase))
                    // ignore weird stuff
                    continue;

                allFoundJobs.Add(sJob);
            }

            SubmitCount = allFoundJobs.Count;

            if (allFoundJobs.Count <= 0) {
                
                isRunning = false;
                wasSuccessful = false;
                isFailed = false;
                
                return;
            }

            //if (allFoundJobs.Count > 1) {
            //    throw new ApplicationException(string.Format("Found {0} Microsoft-HPC-jobs with matching criteria (project is '{1}', name is '{2}', user name is '{3}'). Unable to determine which one correlates to given meta-schedule-job.",allFoundJobs.Count, PrjName, myJob.Name, m_Username));
            //}

            //allFoundJobs[0].SubmitTime 
            SchedulerJob JD = allFoundJobs.ElementAtMax(MsHpcJob => MsHpcJob.SubmitTime);
            
            ISchedulerCollection tasks = JD.GetTaskList(null, null, false);
            foreach (ISchedulerTask t in tasks) {
                DeployDir = t.WorkDirectory;
            }

            switch (JD.State) {
                case JobState.Configuring:
                case JobState.Submitted:
                case JobState.Validating:
                case JobState.ExternalValidation:
                case JobState.Queued:
                isRunning = false;
                wasSuccessful = false;
                isFailed = false;
                break;

                case JobState.Running:
                case JobState.Finishing:
                isRunning = true;
                wasSuccessful = false;
                isFailed = false;
                break;

                case JobState.Finished:
                isRunning = false;
                wasSuccessful = true;
                isFailed = false;
                break;

                case JobState.Failed:
                case JobState.Canceled:
                case JobState.Canceling:
                wasSuccessful = false;
                isFailed = true;
                isRunning = false;
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
        public override object Submit(Job myJob) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            
            ISchedulerJob job = null;
            ISchedulerTask task = null;

            // Create a job and add a task to the job.
            job = m_scheduler.CreateJob();
            job.Name = myJob.Name;
            job.Project = PrjName;
            job.MaximumNumberOfCores = myJob.NumberOfMPIProcs;
            job.MinimumNumberOfCores = myJob.NumberOfMPIProcs;

            job.UserName = m_Username;
            
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

            if(m_ComputeNodes != null) {
                foreach(string node in m_ComputeNodes)
                    job.RequestedNodes.Add(node);
            }


            job.AddTask(task);


            // Start the job.
            m_scheduler.SubmitJob(job, m_Password != null ? m_Username : null, m_Password);

            return job.Id;
        }
    }
}
