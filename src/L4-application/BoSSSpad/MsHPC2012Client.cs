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
using ilPSP.Tracing;
using Microsoft.Hpc.Scheduler;
using Microsoft.Hpc.Scheduler.Properties;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation which uses a Microsoft HPC 2012 server (or later).
    /// </summary>
    [DataContract]
    [Serializable]
    public class MsHPC2012Client : BatchProcessorClient {
        /// <summary>
        /// Empty Constructor for de-serialization
        /// </summary>
        private MsHPC2012Client() : base() {
            //Console.WriteLine("MsHPC2012Client: empty ctor");

            if ( System.Environment.OSVersion.Platform != PlatformID.Win32NT ) {
                throw new NotSupportedException($"The {typeof(MsHPC2012Client).Name} is only supported on MS Windows, but your current platform seems to be {System.Environment.OSVersion.Platform}.");
            }

            base.RuntimeLocation = "win\\amd64";

            //m_AdditionalEnvironmentVars.Add("OMP_PROC_BIND", "spread");
        }


        /// <summary>
        /// Since this is specific for MS Windows systems, it defaults to `win\amd64`
        /// </summary>
        public override string RuntimeLocation {
            get {
                if ( base.RuntimeLocation != null )
                    return base.RuntimeLocation;
                else
                    return "win\\amd64";
            }
            set => base.RuntimeLocation = value;
        }

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
        /// <param name="ComputeNodes">
        /// </param>
        /// <param name="DeployRuntime">
        /// See <see cref="BatchProcessorClient.DeployRuntime"/>.
        /// </param>
        public MsHPC2012Client(string DeploymentBaseDirectory, string ServerName, string Username = null, string[] ComputeNodes = null, bool DeployRuntime = true) : base() {
            if ( System.Environment.OSVersion.Platform != PlatformID.Win32NT ) {
                throw new NotSupportedException($"The {typeof(MsHPC2012Client).Name} is only supported on MS Windows, but your current platform seems to be {System.Environment.OSVersion.Platform}.");
            }

            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            base.DeployRuntime = DeployRuntime;
            base.RuntimeLocation = "win\\amd64";

            this.Username = Username;
            this.ComputeNodes = ComputeNodes;
            this.ServerName = ServerName;

            if ( !Directory.Exists(base.DeploymentBaseDirectory) )
                Directory.CreateDirectory(base.DeploymentBaseDirectory);
#pragma warning disable CA1416
            if ( this.Username == null )
                this.Username = System.Security.Principal.WindowsIdentity.GetCurrent().Name;
#pragma warning restore CA1461

        }

        IScheduler GetInstance() {
            if ( _Scheduler == null ) {
                _Scheduler = new Scheduler();
                _Scheduler.Connect(ServerName);
            }
            return _Scheduler;
        }
        //[NonSerialized]
        //IScheduler m__scheduler;

        /// <summary>
        /// Active Directory user name used on HPC cluster
        /// </summary>
        [DataMember]
        public string Username;

        /// <summary>
        /// Active Directory password used on HPC cluster
        /// </summary>
        [DataMember]
        public string Password;


        /// <summary>
        /// Additional number of cores (for all jobs with more than one MPI rank) which are allocated for 'service', independent of the MPI Size.
        /// </summary>
        [DataMember]
        public int NumOfAdditionalServiceCores = 0;

        /// <summary>
        /// Additional number of cores (for all jobs with only one MPI rank) which are allocated for 'service', independent of the MPI Size;
        /// <see cref="NumOfAdditionalServiceCores"/>.
        /// </summary>
        [DataMember]
        public int NumOfAdditionalServiceCoresMPISerial = 0;


        /// <summary>
        /// Additional number of cores which are allocated for 'service';
        /// <see cref="NumOfAdditionalServiceCores"/>.
        /// </summary>
        [DataMember]
        public int NumOfServiceCoresPerMPIprocess = 0;


        /// <summary>
        /// Active directory computer name of head node
        /// </summary>
        [DataMember]
        public string ServerName;

        /// <summary>
        /// optional: a list of compute node on which some job should run
        /// </summary>
        [DataMember]
        public string[] ComputeNodes;

        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public JobPriority DefaultJobPriority = JobPriority.Normal;

        /// <summary>
        /// Jobs are forced to run on a single node.
        /// </summary>
        [DataMember]
        public bool SingleNode = true;

        /// <summary>
        /// HPC Scheduler Object from Microsoft.Hpc
        /// </summary>
        [NonSerialized]
        protected IScheduler _Scheduler;

        ///// <summary>
        ///// Access to the Microsoft HPC job scheduler interface.
        ///// </summary>
        //IScheduler Scheduler {
        //    get {
        //        if (m__scheduler == null) {
        //            m__scheduler = new Scheduler();
        //            m__scheduler.Connect(ServerName);
        //        }
        //        return m__scheduler;
        //    }
        //}


        /// <summary>
        /// Job status.
        /// </summary>
        public override (BoSSSpad.JobStatus, int? ExitCode) EvaluateStatus(string idToken, object optInfo, string DeployDir) {
            using ( var tr = new FuncTrace() ) {
                int id = int.Parse(idToken);
                var (state, exitCode) = GetStatus(id);

                // TODO: TS Hier rennt der in JobState Failed rein beim Testrunner, das muss gedebuggt werden
                switch ( state ) {
                    case JobState.Configuring:
                    case JobState.Submitted:
                    case JobState.Validating:
                    case JobState.ExternalValidation:
                    case JobState.Queued:
                        return (JobStatus.PendingInExecutionQueue, null);

                    case JobState.Running:
                    case JobState.Finishing:
                    case JobState.Canceling:
                        return (JobStatus.InProgress, null);

                    case JobState.Finished:
                        return (JobStatus.FinishedSuccessful, exitCode);

                    case JobState.Failed:
                    case JobState.Canceled:
                        Console.WriteLine($" ------------ MSHPC FailedOrCanceled; original " + state);
                        return (JobStatus.FailedOrCanceled, exitCode);

                    case JobState.All:
                        return (JobStatus.Unknown, exitCode);

                    default:
                        throw new NotImplementedException("Unknown job state: " + state);
                }
            }
        }

        /// <summary>
        /// get MS HPC status (<see cref="JobState"/>) from Scheduler Job ID;
        /// will be translated to <see cref="BoSSSpad.JobStatus"/> by some other method.
        /// </summary>
        public (JobState s, int? exitCode) GetStatus(int id) {
            using ( var tr = new FuncTrace() ) {
                tr.Info($"Trying to get status for job {id} from scheduler {this.ServerName}");

                if ( this.ServerName.IsEmptyOrWhite() )
                    throw new IOException("'ServerName' for MS HPC scheduler is empty or white");

                var job = _Scheduler.OpenJob(id);
                job.Refresh();
                JobState state = job.State;
                int exitCode = int.MinValue;
                var tasks = job.GetTaskList(null, null, false);
                foreach ( ISchedulerTask task in tasks ) {
                    exitCode = task.ExitCode;
                }

                return (state, exitCode);
            }
        }

        /// <summary>
        /// Path to standard error file.
        /// </summary>
        public override string GetStderrFile(string idToken, string DeployDir) {
            if ( idToken.IsEmptyOrWhite() || DeployDir.IsEmptyOrWhite() )
                return null;
            string fp = Path.Combine(DeployDir, "stderr.txt");
            return fp;
        }
        /// <summary>
        /// Path to standard output file.
        /// </summary>
        public override string GetStdoutFile(string idToken, string DeployDir) {
            if ( idToken.IsEmptyOrWhite() || DeployDir.IsEmptyOrWhite() )
                return null;
            string fp = Path.Combine(DeployDir, "stdout.txt");
            return fp;

        }

        private int GetNumberOfCoresForJobDescription(Job description) {
            int MPISz = description.NumberOfMPIProcs;
            int NumberOfCores = MPISz * description.NumberOfThreads + MPISz * this.NumOfServiceCoresPerMPIprocess + (MPISz > 1 ? this.NumOfAdditionalServiceCores : this.NumOfAdditionalServiceCoresMPISerial);
            return NumberOfCores;
        }

        private ISchedulerTask CreateTaskFromJobDescription(Job description, string DeploymentDirectory, ISchedulerJob job) {
            int NumberOfCores = GetNumberOfCoresForJobDescription(description);
            var task = job.CreateTask();
            task.MaximumNumberOfCores = NumberOfCores;
            task.MinimumNumberOfCores = NumberOfCores;
            task.WorkDirectory = DeploymentDirectory;

            using ( var str = new StringWriter() ) {
                str.Write("mpiexec ");
                if ( !base.DotnetRuntime.IsEmptyOrWhite() )
                    str.Write(base.DotnetRuntime + " ");
                str.Write(description.EntryAssemblyName);
                foreach ( string arg in description.CommandLineArguments ) {
                    str.Write(" ");
                    str.Write(arg);
                }

                task.CommandLine = str.ToString();
            }
            foreach ( var kv in description.EnvironmentVars ) {
                string name = kv.Key;
                string valu = kv.Value;
                task.SetEnvironmentVariable(name, valu);
            }

            task.StdOutFilePath = Path.Combine(DeploymentDirectory, "stdout.txt");
            task.StdErrFilePath = Path.Combine(DeploymentDirectory, "stderr.txt");

            if ( ComputeNodes != null ) {
                foreach ( string node in ComputeNodes )
                    job.RequestedNodes.Add(node);
            }
            return task;
        }

        /// <summary>
        /// Submits the job to the Microsoft HPC server.
        /// </summary>
        public override (string id, object optJobObj) Submit(Job description, string DeploymentDirectory) {
            using ( new FuncTrace() ) {
                int NumberOfCores = GetNumberOfCoresForJobDescription(description);
                string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
                IScheduler Scheduler = GetInstance();
                var job = Scheduler.CreateJob();
                job.Name = description.Name;
                job.Project = PrjName;
                job.MaximumNumberOfCores = NumberOfCores;
                job.MinimumNumberOfCores = NumberOfCores;
                job.SingleNode = this.SingleNode;
                job.Priority = this.DefaultJobPriority;
                job.UserName = Username;

                var task = CreateTaskFromJobDescription(description, DeploymentDirectory, job);

                job.AddTask(task);

                // Start the job.
                Scheduler.SubmitJob(job, Username != null ? Username : null, Password);
                job.Refresh();
                Console.WriteLine($"State of current Job :{job.State}");
                return (job.Id.ToString(), null);
            }
        }

        /// <summary>
        /// Cancels the job with the given id
        /// </summary>
        /// <param name="Id">The identifier for the job</param>
        /// <param name="message">The reason the job was cancelled</param>
        public override void Cancel(string Id, string message) {
            int id = int.Parse(Id);
            IScheduler Scheduler = GetInstance();
            Scheduler.CancelJob(id, message);
        }
        /// <summary>
        /// publicly available wrapper, used for extension methods to <see cref="MsHPC2012Client"/>
        /// </summary>
        /// <returns></returns>
        public ISchedulerJob SubmitJobs(IEnumerable<Job> jobs, string DeploymentDirectory, int NumberOfCores = 64) {
            string PrjName = BoSSSshell.WorkflowMgm.CurrentProject;
            IScheduler Scheduler = GetInstance();
            var job = Scheduler.CreateJob();
            job.Name = PrjName;
            job.Project = PrjName;
            job.MaximumNumberOfCores = NumberOfCores;
            job.MinimumNumberOfCores = jobs.Select(e => e.NumberOfMPIProcs).Max();
            job.SingleNode = this.SingleNode;
            job.Priority = this.DefaultJobPriority;

            job.UserName = Username;

            foreach ( var description in jobs ) {
                var task = CreateTaskFromJobDescription(description, DeploymentDirectory, job);
                job.AddTask(task);
            }

            Scheduler.SubmitJob(job, Username != null ? Username : null, Password);
            return job;
        }

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            string NameString = "";
            if ( !base.Name.IsEmptyOrWhite() )
                NameString = " " + base.Name + " ";

            return $"MS HPC client {NameString}@{this.ServerName}, @{this.DeploymentBaseDirectory}";
        }

        [NonSerialized]
        Dictionary<string, string> m_AdditionalEnvironmentVars = new Dictionary<string, string>();

        /// <summary>
        /// Additional environment variables for the process. 
        /// </summary>
        [DataMember]
        public IDictionary<string, string> AdditionalEnvironmentVars {
            get {
                return m_AdditionalEnvironmentVars;
            }
        }
    }
}

