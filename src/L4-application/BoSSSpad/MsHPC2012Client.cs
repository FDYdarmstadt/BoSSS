﻿/* =======================================================================
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
using System.IO;
using BoSSS.Platform;
using ilPSP;
using System.Runtime.Serialization;
using ilPSP.Tracing;
using System.Diagnostics;
using System.Text;
using System.Threading;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Defines the priorities that you can specify for a job.
    /// </summary>
    public enum JobPriority
    {

        /// <summary>
        /// 
        /// </summary>
        Lowest = 0,

        /// <summary>
        /// 
        /// </summary>
        BelowNormal = 1,

        /// <summary>
        /// 
        /// </summary>
        Normal = 2,

        /// <summary>
        /// 
        /// </summary>
        AboveNormal = 3,

        /// <summary>
        /// 
        /// </summary>
        Highest = 4
    }




    /// <summary>
    /// A <see cref="BatchProcessorClient"/>-implementation which uses a Microsoft HPC 2012 server (or later).
    /// </summary>
    [DataContract]
    [Serializable]
    public class MsHPC2012Client : BatchProcessorClient {

        //
        // Summary:
        //     Defines the state of the job.
        enum JobState
        {
            //
            // Summary:
            //     The job is being configured. The application called the Microsoft.Hpc.Scheduler.IScheduler.CreateJob
            //     method to create the job but has not called the Microsoft.Hpc.Scheduler.IScheduler.AddJob(Microsoft.Hpc.Scheduler.ISchedulerJob)
            //     or Microsoft.Hpc.Scheduler.IScheduler.SubmitJob(Microsoft.Hpc.Scheduler.ISchedulerJob,System.String,System.String)
            //     method to add the job to the scheduler or submit the job to the scheduling queue.
            //     This enumeration member represents a value of 1.
            Configuring = 1,
            //
            // Summary:
            //     The job was submitted to the scheduling queue (see Microsoft.Hpc.Scheduler.IScheduler.SubmitJob(Microsoft.Hpc.Scheduler.ISchedulerJob,System.String,System.String)).
            //     This enumeration member represents a value of 2.
            Submitted = 2,
            //
            // Summary:
            //     The server is determining if the job can run. This enumeration member represents
            //     a value of 4.
            Validating = 4,
            //
            // Summary:
            //     A submission filter is determining if the job can run. For details, see the SubmissionFilterProgram
            //     cluster parameter in the Remarks section of Microsoft.Hpc.Scheduler.IScheduler.SetClusterParameter(System.String,System.String).
            //     This enumeration member represents a value of 8.
            ExternalValidation = 8,
            //
            // Summary:
            //     The job passed validation and was added to the scheduling queue. This enumeration
            //     member represents a value of 16.
            Queued = 16,
            //
            // Summary:
            //     The job is running. This enumeration member represents a value of 32.
            Running = 32,
            //
            // Summary:
            //     The server is cleaning up the resources that were allocated to the job. This
            //     enumeration member represents a value of 64.
            Finishing = 64,
            //
            // Summary:
            //     The job successfully finished (all the tasks in the job finished successfully).
            //     This enumeration member represents a value of 128.
            Finished = 128,
            //
            // Summary:
            //     One or more of the tasks in the job failed or a system error occurred on the
            //     compute node. To get a description of the error, access the Microsoft.Hpc.Scheduler.ISchedulerJob.ErrorMessage
            //     property. This enumeration member represents a value of 256.
            Failed = 256,
            //
            // Summary:
            //     The job was canceled (see Microsoft.Hpc.Scheduler.IScheduler.CancelJob(System.Int32,System.String)).
            //     If the caller provided the reason for canceling the job, then the Microsoft.Hpc.Scheduler.ISchedulerJob.ErrorMessage
            //     property will contain the reason. This enumeration member represents a value
            //     of 512.
            Canceled = 512,
            //
            // Summary:
            //     The job is being canceled. This enumeration member represents a value of 1024.
            Canceling = 1024,
            //
            // Summary:
            //     A mask used to indicate all states. This enumeration member represents a value
            //     of 2047.
            All = 2047
        }


        /// <summary>
        /// Empty Constructor for de-serialization
        /// </summary>
        private MsHPC2012Client() : base() {
            //Console.WriteLine("MsHPC2012Client: empty ctor");

            if(System.Environment.OSVersion.Platform != PlatformID.Win32NT) {
                throw new NotSupportedException($"The {typeof(MsHPC2012Client).Name} is only supported on MS Windows, but your current platform seems to be {System.Environment.OSVersion.Platform}.");
            }
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
        /// <param name="Password">
        /// Password for user <paramref name="Username"/>, can be null if the user is the local user.
        /// </param>
        /// <param name="ComputeNodes">
        /// </param>
        /// <param name="DeployRuntime">
        /// See <see cref="BatchProcessorClient.DeployRuntime"/>.
        /// </param>
        public MsHPC2012Client(string DeploymentBaseDirectory, string ServerName, string Username = null, string Password = null, string[] ComputeNodes = null, bool DeployRuntime = true) : base() {
            if(System.Environment.OSVersion.Platform != PlatformID.Win32NT) {
                throw new NotSupportedException($"The {typeof(MsHPC2012Client).Name} is only supported on MS Windows, but your current platform seems to be {System.Environment.OSVersion.Platform}.");
            }

            
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            base.DeployRuntime = DeployRuntime;


            this.Username = Username;
            this.Password = Password;
            this.ComputeNodes = ComputeNodes;
            this.ServerName = ServerName;

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);
#pragma warning disable CA1416
            if (this.Username == null)
                this.Username = System.Security.Principal.WindowsIdentity.GetCurrent().Name;
#pragma warning restore CA1461

        }

        //[NonSerialized]
        //IScheduler m__scheduler;

        /// <summary>
        /// Active Directory user name used on HPC cluster
        /// </summary>
        [DataMember]
        public string Username;

        /// <summary>
        /// Unsafely stored password
        /// </summary>
        [DataMember]
        public string Password;

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
            using(var tr = new FuncTrace()) {
                
                int id = int.Parse(idToken);
                var intStatus = GetStatus(id);

                switch (intStatus.s)
                {
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
                        //var retCode = (intStatus.exitCode == 0 ? JobStatus.FinishedSuccessful : JobStatus.FailedOrCanceled, intStatus.exitCode);
                        //if (retCode.Item1 != JobStatus.FinishedSuccessful)
                        //    Console.WriteLine($" ------------ MSHPC FailedOrCanceled; original " + JDstate + ", Exit code = " + ExitCode);
                        if (intStatus.exitCode != null && intStatus.exitCode == 0)
                            return (JobStatus.FinishedSuccessful, 0);
                        else
                            return (JobStatus.FailedOrCanceled, intStatus.exitCode);

                    case JobState.Failed:
                    case JobState.Canceled:
                        Console.WriteLine($" ------------ MSHPC FailedOrCanceled; original " + intStatus.s);
                        return (JobStatus.FailedOrCanceled, intStatus.exitCode);

                    default:
                        throw new NotImplementedException("Unknown job state: " + intStatus.s);
                }



                /*
                 * old code, using the .NET API
                 * 

                ISchedulerJob JD;
                //if (optInfo != null && optInfo is ISchedulerJob _JD) {
                //    JD = _JD;
                //} else {
                using(new BlockTrace("Scheduler.OpenJob", tr)) {
                    JD = Scheduler.OpenJob(id);
                }

               
                int ExitCode = int.MinValue;
                using(new BlockTrace("TASK_FILTERING", tr)) {
                    ISchedulerCollection tasks = JD.GetTaskList(null, null, false);
                    foreach(ISchedulerTask t in tasks) {
                        DeployDir = t.WorkDirectory;
                        ExitCode = t.ExitCode;
                    }
                }

                using(new BlockTrace("STATE_EVAL", tr)) {
                    var JDstate = JD.State;

                    switch(JDstate) {
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
                        var retCode = (ExitCode == 0 ? JobStatus.FinishedSuccessful : JobStatus.FailedOrCanceled, ExitCode);
                        if(retCode.Item1 != JobStatus.FinishedSuccessful)
                            Console.WriteLine($" ------------ MSHPC FailedOrCanceled; original " + JDstate + ", Exit code = " + ExitCode);
                        return retCode;

                        case JobState.Failed:
                        case JobState.Canceled:
                        Console.WriteLine($" ------------ MSHPC FailedOrCanceled; original " + JDstate);
                        return (JobStatus.FailedOrCanceled, ExitCode);

                        default:
                        throw new NotImplementedException("Unknown job state: " + JD.State);
                    }
                }

                */




            }
        }

        (JobState s, int? exitCode) GetStatus(int id) {
            // Get job status
            // ==============
            JobState state = JobState.All;
            {
                var args = $"view {id}";
                var Res = ExecuteProcess("job.exe", args, 60000);



                bool bfound = false;
                using(var StandardOutput = new StringReader(Res.stdOut)) {
                    for(string line = StandardOutput.ReadLine(); line != null; line = StandardOutput.ReadLine()) {
                        if(line.StartsWith("state", StringComparison.InvariantCultureIgnoreCase)) {
                            int iRes = line.IndexOf(':');
                            string RestLine = line.Substring(iRes + 1);
                            state = Enum.Parse<JobState>(RestLine);
                            bfound = true;
                            break;
                        } else {
                            continue;
                        }


                        //var parts = line.Split(new char[] { ' ', ' ', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                        //if (parts.Length < 3)
                        //    continue;
                        //if (parts[0].Equals("state", StringComparison.InvariantCultureIgnoreCase))
                        //{
                        //    state = Enum.Parse<JobState>(parts[2]);
                        //    break;
                        //}
                        //else
                        //{
                        //    continue;
                        //}
                    }
                }

                if(!bfound) {
                    throw new IOException("Unable to evalueate status of job " + id + System.Environment.NewLine + Res.stdOut + System.Environment.NewLine + Res.stdErr);
                }
            }

            // Get Exit code from task
            // =======================
            //
            // Note: in our simoplified interface, we only have one task per job

            int? exitcode = null;
            if(state == JobState.Canceled || state == JobState.Failed || state == JobState.Finished) {
                var args2 = $"listtasks {id}";
                var Res2 = ExecuteProcess("job.exe", args2, 60000);


                using(var StandardOutput = new StringReader(Res2.stdOut)) {
                    for(string line = StandardOutput.ReadLine(); line != null; line = StandardOutput.ReadLine()) {
                        if(line.StartsWith("exit code", StringComparison.InvariantCultureIgnoreCase)) {
                            int iRes = line.IndexOf(':');
                            string RestLine = line.Substring(iRes + 1);
                            exitcode = int.Parse(RestLine);
                            break;
                        } else {
                            continue;
                        }


                        //var parts = line.Split(new char[] { ' ', ' ', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                        //if (parts.Length < 3)
                        //    continue;
                        //if (parts[0].Equals("state", StringComparison.InvariantCultureIgnoreCase))
                        //{
                        //    state = Enum.Parse<JobState>(parts[2]);
                        //    break;
                        //}
                        //else
                        //{
                        //    continue;
                        //}
                    }
                }

            }

            return (state, exitcode);
        }


        /// <summary>
        /// Path to standard error file.
        /// </summary>
        public override string GetStderrFile(string idToken, string DeployDir) {
            string fp = Path.Combine(DeployDir, "stderr.txt");
            return fp;
        }
        /// <summary>
        /// Path to standard output file.
        /// </summary>
        public override string GetStdoutFile(string idToken, string DeployDir) {
            string fp = Path.Combine(DeployDir, "stdout.txt");
            return fp;
            
        }

        /// <summary>
        /// Submits the job to the Microsoft HPC server.
        /// </summary>
        public override (string id, object optJobObj) Submit(Job myJob, string DeploymentDirectory) {
            using (new FuncTrace()) {


                // write XML file
                // ==============
                string xmlCode = this.WriteJobXML(myJob, DeploymentDirectory);
                string xmlFilePath;
                {
                    xmlCode = this.WriteJobXML(myJob, DeploymentDirectory);
                    xmlFilePath = Path.Combine(DeploymentDirectory, "job.xml");
                    File.WriteAllText(xmlFilePath, xmlCode);
                }



                // submitt the job
                // ===============
                int id = -1;
                {
                    //var ret = CallJobCmd($"new /jobname:{JobName} /projectname:{PrjName} /scheduler:{server}  /numcores:{NumberOfCores}-{NumberOfCores} /exclusive:{SingleNode.ToString().ToLowerInvariant()} /priority:{Priority} /requestednodes:{Nodes}");

                    string pass = this.Password != null ? ("/password:" + this.Password) : "";
                    var ret = ExecuteProcess("job.exe", $"submit /jobfile:\"{xmlFilePath}\" /scheduler:{this.ServerName} /user:{this.Username} {pass}", 60000);

                    var parts = ret.stdOut.Split(new char[] { ' ', '.', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                    for (int i = 0; i < parts.Length; i++)
                    {
                        if (parts[i].Equals("id:", StringComparison.InvariantCultureIgnoreCase))
                            id = int.Parse(parts[i + 1]);

                    }
                    //Console.WriteLine(ret);
                }

                return (id.ToString(), null);

                /*
                ISchedulerJob MsHpcJob = null;
                ISchedulerTask task = null;

                // Create a job and add a task to the job.
                MsHpcJob = Scheduler.CreateJob();


                
                MsHpcJob.Name = myJob.Name;
                MsHpcJob.Project = PrjName;
                MsHpcJob.MaximumNumberOfCores = myJob.NumberOfMPIProcs;
                MsHpcJob.MinimumNumberOfCores = myJob.NumberOfMPIProcs;
                MsHpcJob.SingleNode = this.SingleNode;
                MsHpcJob.Priority = this.DefaultJobPriority;

                MsHpcJob.UserName = Username;

                task = MsHpcJob.CreateTask();
                task.MaximumNumberOfCores = myJob.NumberOfMPIProcs;
                task.MinimumNumberOfCores = myJob.NumberOfMPIProcs;
                
                task.WorkDirectory = DeploymentDirectory;

                using (var str = new StringWriter()) {
                    str.Write("mpiexec ");
                    if(!base.DotnetRuntime.IsEmptyOrWhite())
                        str.Write(base.DotnetRuntime + " ");
                    str.Write(Path.GetFileName(myJob.EntryAssembly.Location));
                    foreach (string arg in myJob.CommandLineArguments) {
                        str.Write(" ");
                        str.Write(arg);
                    }

                    task.CommandLine = str.ToString();
                }
                foreach (var kv in myJob.EnvironmentVars) {
                    string name = kv.Key;
                    string valu = kv.Value;
                    task.SetEnvironmentVariable(name, valu);
                }

                task.StdOutFilePath = Path.Combine(DeploymentDirectory, "stdout.txt");
                task.StdErrFilePath = Path.Combine(DeploymentDirectory, "stderr.txt");

                if (ComputeNodes != null) {
                    foreach (string node in ComputeNodes)
                        MsHpcJob.RequestedNodes.Add(node);
                }


                MsHpcJob.AddTask(task);

                // Start the job.
                Scheduler.SubmitJob(MsHpcJob, Username != null ? Username : null, Password);
                */

            }
        }

        /// <summary>
        /// 
        /// </summary>
        public override string ToString() {
            string NameString = "";
            if(!base.Name.IsEmptyOrWhite())
                NameString = " " + base.Name + " ";

            return $"MS HPC client {NameString}@{this.ServerName}, @{this.DeploymentBaseDirectory}";
        }

        /// <summary>
        /// 
        /// </summary>
        private (int id, JobState state)[] ListJobs()
        {
           

            string user = this.Username;
            string server = this.ServerName;

            var args = $"list /user:{user} /scheduler:{server} /format:list /state:All";
            var Res = ExecuteProcess("job.exe", args, 60000);
                        

            var states = new List<(int id, JobState state)>();
            using (var StandardOutput = new StringReader(Res.stdOut))
            {
                for (string line = StandardOutput.ReadLine(); line != null; line = StandardOutput.ReadLine())
                {
                    //Console.WriteLine("raw: " + line);

                    int id = -1;

                    //var parts = line.Split(new char[] { ' ', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                    //if (parts.Length < 3)
                    //    continue;
                    if (line.StartsWith("id", StringComparison.InvariantCultureIgnoreCase))
                    {
                        int iRes = line.IndexOf(':');
                        string RestLine = line.Substring(iRes + 1);
                        id = int.Parse(RestLine);
                    }
                    else
                    {
                        continue;
                    }

                    string st_string = null;
                    for (line = StandardOutput.ReadLine(); line != null; line = StandardOutput.ReadLine())
                    {
                        //var _parts = line.Split(new char[] { ' ', '\n', '\t', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                        //if (_parts.Length <= 0)
                        //    break;
                        //if (_parts.Length < 2)
                        //    continue;
                        //if (_parts[0].Equals("state", StringComparison.InvariantCultureIgnoreCase))
                        //    st_string = _parts[2];

                        if (line.StartsWith("state", StringComparison.InvariantCultureIgnoreCase))
                        {
                            int iRes = line.IndexOf(':');
                            string RestLine = line.Substring(iRes + 1);
                            st_string = RestLine;
                        }

                    }
                    if (st_string == null)
                        throw new IOException("unable to parse job list.");
                    JobState st = Enum.Parse<JobState>(st_string);


                    states.Add((id, st));

                }

                //foreach (var t in states)
                //    Console.WriteLine("parsed: " + t);
                return states.ToArray();
            }
        }


        string WriteJobXML(Job myJob, string DeploymentDirectory)
        {
            
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            string JobName = myJob.Name;


            //job modify 190848 /numcores:1 - 1
            int NumberOfCores = myJob.NumberOfMPIProcs;
            bool SingleNode = this.SingleNode;
            string UserName = this.Username;
            var Priority = this.DefaultJobPriority;

            string CommandLine;
            using (var str = new StringWriter())
            {
                str.Write("mpiexec ");
                if (!base.DotnetRuntime.IsEmptyOrWhite())
                    str.Write(base.DotnetRuntime + " ");
                str.Write(Path.GetFileName(myJob.EntryAssembly.Location));
                foreach (string arg in myJob.CommandLineArguments)
                {
                    str.Write(" ");
                    str.Write(arg);
                }

                CommandLine = str.ToString();
            }

            string user = this.Username;

            string WorkDirectory = DeploymentDirectory;
            string StdOutFilePath = Path.Combine(DeploymentDirectory, "stdout.txt");
            string StdErrFilePath = Path.Combine(DeploymentDirectory, "stderr.txt");


            string Nodes = "";
            if (this.ComputeNodes != null)
            {
                foreach (string node in ComputeNodes)
                {
                    if (Nodes.Length > 0)
                        Nodes = ",";
                    Nodes += node;
                }
            }
            

            bool exclusive = false;


            using (var stw = new StringWriter())
            {
                stw.WriteLine($"<?xml version=\"1.0\" encoding=\"utf-8\"?>");
                stw.WriteLine($"<Job Version=\"3.000\" ");
                stw.WriteLine($"     Name=\"{JobName}\" ");
                stw.WriteLine($"	 IsExclusive=\"{exclusive.ToString().ToLowerInvariant()}\" ");
                stw.WriteLine($"	 UnitType=\"Core\" ");
                stw.WriteLine($"	 Owner=\"{user}\" ");
                stw.WriteLine($"	 UserName=\"{user}\" ");
                stw.WriteLine($"	 Project=\"{PrjName}\" ");
                stw.WriteLine($"	 JobType=\"Batch\" ");
                stw.WriteLine($"	 SingleNode = \"{SingleNode.ToString().ToLowerInvariant()}\" ");
                stw.WriteLine($"	 JobTemplate=\"Default\" ");
                stw.WriteLine($"	 Priority=\"{Priority}\" ");
                if (this.ComputeNodes != null)
                    stw.WriteLine($"	 RequestedNodes=\"{Nodes}\" ");
                stw.WriteLine($"	 AutoCalculateMax=\"false\" ");
                stw.WriteLine($"	 AutoCalculateMin=\"false\" ");
                stw.WriteLine($"	 MinCores=\"{NumberOfCores}\" ");
                stw.WriteLine($"	 MaxCores=\"{NumberOfCores}\" ");
                stw.WriteLine($"	 xmlns=\"http://schemas.microsoft.com/HPCS2008R2/scheduler/\">");
                stw.WriteLine($"    <Dependencies />");
                stw.WriteLine($"    <Tasks>");
                stw.WriteLine($"        <Task Version=\"3.000\" ");
                stw.WriteLine($"			  UnitType=\"Core\" ");
                stw.WriteLine($"			  WorkDirectory=\"{WorkDirectory}\" ");
                stw.WriteLine($"			  NiceId=\"1\" ");
                stw.WriteLine($"			  CommandLine=\"{CommandLine}\" ");
                stw.WriteLine($"			  StdOutFilePath=\"{StdOutFilePath}\" ");
                stw.WriteLine($"			  StdErrFilePath=\"{StdErrFilePath}\" ");
                stw.WriteLine($"			  IsExclusive=\"{exclusive.ToString().ToLowerInvariant()}\" ");
                stw.WriteLine($"			  MinCores=\"{NumberOfCores}\" ");
                stw.WriteLine($"			  MaxCores=\"{NumberOfCores}\" ");
                stw.WriteLine($"			  Type=\"Basic\">");
                if (myJob.EnvironmentVars.Count > 0)
                {
                    stw.WriteLine($"            <EnvironmentVariables>");
                    foreach (var kv in myJob.EnvironmentVars)
                    {
                        stw.WriteLine($"                <Variable>");
                        stw.WriteLine($"                    <Name>{kv.Key}</Name>");
                        stw.WriteLine($"                    <Value>{kv.Value}</Value>");
                        stw.WriteLine($"                </Variable>");
                    }
                    stw.WriteLine($"			</EnvironmentVariables>");
                }
                stw.WriteLine($"        </Task>");
                stw.WriteLine($"    </Tasks>");
                stw.WriteLine($"</Job>");


                return stw.ToString();
            }

        }

        /// <summary>
        /// Synchronous wrapper around process execution, 
        /// see https://stackoverflow.com/questions/139593/processstartinfo-hanging-on-waitforexit-why
        /// </summary>
        (int exitcode, string stdOut, string stdErr) ExecuteProcess(string filename, string arguments, int timeout)
        {
            using (Process process = new Process())
            {
                StringBuilder output = new StringBuilder();
                StringBuilder error = new StringBuilder();

                process.StartInfo.FileName = filename;
                process.StartInfo.Arguments = arguments;
                process.StartInfo.UseShellExecute = false;
                process.StartInfo.RedirectStandardOutput = true;
                process.StartInfo.RedirectStandardError = true;

                using (AutoResetEvent outputWaitHandle = new AutoResetEvent(false))
                using (AutoResetEvent errorWaitHandle = new AutoResetEvent(false))
                {
                    process.OutputDataReceived += (sender, e) =>
                    {
                        if (e.Data == null)
                        {
                            outputWaitHandle.Set();
                        }
                        else
                        {
                            output.AppendLine(e.Data);
                        }
                    };
                    process.ErrorDataReceived += (sender, e) =>
                    {
                        if (e.Data == null)
                        {
                            errorWaitHandle.Set();
                        }
                        else
                        {
                            error.AppendLine(e.Data);
                        }
                    };

                    process.Start();

                    process.BeginOutputReadLine();
                    process.BeginErrorReadLine();

                    if (process.WaitForExit(timeout) &&
                        outputWaitHandle.WaitOne(timeout) &&
                        errorWaitHandle.WaitOne(timeout))
                    {
                        if (process.ExitCode != 0)
                        {
                            string modArguments = arguments;
                            modArguments = modArguments.Replace(this.Password, "***"); // make sure we don't send the password to stdout or some other log
                            throw new IOException(filename + " " + modArguments + " exited with code " + process.ExitCode + System.Environment.NewLine + output.ToString() + System.Environment.NewLine + error.ToString());
                        }


                        // Process completed. Check process.ExitCode here.
                        return (process.ExitCode, output.ToString(), error.ToString());
                    }
                    else
                    {
                        // Timed out.
                        string modArguments = arguments;
                        modArguments = modArguments.Replace(this.Password, "***"); // make sure we don't send the password to stdout or some other log
                        throw new IOException("timeout waiting for " + filename + " " + modArguments);
                    }
                }

            }
        }

    }
}
