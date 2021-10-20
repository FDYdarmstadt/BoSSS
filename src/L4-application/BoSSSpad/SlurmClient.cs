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
//using Renci.SshNet;
using System.IO;
using System.Runtime.Serialization;
using ilPSP;
using System.Diagnostics;
using ilPSP.Tracing;
using System.Linq;
//using System.Runtime.InteropServices.WindowsRuntime;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A <see cref="BatchProcessorClient"/> implementation for slurm systems on unix based hpc platforms
    /// </summary>
    [DataContract]
    public class SlurmClient : BatchProcessorClient {

        /// <summary>
        /// Username on the SSH server to connect to.
        /// </summary>
        [DataMember]
        public string Username {
            get;
            private set;
        }

        /// <summary>
        /// Non-recommended SSH password for authentication.
        /// This is not encrypted, the <see cref="PrivateKeyFilePath"/>
        /// </summary>
        [DataMember]
        public string Password {
            get;
            private set;
        }

        /// <summary>
        /// Server name or IP address
        /// </summary>
        [DataMember]
        public string ServerName {
            get;
            set;
        }

        /// <summary>
        /// Preferred SSH authentication method: path to private key file on local system
        /// </summary>
        [DataMember]
        public string PrivateKeyFilePath {
            get;
            set;
        }

        /// <summary>
        /// Additional lines to be added to the auto-generated batch.sh file.
        /// </summary>
        [DataMember]
        public string[] AdditionalBatchCommands {
            get;
            set;
        }

        /// <summary>
        /// Base directory where the executables should be deployed,
        /// i.e. the same location as <see cref="BatchProcessorClient.DeploymentBaseDirectory"/>,
        /// but in the file system of the remote computer on which Slurm is running.
        ///
        /// Example:
        ///  - <see cref="BatchProcessorClient.DeploymentBaseDirectory"/> is set to <tt>C:\serverSSFFSmount\jobdeploy</tt>
        ///  - <see cref="DeploymentBaseDirectoryAtRemote"/> is set to <tt>/home/linuxuser/jobdeploy</tt>
        /// </summary>
        [DataMember]
        public string DeploymentBaseDirectoryAtRemote {
            get;
            protected set;
        }

        string DeploymentDirectoryAtRemote(Job myJob, string DeploymentDirectory) {
            if (!DeploymentBaseDirectoryAtRemote.StartsWith("/")) {
                throw new IOException($"Deployment remote base directory for {this.ToString()} must be rooted/absolute, but '{DeploymentBaseDirectoryAtRemote}' is not.");
            }

            string RelativeDeploymentDirectory = DeploymentDirectory
                .Split(new[] { Path.DirectorySeparatorChar, Path.AltDirectorySeparatorChar }, StringSplitOptions.RemoveEmptyEntries)
                .Last()
                .TrimEnd('/', '\\');
            var tmp = DeploymentBaseDirectoryAtRemote;

            return tmp + "/" + RelativeDeploymentDirectory;
        }


        [NonSerialized]
        SshClient m_SSHConnection;

        SshClient SSHConnection {
            get {
                if (m_SSHConnection == null || m_SSHConnection.IsConnected == false) {
                    // SSHConnection = new SshClient(m_ServerName, m_Username, m_Password);
                    if (PrivateKeyFilePath != null) {
                        var pkf = new PrivateKeyFile(PrivateKeyFilePath);
                        m_SSHConnection = new SshClient(ServerName, Username, pkf);
                    } else if (Password != null) {
                        m_SSHConnection = new SshClient(ServerName, Username, Password);
                    } else if (Password == null) {
                        Console.WriteLine();
                        Console.WriteLine("Please enter your password...");
                        Password = ReadPassword();
                        m_SSHConnection = new SshClient(ServerName, Username, Password);
                    } else {
                        throw new NotSupportedException("Unable to initiate SSH connection -- either a password or private key file is required.");
                    }

                    //m_SSHConnection.Connect();
                }

                return m_SSHConnection;
            }
        }


        /// <summary>
        /// Empty constructor for de-serialization
        /// </summary>
        private SlurmClient() : base() {
        }

        /// <summary>
        /// runs an ls command
        /// </summary>
        public void TestSSH() {
            Console.WriteLine($"Performing test for ssh connection of {this.ToString()} ...");
            var output = SSHConnection.RunCommand("ls", verbose:true);
            //Console.WriteLine(output);
            Console.WriteLine($"Test finished.");
        }

        /// <summary>
        /// Client for submitting jobs directly from the BoSSSpad to slurm systems
        /// </summary>
        public SlurmClient(string DeploymentBaseDirectory, string ServerName, string Username, string PrivateKeyFilePath = null, bool AskForPassword = true) : base() {
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            this.Username = Username;
            this.ServerName = ServerName;
            this.PrivateKeyFilePath = PrivateKeyFilePath;

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            if (AskForPassword) {
                Console.WriteLine();
                Console.WriteLine("Please enter your password...");
                Password = ReadPassword();
                Console.WriteLine("Connecting to " + ServerName + "...");
                Console.WriteLine();
            }

            // test ssh connection
            var ssh = this.SSHConnection;
        }

        /// <summary>
        /// The number of the project where the job shall be executed (see HHLR-Antrag or csum, csreport)
        /// </summary>
        [DataMember]
        public string SlurmAccount {
            set;
            get;
        }

        /// <summary>
        /// If set, SLURM may send email notifications for the current job
        /// </summary>
        [DataMember]
        public string Email {
            set;
            get;
        }

        /// <summary>
        /// .
        /// </summary>
        public override (BoSSSpad.JobStatus,int? ExitCode) EvaluateStatus(string idToken, object optInfo, string DeployDir) { 
        //public override void EvaluateStatus(string idToken, object optInfo, string DeployDir, out bool isRunning, out bool isTerminated, out int ExitCode) {
            using (var tr = new FuncTrace()) {
                //string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
                //DeployDir = null;
                //isRunning = false;
                //wasSuccessful = false;
                //isFailed = false;
                //SubmitCount = 0;

                using (new BlockTrace("FILE_CHECK", tr)) {
                    string exitFile = Path.Combine(DeployDir, "exit.txt");
                    if (File.Exists(exitFile)) {

                        int ExitCode;
                        try {
                            ExitCode = int.Parse(File.ReadAllText(exitFile).Trim());
                        } catch (Exception) {
                            ExitCode = int.MinValue;
                        }
                        return (ExitCode == 0 ? JobStatus.FinishedSuccessful : JobStatus.FailedOrCanceled, ExitCode);
                    }

                    string runningFile = Path.Combine(DeployDir, "isrunning.txt");
                    if (File.Exists(runningFile)) {
                        // no decicion yet;
                        // e.g. assume that slurm terminated the Job after 24 hours => maybe 'isrunning.txt' is not deleted and 'exit.txt' does not exist

                    } else {
                        // no 'isrunning.txt'-token and no 'exit.txt' token: job should be pending in queue

                        return (JobStatus.PendingInExecutionQueue, null);

                        //isRunning = false;
                        //isTerminated = false;
                        //ExitCode = 0;
                        //return;
                    }
                }

                string JobID = idToken;

                using (new BlockTrace("SSH_SLURM_CHECK", tr)) {
                    //using (var output = SSHConnection.RunCommand("squeue -j " + JobID + " -o %T")) {
                    string output = SSHConnection.RunCommand("squeue -j " + JobID + " -o %T").stdout;
                    using(var Reader = new StringReader(output)) {

                        string line = Reader.ReadLine();
                        while(line != null && !line.Equals("state", StringComparison.InvariantCultureIgnoreCase))
                            line = Reader.ReadLine();

                        if(line == null || !line.Equals("state", StringComparison.InvariantCultureIgnoreCase))
                            return (JobStatus.Unknown, null);

                        string jobstatus = Reader.ReadLine();
                        if(jobstatus == null)
                            return (JobStatus.Unknown, null);

                        switch(jobstatus.ToUpperInvariant()) {
                            case "PENDING":
                            return (JobStatus.PendingInExecutionQueue, null);

                            case "RUNNING":
                            case "COMPLETING":
                            return (JobStatus.InProgress, null);

                            case "SUSPENDED":
                            case "STOPPED":
                            case "PREEMPTED":
                            case "FAILED":
                            return (JobStatus.FailedOrCanceled, int.MinValue);

                            case "":
                            case "COMPLETED":
                            // completed, but 'exit.txt' does not exist, something is shady here
                            return (JobStatus.FailedOrCanceled, -1);

                            default:
                            return (JobStatus.Unknown, null);
                        }
                        //}
                    }
                }
            }
        }


        /// <summary>
        /// Returns path to text-file for standard error stream
        /// </summary>
        public override string GetStderrFile(string idToken, string DeployDir) {
            string fp = Path.Combine(DeployDir, "stderr.txt");
            return fp;
        }

        /// <summary>
        /// Returns path to text-file for standard output stream
        /// </summary>
        public override string GetStdoutFile(string idToken, string DeployDir) {
            string fp = Path.Combine(DeployDir, "stdout.txt");
            return fp;
        }


        //void VerifyDatabases() {
        //    foreach (var db in this.AllowedDatabases) {
        //        if (db.AlternateDbPaths.Length <= 0) {
        //            throw new IOException("Missing 'AlternatePaths.txt' in database -- required for sshfs-mounted remote databases.");
        //        }
        //    }
        //}

        /// <summary>
        /// Sets debug flag for Mono
        /// </summary>
        public bool MonoDebug = false;

        /// <summary>
        ///
        /// </summary>
        public override (string id, object optJobObj) Submit(Job myJob, string DeploymentDirectory) {
            using (new FuncTrace()) {
                //VerifyDatabases();


                // load users .bashrc with all dependencies
                buildSlurmScript(myJob, new string[] { "source " + "/home/" + Username + "/.bashrc" }, DeploymentDirectory);

                string jobId = SSHConnection.SubmitJob(DeploymentDirectoryAtRemote(myJob, DeploymentDirectory));
                if(jobId.IsEmptyOrWhite())
                    throw new ApplicationException("missing job id return value from slurm command.");

                return (jobId, null);
            }
        }

        /// <summary>
        /// build batch script with all necessary parameters
        /// </summary>
        void buildSlurmScript(Job myJob, string[] moduleLoad, string DeploymentDirectory) {

            //string jobpath_win = "\\home\\" + Username + myJob.DeploymentDirectory.Substring(2);
            //string jobpath_unix = jobpath_win.Replace("\\", "/");
            string jobpath_unix = DeploymentDirectoryAtRemote(myJob, DeploymentDirectory);

            string jobname = myJob.Name;
            string executiontime = myJob.ExecutionTime;
            int MPIcores = myJob.NumberOfMPIProcs;
            string userName = Username;
            string startupstring;
            string quote = "\"";
            string slurmAccount = this.SlurmAccount;
            //string memPerCPU = "5000";
            //if (myJob.MemPerCPU != null) {
            //    memPerCPU = myJob.MemPerCPU;
            //} else {
            //    memPerCPU = "5000";
            //}
            string email = Email;

            using (var str = new StringWriter()) {
                str.Write($"mpiexec {base.DotnetRuntime} ");
                if (MonoDebug) { 
                    str.Write("-v --debug "); 
                }
                str.Write(jobpath_unix + "/" + Path.GetFileName(myJob.EntryAssembly.Location));
                str.Write(" ");
                str.Write(myJob.EnvironmentVars["BOSSS_ARG_" + 0]);
                str.Write(" ");

                // How the controlfile is handled (serialized or compiled at runtime)
                if (myJob.EnvironmentVars["BOSSS_ARG_1"].Equals("control.obj")) {
                    str.Write(jobpath_unix + "/" + myJob.EnvironmentVars["BOSSS_ARG_1"]);
                } else {
                    str.Write(quote + myJob.EnvironmentVars["BOSSS_ARG_" + 1] + quote);
                }

                startupstring = str.ToString();
            }

            string path = Path.Combine(DeploymentDirectory, "batch.sh");

            using (StreamWriter sw = File.CreateText(path)) {
                sw.NewLine = "\n"; // Unix file endings

                sw.WriteLine("#!/bin/sh");
                sw.WriteLine("#SBATCH -J " + jobname);
                if (slurmAccount != null) {
                    sw.WriteLine("#SBATCH -A " + slurmAccount);
                }
                sw.WriteLine("#SBATCH -o " + jobpath_unix + "/stdout.txt");
                sw.WriteLine("#SBATCH -e " + jobpath_unix + "/stderr.txt");
                sw.WriteLine("#SBATCH -t " + executiontime);
                //sw.WriteLine("#SBATCH --mem-per-cpu=" + myJob.MemPerCPU);
                if (myJob.UseComputeNodesExclusive) {
                    sw.WriteLine("#SBATCH --exclusive");
                }

                sw.WriteLine("#SBATCH -n " + MPIcores);
                if (!email.IsEmptyOrWhite()) {
                    sw.WriteLine("#SBATCH --mail-user=" + email);
                    sw.WriteLine("#SBATCH --mail-type=ALL");
                }
                foreach (var cmd in this.AdditionalBatchCommands ?? Enumerable.Empty<string>()) {
                    sw.WriteLine(cmd);
                }
                //sw.WriteLine("#SBATCH --ntasks-per-node 1");    // Only start one MPI-process per node

                // Load modules
                foreach (string arg in moduleLoad) {
                    sw.WriteLine(arg);
                }

                // Set startupstring
                string RunningToken = DeploymentDirectoryAtRemote(myJob, DeploymentDirectory) + "/isrunning.txt";
                sw.WriteLine($"touch '{RunningToken}'");
                sw.WriteLine("cd " + DeploymentDirectoryAtRemote(myJob, DeploymentDirectory)); // this ensures that any files written out (e.g. .plt-files) are placed in the deployment directory rather than ~
                sw.WriteLine(startupstring);
                sw.WriteLine("echo $? > '" + DeploymentDirectoryAtRemote(myJob, DeploymentDirectory) + "/exit.txt'");
                sw.WriteLine($"rm '{RunningToken}'");
                if (this.DotnetRuntime == "mono") {
                    sw.WriteLine("echo delete mono-crash-dumps, if there are any...");
                    sw.WriteLine($"rm core.*");
                    sw.WriteLine($"rm mono_crash.*");
                    sw.WriteLine($"rm mono_crash.mem.*");
                }
            }

        }

        /// <summary>
        /// Read in log in password for HPC computing system
        /// </summary>
        /// <returns></returns>
        public static string ReadPassword() {
            string password = "";
            ConsoleKeyInfo info = Console.ReadKey(true);
            while (info.Key != ConsoleKey.Enter) {
                if (info.Key != ConsoleKey.Backspace) {
                    Console.Write("*");
                    password += info.KeyChar;
                } else if (info.Key == ConsoleKey.Backspace) {
                    if (!string.IsNullOrEmpty(password)) {
                        // remove one character from the list of password characters
                        password = password.Substring(0, password.Length - 1);
                        // get the location of the cursor
                        int pos = Console.CursorLeft;
                        // move the cursor to the left by one character
                        Console.SetCursorPosition(pos - 1, Console.CursorTop);
                        // replace it with space
                        Console.Write(" ");
                        // move the cursor to the left by one character again
                        Console.SetCursorPosition(pos - 1, Console.CursorTop);
                    }
                }
                info = Console.ReadKey(true);
            }
            // add a new line because user pressed enter at the end of their password
            Console.WriteLine();
            return password;
        }

        /// <summary>
        ///
        /// </summary>
        public override string ToString() {

            string NameString = "";
            if(!base.Name.IsEmptyOrWhite())
                NameString = " " + base.Name + " ";

            return "SlurmClient" + NameString + ": " + Username + "@" + ServerName + ", Slurm account: " + (SlurmAccount ?? "NONE");
        }

    }
}
