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
using System.Threading.Tasks;
using Renci.SshNet;
using System.IO;
using System.Security;
using System.Runtime.Serialization;
using ilPSP;

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
        /// Non-recommended SSH password for authentification.
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

        SshClient m_SSHConnection;

        SshClient SSHConnection {
            get {
                if(m_SSHConnection == null || m_SSHConnection.IsConnected == false) {
                    // SSHConnection = new SshClient(m_ServerName, m_Username, m_Password);
                    if(PrivateKeyFilePath != null) {
                        var pkf = new PrivateKeyFile(PrivateKeyFilePath);
                        m_SSHConnection = new SshClient(ServerName, Username, pkf);
                    } else if(Password == null) {
                        m_SSHConnection = new SshClient(ServerName, Username, Password);
                    } else {
                        throw new NotSupportedException("Unable to initiate SSH connection -- either a password or private key file is required.");
                    }

                    m_SSHConnection.Connect();
                }

                return m_SSHConnection;
            }
        }



        /*
        /// <summary>
        /// Configuration options specific to the <see cref="SlurmClient"/>
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
            public string PrivateKeyFilePath;

            /// <summary>
            /// 
            /// </summary>
            public string SlurmAccount;

            /// <summary>
            /// 
            /// </summary>
            public string Email;


            /// <summary>
            /// %
            /// </summary>
            public override BatchProcessorClient Instance() {
                var r = new SlurmClient(
                    base.DeploymentBaseDirectory,
                    ServerName,
                    Username,
                    PrivateKeyFilePath);
                r.SlurmAccount = SlurmAccount;
                r.Email = Email;

                return r;
            }
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override BatchProcessorClient.Config GetConfig() {
            return new SlurmClient.Config() {
                DeploymentBaseDirectory = this.DeploymentBaseDirectory,
                DeployRuntime = this.DeployRuntime,
                PrivateKeyFilePath = this.m_PrivateKeyFilePath,
                ServerName = this.m_ServerName,
                Username = this.m_Username,
                Email = this.Email,
                SlurmAccount = this.SlurmAccount
            };
        }
        */

        /// <summary>
        /// Empty constructor for de-serialization
        /// </summary>
        private SlurmClient() {
        }

        /// <summary>
        /// runs an ls command 
        /// </summary>
        public void TestSSH() {
            var output = SSHConnection.RunCommand("ls");
            Console.WriteLine(output.Result);
        }

        /// <summary>
        /// Client for submitting jobs directly from the BoSSSpad to slurm systems
        /// </summary>
        public SlurmClient(string DeploymentBaseDirectory, string ServerName, string Username, string PrivateKeyFilePath, bool AskForPassword = true) {
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            this.Username = Username;
            this.ServerName = ServerName;
            this.PrivateKeyFilePath = PrivateKeyFilePath;

            if(!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            if(AskForPassword) {
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
        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            DeployDir = null;
            isRunning = false;
            wasSuccessful = false;
            isFailed = false;

            SshCommand output;

            if (myJob.EnvironmentVars.ContainsKey("JobID")) {
                output = SSHConnection.RunCommand("squeue -j " + myJob.EnvironmentVars["JobID"] + " -o %T");
                int startindex = output.Result.IndexOf("\n");
                int endindex = output.Result.IndexOf("\n", startindex + 1);
                string jobstatus;
                if (startindex == -1 || endindex == -1) {
                    jobstatus = "";
                } else {
                    jobstatus = output.Result.Substring(startindex + 1, (endindex - startindex) - 1);
                }

                switch (jobstatus) {
                    case "RUNNING":
                    case "PENDING":
                    case "COMPLETING":
                        isRunning = true;
                        break;

                    case "":
                        wasSuccessful = true;
                        break;

                    case "FAILED":
                        isFailed = true;
                        break;

                    default:
                        throw new NotImplementedException("Unknown job state: " + jobstatus);
                }
            }


            SubmitCount = 0;
        }

        /// <summary>
        /// Returns path to text-file for standard error stream
        /// </summary>
        public override string GetStderrFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stderr.txt");
            return fp;
        }

        /// <summary>
        /// Returns path to text-file for standard output stream
        /// </summary>
        public override string GetStdoutFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stdout.txt");
            return fp;
        }

        public override object Submit(Job myJob) {

            // load users .bashrc with all dependencies
            buildSlurmScript(myJob, new string[] { "source " + "/home/" + Username + "/.bashrc" });

            string path = "\\home\\" + Username + myJob.DeploymentDirectory.Substring(2);
            // Converting script to unix format
            string convertCmd = " dos2unix " + path + "\\batch.sh";

            // Submitting script to sbatch system
            string sbatchCmd = " sbatch " + path + "\\batch.sh";

            // Convert from Windows to Unix and submit job
            Console.WriteLine();
            var result1 = SSHConnection.RunCommand(convertCmd.Replace("\\", "/"));
            var result2 = SSHConnection.RunCommand(sbatchCmd.Replace("\\", "/"));

            // Otherwise it didn´t work because uploading speed at some clusters is too slow
            if (result1.Error == "" || result2.Result == "") {
                Console.Write("Waiting for file transfer to finish");
                while (result1.Error == "" || result2.Result == "") {
                    Console.Write(".");
                    System.Threading.Thread.Sleep(10000);
                    result1 = SSHConnection.RunCommand(convertCmd.Replace("\\", "/"));
                    result2 = SSHConnection.RunCommand(sbatchCmd.Replace("\\", "/"));
                }
                Console.WriteLine();
            }

            Console.WriteLine(result1.Error);
            Console.WriteLine(result2.Result);

            // Hardcoded extract of JobID
            myJob.EnvironmentVars.Add("JobID", result2.Result.Substring(20, 7));

            return null;
        }

        /// <summary>
        /// build batch script with all necessary parameters
        /// </summary>
        /// <param name="myJob"></param>
        /// <param name="moduleLoad"></param>
        public void buildSlurmScript(Job myJob, string[] moduleLoad) {

            string jobpath_win = "\\home\\" + Username + myJob.DeploymentDirectory.Substring(2);

            string jobpath_unix = jobpath_win.Replace("\\", "/");

            string jobname = myJob.Name;
            string executiontime = myJob.ExecutionTime;
            int MPIcores = myJob.NumberOfMPIProcs;
            string userName = Username;
            string startupstring;
            string quote = "\"";
            string HHLR_project = this.SlurmAccount;
            string memPerCPU;
            if (myJob.MemPerCPU != null) {
                memPerCPU = myJob.MemPerCPU;
            } else {
                memPerCPU = "5000";
            }
            string email = Email;

            using (var str = new StringWriter()) {
                str.Write("mpiexec mono ");
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

            string path = myJob.DeploymentDirectory + "\\batch.sh";

            using (StreamWriter sw = File.CreateText(path)) {
                sw.WriteLine("#!/bin/sh");
                sw.WriteLine("#SBATCH -J " + jobname);
                if (HHLR_project != null) {
                    sw.WriteLine("#SBATCH -A " + HHLR_project);
                }
                sw.WriteLine("#SBATCH -o " + jobpath_unix + "/stdout.txt");
                sw.WriteLine("#SBATCH -e " + jobpath_unix + "/stderr.txt");
                sw.WriteLine("#SBATCH -t " + executiontime);
                sw.WriteLine("#SBATCH --mem-per-cpu=" + memPerCPU);
                if (myJob.UseComputeNodesExclusive) {
                    sw.WriteLine("#SBATCH --exclusive");
                }

                sw.WriteLine("#SBATCH -n " + MPIcores);
                if (!email.IsEmptyOrWhite()) {
                    sw.WriteLine("#SBATCH --mail-user=" + email);
                    sw.WriteLine("#SBATCH --mail-type=ALL");
                }
                sw.WriteLine("#SBATCH -C avx");
                //sw.WriteLine("#SBATCH --ntasks-per-node 1");    // Only start one MPI-process per node

                // Load modules
                foreach (string arg in moduleLoad) {
                    sw.WriteLine(arg);
                }

                // Set startupstring
                sw.WriteLine(startupstring);
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
            return "SlurmClient: " + Username + "@" + ServerName + ", Slurm account: " + (SlurmAccount ?? "NONE");
        }

    }
}
