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

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A <see cref="BatchProcessorClient"/> implementation for slurm systems on unix based hpc platforms
    /// </summary>
    public class SlurmClient : BatchProcessorClient {

        string m_Username;
        string m_Password;
        string m_ServerName;
        SshClient SSHConnection;     

        /// <summary>
        /// Client for submitting jobs directly from the BoSSSpad to slurm systems
        /// </summary>
        /// <param name="DeploymentBaseDirectory"></param>
        /// <param name="ServerName"></param>
        /// <param name="Username"></param>
        public SlurmClient(string DeploymentBaseDirectory, string ServerName, string Username = null) {
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            m_Username = Username;
            m_ServerName = ServerName;

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            Console.WriteLine();
            Console.WriteLine("Please enter your password...");
            m_Password = ReadPassword();
            Console.WriteLine("Connecting to "+ServerName+"...");
            Console.WriteLine();

            SSHConnection = new SshClient(m_ServerName, m_Username, m_Password);


            SSHConnection.Connect();
        }

        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            DeployDir = null;
            isRunning = false;
            wasSuccessful = false;
            isFailed = false;

            SshCommand output;

            if (myJob.EnvironmentVars.ContainsKey("JobID")) {
                output = SSHConnection.RunCommand("squeue -j " + myJob.EnvironmentVars["JobID"] + " -o %T");
                int startindex = output.Result.IndexOf("\n") ;
                int endindex = output.Result.IndexOf("\n",startindex+1);
                string jobstatus;
                if (startindex ==-1 || endindex==-1) {
                    jobstatus = "";
                }
                else {
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

        public override string GetStderrFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stderr.txt");
            if (File.Exists(fp)) {
                return fp;
            }
            else {
                return null;
            }

        }

        public override string GetStdoutFile(Job myJob) {
            string fp = Path.Combine(myJob.DeploymentDirectory, "stdout.txt");
            if (File.Exists(fp)) {
                return fp;
            }
            else {
                return null;
            }
        }

        public override object Submit(Job myJob) {

            // load users .bashrc with all dependencies
            buildSlurmScript(myJob, new string[] { "source " + "/home/" + m_Username + "/.bashrc"});

            string path = "\\home\\" + m_Username + myJob.DeploymentDirectory.Substring(2);
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

            string jobpath_win = "\\home\\" + m_Username + myJob.DeploymentDirectory.Substring(2);

            string jobpath_unix = jobpath_win.Replace("\\", "/");

            string jobname = myJob.Name;
            string executiontime = myJob.ExecutionTime;
            int MPIcores = myJob.NumberOfMPIProcs;
            string usermail = m_Username;
            string startupstring;
            string quote = "\"";

            using (var str = new StringWriter()) {
                str.Write("mpiexec mono ");
                str.Write(jobpath_unix+"/"+Path.GetFileName(myJob.EntryAssembly.Location));              
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
                sw.WriteLine("#SBATCH -o "+jobpath_unix + "/stdout.txt"); 
                sw.WriteLine("#SBATCH -e " + jobpath_unix + "/stderr.txt");
                sw.WriteLine("#SBATCH -t " + executiontime);
                sw.WriteLine("#SBATCH --mem-per-cpu=5000");
                sw.WriteLine("#SBATCH -n " + MPIcores);
                //sw.WriteLine("#SBATCH --mail-user= " + usermail);
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
                }
                else if (info.Key == ConsoleKey.Backspace) {
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

    }
}
