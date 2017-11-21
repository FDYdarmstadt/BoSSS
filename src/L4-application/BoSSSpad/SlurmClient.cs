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
using System.Text;
using System.Threading.Tasks;
using Renci.SshNet;
using System.IO;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// A <see cref="BatchProcessorClient"/> implementation for slurm systems on unix based hpc platforms
    /// </summary>
    public class SlurmClient : BatchProcessorClient {

        string m_Username;
        string m_Password;
        string m_ServerName;
        SshClient SSHConnection;

        public SlurmClient(string DeploymentBaseDirectory, string ServerName, string Username = null, string Password = null) {
            base.DeploymentBaseDirectory = DeploymentBaseDirectory;
            m_Username = Username;
            m_Password = Password;
            m_ServerName = ServerName;

            if (!Directory.Exists(base.DeploymentBaseDirectory))
                Directory.CreateDirectory(base.DeploymentBaseDirectory);

            SSHConnection = new SshClient(m_ServerName, m_Username, m_Password);

            SSHConnection.Connect();
        }

        public override void EvaluateStatus(Job myJob, out int SubmitCount, out bool isRunning, out bool wasSuccessful, out bool isFailed, out string DeployDir) {
            string PrjName = InteractiveShell.WorkflowMgm.CurrentProject;
            DeployDir = null;
            isRunning = false;
            wasSuccessful = false;
            isFailed = false;



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

            buildSlurmScript(myJob, new string[] { "module load gcc", "module load openmpi/hcc/1.6.5", "module load acml" });

            string path = "\\home\\" + m_Username + myJob.DeploymentDirectory.Substring(2);
            string convertCmd = " dos2unix " + path + "\\batch.sh";
            string sbatchCmd = " sbatch "+path + "\\batch.sh";

            // Otherwise it didn´t work
            System.Threading.Thread.Sleep(3000);

            // Convert from Windows to Unix and submit job
            Console.WriteLine();
            var result1 = SSHConnection.RunCommand(convertCmd.Replace("\\", "/"));
            Console.WriteLine(result1.Error);
            var result2 = SSHConnection.RunCommand(sbatchCmd.Replace("\\", "/"));
            Console.WriteLine(result2.Result);

            // Hardcoded extract of JobID
            myJob.EnvironmentVars.Add("JobID", result2.Result.Substring(20,7));

            return null;
        }

        public void buildSlurmScript(Job myJob, string[] moduleLoad) {

            string jobname = myJob.Name;
            string stdoutpath = myJob.Stdout;
            string stderrpath = myJob.Stderr;
            string executiontime = "03:00:00";
            int MPIcores = myJob.NumberOfMPIProcs;
            string usermail = m_Username;
            //string solverdirection = myJob.Solver;
            string startupstring;
            string quote = "\"";

            using (var str = new StringWriter()) {
                str.Write("mpiexec mono ");
                str.Write(Path.GetFileName(myJob.EntryAssembly.Location));
                //for (int i = 0; i<myJob.EnvironmentVars.Count;i++) {
                    str.Write(" ");
                    str.Write(myJob.EnvironmentVars["BOSSS_ARG_"+0]);
                str.Write(" ");
                
                str.Write(quote + myJob.EnvironmentVars["BOSSS_ARG_" +1]+ quote);
                //}
                startupstring = str.ToString();
            }

            string path = myJob.DeploymentDirectory + "\\batch.sh";

            using (StreamWriter sw = File.CreateText(path)) {
                sw.WriteLine("#!/bin/sh");
                sw.WriteLine("#SBATCH -J " + jobname);
                sw.WriteLine("#SBATCH -o " + stdoutpath);
                sw.WriteLine("#SBATCH -e " + stderrpath);
                sw.WriteLine("#SBATCH -t " + executiontime);
                sw.WriteLine("#SBATCH --mem-per-cpu=1750");
                sw.WriteLine("#SBATCH -n " + MPIcores);
                //sw.WriteLine("#SBATCH --mail-user= " + usermail);
                sw.WriteLine("#SBATCH -C avx");

                // Load modules
                foreach (string arg in moduleLoad) {
                    sw.WriteLine(arg);
                }

                // Set startupstring
                sw.WriteLine(startupstring);
            }





        }

    }
}
