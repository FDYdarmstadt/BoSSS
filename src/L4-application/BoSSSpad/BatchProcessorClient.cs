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

using BoSSS.Foundation.IO;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {
    /// <summary>
    /// Abstraction for all kind of batch systems (aka. job managers).
    /// </summary>
    [DataContract]
    abstract public class BatchProcessorClient {


        /// <summary>
        /// Base directory where the executables should be deployed,
        /// accessible from the local machine (e.g. a mounted path if the batch processor deploys on another computer system)
        /// </summary>
        [DataMember]
        public string DeploymentBaseDirectory {
            get;
            protected set;
        }

        /// <summary>
        /// True, if it is required to deploy the runtime on the destination system,
        /// in order to execute a job.
        /// </summary>
        [DataMember]
        public bool DeployRuntime {
            get;
            set;
        }

        /// <summary>
        /// If not null, specifies paths to databases which are accessible to the computer system 
        /// on which this batch processor submits its jobs.
        /// This triggers data synchronization on job submission, if e.g. grid or restart timestep
        /// are in some other database.
        /// </summary>
        [DataMember]
        public string[] AllowedDatabasesPaths;


        List<IDatabaseInfo> m_AllowedDatabases;

        /// <summary>
        /// See <see cref="AllowedDatabasesPaths"/>
        /// </summary>
        public IReadOnlyList<IDatabaseInfo> AllowedDatabases {
            get {
                if (m_AllowedDatabases == null) {
                    m_AllowedDatabases = new List<IDatabaseInfo>();
                    if (AllowedDatabasesPaths != null) {

                        // fill up with empty entries
                        m_AllowedDatabases.AddRange(AllowedDatabasesPaths.Select(path => default(IDatabaseInfo)));

                        // pass 1: search in already-known databases
                        for (int i = 0; i < AllowedDatabasesPaths.Length; i++) {
                            foreach (var db in InteractiveShell.databases) {
                                if (db.PathMatch(AllowedDatabasesPaths[i])) {
                                    // bingo
                                    m_AllowedDatabases[i] = db;
                                    break;
                                }
                            }
                        }

                        // pass 2: try to open other databases
                        for (int i = 0; i < AllowedDatabasesPaths.Length; i++) {
                            if (m_AllowedDatabases[i] == null) {
                                try {
                                    m_AllowedDatabases[i] = new DatabaseInfo(AllowedDatabasesPaths[i]);
                                } catch (Exception e) {
                                    Console.Error.WriteLine($"Unable to open 'allowed database' for {this.ToString()} at path {AllowedDatabasesPaths[i]}. Check configuration file 'BatchProcessorConfig.json'. ({e.GetType().Name} : {e.Message})");
                                    Console.Error.WriteLine($"{this.ToString()} will continue to work, but database syncronization on job submission might not work correctly.");
                                }
                            }
                        }

                        // remove empty entries
                        for (int i = 0; i < m_AllowedDatabases.Count; i++) {
                            if (m_AllowedDatabases[i] == null) {
                                m_AllowedDatabases.RemoveAt(i);
                                i--;
                            }
                        }
                    }
                }

                return m_AllowedDatabases.AsReadOnly();
            }
        }

        string JobDirectoryBaseName(Job myJob) {
            string Exe = Path.GetFileNameWithoutExtension(myJob.EntryAssembly.Location);
            string Proj = InteractiveShell.WorkflowMgm.CurrentProject;
            string Sess = myJob.Name;

            return Proj
                //+ "-" + Sess 
                + "-" + Exe;
        }


        /// <summary>
        /// Returns the directory where the assemblies for <paramref name="myJob"/> 
        /// are deployed if <paramref name="myJob"/> is assigned to this batch processor.
        /// </summary>
        virtual public string GetNewDeploymentDir(Job myJob) {
            if (!Path.IsPathRooted(DeploymentBaseDirectory))
                throw new IOException($"Deployment base directory for {this.ToString()} must be rooted/absolute, but '{DeploymentBaseDirectory}' is not.");

            string ShortName = JobDirectoryBaseName(myJob);
            string DeployDir;
            int Counter = 0;
            do {
                string Suffix = Counter > 0 ? "-" + Counter : "";
                string DateNtime = DateTime.Now.ToString("yyyyMMMdd_HHmmss");
                DeployDir = Path.Combine(DeploymentBaseDirectory, ShortName + DateNtime + Suffix);
                Counter++;
            } while (Directory.Exists(DeployDir) == true);

            return DeployDir;
        }

        /// <summary>
        /// All deployment directories which potentially could match the job on the current batch processor.
        /// </summary>
        virtual public DirectoryInfo[] GetAllExistingDeployDirectories(Job myJob) {
            if (!Path.IsPathRooted(DeploymentBaseDirectory))
                throw new IOException($"Deployment base directory for {this.ToString()} must be rooted/absolute, but '{DeploymentBaseDirectory}' is not.");

            var jobControl = myJob.GetControl();
            if (jobControl == null)
                return null;

            // find all deployment directories relevant for project & job
            // ==========================================================
            string ShortName = JobDirectoryBaseName(myJob);

            var AllDirs = Directory.GetDirectories(DeploymentBaseDirectory, ShortName + "*");


            // filter appropriate ones 
            // =======================

            var filtDirs = new List<DirectoryInfo>();
            foreach (string dir in AllDirs) {
                string ControlObj = Path.Combine(dir, "control.obj");
                if (File.Exists(ControlObj)) {
                    var ctrl = BoSSS.Solution.Control.AppControl.Deserialize(File.ReadAllText(ControlObj));
                    if (InteractiveShell.WorkflowMgm.JobAppControlCorrelation(myJob, ctrl)) {
                        filtDirs.Add(new DirectoryInfo(dir));
                        continue;
                    }
                }

                string ControlScript = Path.Combine(dir, "control.cs");
                if (File.Exists(ControlScript)) {
                    int control_index = 0;
                    int i = myJob.CommandLineArguments.IndexWhere(arg => arg == "--pstudy_case");
                    if (i >= 0) {
                        control_index = int.Parse(myJob.CommandLineArguments[i + 1]);
                    }

                    var ctrl = BoSSS.Solution.Control.AppControl.FromFile(ControlScript, jobControl.GetType(), control_index);
                    if (InteractiveShell.WorkflowMgm.JobAppControlCorrelation(myJob, ctrl)) {
                        filtDirs.Add(new DirectoryInfo(dir));
                        continue;
                    }
                }
            }

            // return
            // ======
            return filtDirs.ToArray();
        }

        /// <summary>
        /// Submits the job to the batch system.
        /// </summary>
        /// <param name="myJob">Job to submit.</param>
        /// <returns>
        /// An optional identifier token (<see cref="Job.BatchProcessorIdentifierToken"/>).
        /// </returns>
        abstract public string Submit(Job myJob);

        /// <summary>
        /// Try to get some information about a job from the job manager.
        /// </summary>
        /// <param name="idToken">Identification within batch processor</param>
        /// <param name="DeployDir"></param>
        /// <param name="isRunning">
        /// True, if <paramref name="myJob"/> is currently running.
        /// </param>
        /// <param name="ExitCode">
        /// if <paramref name="isTerminated"/> is true, the exit code of the application.
        /// </param>
        /// <param name="isTerminated">
        /// True, if the application has exited
        /// </param>
        public abstract void EvaluateStatus(string idToken, string DeployDir, out bool isRunning, out bool isTerminated, out int ExitCode);

        /// <summary>
        /// Path to standard output file, if present - otherwise null.
        /// </summary>
        public abstract string GetStdoutFile(Job myJob);

        /// <summary>
        /// Path to standard error file, if present - otherwise null.
        /// </summary>
        public abstract string GetStderrFile(Job myJob);

        static bool IsNotSystemAssembly(Assembly Ass, string MainAssemblyDir){
            PlatformID CurrentSys = System.Environment.OSVersion.Platform;
            switch(CurrentSys)
            {
                case PlatformID.Unix:
                    { return Path.GetFullPath(Ass.Location).StartsWith("/home/"); }
                case PlatformID.Win32S:
                case PlatformID.Win32Windows:
                default:
                    {
                        return Path.GetDirectoryName(Ass.Location).Equals(MainAssemblyDir)
                            || Path.GetFileName(Ass.Location).StartsWith("BoSSS")
                            || Path.GetFileName(Ass.Location).StartsWith("ilPSP")
                            || !Ass.GlobalAssemblyCache;
                    }
            }
        }

        /// <summary>
        /// Copies the executable files to the <see cref="DeploymentBaseDirectory"/>, 
        /// but does not submit the job.
        /// </summary>
        virtual public void DeployExecuteables(Job myJob, IEnumerable<Tuple<byte[], string>> AdditionalFiles) {

            Console.WriteLine("Deploying executables and additional files ...");

            // Collect files
            List<string> files = new List<string>();
            {
                //string SystemPath = Path.GetDirectoryName(typeof(object).Assembly.Location);
                string MainAssemblyDir = Path.GetDirectoryName(myJob.EntryAssembly.Location);
                foreach (var a in myJob.AllDependentAssemblies) {
                    if (IsNotSystemAssembly(a, MainAssemblyDir)) {
                        files.Add(a.Location);
                    }
                }
            }

            {
                // test for really strange errors
                for(int i = 0; i < files.Count; i++) {
                    for(int j = i+1; j < files.Count; j++) {
                        if (Path.GetFileName(files[i]).Equals(Path.GetFileName(files[j])))
                            throw new ApplicationException("strange internal error");
                    }
                }
            }

            // create deployment directory.
            string DeployDir = myJob.DeploymentDirectory;
            CreateDirectoryWR(DeployDir);

            Console.WriteLine("Deployment directory: " + DeployDir);
            string OriginDir = null;

            // copy files
            foreach (var fOrg in files) {
                string fNmn = Path.GetFileName(fOrg);
                if(fNmn.Equals("mscorlib.dll"))
                    throw new ApplicationException("internal error - something went wrong during filtering.");

                string fTarget = Path.Combine(DeployDir, fNmn);
                
                CopyFileWR(fOrg, fTarget);
                if(OriginDir == null || !OriginDir.Equals(Path.GetDirectoryName(fOrg))) {
                    OriginDir = Path.GetDirectoryName(fOrg);
                }
            }
            Console.WriteLine("copied " + files.Count + " files.");

            // additional files
            if (AdditionalFiles != null) {
                foreach (var t in AdditionalFiles) {
                    string fTarget = Path.Combine(DeployDir, t.Item2);
                    byte[] Content = t.Item1;
                    WriteFileWR(fTarget, Content);
                    Console.WriteLine("   written file: " + t.Item2);
                }
            }

            // deploy runtime
            if (DeployRuntime) {
                string BosssInstall = BoSSS.Foundation.IO.Utils.GetBoSSSInstallDir();
                string BosssBinNative = Path.Combine(BosssInstall, "bin", Path.Combine("native", "win"));
                CopyDirectoryRec(BosssBinNative, DeployDir, "amd64");
                Console.WriteLine("   copied 'amd64' runtime.");
            }

            // finally
            Console.WriteLine("deployment finished.");

            // test
            TestWR(myJob);
        }

        private void TestWR(Job myJob) {

            Exception OP(int iTry) {
                if(myJob.GetControl() != null) {
                    var directories = this.GetAllExistingDeployDirectories(myJob);
                    if(directories == null || directories.Length <= 0) {
                        return new IOException("Job is assigned to batch processor, but no deployment directory can be found.");
                    }
                }
                return null;
            }

            RetryIOop(OP, "testing of job deployment", false);
        }


        /// <summary>
        /// Generic IO operation with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        internal static void RetryIOop(Func<int,Exception> op, string Message, bool SurpressException) {
            int MaxTry = 10;
            Random rnd = null;
            Exception latest = null;
            for(int iTry = 0; iTry < MaxTry; iTry++) {
                // on network file-systems, there seem to be some rare hiccups, sometimes:
                // hundreds of files copied successfully, suddenly an IOException: file already exists.
                // File indeed exists, but is empty -- makes no sense, since deploy directory is freshly created.

                try {

                    latest = op(iTry);
                    if(latest != null)
                        break;
                    else {
                        if(iTry > 0)
                            Console.WriteLine("success.");
                        return;
                    }
                } catch(IOException e) {
                    Console.Error.WriteLine(e.GetType().Name + " " + Message +  " : " + e.Message);
                    
                    Console.WriteLine($"Retrying {iTry + 1} of {MaxTry} (waiting for some time before) ...");
                    if(rnd == null)
                        rnd = new Random();
                    latest = e;
                    int iwait = rnd.Next(77 * 1000);
                    System.Threading.Thread.Sleep(iwait); // sleep for at most 77 seconds...
                }
            }

            if(SurpressException == false && latest != null)
                throw latest;
        }


        /// <summary>
        /// File write with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        private static void WriteFileWR(string fTarget, byte[] Content) {

            Exception OP(int iTry) {
                if(iTry == 0 && File.Exists(fTarget)) {
                    return new IOException("File '" + fTarget + "' already exists - wont over write.");
                }
                File.WriteAllBytes(fTarget, Content);
                return null;
            }

            RetryIOop(OP, "writing file '" + fTarget + "'", false);

            /*
            int MaxTry = 10;
            Random rnd = null;
            Exception latest = null;
            for(int iTry = 0; iTry < MaxTry; iTry++) {
                // on network file-systems, there seem to be some rare hiccups, sometimes:
                // hundreds of files copied successfully, suddenly an IOException: file already exists.
                // File indeed exists, but is empty -- makes no sense, since deploy directory is freshly created.

                try {
                    if(iTry == 0 && File.Exists(fTarget)) {
                        latest = new IOException("File '" + fTarget + "' already exists - wont over write.");
                        break;
                    }
                    File.WriteAllBytes(fTarget, Content);

                    if(iTry > 0) {
                        Console.WriteLine("success.");
                    }
                    return;
                } catch(IOException e) {
                    Console.Error.WriteLine(e.GetType().Name + " during writing of file '" + fTarget + "' : " + e.Message);
                    Console.WriteLine($"Retrying {iTry + 1} of {MaxTry} (waiting for some time before) ...");
                    if(rnd == null)
                        rnd = new Random();
                    latest = e;
                    int iwait = rnd.Next(77 * 1000);
                    System.Threading.Thread.Sleep(iwait); // sleep for at most 77 seconds...
                }
            }

            if(latest != null)
                throw latest;
            */
        }

        /// <summary>
        /// File copy with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        private static void CopyFileWR(string fOrg, string fTarget, bool SurpressException = false) {
            
            Exception op(int iTry) {
                File.Copy(fOrg, fTarget, iTry > 0);
                return null;
            }

            RetryIOop(op, " copy of file '" + fOrg + "' --> '" + fTarget + "'", SurpressException);

            

            /*
            int MaxTry = 10;
            Random rnd = null;
            Exception latest = null;
            for(int iTry = 0; iTry < MaxTry; iTry++) {
                // on network file-systems, there seem to be some rare hiccups, sometimes:
                // hundreds of files copied successfully, suddenly an IOException: file already exists.
                // File indeed exists, but is empty -- makes no sense, since deploy directory is freshly created.

                try {
                    File.Copy(fOrg, fTarget, iTry > 0);
                    if(iTry > 0) {
                        Console.WriteLine("success.");
                    }
                    return;
                } catch(IOException e) {
                    Console.Error.WriteLine(e.GetType().Name + " during copy of file '" + fOrg + "' --> '" + fTarget + "' : " + e.Message);
                    Console.WriteLine($"Retrying {iTry + 1} of {MaxTry} (waiting for some time before) ...");
                    if(rnd == null)
                        rnd = new Random();
                    latest = e;
                    int iwait = rnd.Next(77 * 1000);
                    System.Threading.Thread.Sleep(iwait); // sleep for at most 77 seconds...
                }
            }

            if(SurpressException == false && latest != null)
                throw latest;
            */
        }

        /// <summary>
        /// File copy with re-try (seems to be necessary when working with network file systems).
        /// </summary>
        private static void CreateDirectoryWR(string dstSubDir, bool SurpressException = false) {
            
            Exception op(int i) {
                if(!Directory.Exists(dstSubDir))
                    Directory.CreateDirectory(dstSubDir);
                return null;
            }

            RetryIOop(op, "creation of directory '" + dstSubDir + "'", SurpressException);

            /*
            int MaxTry = 10;
            Random rnd = null;
            Exception latest = null;
            for(int iTry = 0; iTry < MaxTry; iTry++) {
                // on network file-systems, there seem to be some rare hiccups, sometimes:
                // hundreds of files copied successfully, suddenly an IOException: file already exists.
                // File indeed exists, but is empty -- makes no sense, since deploy directory is freshly created.

                try {
                    if (!Directory.Exists(dstSubDir))
                        Directory.CreateDirectory(dstSubDir);

                    if(iTry > 0) {
                        Console.WriteLine("success.");
                    }
                    return;
                } catch(IOException e) {
                    Console.Error.WriteLine(e.GetType().Name + " during creation of directory '" + dstSubDir + "' : " + e.Message);
                    Console.WriteLine($"Retrying {iTry + 1} of {MaxTry} (waiting for some time before) ...");
                    if(rnd == null)
                        rnd = new Random();
                    latest = e;
                    int iwait = rnd.Next(77 * 1000);
                    System.Threading.Thread.Sleep(iwait); // sleep for at most 77 seconds...
                }
            }

            if(SurpressException == false && latest != null)
                throw latest;
            */
        }

        /// <summary>
        /// (tries to) do a recursive copy of a directory
        /// </summary>
        /// <param name="srcDir"></param>
        /// <param name="dstDir"></param>
        /// <param name="filter">search pattern/filter</param>
        public static void CopyDirectoryRec(string srcDir, string dstDir, string filter) {
            string[] srcFiles = Directory.GetFiles(srcDir);

            foreach (string srcFile in srcFiles) {
                //TryCopy(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)));
                CopyFileWR(srcFile, Path.Combine(dstDir, Path.GetFileName(srcFile)), SurpressException:true);
            }

            string[] subDirs;
            if (filter == null)
                subDirs = Directory.GetDirectories(srcDir);
            else
                subDirs = Directory.GetDirectories(srcDir, filter);
            foreach (string srcAbsDir in subDirs) {
                string srcRelDir = Path.GetFileName(srcAbsDir);
                string dstSubDir = Path.Combine(dstDir, srcRelDir);
                CreateDirectoryWR(dstSubDir);
                CopyDirectoryRec(srcAbsDir, dstSubDir, null);
            }
        }

    }
}

